import subprocess
import os
from tempfile import TemporaryDirectory
from concurrent.futures import ThreadPoolExecutor, as_completed
from argparse import ArgumentParser
import json
import uuid
from collections import defaultdict
from tqdm import tqdm

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from SALib.sample import sobol as ssobol
from SALib.analyze import sobol as asobol
from sklearn.metrics import mean_squared_error, r2_score, root_mean_squared_error, mean_absolute_percentage_error, median_absolute_error

def run_single_model(
        X: np.ndarray, 
        in_names: list, 
        out_names: list,
        proc_num: int, 
        tmp_dir: str,
        params: pd.DataFrame,
        POP_NUM: int,
        CONFIG_FILE: str,
        MODEL_DIR: str
    ) -> list[np.ndarray]:
    
    TMP_DIR = f"{tmp_dir}/{proc_num}"
    TMP_PARAM_FILE = f"{TMP_DIR}/params.csv"

    os.makedirs(TMP_DIR, exist_ok=True)

    for i, name in enumerate(in_names):
        params.at[POP_NUM - 1, name] = X[i]

    params.to_csv(TMP_PARAM_FILE, index=False)

    out = subprocess.DEVNULL

    p = subprocess.run(
        [
            "./run",
            TMP_PARAM_FILE,
            CONFIG_FILE,
            str(POP_NUM),
            TMP_DIR
        ],
        cwd=MODEL_DIR,
        stdout=out,
        stderr=out
    )

    if p.returncode != 0:
        raise RuntimeError(f"Subprocess {i} failed with exit code {p.returncode}")
    
    species = params.at[POP_NUM - 1, 'i_sp']
    region = params.at[POP_NUM - 1, 'i_region']
    site = params.at[POP_NUM - 1, 'i_site']

    output_file = os.path.join(TMP_DIR, f"timesteps_output_{species}_{region}_{site}.csv")
    if not os.path.exists(output_file):
        raise FileNotFoundError(f"Expected output file not found: {output_file}")
    
    output_file = pd.read_csv(output_file)
    
    out = output_file[out_names].to_numpy(dtype=float)  # T x Y_D

    return out

def wrapped_garisom(
        X: np.ndarray, 
        in_names: list,
        out_names: list, 
        MAX_WORKERS: int,
        MODEL_DIR: str,
        params: pd.DataFrame,
        POP_NUM: int,
        CONFIG_FILE: str,
        T: int
    ) -> np.ndarray:

    N, D = X.shape
    Y_D = len(out_names)
    print(f"{N} samples")
    res = np.empty((T, N, Y_D))

    with TemporaryDirectory() as tmp:
        print('created temporary directory', tmp)
        pbar = tqdm(total=N)

        with ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
            futures = {
                executor.submit(
                    run_single_model, 
                    X[i, :], 
                    in_names,
                    out_names, 
                    i, 
                    tmp,
                    params,
                    POP_NUM,
                    CONFIG_FILE,
                    MODEL_DIR
                ): i for i in range(N)
            }

            for future in as_completed(futures):
                pbar.update(1)
                idx = futures[future]
                try:
                    out = future.result()
                    res[:, idx, :] = out
                except Exception as e:
                    print(f"Subprocess for index {idx} failed: {e}") 
                    res[:, idx, :] = np.nan

        pbar.close()

    return res

def plot(
        first_order_indices, # shape: (Y_D, D, T)
        second_order_indices,
        outputs,             # shape: (T, N, Y_D)
        param_values, 
        plt_dir,
        res_dir,
        T,
        problem,
        POP_NUM
    ):

    # Plot settings
    start_day = problem['plot_settings']['start_day']
    end_day = problem['plot_settings']['end_day']
    average = problem['plot_settings']['average'] # average over measurement periods
    metric = problem['plot_settings']['metric']

    # Measurement periods in the day (rounded-down and up)
    time_offsets = {
        "GW": {
            "am": (7, 9),
            "pm": (15, 17)
        },
        "E-MD": {
            "am": (7, 9),
            "pm": (15, 17)
        },
        "K-plant": {
            "am": (7, 9),
            "pm": (15, 17)
        },
        "P-PD": {
            "am": (3, 5),
            "pm": (3, 5)
        },
        "P-MD": {
            "am": (13, 15),
            "pm": (13, 15)
        }
    }
        
    ground = None
    match POP_NUM:
        case 1:
            ground = pd.read_csv("data/ccr_hourly_data.csv")
        case 2:
            ground = pd.read_csv("data/jla_hourly_data.csv")
        case 3:
            ground = pd.read_csv("data/tsz_hourly_data.csv")
        case 4:
            ground = pd.read_csv("data/nrv_hourly_data.csv")
        case _:
            raise Exception("Incorrect POP_NUM!")
    
    time = np.arange(T)

    for idx in range(0, len(problem['outputs']), 2):
        fig, axes = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

        for i, name in enumerate(problem['names']):
            axes[0].scatter(time, first_order_indices[idx, i, :], label=f'{name}')
            if idx + 1 < len(problem['outputs']):
                axes[1].scatter(time, first_order_indices[idx + 1, i, :], label=f'{name}')

        axes[0].set_title(f'First-Order Sobol Indices for {problem["outputs"][idx]}')
        if idx + 1 < len(problem['outputs']):
            axes[1].set_title(f'First-Order Sobol Indices for {problem["outputs"][idx + 1]}')

        for ax in axes:
            ax.set_ylabel('Sobol Index')
            ax.legend()
            ax.grid(True)

        axes[1].set_xlabel('Time')
        plt.tight_layout()

        if idx + 1 < len(problem["outputs"]):
            plt.savefig(f"{plt_dir}/first-order_{problem['outputs'][idx]}_{problem['outputs'][idx + 1]}.png")
        else:
            plt.savefig(f"{plt_dir}/first-order_{problem['outputs'][idx]}.png")

    for idx, output_name in enumerate(problem['outputs']):
        fig, axes = plt.subplots(len(problem['names']), 1, figsize=(10, 8), sharex=True)

        for i, name in enumerate(problem['names']):
            plt.figure(figsize=(10, 6))
            plt.scatter(param_values[:, i], outputs[:, :, idx].mean(axis=0), label=f'{output_name}')
            plt.title(f'Mean Output ({output_name}) vs {name}')
            plt.xlabel(name)
            plt.ylabel('Mean Output')
            plt.legend()
            plt.grid(True)

            plt.tight_layout()
            plt.savefig(f"{plt_dir}/mean_output_{output_name}_vs_{name}.png")
            plt.close(fig)

    def cmp_pred_to_ground_metrics(n_ground, n_pred):

        fits = defaultdict(list)

        for ground, pred in zip(n_ground, n_pred):
            
            mse = mean_squared_error(ground, pred)
            rmse = root_mean_squared_error(ground, pred)
            mape = mean_absolute_percentage_error(ground, pred)
            made = median_absolute_error(ground, pred)
            r2 = r2_score(ground, pred)

            fits['mse'].append(mse)
            fits['rmse'].append(rmse)
            fits['mape'].append(mape)
            fits['made'].append(made)
            fits['r2'].append(r2)

        return fits

    errors = {}
    for idx, output_name in enumerate(problem['outputs']):

        if output_name not in time_offsets:
            raise ValueError(f"Unknown output name: {output_name}")

        # Filter ground data based on julian-day and drop NaN values
        col_ground = ground[ground['julian-day'].between(start_day, end_day)][output_name].dropna()

        # Align predictions with the filtered ground data
        col_pred = outputs[:, :, idx]  # (N, T)
        col_pred = pd.DataFrame(col_pred)
        pred_values = col_pred.loc[col_ground.index].T.to_numpy()

        # Get N copies of ground values
        ground_values = np.array([col_ground.copy().to_numpy() for _ in range(outputs[:, :, idx].shape[1])])

        errors[output_name] = cmp_pred_to_ground_metrics(ground_values, pred_values)

    # Save errors
    with open(os.path.join(res_dir, "errors.json"), "w") as f:
        json.dump(errors, f, indent=4)

    # Plot and save the errors for each output in errors
    for output_name, error_values in errors.items():
        for i, param in enumerate(problem['names']):
            plt.figure(figsize=(10, 6))
            plt.scatter(param_values[:, i], error_values[metric], marker='o', label=f'{metric}')
            plt.title(f'Error Metrics Across Samples for {output_name}')
            plt.xlabel(f'{param}')
            # plt.ylabel('Mean Absolute Percent Error (MAPE)')
            plt.grid(True)
            plt.legend()
            plt.tight_layout()
            plt.savefig(f"{plt_dir}/error_{metric}_{param}_{output_name}.png")
            plt.close()

def main():
    parser = ArgumentParser()
    parser.add_argument("--input", "-i", help="File path to sensitivity analysis input file.", required=True, type=str)
    parser.add_argument("--output", "-o", help="Directory path for output files.", default=".", type=str)
    parser.add_argument("--model", "-m", help="Directory path to model executable.",required=True, type=str)
    parser.add_argument("--workers", "-w", help="Number of workers to run SA with.", default=4, type=int)
    parser.add_argument("--samples", "-s", help="Number of samples to run SA with, must be a multiple of two.", default=2**4, type=int)
    parser.add_argument("--pop", "-p", help="Population number to run for model.", default=1, type=int)

    args = parser.parse_args()

    PARAM_FILE = "./DBG/parameters.csv"
    CONFIG_FILE = "../03_test_data/fremont-poplar-data/DBG/configuration.csv"
    DATA_FILE = "./DBG/dataset.csv"
    POP_NUM = args.pop
    MAX_WORKERS = args.workers
    SAMPLES = args.samples
    MODEL_DIR = os.path.join(args.model)
    OUT_DIR = os.path.join(args.output)
    RAND_DIR = str(uuid.uuid4())
    PLT_DIR = os.path.join(OUT_DIR, "plots", RAND_DIR)
    RES_DIR = os.path.join(OUT_DIR, "sa_results", RAND_DIR)

    print("Plots will be saved in: ", PLT_DIR)
    print("Results will be saved in: ", RES_DIR)

    os.makedirs(PLT_DIR, exist_ok=True)
    os.makedirs(RES_DIR, exist_ok=True)

    params = pd.read_csv(PARAM_FILE)
    dataset = pd.read_csv(DATA_FILE)

    T = dataset.shape[0]

    # Loads 'D' parameters for SA
    problem = None
    with open(args.input, "+r") as f:
        problem = json.load(f)

    # Save the problem definition to a JSON file in the results directory
    problem_file = os.path.join(RES_DIR, "problem.json")
    with open(problem_file, "w") as f:
        json.dump(problem, f, indent=4)

    D = problem['num_vars']
    Y_D = len(problem['outputs'])

    param_values = ssobol.sample(problem, N=SAMPLES) # shape : (N, D)

    np.save(f"{RES_DIR}/sample.npy", param_values)

    print("Running model with samples.")
    Y = wrapped_garisom(    # returns shape: (T, N, Y_D)
        param_values, 
        problem['names'],
        problem['outputs'],
        MAX_WORKERS,
        MODEL_DIR,
        params,
        POP_NUM,
        CONFIG_FILE,
        T
    )

    np.save(f"{RES_DIR}/model_output.npy", Y)

    first_order_indices = np.empty((Y_D, D, T))
    second_order_indices = np.empty((Y_D, D, T))

    print("Analyzing indices.")
    for t in range(T):
        for i, out in enumerate(problem['outputs']):
            si = asobol.analyze(problem, Y[t, :, i], print_to_console=False)
            first_order_indices[i, :, t] = np.clip(si['S1'], 0, 1)
            # second_order_indices[i, :, t] = np.clip(si['S2'], 0, 1)

    np.save(f"{RES_DIR}/first_order_indices.npy", first_order_indices)

    print("Creating plots.")
    plot(
        first_order_indices,
        second_order_indices,
        Y,
        param_values,
        PLT_DIR,
        RES_DIR,
        T,
        problem,
        POP_NUM
    )

    print("Finished.")

if __name__ == "__main__":
    main()