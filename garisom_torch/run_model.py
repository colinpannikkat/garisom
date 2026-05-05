#!/usr/bin/env python3
"""
GARISOM-Torch: Main Entry Point

Example script demonstrating how to run the PyTorch implementation
of the GARISOM (Carbon Gain vs Hydraulic Risk Stomatal Optimization Model).

This script shows:
    1. Loading model parameters from CSV
    2. Loading climate forcing data
    3. Running the model for multiple timesteps
    4. Demonstrating gradient computation (differentiability)
    5. Saving results

Usage:
    python run_model.py --params ../03_test_data/parameters.csv \
                        --climate ../03_test_data/climate_forcing_files/dataset_short.csv \
                        --output results.csv
"""

import argparse
import sys
from pathlib import Path

import torch
import numpy as np

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from garisom_torch import (
    Plant,
    load_parameters,
    load_climate_data,
    load_configuration,
    save_results,
    tensor,
)


def run_model(params_path: str,
              climate_path: str,
              output_path: str = None,
              num_layers: int = 5,
              max_timesteps: int = None,
              compute_gradients: bool = False):
    """
    Run the GARISOM-Torch model.
    
    Args:
        params_path: Path to parameters CSV file
        climate_path: Path to climate forcing data CSV
        output_path: Path for output CSV (optional)
        num_layers: Number of soil layers
        max_timesteps: Maximum timesteps to run (None = all)
        compute_gradients: Whether to compute gradients (for demonstration)
    """
    print("=" * 60)
    print("GARISOM-Torch: Differentiable Plant Physiology Model")
    print("=" * 60)
    
    # Load parameters
    print(f"\nLoading parameters from: {params_path}")
    params = load_parameters(params_path)
    
    print(f"  Species: {params.get('species', 'Unknown')}")
    print(f"  Site: {params.get('site', 'Unknown')}")
    print(f"  LAI: {params.get('lai', 'N/A')}")
    print(f"  Height: {params.get('height', 'N/A')} m")
    print(f"  Texture: {params.get('texture', 'loam')}")
    
    # Get soil texture
    texture = params.pop('texture', 'loam')
    
    # Remove non-numeric parameters
    for key in ['species', 'site']:
        params.pop(key, None)
    
    # Create model
    print(f"\nInitializing model with {num_layers} soil layers...")
    model = Plant(num_layers=num_layers, texture=texture)
    
    # Set parameters
    model.set_parameters(params)
    print(f"  Model: {model}")
    
    # Load climate data
    print(f"\nLoading climate data from: {climate_path}")
    climate = load_climate_data(climate_path)
    print(f"  Timesteps: {len(climate)}")
    print(f"  Year range: {climate.year[0].item()} - {climate.year[-1].item()}")
    
    # Limit timesteps if requested
    n_steps = len(climate) if max_timesteps is None else min(max_timesteps, len(climate))
    print(f"  Running {n_steps} timesteps")
    
    # Initialize results storage
    results = {
        'transpiration': [],
        'assimilation': [],
        'conductance': [],
        'predawn': [],
        'midday': [],
        'leaf_temp': [],
        'ppfd_sun': [],
        'is_active': [],
    }
    
    # Initialize soil water content (at field capacity)
    soil_water = torch.ones(num_layers + 1) * 0.3  # ~30% VWC
    
    # Run model
    print("\nRunning simulation...")
    
    # If computing gradients, make certain parameters require gradients
    if compute_gradients:
        # Example: make Vcmax trainable
        model.params.v_max25.requires_grad = True
        model.params.j_max25.requires_grad = True
        print("  Gradient computation enabled for v_max25, j_max25")
    
    total_transpiration = tensor(0.0)
    total_assimilation = tensor(0.0)
    active_timesteps = 0
    
    for t in range(n_steps):
        # Get forcing for this timestep
        forcing = climate[t]
        
        # Run model forward pass
        output = model(
            solar=forcing['solar'],
            vpd=forcing['vpd'],
            tair=forcing['tair'],
            wind=forcing['wind'],
            julian_day=forcing['julian_day'],
            hour=forcing['hour'],
            soil_water=soil_water,
        )
        
        # Store results
        for key in results:
            if key in output:
                val = output[key]
                if isinstance(val, torch.Tensor):
                    results[key].append(val.detach().item())
                else:
                    results[key].append(float(val))
        
        # Update soil water (simple extraction)
        if output['is_active']:
            model.update_soil_water(output['transpiration'], timestep_hours=1.0)
            total_transpiration = total_transpiration + output['transpiration']
            total_assimilation = total_assimilation + output['assimilation']
            active_timesteps += 1
        
        # Progress indicator
        if (t + 1) % 100 == 0 or t == n_steps - 1:
            print(f"  Step {t+1}/{n_steps} | "
                  f"E={output['transpiration'].item():.3f} | "
                  f"A={output['assimilation'].item():.2f} | "
                  f"Active: {active_timesteps}")
    
    print("\n" + "-" * 60)
    print("Simulation Complete")
    print("-" * 60)
    print(f"  Total active timesteps: {active_timesteps}")
    print(f"  Total transpiration: {total_transpiration.item():.2f} kg/m² BA")
    print(f"  Total assimilation: {total_assimilation.item():.2f} µmol/m² leaf")
    
    # Demonstrate gradient computation
    if compute_gradients:
        print("\n" + "-" * 60)
        print("Gradient Demonstration")
        print("-" * 60)
        
        # Create a loss (e.g., negative total assimilation)
        loss = -total_assimilation
        
        # Backward pass
        loss.backward()
        
        print(f"  ∂Loss/∂v_max25 = {model.params.v_max25.grad.item():.6f}")
        print(f"  ∂Loss/∂j_max25 = {model.params.j_max25.grad.item():.6f}")
        print("\n  These gradients show how changing Vcmax and Jmax")
        print("  affects total carbon assimilation!")
    
    # Save results
    if output_path:
        # Convert lists to tensors
        results_tensors = {k: torch.tensor(v) for k, v in results.items()}
        save_results(results_tensors, output_path, climate)
    
    return model, results


def demo_optimization():
    """
    Demonstrate parameter optimization using PyTorch autograd.
    """
    print("\n" + "=" * 60)
    print("Parameter Optimization Demo")
    print("=" * 60)
    
    # Create a simple model
    model = Plant(num_layers=3, texture='loam')
    
    # Make some parameters trainable
    model.params.v_max25 = torch.nn.Parameter(tensor(50.0))
    model.params.j_max25 = torch.nn.Parameter(tensor(100.0))
    
    # Create optimizer
    optimizer = torch.optim.Adam([model.params.v_max25, model.params.j_max25], lr=1.0)
    
    # Target assimilation rate
    target_a = 15.0  # µmol/m²/s
    
    # Simple test conditions
    solar = tensor(800.0)
    vpd = tensor(1.5)
    tair = tensor(25.0)
    wind = tensor(2.0)
    julian_day = tensor(180)
    hour = tensor(12.0)
    
    print(f"\nTarget assimilation: {target_a} µmol/m²/s")
    print(f"Initial v_max25: {model.params.v_max25.item():.1f}")
    print(f"Initial j_max25: {model.params.j_max25.item():.1f}")
    
    # Optimization loop
    print("\nOptimizing...")
    for i in range(20):
        optimizer.zero_grad()
        
        # Forward pass
        output = model(
            solar=solar, vpd=vpd, tair=tair, wind=wind,
            julian_day=julian_day, hour=hour
        )
        
        # Loss: squared error from target
        loss = (output['assimilation'] - target_a) ** 2
        
        # Backward pass
        loss.backward()
        
        # Update parameters
        optimizer.step()
        
        # Clamp to reasonable ranges
        with torch.no_grad():
            model.params.v_max25.clamp_(10, 200)
            model.params.j_max25.clamp_(20, 400)
        
        if (i + 1) % 5 == 0:
            print(f"  Step {i+1}: A={output['assimilation'].item():.2f}, "
                  f"v_max25={model.params.v_max25.item():.1f}, "
                  f"j_max25={model.params.j_max25.item():.1f}, "
                  f"loss={loss.item():.4f}")
    
    print(f"\nFinal v_max25: {model.params.v_max25.item():.1f}")
    print(f"Final j_max25: {model.params.j_max25.item():.1f}")
    print(f"Final assimilation: {output['assimilation'].item():.2f} µmol/m²/s")


def main():
    parser = argparse.ArgumentParser(
        description='GARISOM-Torch: Differentiable Plant Physiology Model'
    )
    
    parser.add_argument(
        '--params', '-p',
        default='../03_test_data/parameters.csv',
        help='Path to parameters CSV file'
    )
    
    parser.add_argument(
        '--climate', '-c',
        default='../03_test_data/climate_forcing_files/dataset_short.csv',
        help='Path to climate forcing data CSV'
    )
    
    parser.add_argument(
        '--output', '-o',
        default=None,
        help='Path for output CSV file'
    )
    
    parser.add_argument(
        '--layers', '-l',
        type=int,
        default=5,
        help='Number of soil layers'
    )
    
    parser.add_argument(
        '--max-steps', '-m',
        type=int,
        default=100,
        help='Maximum timesteps to run'
    )
    
    parser.add_argument(
        '--gradients', '-g',
        action='store_true',
        help='Demonstrate gradient computation'
    )
    
    parser.add_argument(
        '--optimize-demo',
        action='store_true',
        help='Run parameter optimization demo'
    )
    
    args = parser.parse_args()
    
    # Check file paths
    params_path = Path(args.params)
    climate_path = Path(args.climate)
    
    if not params_path.exists():
        print(f"Error: Parameters file not found: {params_path}")
        sys.exit(1)
    
    if not climate_path.exists():
        print(f"Error: Climate file not found: {climate_path}")
        sys.exit(1)
    
    # Run model
    run_model(
        params_path=str(params_path),
        climate_path=str(climate_path),
        output_path=args.output,
        num_layers=args.layers,
        max_timesteps=args.max_steps,
        compute_gradients=args.gradients,
    )
    
    # Optionally run optimization demo
    if args.optimize_demo:
        demo_optimization()


if __name__ == '__main__':
    main()
