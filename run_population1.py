"""
Run GARISOM-Torch with Fremont Poplar data for Population 1 (CCR-COL).
Collects leaf temperature gradients at each timestep and plots outputs.
"""

import torch
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

# Add the package to path
import sys
sys.path.insert(0, str(Path(__file__).parent))

from garisom_torch import Plant, tensor
from garisom_torch.constants import get_device, get_dtype

# Set up paths - use DBG folder
data_dir = Path(__file__).parent / "03_test_data" / "fremont-poplar-data" / "DBG"
params_file = data_dir / "parameters.csv"
climate_file = data_dir / "dataset.csv"

print("="*70)
print("GARISOM-Torch: Running Population 1 (CCR-COL)")
print("="*70)

# Load parameters for population 1
params_df = pd.read_csv(params_file)
pop1 = params_df.iloc[0]  # First row = population 1

print(f"\nSpecies: {pop1['i_sp']}")
print(f"Region: {pop1['i_region']}")
print(f"Site: {pop1['i_site']}")
print(f"LAI: {pop1['i_leafAreaIndex']}")
print(f"Height: {pop1['i_height']} m")
print(f"Vmax25: {pop1['i_vmax25']} µmol/m²/s")
print(f"Jmax25: {pop1['i_jmax25']} µmol/m²/s")

# Load climate data
climate_df = pd.read_csv(climate_file)
print(f"\nClimate data: {len(climate_df)} timesteps")
print(f"Days: {climate_df['Day'].min()} to {climate_df['Day'].max()}")

# Create model with population 1 parameters
num_layers = int(pop1['i_layers'])
texture = pop1['i_texture'] if pd.notna(pop1['i_texture']) else 'loam'

model = Plant(num_layers=max(1, num_layers), texture='loam')  # Default to loam

# Set parameters from CSV
model.params.lai = torch.nn.Parameter(tensor(float(pop1['i_leafAreaIndex'])), requires_grad=False)
model.params.height = torch.nn.Parameter(tensor(float(pop1['i_height'])), requires_grad=False)
model.params.lat = torch.nn.Parameter(tensor(float(pop1['i_latitude']) * np.pi / 180), requires_grad=False)
model.params.leaf_width = torch.nn.Parameter(tensor(float(pop1['i_leafWidth'])), requires_grad=False)

# Make Vmax25 trainable to compute gradients
model.params.v_max25 = torch.nn.Parameter(tensor(float(pop1['i_vmax25'])))
model.params.j_max25 = torch.nn.Parameter(tensor(float(pop1['i_jmax25'])))

# Set Weibull parameters if available
if pd.notna(pop1['i_cl']) and pd.notna(pop1['i_bl']):
    model.xylem.leaf.b = torch.nn.Parameter(tensor(float(pop1['i_bl'])), requires_grad=False)
    model.xylem.leaf.c = torch.nn.Parameter(tensor(float(pop1['i_cl'])), requires_grad=False)

if pd.notna(pop1['i_cs']) and pd.notna(pop1['i_bs']):
    model.xylem.stem.b = torch.nn.Parameter(tensor(float(pop1['i_bs'])), requires_grad=False)
    model.xylem.stem.c = torch.nn.Parameter(tensor(float(pop1['i_cs'])), requires_grad=False)

print(f"\nModel: {model}")

# Storage for outputs
results = {
    'day': [],
    'hour': [],
    'datetime': [],
    'solar': [],
    'tair': [],
    'vpd': [],
    'transpiration': [],
    'assimilation': [],
    'conductance': [],
    'predawn': [],
    'midday': [],
    'leaf_temp': [],
    'leaf_temp_grad_vmax': [],  # Gradient of leaf temp wrt vmax25
}

print("\nRunning simulation (batched)...")

# Convert climate data to tensors for batched processing
n_timesteps = len(climate_df)
solar_all = tensor(climate_df['Solar_Wm2'].values)
vpd_all = tensor(climate_df['D_kPa'].values)
tair_all = tensor(climate_df['Tair_C'].values)
wind_all = tensor(np.maximum(0.5, climate_df['Wind_ms.1'].values))
julian_day_all = tensor(climate_df['Day'].values)
hour_all = tensor(climate_df['Hour'].values)

# Process in batches for speed, but model needs per-timestep currently
# For now, vectorize the storage and run sequentially but efficiently
print(f"  Processing {n_timesteps} timesteps...")

# Pre-allocate numpy arrays for speed
transpiration_arr = np.zeros(n_timesteps)
assimilation_arr = np.zeros(n_timesteps)
conductance_arr = np.zeros(n_timesteps)
predawn_arr = np.zeros(n_timesteps)
midday_arr = np.zeros(n_timesteps)
leaf_temp_arr = np.zeros(n_timesteps)
leaf_temp_grad_arr = np.zeros(n_timesteps)

import time
start_time = time.time()

# Run with no_grad for speed, compute gradients only at select timesteps
grad_timesteps = list(range(12, n_timesteps, 24))  # Noon each day

with torch.no_grad():
    for idx in range(n_timesteps):
        output = model(
            solar=solar_all[idx],
            vpd=vpd_all[idx],
            tair=tair_all[idx],
            wind=wind_all[idx],
            julian_day=julian_day_all[idx],
            hour=hour_all[idx],
        )
        
        transpiration_arr[idx] = output['transpiration'].item()
        assimilation_arr[idx] = output['assimilation'].item()
        conductance_arr[idx] = output['conductance'].item()
        predawn_arr[idx] = output['predawn'].item()
        midday_arr[idx] = output['midday'].item()
        leaf_temp_arr[idx] = output['leaf_temp'].item()
        
        if idx % 500 == 0:
            elapsed = time.time() - start_time
            rate = (idx + 1) / elapsed if elapsed > 0 else 0
            print(f"  Timestep {idx}/{n_timesteps} ({rate:.1f} steps/sec)")

# Now compute gradients at select timesteps
print(f"\nComputing gradients at {len(grad_timesteps)} timesteps...")
for idx in grad_timesteps:
    model.zero_grad()
    output = model(
        solar=solar_all[idx],
        vpd=vpd_all[idx],
        tair=tair_all[idx],
        wind=wind_all[idx],
        julian_day=julian_day_all[idx],
        hour=hour_all[idx],
    )
    if output['leaf_temp'].requires_grad:
        output['leaf_temp'].backward()
        if model.params.v_max25.grad is not None:
            leaf_temp_grad_arr[idx] = model.params.v_max25.grad.item()

elapsed = time.time() - start_time
print(f"\nSimulation complete: {n_timesteps} timesteps in {elapsed:.1f}s ({n_timesteps/elapsed:.1f} steps/sec)")

# Build results DataFrame - handle potential NaNs in input
results_df = pd.DataFrame({
    'day': climate_df['Day'].fillna(0).values.astype(int),
    'hour': climate_df['Hour'].fillna(0).values.astype(int),
    'solar': climate_df['Solar_Wm2'].values,
    'tair': climate_df['Tair_C'].values,
    'vpd': climate_df['D_kPa'].values,
    'transpiration': transpiration_arr,
    'assimilation': assimilation_arr,
    'conductance': conductance_arr,
    'predawn': predawn_arr,
    'midday': midday_arr,
    'leaf_temp': leaf_temp_arr,
    'leaf_temp_grad_vmax': leaf_temp_grad_arr,
})
results_df['datetime'] = results_df['day'].astype(str) + ':' + results_df['hour'].apply(lambda x: f'{x:02d}')

# Summary statistics
print("\n" + "="*70)
print("Summary Statistics")
print("="*70)
print(f"Transpiration: {results_df['transpiration'].mean():.4f} ± {results_df['transpiration'].std():.4f} kg/hr/m² BA")
print(f"Assimilation: {results_df['assimilation'].mean():.4f} ± {results_df['assimilation'].std():.4f} µmol/m²/s")
print(f"Leaf temp: {results_df['leaf_temp'].mean():.2f} ± {results_df['leaf_temp'].std():.2f} °C")
print(f"Air temp: {results_df['tair'].mean():.2f} ± {results_df['tair'].std():.2f} °C")

# Count non-zero gradients
nonzero_grads = (results_df['leaf_temp_grad_vmax'] != 0).sum()
print(f"\nLeaf temp gradients computed: {nonzero_grads}/{len(results_df)} timesteps")

# Create figure with subplots
fig, axes = plt.subplots(4, 2, figsize=(14, 16))
fig.suptitle(f"GARISOM-Torch: Population 1 (POFR CCR-COL)\nVmax25={pop1['i_vmax25']:.1f}, LAI={pop1['i_leafAreaIndex']:.1f}", 
             fontsize=14, fontweight='bold')

# Create x-axis as continuous hours
x = np.arange(len(results_df))
x_days = results_df['day'].values

# 1. Solar radiation and air temperature
ax1 = axes[0, 0]
ax1.plot(x, results_df['solar'], 'orange', label='Solar (W/m²)', alpha=0.7)
ax1.set_ylabel('Solar Radiation (W/m²)', color='orange')
ax1.tick_params(axis='y', labelcolor='orange')
ax1b = ax1.twinx()
ax1b.plot(x, results_df['tair'], 'red', label='Tair (°C)', alpha=0.7)
ax1b.set_ylabel('Air Temperature (°C)', color='red')
ax1b.tick_params(axis='y', labelcolor='red')
ax1.set_title('Environmental Forcing')
ax1.set_xlabel('Hour')

# 2. VPD
ax2 = axes[0, 1]
ax2.plot(x, results_df['vpd'], 'purple', alpha=0.7)
ax2.set_ylabel('VPD (kPa)')
ax2.set_xlabel('Hour')
ax2.set_title('Vapor Pressure Deficit')

# 3. Transpiration
ax3 = axes[1, 0]
ax3.plot(x, results_df['transpiration'], 'blue', alpha=0.7)
ax3.set_ylabel('Transpiration (kg/hr/m² BA)')
ax3.set_xlabel('Hour')
ax3.set_title('Transpiration')
ax3.axhline(y=0, color='gray', linestyle='--', alpha=0.3)

# 4. Assimilation
ax4 = axes[1, 1]
ax4.plot(x, results_df['assimilation'], 'green', alpha=0.7)
ax4.set_ylabel('Assimilation (µmol/m²/s)')
ax4.set_xlabel('Hour')
ax4.set_title('Carbon Assimilation')
ax4.axhline(y=0, color='gray', linestyle='--', alpha=0.3)

# 5. Leaf temperature vs air temperature
ax5 = axes[2, 0]
ax5.plot(x, results_df['tair'], 'red', label='Air temp', alpha=0.5)
ax5.plot(x, results_df['leaf_temp'], 'darkgreen', label='Leaf temp', alpha=0.7)
ax5.set_ylabel('Temperature (°C)')
ax5.set_xlabel('Hour')
ax5.set_title('Leaf vs Air Temperature')
ax5.legend()

# 6. Leaf temperature gradient wrt Vmax25
ax6 = axes[2, 1]
ax6.plot(x, results_df['leaf_temp_grad_vmax'], 'brown', alpha=0.7)
ax6.set_ylabel('∂Tleaf/∂Vmax25')
ax6.set_xlabel('Hour')
ax6.set_title('Leaf Temp Gradient w.r.t. Vmax25')
ax6.axhline(y=0, color='gray', linestyle='--', alpha=0.3)

# 7. Conductance
ax7 = axes[3, 0]
ax7.plot(x, results_df['conductance'], 'cyan', alpha=0.7)
ax7.set_ylabel('Conductance (mmol/m²/s)')
ax7.set_xlabel('Hour')
ax7.set_title('Stomatal Conductance')

# 8. Water potential
ax8 = axes[3, 1]
ax8.plot(x, results_df['predawn'], 'navy', label='Predawn', alpha=0.7)
ax8.plot(x, results_df['midday'], 'darkred', label='Midday', alpha=0.7)
ax8.set_ylabel('Water Potential (MPa)')
ax8.set_xlabel('Hour')
ax8.set_title('Plant Water Potential')
ax8.legend()

plt.tight_layout()
plt.savefig('population1_output.png', dpi=150, bbox_inches='tight')
print(f"\nPlot saved to: population1_output.png")

# Show plot
plt.show()

# Save results to CSV
results_df.to_csv('population1_results.csv', index=False)
print(f"Results saved to: population1_results.csv")

print("\n" + "="*70)
print("Done!")
print("="*70)
