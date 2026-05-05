"""
Compare GARISOM-Torch outputs with original C++ GARISOM model.

Key differences to note:
1. C++ uses positive water potentials, PyTorch uses negative (suction convention)
2. C++ has full water balance (rainfall, drainage), PyTorch has simplified water
3. C++ uses discrete optimization, PyTorch uses differentiable soft-argmax
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

print("="*70)
print("GARISOM: PyTorch vs C++ vs Ground Truth Comparison")
print("="*70)

# Load PyTorch results
torch_results = pd.read_csv('population1_results.csv')
print(f"\nPyTorch results: {len(torch_results)} timesteps")

# Load C++ results
cpp_results = pd.read_csv('02_program_code/timesteps_output_POFR_CCR-COL_DBG.csv')
print(f"C++ results: {len(cpp_results)} timesteps")

# Load ground truth data
ground_truth = pd.read_csv('03_test_data/fremont-poplar-data/data/ground/ccr_hourly_data.csv')
print(f"Ground truth (hourly): {len(ground_truth)} timesteps")

# Load leaf temperature ground truth (more complete)
leaftemp_gt = pd.read_csv('03_test_data/fremont-poplar-data/data/ground/ccr_leaftemp.csv')
# Average across loggers for each timestep
leaftemp_gt_avg = leaftemp_gt.groupby(['year', 'julian-day', 'standard-time']).agg({
    'leaftemp': 'mean',
    'Tair_C': 'mean'
}).reset_index()
leaftemp_gt_avg['day_hour'] = leaftemp_gt_avg['julian-day'].astype(str) + '_' + leaftemp_gt_avg['standard-time'].astype(str)
print(f"Ground truth leaftemp: {len(leaftemp_gt_avg)} unique timesteps (from {leaftemp_gt['logger.id'].nunique()} loggers)")

# Count non-empty ground truth measurements
gt_gw_count = ground_truth['GW'].notna().sum()
gt_emd_count = ground_truth['E-MD'].notna().sum()
gt_leaftemp_count = ground_truth['leaftemp'].notna().sum()
print(f"  - GW measurements: {gt_gw_count}")
print(f"  - E-MD measurements: {gt_emd_count}")
print(f"  - sparse leaftemp: {gt_leaftemp_count}")
print(f"  - continuous leaftemp: {len(leaftemp_gt_avg)}")

# NOTE: C++ uses positive water potentials, PyTorch uses negative
print("\nSign convention: C++ positive pressures, PyTorch negative (suction)")

# Align by day and hour
torch_results['day_hour'] = torch_results['day'].astype(str) + '_' + torch_results['hour'].astype(str)
cpp_results['day_hour'] = cpp_results['julian-day'].astype(str) + '_' + cpp_results['standard-time'].astype(str)
ground_truth['day_hour'] = ground_truth['julian-day'].astype(str) + '_' + ground_truth['standard-time'].astype(str)

# Rename C++ columns before merge to be explicit
cpp_results = cpp_results.rename(columns={
    'leaftemp': 'leaftemp_cpp',
    'E-MD': 'E-MD_cpp',
    'P-MD': 'P-MD_cpp',
    'P-PD': 'P-PD_cpp',
    'GW': 'GW_cpp'
})

# Merge PyTorch and C++ results
merged = pd.merge(
    torch_results, 
    cpp_results, 
    on='day_hour', 
    suffixes=('_torch', '_cpp'),
    how='inner'
)

# Merge with ground truth (sparse measurements)
merged = pd.merge(
    merged,
    ground_truth[['day_hour', 'GW', 'E-MD', 'P-PD', 'P-MD', 'leaftemp', 'leaf-air-temp-diff']],
    on='day_hour',
    how='left',
)
# Rename ground truth columns to avoid confusion
merged = merged.rename(columns={
    'GW': 'GW_gt', 
    'E-MD': 'E-MD_gt', 
    'P-PD': 'P-PD_gt', 
    'P-MD': 'P-MD_gt',
    'leaftemp': 'leaftemp_sparse_gt',
    'leaf-air-temp-diff': 'leaf-air-temp-diff_gt'
})

# Merge with continuous leaf temperature ground truth
merged = pd.merge(
    merged,
    leaftemp_gt_avg[['day_hour', 'leaftemp', 'Tair_C']],
    on='day_hour',
    how='left',
    suffixes=('', '_leaftemp_gt')
)
merged = merged.rename(columns={
    'leaftemp': 'leaftemp_gt',
    'Tair_C': 'Tair_gt'
})

print(f"Matched timesteps: {len(merged)}")
print(f"Timesteps with leaf temp ground truth: {merged['leaftemp_gt'].notna().sum()}")

# Convert PyTorch predawn to positive for comparison (flip sign)
merged['predawn_pos'] = -merged['predawn']  # Make positive for comparison
merged['midday_pos'] = -merged['midday']    # Make positive for comparison

# Key variables to compare
print("\n" + "="*70)
print("Variable Comparison (Mean ± Std)")
print("="*70)

comparisons = [
    ('leaf_temp', 'leaftemp_cpp', 'leaftemp_gt', 'Leaf Temperature (°C)'),
    ('assimilation', 'Anet-la', None, 'Assimilation (µmol/m²/s)'),
    ('transpiration', 'E-MD_cpp', 'E-MD_gt', 'Transpiration'),
    ('predawn_pos', 'P-PD', None, 'Predawn Pressure (MPa, sign-corrected)'),
    ('midday_pos', 'P-MD_cpp', None, 'Midday Pressure (MPa, sign-corrected)'),
]

for item in comparisons:
    if len(item) == 4:
        torch_col, cpp_col, gt_col, label = item
    else:
        torch_col, cpp_col, label = item
        gt_col = None
        
    if torch_col in merged.columns and cpp_col in merged.columns:
        torch_vals = merged[torch_col]
        cpp_vals = merged[cpp_col]
        
        print(f"\n{label}:")
        print(f"  PyTorch: {torch_vals.mean():.4f} ± {torch_vals.std():.4f}")
        print(f"  C++:     {cpp_vals.mean():.4f} ± {cpp_vals.std():.4f}")
        
        # Ground truth if available
        if gt_col and gt_col in merged.columns:
            gt_vals = merged[gt_col].dropna()
            if len(gt_vals) > 0:
                print(f"  Ground:  {gt_vals.mean():.4f} ± {gt_vals.std():.4f} (n={len(gt_vals)})")
        
        # Correlation PyTorch vs C++
        valid_mask = ~torch_vals.isna() & ~cpp_vals.isna() & (cpp_vals != 0)
        if valid_mask.sum() > 10:
            corr = np.corrcoef(torch_vals[valid_mask], cpp_vals[valid_mask])[0, 1]
            print(f"  Corr (Torch vs C++): {corr:.4f}")
        
        # RMSE
        rmse = np.sqrt(np.mean((torch_vals - cpp_vals)**2))
        print(f"  RMSE (Torch vs C++): {rmse:.4f}")

# Create comparison plots
fig, axes = plt.subplots(3, 2, figsize=(14, 12))
fig.suptitle("GARISOM: PyTorch vs C++ vs Ground Truth\nPopulation 1 (POFR CCR-COL)", 
             fontsize=14, fontweight='bold')

x = np.arange(len(merged))

# 1. Leaf Temperature with ground truth
ax1 = axes[0, 0]
ax1.plot(x, merged['leaftemp_cpp'], 'r-', label='C++', alpha=0.7, linewidth=1)
ax1.plot(x, merged['leaf_temp'], 'b--', label='PyTorch', alpha=0.7, linewidth=1)
if 'leaftemp_gt' in merged.columns:
    gt_mask = merged['leaftemp_gt'].notna()
    ax1.scatter(x[gt_mask], merged.loc[gt_mask, 'leaftemp_gt'], 
                c='green', s=5, alpha=0.5, label='Ground Truth', zorder=5)
ax1.set_ylabel('Leaf Temperature (°C)')
ax1.set_xlabel('Timestep')
ax1.set_title('Leaf Temperature')
ax1.legend()

# 2. Assimilation
ax2 = axes[0, 1]
ax2.plot(x, merged['Anet-la'], 'r-', label='C++', alpha=0.7, linewidth=1)
ax2.plot(x, merged['assimilation'], 'b--', label='PyTorch', alpha=0.7, linewidth=1)
ax2.set_ylabel('Assimilation (µmol/m²/s)')
ax2.set_xlabel('Timestep')
ax2.set_title('Net Assimilation')
ax2.legend()

# 3. Scatter plot - Leaf temp (PyTorch vs Ground Truth)
ax3 = axes[1, 0]
if 'leaftemp_gt' in merged.columns:
    valid = merged['leaftemp_gt'].notna()
    ax3.scatter(merged.loc[valid, 'leaftemp_gt'], merged.loc[valid, 'leaf_temp'], 
                alpha=0.3, s=10, c='blue', label='PyTorch')
    ax3.scatter(merged.loc[valid, 'leaftemp_gt'], merged.loc[valid, 'leaftemp_cpp'], 
                alpha=0.3, s=10, c='red', label='C++')
    ax3.plot([20, 50], [20, 50], 'k--', label='1:1 line')
    ax3.set_xlabel('Ground Truth Leaf Temp (°C)')
    ax3.set_ylabel('Model Leaf Temp (°C)')
    ax3.set_title('Leaf Temperature: Model vs Ground Truth')
    ax3.legend()
else:
    valid = merged['leaftemp_cpp'] != 0
    ax3.scatter(merged.loc[valid, 'leaftemp_cpp'], merged.loc[valid, 'leaf_temp'], alpha=0.3, s=10)
    ax3.plot([20, 50], [20, 50], 'k--', label='1:1 line')
    ax3.set_xlabel('C++ Leaf Temp (°C)')
    ax3.set_ylabel('PyTorch Leaf Temp (°C)')
    ax3.set_title('Leaf Temperature: 1:1 Comparison')
    ax3.legend()

# 4. Scatter plot - Assimilation
ax4 = axes[1, 1]
valid = merged['Anet-la'] != 0
ax4.scatter(merged.loc[valid, 'Anet-la'], merged.loc[valid, 'assimilation'], alpha=0.3, s=10)
ax4.plot([-5, 15], [-5, 15], 'k--', label='1:1 line')
ax4.set_xlabel('C++ Assimilation (µmol/m²/s)')
ax4.set_ylabel('PyTorch Assimilation (µmol/m²/s)')
ax4.set_title('Assimilation: 1:1 Comparison')
ax4.legend()

# 5. Water Potential (sign-corrected)
ax5 = axes[2, 0]
ax5.plot(x, merged['P-PD_cpp'], 'r-', label='C++ Predawn', alpha=0.7, linewidth=1)
ax5.plot(x, merged['predawn_pos'], 'b--', label='PyTorch Predawn', alpha=0.7, linewidth=1)
ax5.plot(x, merged['P-MD_cpp'], 'darkred', label='C++ Midday', alpha=0.5, linewidth=1)
ax5.plot(x, merged['midday_pos'], 'navy', linestyle='--', label='PyTorch Midday', alpha=0.5, linewidth=1)
ax5.set_ylabel('Water Potential (MPa)')
ax5.set_xlabel('Timestep')
ax5.set_title('Plant Water Potential (sign-corrected)')
ax5.legend(fontsize=8)

# 6. Residuals histogram - Model vs Ground Truth
ax6 = axes[2, 1]
if 'leaftemp_gt' in merged.columns:
    valid = merged['leaftemp_gt'].notna()
    residuals_torch = merged.loc[valid, 'leaf_temp'] - merged.loc[valid, 'leaftemp_gt']
    residuals_cpp = merged.loc[valid, 'leaftemp_cpp'] - merged.loc[valid, 'leaftemp_gt']
    ax6.hist(residuals_torch, bins=50, alpha=0.5, label=f'PyTorch (bias={residuals_torch.mean():.1f}°C)', color='blue')
    ax6.hist(residuals_cpp, bins=50, alpha=0.5, label=f'C++ (bias={residuals_cpp.mean():.1f}°C)', color='red')
    ax6.axvline(x=0, color='k', linestyle='--')
    ax6.set_xlabel('Leaf Temp Residual (Model - Ground Truth)')
    ax6.set_ylabel('Count')
    ax6.set_title('Leaf Temperature Residuals vs Ground Truth')
    ax6.legend()
else:
    valid = (merged['leaftemp_cpp'] != 0) & (merged['leaf_temp'] != 0)
    residuals = merged.loc[valid, 'leaf_temp'] - merged.loc[valid, 'leaftemp_cpp']
    ax6.hist(residuals, bins=50, edgecolor='black', alpha=0.7)
    ax6.axvline(x=0, color='r', linestyle='--', label='Zero')
    ax6.axvline(x=residuals.mean(), color='g', linestyle='-', label=f'Mean={residuals.mean():.2f}')
    ax6.set_xlabel('Leaf Temp Residual (PyTorch - C++)')
    ax6.set_ylabel('Count')
    ax6.set_title('Leaf Temperature Residuals')
    ax6.legend()

plt.tight_layout()
plt.savefig('model_comparison.png', dpi=150, bbox_inches='tight')
print(f"\nComparison plot saved to: model_comparison.png")

plt.show()

# Detailed stats for daytime only (when assimilation > 0)
print("\n" + "="*70)
print("Daytime-only Comparison (Assimilation > 0)")
print("="*70)

daytime = merged[merged['Anet-la'] > 0].copy()
print(f"Daytime timesteps: {len(daytime)}")

if len(daytime) > 0:
    print(f"\nLeaf Temperature:")
    print(f"  C++ mean:     {daytime['leaftemp_cpp'].mean():.2f} °C")
    print(f"  PyTorch mean: {daytime['leaf_temp'].mean():.2f} °C")
    print(f"  Bias (Torch-C++): {(daytime['leaf_temp'] - daytime['leaftemp_cpp']).mean():+.2f} °C")
    
    print(f"\nAssimilation:")
    print(f"  C++ mean:     {daytime['Anet-la'].mean():.2f} µmol/m²/s")
    print(f"  PyTorch mean: {daytime['assimilation'].mean():.2f} µmol/m²/s")
    print(f"  Bias (Torch-C++): {(daytime['assimilation'] - daytime['Anet-la']).mean():+.2f} µmol/m²/s")
    
    print(f"\nTranspiration:")
    print(f"  C++ mean:     {daytime['E-MD_cpp'].mean():.2f}")
    print(f"  PyTorch mean: {daytime['transpiration'].mean():.2f}")
    print(f"  Bias (Torch-C++): {(daytime['transpiration'] - daytime['E-MD_cpp']).mean():+.2f}")
    
    print(f"\nMidday Water Potential (sign-corrected):")
    print(f"  C++ mean:     {daytime['P-MD_cpp'].mean():.3f} MPa")
    print(f"  PyTorch mean: {daytime['midday_pos'].mean():.3f} MPa")
    print(f"  Bias (Torch-C++): {(daytime['midday_pos'] - daytime['P-MD_cpp']).mean():+.3f} MPa")

# Ground truth comparison
print("\n" + "="*70)
print("Ground Truth Comparison (Leaf Temperature)")
print("="*70)

if 'leaftemp_gt' in merged.columns:
    gt_valid = merged['leaftemp_gt'].notna()
    print(f"Timesteps with ground truth: {gt_valid.sum()}")
    
    if gt_valid.sum() > 0:
        gt_data = merged[gt_valid].copy()
        
        # PyTorch vs Ground Truth
        torch_rmse = np.sqrt(np.mean((gt_data['leaf_temp'] - gt_data['leaftemp_gt'])**2))
        torch_bias = (gt_data['leaf_temp'] - gt_data['leaftemp_gt']).mean()
        torch_corr = np.corrcoef(gt_data['leaf_temp'], gt_data['leaftemp_gt'])[0, 1]
        
        # C++ vs Ground Truth
        cpp_rmse = np.sqrt(np.mean((gt_data['leaftemp_cpp'] - gt_data['leaftemp_gt'])**2))
        cpp_bias = (gt_data['leaftemp_cpp'] - gt_data['leaftemp_gt']).mean()
        cpp_corr = np.corrcoef(gt_data['leaftemp_cpp'], gt_data['leaftemp_gt'])[0, 1]
        
        print(f"\nGround Truth mean: {gt_data['leaftemp_gt'].mean():.2f} ± {gt_data['leaftemp_gt'].std():.2f} °C")
        print(f"\nPyTorch vs Ground Truth:")
        print(f"  Mean: {gt_data['leaf_temp'].mean():.2f} °C")
        print(f"  Bias: {torch_bias:+.2f} °C")
        print(f"  RMSE: {torch_rmse:.2f} °C")
        print(f"  Corr: {torch_corr:.4f}")
        
        print(f"\nC++ vs Ground Truth:")
        print(f"  Mean: {gt_data['leaftemp_cpp'].mean():.2f} °C")
        print(f"  Bias: {cpp_bias:+.2f} °C")
        print(f"  RMSE: {cpp_rmse:.2f} °C")
        print(f"  Corr: {cpp_corr:.4f}")

# Summary of key issues
print("\n" + "="*70)
print("KEY DIFFERENCES / ISSUES")
print("="*70)
print("""
1. WATER BALANCE: PyTorch doesn't update soil water between timesteps.
   - Predawn stays constant near field capacity (~0 MPa)
   - C++ has dynamic water balance with rainfall/drainage

2. TRANSPIRATION: PyTorch transpiration weakly correlates with C++
   - May be due to different optimization (soft-argmax vs discrete)
   - Different water potentials affect supply curve

3. LEAF TEMPERATURE: PyTorch runs ~7°C hotter
   - Energy balance may differ
   - Transpiration cooling is reduced (lower E)

4. ASSIMILATION: PyTorch is lower than C++
   - Follows from lower transpiration (stomatal limitation)
   - Higher leaf temp may reduce photosynthesis
""")

print("="*70)
print("Done!")
print("="*70)
