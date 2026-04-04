#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import sys

project_root = '/home/ksenya/sobkina_ke_sle'
build_dir = os.path.join(project_root, 'build')
convergence_csv = os.path.join(build_dir, 'convergence2.csv')
timing_csv = os.path.join(build_dir, 'timing2.csv')
output_convergence = os.path.join(project_root, 'convergence_v2.png')
output_time = os.path.join(project_root, 'convergence_time_v2.png')

plt.rcParams['figure.figsize'] = (10, 6)
plt.rcParams['font.size'] = 11

methods = {
    'jacobi_error': {'color': '#1f77b4', 'label': 'Jacobi'},
    'gaussseidel_error': {'color': '#ff7f0e', 'label': 'Gauss-Seidel'},
    'symmetricgs_error': {'color': '#2ca02c', 'label': 'Symmetric GS'},
    'chebyshevsym_error': {'color': '#d62728', 'label': 'Chebyshev Sym'}
}

if not os.path.exists(convergence_csv):
    print(f"Error: file not found: {convergence_csv}")
    sys.exit(1)

df = pd.read_csv(convergence_csv)
iter_col = df.columns[0]

plt.figure()
for col_name, props in methods.items():
    if col_name in df.columns:
        data = df[[iter_col, col_name]].dropna()
        if len(data) > 1:
            plt.semilogy(
                data[iter_col],
                data[col_name],
                label=props['label'],
                color=props['color'],
                linewidth=2,
                marker='.',
                markersize=3
            )

plt.xlabel('N iter', fontsize=12)
plt.ylabel('Error', fontsize=12)
plt.title('Convergence Comparison (Symmetric Matrix)', fontsize=14, fontweight='bold')
plt.legend(loc='upper right', fontsize=10)
plt.grid(True, which='both', alpha=0.6, linestyle='--')
plt.tight_layout()
plt.savefig(output_convergence, dpi=300, bbox_inches='tight')
plt.close()
print(f"Saved: {output_convergence}")

if os.path.exists(timing_csv):
    df_time = pd.read_csv(timing_csv)

    plt.figure()
    for col_name, props in methods.items():
        if col_name not in df.columns:
            continue

        method_key = col_name.replace('_error', '')
        time_data = df_time[df_time['method'] == method_key]

        if time_data.empty:
            continue

        total_time = time_data['time_ms'].values[0]
        conv_data = df[[iter_col, col_name]].dropna()

        if len(conv_data) < 2 or total_time <= 0:
            continue

        times = np.linspace(0, total_time, len(conv_data))

        plt.semilogy(
            times,
            conv_data[col_name].values,
            label=props['label'],
            color=props['color'],
            linewidth=2
        )

    plt.xlabel('Time (ms)', fontsize=12)
    plt.ylabel('Error', fontsize=12)
    plt.title('Convergence vs Time', fontsize=14, fontweight='bold')
    plt.legend(loc='upper right', fontsize=10)
    plt.grid(True, which='both', alpha=0.6, linestyle='--')
    plt.tight_layout()
    plt.savefig(output_time, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_time}")

print("Done!")