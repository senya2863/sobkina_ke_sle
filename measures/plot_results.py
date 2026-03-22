#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import sys

# Проверяем, что запускаем из корня проекта
if not os.path.exists('measures/convergence.csv'):
    print("❌ Ошибка: запустите скрипт из корня проекта (/home/ksenya/sobkina_ke_sle/)")
    print("   Или проверьте, что build/convergence.csv существует")
    sys.exit(1)

# === Чтение данных ===
df = pd.read_csv('/home/ksenya/sobkina_ke_sle/measures/convergence.csv')
timing_df = pd.read_csv('/home/ksenya/sobkina_ke_sle/measures/timing.csv')

# === График 1: Ошибка от числа итераций ===
plt.figure(figsize=(10, 6))
plt.semilogy(df['Iter'], df['MPI_error'], '--', label='MPI', linewidth=2)
plt.semilogy(df['Iter'], df['Jacobi_error'], '--', label='Jacobi', linewidth=2)
plt.semilogy(df['Iter'], df['GaussSeidel_error'], '--', label='Gauss-Seidel', linewidth=2)

if 'Chebyshev_error' in df.columns:
    plt.semilogy(df['Iter'], df['Chebyshev_error'], '--', label='Chebyshev', linewidth=2)

plt.xlabel('N iter', fontsize=12)
plt.ylabel('Error', fontsize=12)
plt.title('Convergence Comparison', fontsize=14)
plt.legend(loc='upper right')
plt.grid(True, which='both', alpha=0.7)
plt.tight_layout()
plt.savefig('convergence_iter.png')
plt.close()
print("✓ Saved: convergence_iter.png")

# === График 2: Ошибка от времени ===
plt.figure(figsize=(10, 6))

# Сопоставление: имя метода в timing.csv → имя колонки в convergence.csv
methods_map = {
    'MPI': 'MPI_error',
    'Jacobi': 'Jacobi_error',
    'GaussSeidel': 'GaussSeidel_error',
    'Chebyshev': 'Chebyshev_error'
}

colors = ['blue', 'orange', 'green', 'red']

for idx, (method_key, col_name) in enumerate(methods_map.items()):
    # Проверяем, есть ли такая колонка с ошибками в нашем CSV
    if col_name not in df.columns:
        continue

    # Ищем время для этого метода в timing.csv
    time_row = timing_df[timing_df['Method'] == method_key]
    if time_row.empty:
        continue

    total_time = time_row['Time_ms'].values[0]

    # Берём только валидные данные (без NaN)
    valid_data = df[[ 'Iter', col_name ]].dropna()
    if len(valid_data) < 2:
        continue

    # Создаём временную шкалу
    times = np.linspace(0, total_time, len(valid_data))
    errors = valid_data[col_name].values

    plt.semilogy(times, errors, label=method_key, color=colors[idx], linewidth=2)

plt.xlabel('Time (ms)', fontsize=12)
plt.ylabel('Error', fontsize=12)
plt.title('Convergence vs Time', fontsize=14)
plt.legend(loc='upper right')
plt.grid(True, which='both', alpha=0.7)
plt.tight_layout()
plt.savefig('convergence_time.png')
plt.close()
print("✓ Saved: convergence_time.png")
