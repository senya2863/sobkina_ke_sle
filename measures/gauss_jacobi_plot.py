import pandas as pd
import matplotlib.pyplot as plt
import os

csv_path = "convergence.csv"

if not os.path.exists(csv_path):
    raise FileNotFoundError(f"CSV не найден: {csv_path}")

data = pd.read_csv(csv_path)

plt.figure(figsize=(8,6))

plt.plot(data["Iter"], data["Jacobi"], label="Jacobi", marker='o')
plt.plot(data["Iter"], data["GaussSeidel"], label="Gauss-Seidel", marker='x')
plt.plot(data["Iter"], data["SimpleIter"], label="Simple iteration", marker='s')

plt.yscale("log")

plt.xlabel("Iteration")
plt.ylabel("Norm of error")
plt.title("Convergence of iterative methods")

plt.grid(True, which="both", linestyle="--", linewidth=0.5)
plt.legend()

plt.savefig("convergence.png", dpi=300)
plt.show()