import pandas as pd
import matplotlib.pyplot as plt

data = pd.read_csv("times.txt", sep=" ")

plt.figure()

plt.plot(data["size"], data["dense"], marker="o", label="Dense")
plt.plot(data["size"], data["csr"], marker="o", label="CSR")

plt.xlabel("Matrix size")
plt.ylabel("Time (seconds)")
plt.title("Dense vs CSR matrix-vector multiplication")

plt.legend()
plt.grid()

plt.savefig("graph.png", dpi=300)
plt.savefig("graph.pdf")

plt.show()