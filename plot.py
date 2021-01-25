import matplotlib.pyplot as plt
import pandas as pd

df = pd.read_csv("data/rate_perturbed_poincare/poincare_20.000000_100_1000000.txt", usecols=(0,1))
df.columns = ["x", "v"]

plt.scatter(df.x, df.v, c="black", s=0.25)
plt.grid()
plt.show()