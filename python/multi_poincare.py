import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

file_path = "data/poincare_multi_13.000000_100_100000.txt"
data = np.loadtxt(file_path)
df = pd.read_csv(file_path)

positions = data[:,0]
velocities = data[:,1]

plt.scatter(positions, velocities, s=0.25)
plt.show()