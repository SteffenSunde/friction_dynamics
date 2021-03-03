import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

file_path = "data/poincare_multi_13.000000_100_2000000.txt"
num_blocks = 100
num_initials = 10

#data = np.loadtxt(file_path)
df = pd.read_csv(file_path)
df.columns = ["pos", "vel"]

num_points = df.shape[0]//num_initials

# positions = data[:,0]
# velocities = data[:,1]

for i in range(num_initials):
    start = i*num_points
    stop = start + num_points
    plt.scatter(df.pos.iloc[start:stop], df.vel.iloc[start:stop], s=0.25)
plt.show()