import matplotlib.pyplot as plt
from math import pi
import numpy as np
import pandas as pd

def main():
    file = "data/history_single.csv"
    amplitude_displacement = 0.01
    frequency = 15

    df = pd.read_csv(file)
    df.columns = ["x", "v", "t"]  # TODO: Colorize according to initial condition

    df["x_belt"] = amplitude_displacement * np.sin(2*pi*frequency*df.t)
    df["v_belt"] = 2*pi*frequency * amplitude_displacement*np.cos(2*pi*frequency*df.t)

    fig, ax = plt.subplots()
    ax.plot(df.x, df.v, label="Block")
    ax.plot(df.x_belt, df.v_belt, label="Belt")
    plt.show()


if __name__ == "__main__":
    main()