import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import matplotlib.ticker as mtick
import pandas as pd
import numpy as np
import scipy.fftpack
from scipy.signal import find_peaks
from math import pi, log10


def main():
    data_file = "data/shear_13.000000Hz_100blocks.csv"
    df = pd.read_csv(data_file)
    df.columns = ["time", "shear"]
    dt = df.time.iloc[1]-df.time.iloc[0]

    num_steps = len(df.time)
    yf = scipy.fftpack.fft(df.shear.values)
    yf = 2.0/num_steps * np.abs(yf)
    xf =  np.linspace(0.0, 1.0//(2.0*dt), num_steps//2)

    fig = plt.figure(figsize=(8,3))
    ax1 = plt.subplot(121)
    ax2 = plt.subplot(122)

    ax1.plot(df.time, df.shear)
    ax1.set_xlabel("Time")
    ax1.set_ylabel("Shear")

    ax2.plot(xf, yf[:num_steps//2])
    xmax, ymax = calc_fft_limits(xf,yf)

    plt.show()

def calc_fft_limits(x, y): # TODO: Fix. rediculux algorithm
    threshold = -1
    num_steps = len(x)
    max_freq = 0
    for i in range(1, len(x[:num_steps//2])): #len(power_spectrum)//2):
        if y[i] > max_freq:
            max_freq = y[i]
            frequency_shift = i
        max_freq = max(max_freq, y[i])
        if y[i] > 0.1 * max_freq:
            threshold = i
        
    threshold_index = int(1.2 * threshold)
    x_max = x[threshold_index]
    y_max = 1.2 * max_freq
    return x_max, y_max


if __name__ == "__main__":
    main()