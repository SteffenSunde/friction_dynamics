import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import pandas as pd
import numpy as np
import scipy.fftpack
from scipy.signal import find_peaks
from math import pi

def main(dof=0):
    file = "data/shear_13.000000Hz_100blocks.csv"
    frequency = 11
    displacement = 0.01
    df = pd.read_csv(file)
    plot_result(df.iloc[:,0].values, df.iloc[:,1].values, df.iloc[:,2].values, df.iloc[:,3].values, df.iloc[:,4].values, frequency, displacement)


def plot_result(time, values1, values2, values3, values4, frequency_peaks=True, meta=""):
    dt = time[1]-time[0]  # Note: assuming equispaced timesteps.
    num_steps = len(time)

    ## Calculate fft of given dof
    yf1 = scipy.fftpack.fft(values1)
    yf1 = 2.0/num_steps * np.abs(yf1)
    xf1 =  np.linspace(0.0, 1.0//(2.0*dt), num_steps//2)

    yf2 = scipy.fftpack.fft(values2)
    yf2 = 2.0/num_steps * np.abs(yf2)
    xf2 =  np.linspace(0.0, 1.0//(2.0*dt), num_steps//2)

    yf3 = scipy.fftpack.fft(values3)
    yf3 = 2.0/num_steps * np.abs(yf3)
    xf3 =  np.linspace(0.0, 1.0//(2.0*dt), num_steps//2)

    yf4 = scipy.fftpack.fft(values4)
    yf4 = 2.0/num_steps * np.abs(yf4)
    xf4 =  np.linspace(0.0, 1.0//(2.0*dt), num_steps//2)

    # Find max for plot limits etc.


    fig = plt.figure(figsize=(10,10))
    fmt = FuncFormatter(lambda x, pos: tickformat(x / 1e-2))
    #ax.xaxis.set_major_formatter(fmt)
    #fig.suptitle(f"Block {block} ({meta})", fontsize=16)

    ax1 = plt.subplot(221)
    ax2 = plt.subplot(222)
    ax3 = plt.subplot(223)
    ax4 = plt.subplot(224)

    #ax3.plot(freq[1:], power_spectrum[1:])
    ax1.plot(xf1, yf1[:num_steps//2])
    xmax, ymax = calc_fft_limits(xf1, yf1)
    ax1.set_xlim((0, xmax))
    ax1.set_ylim((0, ymax))
    #ax1.ticklabel_format(axis='y',style='sci')
    #ax1.tick_params(axis='x', which='major')
    #ax1.yaxis.set_major_formatter(fmt)
    ax1.set_xlabel("Frequency [Hz]")
    #ax1.set_ylabel(r"Slip [($\times 10^{-2}$ mm)]")
    ax1.set_ylabel("Slip [mm]")
    ax1.set_title("Trailing end")

    ax2.plot(xf2, yf2[:num_steps//2])
    xmax, ymax = calc_fft_limits(xf2, yf2)
    ax2.set_xlim((0, xmax))
    ax2.set_ylim((0, ymax))
    #ax2.tick_params(axis='x', which='major')
    #ax2.yaxis.set_major_formatter(fmt)
    ax2.set_xlabel("Frequency [Hz]")
    #ax2.set_ylabel(r"Slip [($\times 10^{-2}$ mm)]")
    ax2.set_ylabel("Slip [mm]")
    ax2.set_title("Trailing contact edge")

    ax3.plot(xf3, yf3[:num_steps//2])
    xmax, ymax = calc_fft_limits(xf3, yf3)
    ax3.set_xlim((0, xmax))
    ax3.set_ylim((0, ymax))
    ax3.set_xlabel("Frequency [Hz]")
    ax3.set_ylabel("Displacement [?]")
    ax3.set_title("Contact center")

    ax4.plot(xf4, yf4[:num_steps//2])
    xmax, ymax = calc_fft_limits(xf4, yf4)
    ax4.set_xlim((0, xmax))
    ax4.set_ylim((0, ymax))
    ax4.set_xlabel("Frequency [Hz]")
    ax4.set_ylabel("Displacement [?]")
    ax4.set_title("Total displacement")

    # if frequency_peaks:
    #     peak_indices, _  = find_peaks(yf[:num_steps//2])
    #     peak_frequencies = xf[peak_indices]
    #     peak_heights = yf[peak_indices]
    #     ax2.scatter(peak_frequencies, peak_heights, s=60, edgecolors='r', facecolors='none')
    #     for i in range(len(peak_frequencies)):
    #         ax2.annotate(" {:.1f}".format(peak_frequencies[i]), (peak_frequencies[i], peak_heights[i]), horizontalalignment="left", rotation=45, size=12)

    #ax2.ticklabel_format(axis='y', style='sci', scilimits=(0,0))

    #fig.savefig(f"hertzstatesine_1dof_{frequency}Hz_{amplitude}_{meta}.png", dpi=400)
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

def tickformat(x):
    if int(x) == float(x):
        return str(int(x))
    else:
        return str(x) 

if __name__ == "__main__":
    main()