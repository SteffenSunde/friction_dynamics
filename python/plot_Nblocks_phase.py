import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import scipy.fftpack
from scipy.signal import find_peaks
from math import pi

def main(dof=50):
    file = "data/mdof/mdof_slip_history_f15.000000_P25.000000_xi0.050000_del1.000000.csv"
    frequency = 15
    displacement = 0.01
    df = pd.read_csv(file)
    plot_result(df.iloc[:,3*dof].values, df.iloc[:,3*dof+1].values, df.iloc[:,-1].values, frequency, displacement)


def plot_result(positions, velocities, time, frequency, amplitude, frequency_peaks=True, meta=""):
    dt = time[1]-time[0]  # Note: assuming equispaced timesteps.
    num_steps = len(time)

    belt_positions = amplitude*np.sin(2*pi*frequency*time)  #data[-2, transient:] #
    belt_velocities = 2*pi*frequency*amplitude*np.cos(2*pi*frequency*time)  #data[-1, transient:]
    slip = positions - belt_positions

    ## Calculate fft of given dof
    yf = scipy.fftpack.fft(positions)
    yf = 2.0/num_steps * np.abs(yf)
    xf =  np.linspace(0.0, 1.0//(2.0*dt), num_steps//2)

    # Find max for plot limits etc.
    max_freq = 0
    threshold = -1
    for i in range(1, len(xf[:num_steps//2])): #len(power_spectrum)//2):
        if yf[i] > max_freq:
            max_freq = yf[i]
            frequency_shift = i
        max_freq = max(max_freq, yf[i])
        if yf[i] > 0.1 * max_freq:
            threshold = i
        
    threshold_index = int(1.2 * threshold)
    x_max = xf[threshold_index]
    y_max = 1.2 * max_freq

    fig = plt.figure(figsize=(10,10))
    #fig.suptitle(f"Block {block} ({meta})", fontsize=16)

    ax1 = plt.subplot(221)
    ax2 = plt.subplot(222)
    ax3 = plt.subplot(223)
    ax4 = plt.subplot(224)

    ax1.plot(time, positions)
    ax1.set_xlabel("Time [s]")
    ax1.set_ylabel("Displacement [mm]")
    ax1.ticklabel_format(axis='y', style='sci', scilimits=(0,0))

    ax2.plot(positions, velocities, label="Block")
    ax2.plot(belt_positions, belt_velocities, label="Belt")
    ax2.set_xlabel("Displacement [mm]")
    ax2.set_ylabel("Velocity [mm/s]")
    ax2.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
    ax2.legend(loc="upper left")

    #ax3.plot(freq[1:], power_spectrum[1:])
    ax3.plot(xf, yf[:num_steps//2])
    ax3.set_xlim((0, x_max))
    ax3.set_ylim((0, y_max))
    ax3.set_xlabel("Frequency [Hz]")
    ax3.set_ylabel("Displacement [mm]")

    if frequency_peaks:
        peak_indices, _  = find_peaks(yf[:num_steps//2])
        peak_frequencies = xf[peak_indices]
        peak_heights = yf[peak_indices]
        ax3.scatter(peak_frequencies, peak_heights, s=60, edgecolors='r', facecolors='none')
        for i in range(len(peak_frequencies)):
            ax3.annotate(" {:.1f}".format(peak_frequencies[i]), (peak_frequencies[i], peak_heights[i]), horizontalalignment="left", rotation=45, size=12)

    ax3.ticklabel_format(axis='y', style='sci', scilimits=(0,0))

    #ax4.plot(slip, velocities)
    #ax4.set_xlabel("Slip [mm]")
    ax4.plot(time, slip)
    ax4.set_xlabel("Time [s]")
    ax4.set_ylabel("Slip [mm]")
    ax4.ticklabel_format(axis='y', style='sci', scilimits=(0,0))

    #meta = {"N": f"{num_blocks}", "transient": "}
    
    #fig.savefig(f"hertzstatesine_1dof_{frequency}Hz_{amplitude}_{meta}.png", dpi=400)
    plt.show() 


if __name__ == "__main__":
    main()