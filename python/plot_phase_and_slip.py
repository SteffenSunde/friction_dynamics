import matplotlib.pyplot as plt
import matplotlib
import scipy.fftpack
from math import pi
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d, Axes3D

def parse_header(file_path):
    arguments = {}
    with open(file_path) as f:
        line = f.readline()
        while (line.startswith('#')):
            wordlist = line[1:].split(',')
            for word in wordlist:
                key, value = word.split(':')
                arguments[key] = value
            line = f.readline()
    return arguments

def main():
    #matplotlib.rcParams.update({'font.size': 14})

    file = "data/history_single.csv"

    meta = parse_header(file)
    figure_title = "".join(["{}: {} ".format(key, value) for key, value in meta.items()])

    save_figure = False
    figure_path = "data/history_single_f10_chaos.pdf"
    amplitude_displacement = float(meta['disp'])
    frequency = float(meta['f'])

    df = pd.read_csv(file, skiprows=1)
    df.columns = ["x", "v", "t"]
    num_steps = df.shape[0]
    dt = df.t[1] - df.t[0]

    df["x_belt"] = amplitude_displacement * np.sin(2*pi*frequency*df.t)
    df["v_belt"] = 2*pi*frequency * amplitude_displacement*np.cos(2*pi*frequency*df.t)
    df["slip"] = -df.x + df.x_belt

    fig = plt.figure(figsize=(12,4))
    
    # yf1 = scipy.fftpack.fft(df.x.values)
    # yf1 = 2.0/num_steps * np.abs(yf1)
    # xf1 =  np.linspace(0.0, 1.0//(2.0*dt), num_steps//2)
    # ax1 = fig.add_subplot(131)
    # ax1.plot(xf1, yf1[:num_steps//2])
    # ax1.set_xlim((0,100))
    # ax1.set_xlabel(r"Frequency [Hz]", fontsize=14)
    # ax1.set_ylabel(r"$x$", fontsize=14)
    # ax1.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

    ax2 = fig.add_subplot(121)
    ax2.plot(df.x, df.v, label="Block", linestyle="-", linewidth=0.75)
    ax2.plot(df.x_belt, df.v_belt, label="Belt", linewidth=0.75)
    ax2.legend()
    ax2.set_xlabel(r"$Displacement [mm]$")
    ax2.set_ylabel(r"$Velocity [mm/s]$")
    ax2.ticklabel_format(style='sci', axis='both', scilimits=(0,0))

    ax3 = fig.add_subplot(122)
    ax3.plot(df.t, df.slip)
    ax3.plot(df.t, df.x_belt)
    ax3.set_xlim((df.t.iloc[0], df.t.iloc[0]+2/frequency))
    ax3.set_xlabel(r"$Time [s]$")
    ax3.set_ylabel(r"$Slip [mm]$")
    #ax3.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
    #ax3.subplots_adjust(left=-0.5, right=0.5, top=0, bottom=0)
    #ax1.set_yticks([-10, 10])
    
    fig.suptitle(figure_title)
    plt.tight_layout()
    if save_figure: fig.savefig(figure_path)
    plt.show()


if __name__ == "__main__":
    main()