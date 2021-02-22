import matplotlib.pyplot as plt
import matplotlib
import scipy.fftpack
from math import pi
import numpy as np
import pandas as pd


# def parse_header(file_path):
#     arguments = {}
#     with open(file_path) as f:
#         line = f.readline()
#         while (line.startswith('#')):
#             wordlist = line.split(',')
#             for word in wordlist:
#                 key, value = word.split(':')
#                 arguments[key] = value
#             line = f.readline()


def main():
    #matplotlib.rcParams.update({'font.size': 14})

    file = "data/history_single.csv"

    #print(parse_header(file))

    save_figure = True
    figure_path = "data/history_single_f10_chaos.pdf"
    amplitude_displacement = 0.01
    frequency = 15

    df = pd.read_csv(file, skiprows=1)
    df.columns = ["x", "v", "t"]
    num_steps = df.shape[0]
    dt = df.t[1] - df.t[0]

    df["x_belt"] = amplitude_displacement * np.sin(2*pi*frequency*df.t)
    df["v_belt"] = 2*pi*frequency * amplitude_displacement*np.cos(2*pi*frequency*df.t)
    df["slip"] = df.x - df.x_belt

    fig = plt.figure(figsize=(12,4))
    
    yf1 = scipy.fftpack.fft(df.x.values)
    yf1 = 2.0/num_steps * np.abs(yf1)
    xf1 =  np.linspace(0.0, 1.0//(2.0*dt), num_steps//2)
    ax1 = fig.add_subplot(131)
    ax1.plot(xf1, yf1[:num_steps//2])
    ax1.set_xlim((0,50))
    ax1.set_xlabel(r"Frequency [Hz]", fontsize=14)
    ax1.set_ylabel(r"$x$", fontsize=14)
    ax1.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

    ax2 = fig.add_subplot(132)
    ax2.plot(df.x, df.v, label="Block", linestyle="-", linewidth=0.75)
    ax2.plot(df.x_belt, df.v_belt, label="Belt", linewidth=0.75)
    ax2.legend()
    ax2.set_xlabel(r"$x$", fontsize=14)
    ax2.set_ylabel(r"$\dot{x}$", fontsize=14)
    ax2.ticklabel_format(style='sci', axis='both', scilimits=(0,0))

    ax3 = fig.add_subplot(133, projection='3d')
    ax3.plot(df.x, df.v, df.x_belt, linewidth=0.75)
    ax3.set_xlabel(r"$x$", fontsize=14)
    ax3.set_ylabel(r"$\dot{x}$", fontsize=14)
    ax3.set_zlabel(r"$x_0$", fontsize=14)
    ax3.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
    ax3.view_init(elev=20, azim=55)
    #ax3.subplots_adjust(left=-0.5, right=0.5, top=0, bottom=0)
    #ax1.set_yticks([-10, 10])
    
    
    plt.tight_layout()
    if save_figure: fig.savefig(figure_path)
    plt.show()


if __name__ == "__main__":
    main()