import matplotlib.pyplot as plt
import matplotlib
import scipy.fftpack
from math import pi
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d, Axes3D
import glob
from matplotlib.ticker import MaxNLocator
from matplotlib.ticker import FormatStrFormatter

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

    #file = "data/mdof/mdof_slip_history_f15.000000_P25.000000_xi0.050000_del1.000000.csv"
    file = "data/mdof/mass_damped/mdof_slip_history_f15.000000_P100.000000_xi0.050000_del1.000000.csv"

    df = pd.read_csv(file, skiprows=1)
    df.columns = ["time", "slip_end", "slip_edge", "slip_center", "trash"]
    meta = parse_header(file)

    save_figure = True
    figure_path = "data/compare_blocks_stick.pdf"
    
    amplitude_displacement = float(meta['disp'])
    frequency = float(meta['f'])

    df["x_belt"] = -amplitude_displacement * np.sin(2*pi*frequency*df.time)
    df["v_belt"] = 2*pi*frequency * amplitude_displacement*np.cos(2*pi*frequency*df.time)

    #fig = plt.figure(figsize=(10,4))
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(9,3), sharey=True)
    
    #ax1 = fig.add_subplot(131, sharey="all")
    ax1.plot(df.time, df.slip_end, label="Block")
    ax1.plot(df.time, df.x_belt, label="Belt")
    ax1.legend()
    #ax1.set_xlabel("Time [s]")
    ax1.set_ylabel("Slip [mm]")
    ax1.set_xlim((df.time.iloc[0], df.time.iloc[0]+1.1/frequency))
    ax1.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
    #ax1.set_xticks([ax1.get_xticks()[1], ax1.get_xticks()[-1]])
    #ax1.set_xticks([33.3, 33.4])
    ax1.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax1.legend(loc="upper left")
    ax1.set_title("Chain end")
    ax1.grid()

    #ax2 = fig.add_subplot(132)
    ax2.plot(df.time, df.slip_edge, label="Blo")
    ax2.plot(df.time, df.x_belt, label="Belt")
    ax2.set_xlim((df.time.iloc[0], df.time.iloc[0]+1.1/frequency))
    #ax2.set_xlabel("Time [s]")
    #ax2.set_ylabel("Slip [mm]")
    #ax2.legend(loc="upper left")
    #ax2.set_xticks([33.3, 33.4])
    ax2.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax2.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    ax2.set_title("Contact edge")
    ax2.grid()

    #ax3 = fig.add_subplot(133)
    ax3.plot(df.time, df.slip_center, label="Center", linestyle="-")
    ax3.plot(df.time, df.x_belt, label="Belt")
    #ax3.set_xlabel("Time [s]")
    #ax3.set_ylabel("Slip [mm]")
    #ax3.set_xticks([ax1.get_xticks()[1], ax1.get_xticks()[-1]])
    ax3.set_xlim((df.time.iloc[0], df.time.iloc[0]+1.1/frequency))
    ax3.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax3.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    ax3.set_title("Contact center")
    ax3.grid()
    #ax3.legend(loc="upper left")
    #ax3.grid()

    fig.text(0.53, 0.01, 'Time [s]', ha='center')

    plt.tight_layout()
    if save_figure: fig.savefig(figure_path)
    plt.show()


if __name__ == "__main__":
    main()
