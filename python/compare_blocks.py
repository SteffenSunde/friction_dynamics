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

    file = "data/mdof/mdof_history_f15.000000_P25.000000_xi0.050000_del1.000000.csv"

    df = pd.read_csv(files, skiprows=1)
    meta = parse_header(file)

    save_figure = False
    figure_path = "data/compare_blocks.pdf"

    
    amplitude_displacement = float(meta['disp'])
    frequency = float(meta['f'])

    df = pd.read_csv(file, skiprows=1)
    df.columns = ["x", "v", "t"]
    df["x_belt"] = amplitude_displacement * np.sin(2*pi*frequency*df.t)
    df["v_belt"] = 2*pi*frequency * amplitude_displacement*np.cos(2*pi*frequency*df_stick.t)
    df["slip"] = - df.x + df.x_belt


    fig = plt.figure(figsize=(8,6))
    
    ax1 = fig.add_subplot(221)
    ax1.plot(df.x, df.v, label="Block")
    ax1.plot(df.x_belt, df.v_belt, label="Belt")
    ax1.legend()
    ax1.set_xlabel("Displacement [mm]")
    ax1.set_ylabel("Velocity [mm/s]")
    ax1.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
    ax1.legend(loc="upper left")
    #ax1.grid()

    ax2 = fig.add_subplot(222)
    ax2.plot(df.t, df.slip, label="Slip")
    ax2.plot(df.t, df.x_belt, label="Belt")
    ax2.set_xlim((df.t.iloc[0], df.t.iloc[0]+2/frequency))
    ax2.set_xlabel("Time [s]")
    ax2.set_ylabel("Slip [mm]")
    ax2.legend(loc="upper left")
    ax2.set_xticks([ax2.get_xticks()[1], ax2.get_xticks()[-1]])
    ax2.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    ax2.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    #ax2.grid()

    ax3 = fig.add_subplot(223)
    ax3.plot(df.x, df.v, label="Block", linestyle="-")
    ax3.plot(df.x_belt, df.v_belt, label="Belt")
    ax3.legend()
    ax3.set_xlabel("Displacement [mm]")
    ax3.set_ylabel("Velocity [mm/s]")
    ax3.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
    ax3.legend(loc="upper left")
    #ax3.grid()

    ax4 = fig.add_subplot(224)
    ax4.plot(df.t, df.slip, label="Slip")
    ax4.plot(df.t, df.x_belt, label="Belt")
    ax4.set_xlim((df.t.iloc[0], df.t.iloc[0]+2/frequency))
    ax4.set_xlabel("Time [s]")
    ax4.set_ylabel("Slip [mm]")
    ax4.legend(loc="upper left")
    ax4.set_xticks([ax4.get_xticks()[1], ax4.get_xticks()[-1]])
    ax4.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    ax4.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    #ax4.grid()
    
    #fig.suptitle(figure_title)
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    if save_figure: fig.savefig(figure_path)
    plt.show()


if __name__ == "__main__":
    main()