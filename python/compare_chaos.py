import matplotlib.pyplot as plt
import matplotlib
import scipy.fftpack
from math import pi
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
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

    figure_title = "Regular vs Chaotic solution"

    file_regular = "data/sdof/newdamp/compare_chaos/history_single_p1200.000000_xi0.000000.csv"
    meta_regular = parse_header(file_regular)

    file_chaos = "data/sdof/newdamp/compare_chaos/history_single_p1200.000000_xi0.050000.csv"
    meta_chaos = parse_header(file_chaos)
    #figure_title = "".join(["{}: {} ".format(key, value) for key, value in meta.items()])

    print(meta_chaos)

    save_figure = False
    figure_path = "data/compare_chaos.pdf"
    for key, val in meta_regular.items():
        if val != meta_chaos[key]:
            print("Difference (regular/chaos) in {}: {} vs {}".format(key, val, meta_chaos[key]))

    amplitude_displacement = float(meta_regular['disp'])
    frequency = float(meta_regular['f'])

    df_regular = pd.read_csv(file_regular, skiprows=1)
    df_regular.columns = ["x", "v", "t"]
    df_regular["x_belt"] = amplitude_displacement * np.sin(2*pi*frequency*df_regular.t)
    df_regular["v_belt"] = 2*pi*frequency * amplitude_displacement*np.cos(2*pi*frequency*df_regular.t)
    df_regular["slip"] = - df_regular.x + df_regular.x_belt

    df_chaos = pd.read_csv(file_chaos, skiprows=1)
    df_chaos.columns = ["x", "v", "t"]
    df_chaos["x_belt"] = amplitude_displacement * np.sin(2*pi*frequency*df_chaos.t)
    df_chaos["v_belt"] = 2*pi*frequency * amplitude_displacement*np.cos(2*pi*frequency*df_chaos.t)
    df_chaos["slip"] = df_chaos.x - df_chaos.x_belt

    fig = plt.figure(figsize=(8,6))
    
    ax1 = fig.add_subplot(221)
    ax1.plot(df_regular.x, df_regular.v, label="Block", linestyle="-", linewidth=0.75)
    ax1.plot(df_chaos.x_belt, df_chaos.v_belt, label="Belt", linewidth=0.75)
    ax1.legend()
    ax1.set_xlabel("Displacement [mm]")
    ax1.set_ylabel("Velocity [mm/s]")
    ax1.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
    ax1.legend(loc="upper left")
    ax1.grid()

    ax2 = fig.add_subplot(222)
    ax2.plot(df_regular.t, df_regular.slip, label="Slip")
    ax2.plot(df_regular.t, df_regular.x_belt, label="Belt")
    ax2.set_xlim((df_regular.t.iloc[0], df_regular.t.iloc[0]+2/frequency))
    ax2.set_xlabel("Time [s]")
    ax2.set_ylabel("Displacement [mm]")
    ax2.legend(loc="upper left")
    ax2.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    ax2.grid()

    ax3 = fig.add_subplot(223)
    ax3.plot(df_chaos.x, df_chaos.v, label="Block", linestyle="-", linewidth=0.75)
    ax3.plot(df_chaos.x_belt, df_chaos.v_belt, label="Belt", linewidth=0.75)
    ax3.legend()
    ax3.set_xlabel("Displacement [mm]")
    ax3.set_ylabel("Velocity [mm/s]")
    ax3.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
    ax3.legend(loc="upper left")
    ax3.grid()

    ax4 = fig.add_subplot(224)
    ax4.plot(df_chaos.t, df_chaos.slip, label="Slip")
    ax4.plot(df_chaos.t, df_chaos.x_belt, label="Belt")
    ax4.set_xlim((df_chaos.t.iloc[0], df_chaos.t.iloc[0]+2/frequency))
    ax4.set_xlabel("Time [s]")
    ax4.set_ylabel("Displacement [mm]")
    ax4.legend(loc="upper left")
    ax4.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    ax4.grid()
    
    
    fig.suptitle(figure_title)
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    if save_figure: fig.savefig(figure_path)
    plt.show()


if __name__ == "__main__":
    main()