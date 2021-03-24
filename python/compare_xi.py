import matplotlib.pyplot as plt
import matplotlib
import scipy.fftpack
from math import pi
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d, Axes3D
import glob
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

    files = glob.glob("data/sdof/newdamp/compare_xi/history*xi*")

    dataframes = []
    metas = []
    for f in files:
        metas.append(parse_header(f))
        amplitude_displacement = float(metas[0]['disp'])
        frequency = float(metas[0]['f'])
        df = pd.read_csv(f, skiprows=1)
        df.columns = ["x", "v", "t"]
        df["x_belt"] = amplitude_displacement * np.sin(2*pi*frequency*df.t)
        df["v_belt"] = 2*pi*frequency * amplitude_displacement*np.cos(2*pi*frequency*df.t)
        df["slip"] = -df.x + df.x_belt
        dataframes.append(df)
        
    # Verify that cases do not have different settings
    for meta in metas:
        for key,val in metas[0].items():
            if meta[key] != val:
                print("Difference in {}: {} vs {}".format(key, val, meta[key]))

    figure_title = r"Comparing Damping ratios"

    save_figure = True
    figure_path = "data/compare_xi.pdf"

    fig = plt.figure(figsize=(8,4))

    ax2 = fig.add_subplot(121)
    ax2.plot(dataframes[0].x_belt, dataframes[0].v_belt, label="Belt")
    for i, df in enumerate(dataframes):
        ax2.plot(df.x, df.v, label=r"$\xi={:.3f}$".format(float(metas[i]["xi"])))
    ax2.legend(loc="lower left")
    ax2.set_xlabel("Displacement [mm]")
    ax2.set_ylabel("Velocity [mm/s]")
    ax2.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
    ax2.set_ylim((-4, 1.05))
    ax2.set_xlim((0, 0.0101))
    ax2.grid()

    ax3 = fig.add_subplot(122)
    ax3.plot(dataframes[0].t, dataframes[0].x_belt, label="Belt")
    for i, df in enumerate(dataframes):
        ax3.plot(df.t, df.slip, label=r"$\xi={:.3f}$".format(float(metas[i]["xi"])))
    ax3.set_xlim((df.t.iloc[0]+0.01, df.t.iloc[0]+0.5/frequency+0.01))
    ax3.set_xlabel("Time [s]")
    ax3.set_ylabel("Slip [mm]")
    ax3.legend(loc="lower left")
    ax3.set_xticks([ax3.get_xticks()[1], ax3.get_xticks()[-1]])
    ax3.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax3.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    ax3.grid()

    #fig.suptitle(figure_title)
    
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    if save_figure: fig.savefig(figure_path)
    plt.show()


if __name__ == "__main__":
    main()