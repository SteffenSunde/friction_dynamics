import matplotlib.pyplot as plt
import matplotlib
import scipy.fftpack
from math import pi
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d, Axes3D
import glob

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

    files = glob.glob("data/trans/history_single*trans*")

    dataframes = []
    metas = []
    for i, f in enumerate(files):
        metas.append(parse_header(f))
        amplitude_displacement = float(metas[0]['disp'])
        frequency = float(metas[0]['f'])
        df = pd.read_csv(f, skiprows=1)
        df.columns = ["x", "v", "t"]
        df["x_belt"] = amplitude_displacement * np.sin(2*pi*frequency*df.t)
        df["v_belt"] = 2*pi*frequency * amplitude_displacement*np.cos(2*pi*frequency*df.t)
        df["slip"] = -df.x + df.x_belt
        metas[i]["trans"] = f.split("\\")[-1].split("trans")[-1][:-4]
        dataframes.append(df)
        
    print(files)

    # Check for differences
    for meta in metas:
        for key,val in metas[0].items():
            if meta[key] != val:
                print("Difference in {}: {} vs {}".format(key, val, meta[key]))

    figure_title = r"Comparing Transient discards"

    save_figure = False
    figure_path = "data/compare_trans.pdf"

    fig = plt.figure(figsize=(12,4))

    ax2 = fig.add_subplot(121)
    
    ax2.plot(df.x_belt, df.v_belt, label="Belt", linewidth=0.75)
    
    for i, df in enumerate(dataframes):
        ax2.plot(df.x, df.v, label="Trans={}".format(metas[i]["trans"]), linestyle="-", linewidth=0.75)
    
    ax2.legend(loc="upper left")
    ax2.set_xlabel("Displacement [mm]")
    ax2.set_ylabel("Velocity [mm/s]")
    ax2.ticklabel_format(style='sci', axis='both', scilimits=(0,0))

    ax3 = fig.add_subplot(122)
    ax3.plot(df.t, df.x_belt, label="Belt")

    for i, df in enumerate(dataframes):
        ax3.plot(df.t, df.slip, label="Trans={}".format(metas[i]["trans"]))

    ax3.set_xlim((df.t.iloc[0], df.t.iloc[0]+2/frequency))
    ax3.set_xlabel("Time [s]")
    ax3.set_ylabel("Displacement [mm]")
    ax3.legend(loc="upper left")

    fig.suptitle(figure_title)
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    if save_figure: fig.savefig(figure_path)
    plt.show()


if __name__ == "__main__":
    main()