import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from math import sin, pi

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
    file = "data/mdof/evolve/fretting_loop_f15.000000_p10.000000_del1.000000_xi0.050000.csv"
    savefig = False
    fig_path = "data/evolved_loop_p10_rate001.pdf"

    df = pd.read_csv(file, skiprows=1)
    df.columns = ["time", "Q", "mu_edge", "mu_half", "mu_center"]
    meta = parse_header(file)
    num_steps = int(meta["num_loops"])
    frequency = float(meta["f"])
    stroke = float(meta["disp"])
    df["d"] = stroke*np.sin(2*pi*frequency*df.time)

    num_points = df.shape[0]
    step_size = int(num_points/num_steps)
    num_cycles = 5
    cycle_size = int(step_size/num_cycles)

    df = df.iloc[:cycle_size*num_cycles*num_steps]
    data = df.Q.values.reshape(num_cycles*num_steps, cycle_size)
    print(data.shape)

    cycle_qs = []
    for i in range(num_steps):
        start = i*num_cycles
        end = (i+1)*num_cycles
        mean = data[start:end, :].mean(axis=0)
        cycle_qs.append(mean)

    # plt.plot(df.time)
    # plt.show()

    fig, (ax1,ax2) = plt.subplots(1, 2, figsize=(8,4))
    for i in range(num_steps):
        time = df.time.iloc[i*step_size+i]
        ax1.plot(df.d.iloc[:cycle_size], cycle_qs[i], label="{:.0f} s".format(time))

    ax1.set_xlabel(r"Displacement [mm]")
    ax1.set_ylabel("Shear force [N]")
    ax1.grid()
    ax1.legend()
    #ax1.ticklabel_format(axis='both', style='sci', scilimits=(0,0))

    ax2.plot(df.time, df.mu_edge, label="Edge")
    ax2.plot(df.time, df.mu_half, label="Quart")
    ax2.plot(df.time, df.mu_center, label="Center")
    ax2.grid()
    ax2.set_xlabel("Time [s]")
    ax2.set_ylabel("COF [-]")
    ax2.legend()

    fig.tight_layout()
    if savefig: fig.savefig(fig_path)
    plt.show()


if __name__ == "__main__":
    main()