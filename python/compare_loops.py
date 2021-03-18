#%% Import dependencies
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from matplotlib import rcParams

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

    file = "data/mdof/mdof_fretting_loop15.000000_P15.000000_xi0.050000.csv"
    fig_path = "data/mdof/fretting_cycle_f15_P15_xi005.pdf"
    savefig = True
    
    df = pd.read_csv(file, skiprows=1)
    df.columns = ["Q", "d"]
    meta = parse_header(file)
    frequency = float(meta["f"])

    num_cycles = int(10.0 * frequency + 0.5)
    chunk_size = int(df.shape[0] / num_cycles)

    data = np.zeros((num_cycles-1, chunk_size))
    for i in range(num_cycles-1):
        data[i,:] = df.Q.iloc[i*chunk_size:(i+1)*chunk_size]

    fig, ax = plt.subplots(figsize=(4,3))
    ax.plot(df.d.iloc[:chunk_size], data[i,:], color="green")
    ax.set_xticks((-0.01, 0.0, 0.01))
    ax.set_xlabel(r"Displacement, $\delta$ [mm]")
    ax.set_ylabel(r"Shear, $Q$ [N]")
    plt.tight_layout()
    ax.grid()
    if savefig: fig.savefig(fig_path)
    plt.show()
    
    return

    ## Trash to follow
    
    plt.plot(x_cycle[:half], upper[:half], label="Upper")
    plt.plot(x_cycle[:half], lower[:half], label="Lower")
    # plt.plot(x_cycle, outer, label="Outer")
    # plt.plot(x_cycle, inner, label="Inner")
    plt.fill_between(x_cycle, inner, alpha=0.5)
    plt.plot(x_cycle, mean, label="Mean")
    plt.legend()
    plt.show()

    outer = np.concatenate((lower[:half], upper[half:]))
    inner = np.concatenate((upper[:half], lower[half:]))
    
    # plt.plot(x_cycle, upper)
    # plt.plot(x_cycle, lower)
    # plt.show()

    fig, ax = plt.subplots(figsize=(5,4))
    plt.plot(x_cycle, mean, label="Mean")
    plt.plot(x_cycle, mean+std, label="+std")
    plt.plot(x_cycle, mean-std, label="-std")
    plt.plot(x_cycle, outer, label="Outer")
    plt.plot(x_cycle, inner, label="Inner")

    ax.fill(x_cycle, outer, alpha=0.5)
    #plt.fill_between(df.d[:chunk_size//2], mean[:half]+std[:half], alpha=0.5, zorder=1000)
    #plt.fill_between(df.d[:chunk_size//2], 0, mean[half:]+std[half:], alpha=0.5, zorder=1000)

    ax.set_xlabel(r"Displacement [mm]")
    ax.set_ylabel("Shear force [N]")
    ax.set_xticks([-0.01, 0.0, 0.01])
    ax.legend()
    #ax.set_yticks([df.Q.min(), 0.0, df.Q.max()])
    ax.ticklabel_format(axis='both', style='sci', scilimits=(0,0))
    ax.grid()

    #fig.savefig("data/fretting_map.pdf")
    plt.show()




#%%
if __name__ == "__main__":
    main()