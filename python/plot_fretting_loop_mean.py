#%% Import dependencies
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from matplotlib import rcParams
#rcParams.update({'font.size': 14})

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


#%% Import data
def main():

    file = "data/mdof/mass_damped/massdamped_mdof_fretting_loop_f15.000000_P10.000000_xi0.050000_del1.000000.csv"
    
    meta = parse_header(file)
    fig_path = "data/mdof/for_thesis/fretting_loop_f{:.0f}_P{:.0f}_xi{:.4f}_del{:.4f}.pdf".format(float(meta["f"]), float(meta["p"]), float(meta["xi"]), float(meta["del"])) 
    x_lower, y_lower = -0.01, 0.01#-0.01,0.01

    savefig = False

    df = pd.read_csv(file, skiprows=1)
    df.columns = ["Q", "d"]

    # Attempt to read proper displacement
    file_displ = "data/mdof/mass_damped/mdof_slip_history_f15.000000_P10.000000_xi0.050000_del1.000000.csv"
    df_displ = pd.read_csv(file_displ, skiprows=1)
    df_displ.columns = ["time", "slip_end", "slip_edge", "slip_center", "displ"]

    
    print(meta)
    frequency = float(meta["f"])
    pressure = float(meta["p"])
    damping = float(meta["xi"])
    slope = float(meta["del"])

    num_cycles = int(2.0 * frequency + 0.5)
    chunk_size = int(df.shape[0] / num_cycles)

    data = np.zeros((num_cycles-1, chunk_size))
    for i in range(num_cycles-1):
        data[i,:] = df.Q.iloc[i*chunk_size:(i+1)*chunk_size]

    std = data.std(0)
    mean = data.mean(0)
    half = len(mean)//2
    quart = half//2

    #x_cycle = df.d.iloc[:chunk_size]
    x_cycle = df_displ.displ[:chunk_size]
    upper = mean+std
    lower = mean-std

    outer = np.concatenate((upper[:quart], lower[quart:half+quart], upper[half+quart:]))
    inner = np.concatenate((lower[:quart], upper[quart:half+quart], lower[half+quart:]))

    filter_inner = []
    x_inner = []
    for i, val in enumerate(x_cycle):
        if val > x_lower and val < y_lower:  # TODO: For now these are adjusted manually! 
            filter_inner.append(inner[i])
            x_inner.append(x_cycle[i])
        else:
            filter_inner.append(0)
            x_inner.append(0)


    fig, ax = plt.subplots(figsize=(5,4))
    ax.plot(x_cycle, mean, label="Mean cycle (Mean SD {:.2f})".format(std.mean()), color="green")
    #ax.plot(x_cycle, outer, label="Outer")
    #ax.plot(x_inner, filter_inner, label="Inner")
    ax.fill_between(x_cycle, outer, filter_inner, label="SD. (Mean {:.2f})".format(std.mean()), alpha=0.25)
    #plt.fill_between(x_cycle, inner, outer, where=filter_inner, alpha=0.5)
    #ax.set_xticks((-0.01, -0.005, 0.0, 0.005, 0.01))
    ax.set_xlabel(r"Displacement, $d$ [mm]")
    ax.set_ylabel(r"Shear, $Q$ [N]")
    
    ax.legend(loc="upper left")
    ax.grid()

    bbox_props = dict(boxstyle="round", fc="w", ec="0.5", alpha=0.8)
    text_props = dict(fontsize=10)
    text = r"$p_0$={:.1f}".format(pressure) \
        + "\n" + r"$\xi$={:.2f}".format(damping) \
        + "\n" + r"$\delta$={:.1f}".format(slope)
    ax.text(0.60*x_cycle.max(), 0.95*mean.min(), text , bbox=bbox_props, fontdict=text_props)
    
    plt.tight_layout()
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