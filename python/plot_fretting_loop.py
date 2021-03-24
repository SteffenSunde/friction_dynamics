#%% Import dependencies
import matplotlib.pyplot as plt
import pandas as pd

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
    file = "data/mdof/stiff/mdof_fretting_loop_f15.000000_P10.000000_xi0.050000_del1.000000.csv"
    df = pd.read_csv(file, skiprows=1)
    df.columns = ["Q", "d"]

    meta = parse_header(file)
    frequency = float(meta["f"])
    pressure = float(meta["p"])
    damping = float(meta["xi"])
    slope = float(meta["del"])

    # Attempt to read proper displacement
    file_displ = "data/mdof/stiff/mdof_slip_history_f15.000000_P10.000000_xi0.050000_del1.000000.csv"
    df_displ = pd.read_csv(file_displ, skiprows=1)
    df_displ.columns = ["time", "slip_end", "slip_edge", "slip_center", "displ"]

    fig, ax = plt.subplots(figsize=(5,4))
    #plt.plot(df.d, df.Q)
    #plt.plot(df_displ.displ, df.Q)  # Is this more strictly correct hysteresis loop?? weird shape nonetheless
    ax.set_xlabel(r"Displacement [mm]")
    ax.set_ylabel("Shear force [N]")
    ax.grid()
    ax.ticklabel_format(axis='both', style='sci', scilimits=(0,0))

    bbox_props = dict(boxstyle="round", fc="w", ec="0.5", alpha=0.8)
    text_props = dict(fontsize=10)
    text = r"$p_0$={:.1f}".format(pressure) \
        + "\n" + r"$\xi$={:.2f}".format(damping) \
        + "\n" + r"$\delta$={:.1f}".format(slope)
    ax.text(0.006, -350, text , bbox=bbox_props, fontdict=text_props)


    #fig.savefig("data/fretting_map.pdf")
    plt.show()


#%%
if __name__ == "__main__":
    main()