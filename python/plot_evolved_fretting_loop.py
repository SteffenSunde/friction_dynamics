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

def main():
    file = "data/mdof/evolve/fretting_loop_f15.000000_del5.000000_xi0.050000.csv"
    df = pd.read_csv(file, skiprows=1)
    df.columns = ["Q", "d", "mu_edge", "mu_center"]

    meta = parse_header(file)
    num_loops = int(meta["num_loops"])
    num_points = df.shape[0]
    loop_size = int(num_points/num_loops)

    plt.plot(df.mu_edge)

    fig, ax = plt.subplots(figsize=(5,4))
    for i in range(num_loops):
        start = i*loop_size
        end = (i+1)*loop_size
        plt.plot(df.d.iloc[start:end], df.Q.iloc[start:end], label="{}".format(i))
    
    ax.set_xlabel(r"Displacement [mm]")
    ax.set_ylabel("Shear force [N]")
    ax.grid()
    ax.legend()
    ax.ticklabel_format(axis='both', style='sci', scilimits=(0,0))

    #fig.savefig("data/fretting_map.pdf")
    plt.show()


if __name__ == "__main__":
    main()