import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

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
    file = "data/poincare_single.csv"
    save_fig = False
    #fig_path = "data/poincare_f10_chaos.png"
    num_colors = 4

    meta = parse_header(file)
    figure_title = "f: {}, del: {}, xi: {}".format(meta["f"], meta["del"], meta["xi"])

    #parse_header(file)
    df = pd.read_csv(file, skiprows=1)
    df.columns = ["x", "v", "t"]  # TODO: Colorize according to initial condition?
    print("Max step: ", np.diff(df.t).max())
    print("Min step: ", np.diff(df.t).min())
    fig, ax = plt.subplots()
    ax.scatter(df.x, df.v, color="black", s=0.25)  #c=range(df.shape[0])
    ax.set_xlabel(r"$x$", fontsize=14)
    ax.set_ylabel(r"$\dot{x}$", fontsize=14)

    fig.suptitle(figure_title)
    plt.tight_layout()
    if save_fig: fig.savefig(fig_path, dpi=500)
    plt.show()


if __name__ == "__main__":
    main()