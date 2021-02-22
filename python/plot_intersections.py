import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

def main():
    file = "data/poincare_single.csv"
    save_fig = False
    #fig_path = "data/poincare_f10_chaos.png"
    num_colors = 4

    #parse_header(file)
    df = pd.read_csv(file, skiprows=1)
    df.columns = ["x", "v", "t"]  # TODO: Colorize according to initial condition
    print(np.diff(df.t).max())
    fig, ax = plt.subplots()
    # ax.scatter(df.x, df.v, c=range(df.shape[0]), s=0.01)
    # ax.set_xlabel(r"$x$", fontsize=14)
    # ax.set_ylabel(r"$\dot{x}$", fontsize=14)
    # for i in range(num_colors):
    #     size = df.shape[0]//num_colors
    #     start = i * size
    #     stop = (i+1)*size
    #     ax.scatter(df.x[start:stop], df.v[start:stop], s=0.01)
    
    #ax.scatter(df.x, df.v, s=0.01, color="black")
    # ax.set_xlim((-0.01, 0.01))
    # ax.set_ylim((0.0, 1.5))
    #ax.set_aspect("equal")
    plt.tight_layout()
    if save_fig: fig.savefig(fig_path, dpi=500)
    plt.show()


if __name__ == "__main__":
    main()