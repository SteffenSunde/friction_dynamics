#%% Import dependencies
import matplotlib.pyplot as plt
import pandas as pd

#%% Import data
def main():
    file = "data/mdof/stiff_fretting_loop_P15.000000.csv"
    df = pd.read_csv(file, skiprows=1)
    df.columns = ["Q", "d"]

    fig, ax = plt.subplots(figsize=(5,4))
    plt.plot(df.d, df.Q)
    ax.set_xlabel(r"Displacement [mm]")
    ax.set_ylabel("Shear force [N]")
    ax.grid()
    ax.ticklabel_format(axis='both', style='sci', scilimits=(0,0))

    #fig.savefig("data/fretting_map.pdf")
    plt.show()




#%%
if __name__ == "__main__":
    main()