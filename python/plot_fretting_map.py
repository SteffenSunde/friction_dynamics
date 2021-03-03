import matplotlib.pyplot as plt
import pandas as pd

def main():
    file = "data/fretting_map.csv"
    df = pd.read_csv(file)
    df.columns = ["Q", "d"]  # TODO: Colorize according to initial condition
    fig, ax = plt.subplots()
    plt.plot(df.d, df.Q)
    plt.show()


if __name__ == "__main__":
    main()