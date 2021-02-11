import matplotlib.pyplot as plt
import pandas as pd

def main():
    file = "data/poincare_263.396127.txt"
    df = pd.read_csv(file)
    df.columns = ["x", "v"]  # TODO: Colorize according to initial condition
    fig, ax = plt.subplots()
    # for i in range(100):
    #     ax.scatter(df.x[101*i:100*i+100], df.v[100*i:100*i+100], s=0.25)
    ax.scatter(df.x, df.v, s=0.25)
    plt.show()


if __name__ == "__main__":
    main()