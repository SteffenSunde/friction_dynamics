import matplotlib.pyplot as plt
import pandas as pd

def main():
    file = "data/rate_perturbed_poincare/poincare_20.000000_100_1000000.txt"
    df = pd.read_csv(file)
    df.columns = ["x", "v"]  # TODO: Colorize according to initial condition
    fig, ax = plt.subplots()
    for i in range(100):
        ax.scatter(df.x[101*i:100*i+100], df.v[100*i:100*i+100], s=0.25)
    plt.show()


if __name__ == "__main__":
    main()