import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

df = pd.read_csv("neutrality_benchmark_rho_new_simplification.txt", sep=" ")


def normalize(g):
    if g.name == "mem":
        return g / g.min()
    return g


df["time"] /= 60.0
df["rho"] = np.log10(df["rho"])


dfn = df.groupby(["N", "rho"]).transform(normalize)

df["mem"] = dfn["mem"]

g = sns.FacetGrid(df, col="N", sharey=False)

g = (
    g.map(
        sns.scatterplot,
        "rho",
        "time",
        data=df,
        style="method",
        hue="simulator",
        alpha=0.8,
    )
    .add_legend()
    .set_axis_labels(r"$log_{10}(4Nr)$", "Run time (minutes)")
)

plt.savefig("neutrality_benchmark_rec_new_simplification.png")

g = sns.FacetGrid(df, col="N", sharey=False)

g = (
    g.map(
        sns.scatterplot,
        "rho",
        "mem",
        data=df,
        style="method",
        hue="simulator",
        alpha=0.8,
    )
    .add_legend()
    .set_axis_labels(r"$log_{10}(4Nr)$", "Relative memory use.")
)

plt.savefig("neutrality_benchmark_rec_new_simplification_memory.png")
