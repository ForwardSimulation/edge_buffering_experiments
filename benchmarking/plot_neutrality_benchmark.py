import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


def normalize(g):
    if g.name == "mem":
        return g / g.min()
    return g


df = pd.read_csv("neutrality_benchmark.txt", sep=" ")
df["time"] /= 60.0


dfn = df.groupby(["N"]).transform(normalize)

df["mem"] = df["mem"] / (1024 ** 2)

f, (ax1, ax2) = plt.subplots(1, 2, sharex=False, sharey=False)

p = sns.scatterplot(
    "N",
    "time",
    data=df[df.N < 2.5e4],
    hue="simulator",
    style="method",
    alpha=0.6,
    size="simulator",
    sizes=[60, 60, 60, 60],
    ax=ax1,
)

p.set(
    xlabel="N (number of diploids)",
    ylabel="Run time (minutes)",
    xticks=[1000, 5000, 10000],
)

p2 = sns.scatterplot(
    "N",
    "time",
    data=df,
    hue="simulator",
    style="method",
    alpha=0.6,
    size="simulator",
    sizes=[60, 60, 60, 60],
    ax=ax2,
)

p2.set(
    xlabel="N (number of diploids)", ylabel=None, xticks=[1000, 10000, 25000],
)

ax1.get_legend().remove()
ax1.legend(ncol=1)
ax2.get_legend().remove()

f.suptitle("WF, non-overlapping, no selection, no recombination\n5N generations")
p.get_figure().savefig("neutrality_benchmark.png", dpi=300)

f, ax1 = plt.subplots(1, 1, sharex=False, sharey=False)

p = sns.scatterplot(
    "N",
    "mem",
    data=df,
    style="method",
    hue="simulator",
    size="simulator",
    sizes=[60, 60, 60, 60],
    alpha=0.7,
)
p.set(
    xlabel="N (number of diploids)",
    ylabel="Peak RAM (GB)",
    xticks=[1000, 10000, 25000],
)

plt.savefig("neutrality_benchmark_memory.png")
