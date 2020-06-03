import pandas as pd
import seaborn as sns

df = pd.read_csv("neutrality_benchmark.txt", sep=" ")

df = df[
    (df.simulator != "slim") & (df.method != "cppsort_par") & (df.method != "buffer")
]
print(df)

df["time"] /= 60.0

p = sns.scatterplot(
    "N",
    "time",
    data=df,
    hue="simulator",
    style="method",
    alpha=0.8,
    size="simulator",
    sizes=[60, 60],
)

p.set(
    xlabel="N (number of diploids)",
    ylabel="Run time (minutes)",
    title="WF, non-overlapping, no selection, no recombination\n5N generations",
)

p.get_figure().savefig("neutrality_benchmark_simpler.png")
