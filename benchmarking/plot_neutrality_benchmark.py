import pandas as pd
import seaborn as sns

df = pd.read_csv("neutrality_benchmark.txt", sep=" ")

df["time"] /= 60.

p = sns.scatterplot("N", "time", data=df, hue="simulator", style="method")
p.set(
    xlabel="N (number of diploids)",
    ylabel="Run time (minutes)",
    title="WF, non-overlapping, no selection, no recombination\n5N generations",
)
p.get_figure().savefig("neutrality_benchmark.png")
