import pandas as pd
import seaborn as sns

df = pd.read_csv("neutrality_time_sorting.txt", sep=" ")

p = sns.scatterplot(
    "N",
    "percent",
    data=df,
    hue="simulator",
    style="method",
    alpha=0.8,
    size="simulator",
    sizes=[60, 60],
)

p.set(
    xlabel="N (number of diploids)",
    ylabel="Percent of run time spent sorting",
    title="WF, non-overlapping, no selection, no recombination\n5N generations",
    ylim=(0.0, 100.0),
)

p.get_figure().savefig("neutrality_time_sorting.png")
