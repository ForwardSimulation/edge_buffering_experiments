import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

df = pd.read_csv("neutrality_benchmark_rho.txt", sep=" ")

df["time"] /= 60.0
df["rho"] = np.log10(df["rho"])

g = sns.FacetGrid(df, col="N", sharey=False)

g = (
    g.map(sns.scatterplot, "rho", "time", data=df, hue="method", alpha=0.8,)
    .add_legend()
    .set_axis_labels(r"$log_{10}(4Nr)$", "Run time (minutes)")
)

plt.savefig("neutrality_benchmark_rec.png")
