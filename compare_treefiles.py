import sys

import numpy as np
import tskit


def pairwise_distance_branch(ts: tskit.TreeSequence, samples: np.array):
    sample_sets = []
    indexes = []
    for i in range(len(samples)):
        sample_sets.append([i])
        for j in range(i + 1, len(samples)):
            indexes.append((i, j))

    div = ts.divergence(sample_sets, indexes=indexes, mode="branch")
    return div


treefile0 = sys.argv[1]
treefile1 = sys.argv[2]


ts0 = tskit.load(treefile0)
ts1 = tskit.load(treefile1)
print(ts0.num_trees, ts1.num_trees)

# for i,j in zip([treefile0, treefile1], [ts0, ts1]):
#     fn = i.replace("trees","svg")
#     print(fn)
#     ts0.first().draw(path=fn, height=5000,width=5000,format="svg")

assert len(ts0.tables.nodes) == len(
    ts1.tables.nodes
), f"{len(ts0.tables.nodes)} {len(ts1.tables.nodes)}"
assert len(ts0.tables.edges) == len(
    ts1.tables.edges
), f"{len(ts0.tables.edges)} {len(ts1.tables.edges)}"

print(
    f"Comparing node times: {np.array_equal(ts0.tables.nodes.time, ts1.tables.nodes.time)}"
)

pd0 = pairwise_distance_branch(ts0, [i for i in ts0.samples()])
pd1 = pairwise_distance_branch(ts1, [i for i in ts1.samples()])
#
print(f"Comparing pairwise distance matrix: {np.array_equal(pd0, pd1)}")
