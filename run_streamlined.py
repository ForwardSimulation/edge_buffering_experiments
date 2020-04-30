import numpy as np

import streamlined_prototype

np.random.seed(666)
# tables, tables2 = streamlined_prototype.wright_fisher(500, 1000, 0.0, 33)
# tables, tables2 = streamlined_prototype.wright_fisher(500, 1000, 0.5, 33)
tables, tables2 = streamlined_prototype.wright_fisher(50, 1000, 0.5, 33)
assert np.array_equal(tables.nodes.time, tables2.nodes.time)
ts = tables.tree_sequence()
ts2 = tables2.tree_sequence()
m = streamlined_prototype.pairwise_distance_branch(ts, [i for i in ts.samples()])
m2 = streamlined_prototype.pairwise_distance_branch(ts2, [i for i in ts2.samples()])
assert np.array_equal(m, m2)
ts.first().draw(format="svg", path="foo.svg", height=5000, width=5000)
ts.first().draw(format="svg", path="foo2.svg", height=5000, width=5000)
