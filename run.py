import numpy as np

import wright_fisher

np.random.seed(333)

pstate = wright_fisher.PopState(100)

pstate = wright_fisher.wright_fisher(2000, 0.9, pstate)

flags = np.zeros(len(pstate.tables.nodes), dtype=np.uint32)
for p in pstate.parents:
    flags[p.n0] = 1
    flags[p.n1] = 1

pstate.tables.nodes.set_columns(
    flags=flags, time=-1.0 * (pstate.tables.nodes.time - pstate.tables.nodes.time.max())
)

pstate.tables.sort()
idmap = pstate.tables.simplify()
samples = np.where(flags == 1)[0]
ts = pstate.tables.tree_sequence()
node_colors = {idmap[i]: "green" for i in samples}
t = next(ts.trees())
t.draw(path="tree.svg", format="svg", height=1000, width=1000, node_colours=node_colors)
