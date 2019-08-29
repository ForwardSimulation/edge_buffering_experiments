import wright_fisher
import numpy as np

np.random.seed(333)

pstate = wright_fisher.PopState(100)

pstate = wright_fisher.wright_fisher_eb(2000, 0.9, pstate)

for o in reversed(pstate.generation_offsets):
    for i in range(*o):
        for j in pstate.buffered_edges[i][0]:
            pstate.tables.edges.add_row(*j)
        for j in pstate.buffered_edges[i][1]:
            pstate.tables.edges.add_row(*j)


# for eb in pstate.buffered_edges:
#     print(eb)

flags = np.zeros(len(pstate.tables.nodes), dtype=np.uint32)
for p in pstate.parents:
    flags[p.n0] = 1
    flags[p.n1] = 1

pstate.tables.nodes.set_columns(
    flags=flags, time=-1.0*(pstate.tables.nodes.time - pstate.tables.nodes.time.max()))

idmap = pstate.tables.simplify()
samples = np.where(flags == 1)[0]
ts = pstate.tables.tree_sequence()
node_colors = {idmap[i]: 'green' for i in samples}
t = next(ts.trees())
t.draw(path="tree_eb.svg", format="svg", height=1000,
       width=1000, node_colours=node_colors)
