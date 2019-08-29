import wright_fisher
import numpy as np
import tskit

np.random.seed(333)

pstate = wright_fisher.PopState(100)

pstate = wright_fisher.wright_fisher_eb(2000, 0.9, pstate)

for o in reversed(pstate.generation_offsets):
    for i in range(*o):
        for j in pstate.buffered_edges[i][0]:
            pstate.tables.edges.add_row(*j)
        for j in pstate.buffered_edges[i][1]:
            pstate.tables.edges.add_row(*j)
        # reset the buffer
        for j in pstate.buffered_edges[i]:
            j.clear()

# Reset offsets
pstate.generation_offsets = [(0, len(pstate.buffered_edges))]

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
t.draw(path="tree_eb_multiple_step1.svg", format="svg", height=1000,
       width=1000, node_colours=node_colors)

# Remap the parental nodes
for p in pstate.parents:
    assert p.n0 != tskit.NULL, "NULL node badness"
    assert p.n1 != tskit.NULL, "NULL node badness"
    p.n0 = idmap[p.n0]
    p.n1 = idmap[p.n1]

# Annoyane arising from recording time forwards:
flags = np.zeros(len(pstate.tables.nodes), dtype=np.uint32)
pstate.tables.nodes.set_columns(
    flags=flags, time=-1.0*(pstate.tables.nodes.time + pstate.tables.nodes.time.max()))

# Simulate again
pstate = wright_fisher.wright_fisher_eb(2000, 0.9, pstate)

# Copy existing tables and do the straigtforward thing
tcopy = tskit.TableCollection(pstate.tables.sequence_length)
flags = np.zeros(len(pstate.tables.nodes), dtype=np.uint32)
for p in pstate.parents:
    flags[p.n0] = 1
    flags[p.n1] = 1

tcopy.nodes.set_columns(
    flags=flags, time=-1.0*(pstate.tables.nodes.time - pstate.tables.nodes.time.max()))

tcopy.edges.set_columns(pstate.tables.edges.left, pstate.tables.edges.right,
                        pstate.tables.edges.parent, pstate.tables.edges.child)

for eb in pstate.buffered_edges:
    for i in eb[0] + eb[1]:
        tcopy.edges.add_row(*i)

tcopy.sort()
idmap = tcopy.simplify()
samples = np.where(flags == 1)[0]
ts = tcopy.tree_sequence()
node_colors = {idmap[i]: 'green' for i in samples}
t = next(ts.trees())
t.draw(path="tree_eb_multiple_step2a.svg", format="svg", height=1000,
       width=1000, node_colours=node_colors)
