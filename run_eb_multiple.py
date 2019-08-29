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
node_labels = {idmap[i]: ts.tables.nodes.time[idmap[i]] for i in samples}
t = next(ts.trees())
t.draw(path="tree_eb_multiple_step1.svg", format="svg", height=1000,
       width=1000, node_colours=node_colors, node_labels=node_labels)

# Remap the parental nodes
for p in pstate.parents:
    assert p.n0 != tskit.NULL, "NULL node badness"
    assert p.n1 != tskit.NULL, "NULL node badness"
    p.n0 = idmap[p.n0]
    p.n1 = idmap[p.n1]
    assert p.n0 < len(pstate.tables.nodes)
    assert p.n1 < len(pstate.tables.nodes)

# Annoyane arising from recording time forwards:
flags = np.zeros(len(pstate.tables.nodes), dtype=np.uint32)
pstate.tables.nodes.set_columns(
    flags=flags, time=-1.*(pstate.tables.nodes.time - pstate.tables.nodes.time.max()))

# Simulate again
old_node_table_len = len(pstate.tables.nodes)
pstate = wright_fisher.wright_fisher_eb(20, 0.9, pstate)

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
num_edges = len(tcopy.edges)
for eb in pstate.buffered_edges:
    for i in eb[0] + eb[1]:
        tcopy.edges.add_row(*i)
tcopy.sort()
tcopy_num_edges_b4_simplify = len(tcopy.edges)
idmap = tcopy.simplify()
samples = np.where(flags == 1)[0]
ts = tcopy.tree_sequence()
node_colors = {idmap[i]: 'green' for i in samples}
node_labels = {idmap[i]: ts.tables.nodes.time[idmap[i]] for i in samples}
t = next(ts.trees())
t.draw(path="tree_eb_multiple_step2a.svg", format="svg", height=1000,
       width=1000, node_colours=node_colors, node_labels=node_labels)

# for i, eb in enumerate(pstate.buffered_edges):
#     p1, p2 = tskit.NULL, tskit.NULL
#     for j in eb[0]:
#         p1 = j[2]
#         break
#     for j in eb[1]:
#         p2 = j[2]
#         break
#     prev1 = None
#     if p1 != tskit.NULL and p1 < old_node_table_len:
#         prev1 = np.where(pstate.tables.edges.parent == p1)[0]
#         if len(prev1) > 0:
#             print(p1, pstate.tables.edges.parent[prev1])
#     prev2 = None
#     if p2 != tskit.NULL and p2 < old_node_table_len:
#         prev2 = np.where(pstate.tables.edges.parent == p2)[0]
#         if len(prev2) > 0:
#             print(p2, pstate.tables.edges.parent[prev2])

pstate.tables.nodes.set_columns(
    flags=flags, time=-1.0*(pstate.tables.nodes.time - pstate.tables.nodes.time.max()))

new_edges = 0
for e in pstate.buffered_edges:
    for j in e[0] + e[1]:
        new_edges += 1

assert len(pstate.tables.edges) + new_edges == tcopy_num_edges_b4_simplify

temp_edges = tskit.EdgeTable()
temp_edges_from_before = tskit.EdgeTable()

new_edges2 = 0
edges_added = 0
for o in reversed(pstate.generation_offsets):
    print(*o)
    for i in range(*o):
        # Get the parent node IDs
        # NOTE: this could be more efficient
        pnodes = [None, None]
        for n in [0, 1]:
            for j in pstate.buffered_edges[i][n]:
                pnodes[n] = j[2]
                break
        if pnodes[0] is not None and pnodes[1] is not None:
            assert pnodes[0] < pnodes[1]
        for n in [0, 1]:
            if pnodes[n] is not None:
                new_edges2 += len(pstate.buffered_edges[i][n])
                if pnodes[n] < old_node_table_len:
                    for k in reversed(pstate.buffered_edges[i][n]):
                        edges_added += 1
                        temp_edges_from_before.add_row(*k)
                    ex = np.flip(
                        np.where(pstate.tables.edges.parent == pnodes[n])[0])
                    if len(ex) > 0:
                        v = len(temp_edges_from_before)
                        temp_edges_from_before.append_columns(
                            pstate.tables.edges.left[ex],
                            pstate.tables.edges.right[ex],
                            pstate.tables.edges.parent[ex],
                            pstate.tables.edges.child[ex])
                        edges_added += (len(temp_edges_from_before)-v)
                        assert len(temp_edges_from_before) - v == len(ex)
                else:
                    for k in pstate.buffered_edges[i][n]:
                        edges_added += 1
                        temp_edges.add_row(*k)
print("result", edges_added, new_edges2, new_edges2 + len(pstate.tables.edges), tcopy_num_edges_b4_simplify,
      len(temp_edges) + len(temp_edges_from_before))
assert new_edges == new_edges2
print(len(pstate.tables.edges), len(temp_edges), new_edges2)

pstate.tables.edges.set_columns(
    temp_edges.left, temp_edges.right, temp_edges.parent, temp_edges.child)
pstate.tables.edges.append_columns(
    temp_edges_from_before.left, temp_edges_from_before.right, np.flip(
        temp_edges_from_before.parent),
    np.flip(temp_edges_from_before.child))

print(len(pstate.tables.edges), len(temp_edges), new_edges2)

# for i in pstate.tables.edges:
#     print(i, pstate.tables.nodes.time[i.parent])
print(tcopy_num_edges_b4_simplify, len(pstate.tables.edges))
pstate.tables.simplify()
ts2 = pstate.tables.tree_sequence()
print(len(pstate.tables.edges), len(temp_edges), new_edges2)
print(len(ts.tables.edges), len(ts2.tables.edges))
