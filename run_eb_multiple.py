import wright_fisher
import numpy as np
import tskit
import sys

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

with open('old_edges.txt', 'w') as f:
    for e in pstate.tables.edges:
        f.write(f'{e} {pstate.tables.nodes.time[e.parent]}\n')
# NOTE: the cleanup steps below should be part of a more general
# API

# Remap the parental nodes
for p in pstate.parents:
    assert p.n0 != tskit.NULL, "NULL node badness"
    assert p.n1 != tskit.NULL, "NULL node badness"
    p.n0 = idmap[p.n0]
    p.n1 = idmap[p.n1]
    assert p.n0 < len(pstate.tables.nodes)
    assert p.n1 < len(pstate.tables.nodes)

# We need to know when the parental nodes first appear in the edge table
pwhere = [None] * len(pstate.parents)
edge_proxy = np.rec.fromarrays((pstate.tables.nodes.time[pstate.tables.edges.parent],
                                pstate.tables.edges.parent),
                               dtype=[('time', np.float64), ('parent', np.int32)])
with open("before.txt", 'w') as f:
    for i, p in enumerate(pstate.parents):
        l0 = np.where(pstate.tables.edges.parent == p.n0)[0]
        l1 = np.where(pstate.tables.edges.parent == p.n1)[0]
        f.write(f"{p.n0} -> {l0}, {p.n1} -> {l1}\n")
        ptime = pstate.tables.nodes.time[p.n0]
        loc0 = len(pstate.tables.edges)
        loc1 = len(pstate.tables.edges)
        is_edge0 = False
        is_edge1 = False
        if len(l0) > 0:
            loc0 = l0[0]
            is_edge0 = True
        if len(l1) > 0:
            loc1 = l1[0]
            is_edge1 = True

        loc0a = edge_proxy.searchsorted(
            np.array((ptime, p.n0), dtype=edge_proxy.dtype))
        loc1a = edge_proxy.searchsorted(
            np.array((ptime, p.n1), dtype=edge_proxy.dtype))
        loc0 = loc0a
        loc1 = loc1a

        print(f"loc: {p.n0} {p.n1}-> {loc0a} {loc1a}")

        pwhere[i] = (p, (ptime, loc0, loc1, is_edge0, is_edge1))
# Relies on Python's sort being a stable sort!
pwhere = sorted(pwhere,
                key=lambda x: (x[1][0], min(x[1][1], x[1][2])))
pstate.parents = [i[0] for i in pwhere]
pwhere = [i[1] for i in pwhere]
with open("after.txt", 'w') as f:
    for i, p in enumerate(pstate.parents):
        p.index = i
        f.write(
            f"{p.n0} {p.n1}-> {pwhere[p.index]}, {pstate.tables.nodes.time[p.n0]}  {pstate.tables.nodes.time[p.n1]}\n")

# Show we can rebuild
E = 0
e = tskit.tables.EdgeTable()
for p, w in zip(pstate.parents, pwhere):
    for i in [0, 1]:
        while E < len(pstate.tables.edges) and E <= w[i+1] and w[i+1] < len(pstate.tables.edges):
            e.add_row(pstate.tables.edges[E].left,
                      pstate.tables.edges[E].right,
                      pstate.tables.edges[E].parent,
                      pstate.tables.edges[E].child)
            E += 1

while E < len(pstate.tables.edges):
    e.add_row(pstate.tables.edges[E].left,
              pstate.tables.edges[E].right,
              pstate.tables.edges[E].parent,
              pstate.tables.edges[E].child)
    E += 1

assert e == pstate.tables.edges


# Sort the parent tracker on birth order and reindex
# pstate.parents = sorted(pstate.parents, key=lambda x: (
#     pstate.tables.nodes.time[x.n0], x.n0))
# for i, p in enumerate(pstate.parents):
#     p.index = i
# assert max([i.index for i in pstate.parents]) < len(pstate.parents)

# Reset the buffer and the index
pstate.buffered_edges = [[[], []] for i in range(len(pstate.parents))]
pstate.next_parent = len(pstate.buffered_edges)
# Reset offsets
pstate.generation_offsets = [(0, len(pstate.buffered_edges))]
# Reset master parent node list
pstate.pnodes = [(i.n0, i.n1) for i in pstate.parents]
# print(len(pstate.parents), len(pstate.buffered_edges),
#       pstate.next_parent, pstate.generation_offsets)

# Annoyance arising from recording time forwards:
flags = np.zeros(len(pstate.tables.nodes), dtype=np.uint32)
ot = np.copy(pstate.tables.nodes.time)
newtime = np.array([pstate.current_generation -
                    i for i in pstate.tables.nodes.time])
pstate.tables.nodes.set_columns(flags=flags, time=newtime)

# Validate that we've not screwed up parent birth times
# NOTE: gotta thing forwards
for i in range(2, len(pstate.tables.edges)):
    e = pstate.tables.edges[i-1]
    ee = pstate.tables.edges[i]
    t = pstate.tables.nodes.time[e.parent]
    tt = pstate.tables.nodes.time[ee.parent]
    if t < tt:
        print(i, e, ee, t, tt)
        raise RuntimeError("messed up parent birth order")

with open("times.txt", "w") as f:
    for i, j in zip(ot, newtime):
        f.write(f"{i} {j}\n")

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
E = 0
# NOTE: issue now is corner case of pre-existing edge
# not in edge table, but then leaving descendants later
for o in reversed(pstate.generation_offsets):
    print("range =", *o)
    for i in range(*o):
        # Fetch the parent node IDs
        pnodes = pstate.pnodes[i]
        if i < len(pwhere):
            pnodes = pstate.pnodes[i]
            where0 = pwhere[i][1]
            where1 = pwhere[i][2]
            print(pnodes)
            # Add in all pre-existing edges ancestral to pnodes[0]
            while E < len(pstate.tables.edges) and E < where0:
                e = pstate.tables.edges[E]
                print(f"0: {pnodes} {e}")
                temp_edges_from_before.add_row(e.left, e.right, e.parent, e.child)
                E += 1

            # Add in all pre-exising edges descending from pnodes[0]
            while E < len(pstate.tables.edges) and pstate.tables.edges.parent[E] == pnodes[0]:
                e = pstate.tables.edges[E]
                print(f"0a: {pnodes} {e}")
                temp_edges_from_before.add_row(e.left, e.right, e.parent, e.child)
                E += 1

            for k in pstate.buffered_edges[i][0]:
                temp_edges_from_before.add_row(*k)
                print(f"I want to add {k} from {pnodes[0]}")

            # Add in all pre-existing edges ancestral to pnodes[1]
            while E < len(pstate.tables.edges) and E < where1:
                e = pstate.tables.edges[E]
                print(f"0b: {pnodes} {e}")
                temp_edges_from_before.add_row(e.left, e.right, e.parent, e.child)
                E += 1

            # Add in all pre-exising edges descending from pnodes[1]
            while E < len(pstate.tables.edges) and pstate.tables.edges.parent[E] == pnodes[1]:
                e = pstate.tables.edges[E]
                print(f"0c: {pnodes} {e}")
                temp_edges_from_before.add_row(e.left, e.right, e.parent, e.child)
                E += 1

            for k in pstate.buffered_edges[i][1]:
                temp_edges_from_before.add_row(*k)
                print(f"I want to add {k} from {pnodes[1]}")
        else:
            for n in [0, 1]:
                for k in pstate.buffered_edges[i][n]:
                    assert k[2] == pnodes[n]
                    edges_added += 1
                    temp_edges.add_row(*k)
                    print(f"4: adding {temp_edges[-1]}")

        # for n in [0, 1]:
        #     # NOTE: the error here is that we do not process
        #     # previous edge table entries w.r.to all in pnodes
        #     # prior to adding in the buffered edges
        #     # In the inner while loops, the correct procedure
        #     # should be:
        #     # 1. Add in all edges whose parents predate
        #     #    any of the two parental nodes.
        #     # 2. Add in all pre-existing edges from parent node 1
        #     # 3. Add in all new edges from parent node 1
        #     # 4. Add in all pre-exising edges from parent node 2
        #     # 5. Add in all new edges from parent node 2
        #     if pnodes[n] < len(pstate.tables.edges):
        #         ptime = pstate.tables.nodes.time[pnodes[0]]
        #         assert ptime == pstate.tables.nodes.time[pnodes[1]]
        #         if pwhere[i][1+n] < len(pstate.tables.edges):
        #             print(
        #                 f"Taking care of {pnodes[n]} {old_node_table_len} {n} {pwhere[i]}")
        #             while E < len(pstate.tables.edges) and E < pwhere[i][1+n]:
        #                 print(
        #                     f"0: adding {pstate.tables.edges[E]} {pstate.tables.nodes.time[pstate.tables.edges[E].parent]}")
        #                 temp_edges_from_before.add_row(
        #                     pstate.tables.edges[E].left,
        #                     pstate.tables.edges[E].right,
        #                     pstate.tables.edges[E].parent,
        #                     pstate.tables.edges[E].child)
        #                 E += 1
        #             while E < len(pstate.tables.edges) and pstate.tables.edges.parent[E] == pnodes[n]:
        #                 print(
        #                     f"1: adding {pstate.tables.edges[E]} {pstate.tables.nodes.time[pstate.tables.edges[E].parent]}")
        #                 temp_edges_from_before.add_row(
        #                     pstate.tables.edges[E].left,
        #                     pstate.tables.edges[E].right,
        #                     pstate.tables.edges[E].parent,
        #                     pstate.tables.edges[E].child)
        #                 E += 1
        #             print("done with pre-existing parent edges")
            # else:

while E < len(pstate.tables.edges):
    print(f"5: adding {pstate.tables.edges[E]}")
    temp_edges_from_before.add_row(
        pstate.tables.edges[E].left,
        pstate.tables.edges[E].right,
        pstate.tables.edges[E].parent,
        pstate.tables.edges[E].child)
    E += 1

for e in temp_edges_from_before:
    z = np.where(temp_edges_from_before.parent == e.parent)[0]
    print(f"edges: {e} -> {np.diff(z)}")

sys.exit(0)

print("result", edges_added, new_edges2, new_edges2 + len(pstate.tables.edges), tcopy_num_edges_b4_simplify,
      len(temp_edges) + len(temp_edges_from_before))
# assert new_edges == new_edges2
print(len(pstate.tables.edges), len(temp_edges), new_edges2)

print(np.where(temp_edges.parent == 44)[0])
print(np.where(temp_edges_from_before.parent == 44)[0])

# Test parent order prior to merging the two tables together.
for i in range(2, len(temp_edges)):
    ei = temp_edges[i-1]
    eip1 = temp_edges[i]
    ti = pstate.tables.nodes.time[ei.parent]
    tip1 = pstate.tables.nodes.time[eip1.parent]
    if ti > tip1:
        print(i, ti, tip1, ei, eip1)
        print("unsorted from new edges")
        sys.exit(0)

for i in range(2, len(temp_edges_from_before)):
    ei = temp_edges_from_before[i-1]
    eip1 = temp_edges_from_before[i]
    ti = pstate.tables.nodes.time[ei.parent]
    tip1 = pstate.tables.nodes.time[eip1.parent]
    if ti > tip1:
        print(i, ti, tip1, ei, eip1)
        print("unsorted from b4")

pstate.tables.edges.set_columns(
    temp_edges.left, temp_edges.right, temp_edges.parent, temp_edges.child)
pstate.tables.edges.append_columns(
    temp_edges_from_before.left, temp_edges_from_before.right,
    temp_edges_from_before.parent,
    temp_edges_from_before.child)

# Test contiguity
up = np.unique(pstate.tables.edges.parent)
for u in up:
    w = np.where(pstate.tables.edges.parent == u)[0]
    d = np.diff(w)
    for wi in d:
        if wi > 1:
            print("//")
            for j in w:
                print("Bad", j, pstate.tables.edges[j])
            break

# Test parent order
for i in range(2, len(pstate.tables.edges)):
    ei = pstate.tables.edges[i-1]
    eip1 = pstate.tables.edges[i]
    ti = pstate.tables.nodes.time[ei.parent]
    tip1 = pstate.tables.nodes.time[eip1.parent]
    if ti > tip1:
        print(i, ti, tip1, ei, eip1)
        print("unsorted")
        sys.exit(0)


print(len(pstate.tables.edges), len(temp_edges), new_edges2)

# for i in pstate.tables.edges:
#     print(i, pstate.tables.nodes.time[i.parent])
print(tcopy_num_edges_b4_simplify, len(pstate.tables.edges))
pstate.tables.simplify()
ts2 = pstate.tables.tree_sequence()
print(len(pstate.tables.edges), len(temp_edges), new_edges2)
print(len(ts.tables.edges), len(ts2.tables.edges))
