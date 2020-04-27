import argparse
import sys

import numpy as np
import tskit


class Parent(object):
    def __init__(self, index, n0, n1):
        self.index = index
        self.n0 = n0
        self.n1 = n1


class PopState(object):
    def __init__(self, N):
        self.parents = [Parent(i, 2 * i, 2 * i + 1) for i in range(N)]
        self.next_parent = N
        self.tables = tskit.TableCollection(1.0)
        self.buffered_edges = [[[], []] for i in range(N)]
        self.pnodes = [(2 * i, 2 * i + 1) for i in range(N)]
        self.generation_offsets = [(0, len(self.buffered_edges))]
        self.current_generation = 0

        # Measure time going forwards.
        # Will reverse later
        for i in range(N):
            self.tables.nodes.add_row(time=0.0)
            self.tables.nodes.add_row(time=0.0)


def wright_fisher(ngens, psurvival, popstate):
    if psurvival >= 1.0 or psurvival < 0:
        raise ValueError("unhelpful survival probability")

    for gen in range(1, ngens + 1):
        # regulation
        dead = []
        parent_list = []
        for i, p in enumerate(popstate.parents):
            if p.index == -1:
                raise RuntimeError("oops, dead!")
            if np.random.uniform() > psurvival:
                parents = np.random.choice(len(popstate.parents), 2)
                # "Mendel"
                p0node = popstate.parents[parents[0]].n0
                i0 = 0
                if np.random.uniform() < 0.5:
                    p0node = popstate.parents[parents[0]].n1
                    i0 = 1
                p1node = popstate.parents[parents[1]].n0
                i1 = 0
                if np.random.uniform() < 0.5:
                    p1node = popstate.parents[parents[1]].n1
                    i1 = 1
                parent_list.append(
                    (
                        p0node,
                        p1node,
                        i0,
                        i1,
                        popstate.parents[parents[0]].index,
                        popstate.parents[parents[1]].index,
                    )
                )
                dead.append(i)
        x = len(popstate.buffered_edges)
        for d, p in zip(dead, parent_list):
            # NOTE: apply "dead" flag here
            # so that we aren't giving invalid
            # indexes in the regulation step above
            popstate.parents[d].index = -1
            n0 = popstate.tables.nodes.add_row(time=popstate.current_generation + gen)
            n1 = popstate.tables.nodes.add_row(time=popstate.current_generation + gen)
            popstate.buffered_edges[p[4]][p[2]].append((0, 1, p[0], n0))
            popstate.buffered_edges[p[5]][p[3]].append((0, 1, p[1], n1))
            popstate.pnodes.append((n0, n1))
            popstate.parents[d] = Parent(popstate.next_parent, n0, n1)
            popstate.next_parent += 1
            popstate.buffered_edges.append([[], []])

        if len(dead) > 0:
            popstate.generation_offsets.append((x, len(popstate.buffered_edges)))
    popstate.current_generation += gen

    return popstate


def remap_alive_parent_nodes(pstate, idmap):
    for p in pstate.parents:
        p.n0 = idmap[p.n0]
        p.n1 = idmap[p.n1]

        for i in [p.n0, p.n1]:
            assert i != tskit.NULL, "alive parent node mapped to NULL"


def index_alive_parent_nodes(alive_parents, nnodes, edges):
    isparent = [False for i in range(nnodes)]
    ischild = [False for i in range(nnodes)]
    isalive = [False for i in range(nnodes)]
    minparent = [tskit.NULL for i in range(nnodes)]
    maxparent = [tskit.NULL for i in range(nnodes)]

    for p in alive_parents:
        for n in [p.n0, p.n1]:
            isalive[n] = True

    for i, e in enumerate(edges):
        if isalive[e.parent] is True:
            isparent[e.parent] = True
            maxparent[e.parent] = i
            if minparent[e.parent] == tskit.NULL:
                minparent[e.parent] = i
        if isalive[e.child] is True:
            ischild[e.child] = True

    return isparent, ischild, minparent, maxparent


def sort_and_index_alive_parents(alive_parents, nodes, minparent):
    times = nodes.time[:]

    def f(l1, l2):
        if l1 == tskit.NULL and l2 == tskit.NULL:
            return np.iinfo(np.uint32).max
        if l1 == tskit.NULL:
            return l2
        elif l2 == tskit.NULL:
            return l1
        return min(l1, l2)

    alive_parents[:] = sorted(
        alive_parents, key=lambda x: (times[x.n0], f(minparent[x.n0], minparent[x.n1]))
    )
    for i, p in enumerate(alive_parents):
        p.index = i


def brute_force_merge_and_simplify(pstate):
    tc = tskit.TableCollection(pstate.tables.sequence_length)

    flags = np.zeros(len(pstate.tables.nodes), dtype=np.uint32)
    for p in pstate.parents:
        flags[p.n0] = 1
        flags[p.n1] = 1
    tc.nodes.set_columns(
        flags=flags,
        time=-1.0 * (pstate.tables.nodes.time - pstate.tables.nodes.time.max()),
    )

    tc.edges.set_columns(
        pstate.tables.edges.left,
        pstate.tables.edges.right,
        pstate.tables.edges.parent,
        pstate.tables.edges.child,
    )
    for eb in pstate.buffered_edges:
        for i in eb[0] + eb[1]:
            tc.edges.add_row(*i)
    tc.sort()
    tc.simplify()
    return tc.tree_sequence()


def reedgeucation(pstate, ischild, minparent, maxparent):
    """
    Great function name, or best function name ever?
    """

    E = 0
    edges_new_births = tskit.EdgeTable()
    edges_previous_births = tskit.EdgeTable()
    for o in reversed(pstate.generation_offsets):
        for i in range(*o):
            # Get parent node IDs
            pnodes = pstate.pnodes[i]
            if i < len(pstate.parents):
                if minparent[pnodes[0]] == tskit.NULL:
                    isparent0 = False
                else:
                    isparent0 = True
                if minparent[pnodes[1]] == tskit.NULL:
                    isparent1 = False
                else:
                    isparent1 = True
                # mn0 = minparent[pnodes[0]]
                mx0 = maxparent[pnodes[0]]
                mn1 = minparent[pnodes[1]]
                mx1 = maxparent[pnodes[1]]
                if isparent0 is True and isparent1 is True:
                    assert mx0 != tskit.NULL
                    assert mx1 != tskit.NULL
                    edges_previous_births.append_columns(
                        pstate.tables.edges.left[E : mx0 + 1],
                        pstate.tables.edges.right[E : mx0 + 1],
                        pstate.tables.edges.parent[E : mx0 + 1],
                        pstate.tables.edges.child[E : mx0 + 1],
                    )
                    E = mx0 + 1
                    for k in pstate.buffered_edges[i][0]:
                        assert k[2] == pnodes[0]
                        edges_previous_births.add_row(*k)
                    edges_previous_births.append_columns(
                        pstate.tables.edges.left[E : mx1 + 1],
                        pstate.tables.edges.right[E : mx1 + 1],
                        pstate.tables.edges.parent[E : mx1 + 1],
                        pstate.tables.edges.child[E : mx1 + 1],
                    )
                    E = mx1 + 1
                    for k in pstate.buffered_edges[i][1]:
                        assert k[2] == pnodes[1]
                        edges_previous_births.add_row(*k)
                elif isparent0 is True:
                    assert mx0 != tskit.NULL
                    assert isparent1 is False
                    edges_previous_births.append_columns(
                        pstate.tables.edges.left[E : mx0 + 1],
                        pstate.tables.edges.right[E : mx0 + 1],
                        pstate.tables.edges.parent[E : mx0 + 1],
                        pstate.tables.edges.child[E : mx0 + 1],
                    )
                    E = mx0 + 1
                    for k in pstate.buffered_edges[i][0]:
                        edges_previous_births.add_row(*k)
                    for k in pstate.buffered_edges[i][1]:
                        edges_previous_births.add_row(*k)
                elif isparent1 is True:
                    assert mn1 != tskit.NULL
                    assert mx1 != tskit.NULL
                    assert isparent0 is False
                    edges_previous_births.append_columns(
                        pstate.tables.edges.left[E:mn1],
                        pstate.tables.edges.right[E:mn1],
                        pstate.tables.edges.parent[E:mn1],
                        pstate.tables.edges.child[E:mn1],
                    )
                    for k in pstate.buffered_edges[i][0]:
                        edges_previous_births.add_row(*k)
                    edges_previous_births.append_columns(
                        pstate.tables.edges.left[mn1 : mx1 + 1],
                        pstate.tables.edges.right[mn1 : mx1 + 1],
                        pstate.tables.edges.parent[mn1 : mx1 + 1],
                        pstate.tables.edges.child[mn1 : mx1 + 1],
                    )
                    for k in pstate.buffered_edges[i][1]:
                        edges_previous_births.add_row(*k)
                    E = mx1 + 1
                else:
                    ptime = pstate.tables.nodes.time[pnodes[0]]
                    if ischild[pnodes[0]] or ischild[pnodes[1]]:
                        while (
                            E < len(pstate.tables.edges)
                            and pstate.tables.nodes.time[pstate.tables.edges.parent[E]]
                            < ptime
                        ):
                            e = pstate.tables.edges[E]
                            edges_previous_births.add_row(
                                e.left, e.right, e.parent, e.child
                            )
                            E += 1

                    for n in [0, 1]:
                        for k in pstate.buffered_edges[i][n]:
                            assert k[2] == pnodes[n], f"{k} {pnodes}"
                            edges_previous_births.add_row(*k)
            else:
                for n in [0, 1]:
                    for k in pstate.buffered_edges[i][n]:
                        assert k[2] == pnodes[n]
                        edges_new_births.add_row(*k)

    while E < len(pstate.tables.edges):
        edges_previous_births.add_row(
            pstate.tables.edges[E].left,
            pstate.tables.edges[E].right,
            pstate.tables.edges[E].parent,
            pstate.tables.edges[E].child,
        )
        E += 1
    pstate.tables.edges.set_columns(
        edges_new_births.left,
        edges_new_births.right,
        edges_new_births.parent,
        edges_new_births.child,
    )
    pstate.tables.edges.append_columns(
        edges_previous_births.left,
        edges_previous_births.right,
        edges_previous_births.parent,
        edges_previous_births.child,
    )


def make_parser():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    optional = parser.add_argument_group("Simulation parameters")
    optional.add_argument(
        "--popsize", "-N", default=100, type=int, help="Diploid population size"
    )
    optional.add_argument("--seed", type=int, default=42, help="RNG seed")
    optional.add_argument(
        "--psurvival",
        "-p",
        type=float,
        default=0.9,
        help="Survival probability per generation",
    )
    optional.add_argument(
        "--burnin", "-b", type=int, default=20, help="Burnin time. Multiple of N"
    )
    return parser


if __name__ == "__main__":
    parser = make_parser()
    args = parser.parse_args(sys.argv[1:])
    np.random.seed(args.seed)

    pstate = PopState(args.popsize)

    pstate = wright_fisher(args.burnin * args.popsize, args.psurvival, pstate)

    # After one bout of simulation, things are easy,
    # an we can simply populate the edge table
    for o in reversed(pstate.generation_offsets):
        for i in range(*o):
            for j in pstate.buffered_edges[i][0]:
                pstate.tables.edges.add_row(*j)
            for j in pstate.buffered_edges[i][1]:
                pstate.tables.edges.add_row(*j)
            # reset the buffer
            for j in pstate.buffered_edges[i]:
                j.clear()

    # Simplify with respect to currently alive individuals
    # 1. Set the sample flag for them
    flags = np.zeros(len(pstate.tables.nodes), dtype=np.uint32)
    for p in pstate.parents:
        for i in [p.n0, p.n1]:
            flags[i] = tskit.NODE_IS_SAMPLE
    # 2. Set the flags + reverse time
    pstate.tables.nodes.set_columns(
        flags=flags,
        time=-1.0 * (pstate.tables.nodes.time - pstate.tables.nodes.time.max()),
    )
    # 3. Simplify
    idmap = pstate.tables.simplify()

    # Post-simplification cleanup required if we are going
    # to simulate more and then simplify again
    # 1. Remap alive parent nodes.  This is normal/expected procedure.
    remap_alive_parent_nodes(pstate, idmap)
    # 2. Figure out if & where each parental node is in the
    #    simplified edge table
    #    Complexity: O(no. edges)
    isparent, ischild, minparent, maxparent = index_alive_parent_nodes(
        pstate.parents, len(pstate.tables.nodes), pstate.tables.edges
    )
    # 3. We need to sort the alive parents, and re-index them
    sort_and_index_alive_parents(pstate.parents, pstate.tables.nodes, minparent)
    for p in pstate.parents:
        w = minparent[p.n0]
        w2 = minparent[p.n1]

    with open("after_sorting_parents_cleaner.txt", "w") as f:
        for p in pstate.parents:
            ptime = pstate.tables.nodes.time[p.n0]
            w = minparent[p.n0]
            w2 = minparent[p.n1]
            f.write(
                f"{p.n0} {p.n1}-> ({ptime} {w} {w2} {isparent[p.n0]} {isparent[p.n1]})\n"
            )

    # sys.exit(0)

    # 4. Reset the buffer and the index
    pstate.buffered_edges = [[[], []] for i in range(len(pstate.parents))]
    # 5. Update the index value of the next birth
    pstate.next_parent = len(pstate.buffered_edges)
    # 6. Reset offsets
    pstate.generation_offsets = [(0, len(pstate.buffered_edges))]
    # 7. Reset master parent node list
    pstate.pnodes = [(i.n0, i.n1) for i in pstate.parents]
    # 8. Reset flags and change time from backwards to forwards
    flags = np.zeros(len(pstate.tables.nodes), dtype=np.uint32)
    newtime = np.array(
        [pstate.current_generation - i for i in pstate.tables.nodes.time]
    )
    pstate.tables.nodes.set_columns(flags=flags, time=newtime)

    # Evolve again for a short bit of time
    pstate = wright_fisher(20, args.psurvival, pstate)

    # For testing, we merge the data "the old way"
    # and return a tree sequence
    ts_classic_method = brute_force_merge_and_simplify(pstate)

    flags = np.zeros(len(pstate.tables.nodes), dtype=np.uint32)
    for p in pstate.parents:
        flags[p.n0] = 1
        flags[p.n1] = 1
    pstate.tables.nodes.set_columns(
        flags=flags,
        time=-1.0 * (pstate.tables.nodes.time - pstate.tables.nodes.time.max()),
    )

    reedgeucation(pstate, ischild, minparent, maxparent)

    idmap = pstate.tables.simplify()
    ts = pstate.tables.tree_sequence()

    # print trees to file
    next(ts_classic_method.trees()).draw(
        path="brute_force.svg", format="svg", height=1000, width=1000
    )
    next(ts.trees()).draw(path="buffered.svg", format="svg", height=1000, width=1000)

    # for e in ts.tables.edges:
    #     print(e)
