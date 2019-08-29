import tskit
import numpy as np


class Parent(object):
    def __init__(self, index, n0, n1):
        self.index = index
        self.n0 = n0
        self.n1 = n1


class PopState(object):
    def __init__(self, N):
        self.parents = [Parent(i, 2*i, 2*i+1) for i in range(N)]
        self.next_parent = N
        self.tables = tskit.TableCollection(1.0)
        self.buffered_edges = [[[], []] for i in range(N)]
        self.generation_offsets = [(0, len(self.buffered_edges))]
        self.current_generation = 0

        # Measure time going forwards.
        # Will reverse later
        for i in range(N):
            self.tables.nodes.add_row(time=0.0)
            self.tables.nodes.add_row(time=0.0)


def wright_fisher(ngens, psurvival, popstate):
    if psurvival >= 1. or psurvival < 0:
        raise ValueError("unhelpful survival probability")

    for gen in range(1, ngens+1):
        # regulation
        dead = []
        parent_list = []
        for i, p in enumerate(popstate.parents):
            if p.index == -1:
                raise RuntimeError("oops, dead!")
            if np.random.uniform() > psurvival:
                p.index = -1
                parents = np.random.choice(len(popstate.parents), 2)
                # "Mendel"
                p0node = popstate.parents[parents[0]].n0
                if np.random.uniform() < 0.5:
                    p0node = popstate.parents[parents[0]].n1
                p1node = popstate.parents[parents[1]].n0
                if np.random.uniform() < 0.5:
                    p1node = popstate.parents[parents[1]].n1
                parent_list.append((p0node, p1node))
                dead.append(i)
        for d, p in zip(dead, parent_list):
            n0 = popstate.tables.nodes.add_row(time=gen)
            n1 = popstate.tables.nodes.add_row(time=gen)
            popstate.tables.edges.add_row(
                left=0, right=1, parent=p[0], child=n0)
            popstate.tables.edges.add_row(
                left=0, right=1, parent=p[1], child=n1)
            popstate.parents[d] = Parent(popstate.next_parent, n0, n1)
            popstate.next_parent += 1

    return popstate


def wright_fisher_eb(ngens, psurvival, popstate):
    if psurvival >= 1. or psurvival < 0:
        raise ValueError("unhelpful survival probability")

    for gen in range(1, ngens+1):
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
                parent_list.append((p0node, p1node, i0, i1,
                                    popstate.parents[parents[0]].index,
                                    popstate.parents[parents[1]].index))
                dead.append(i)
        x = len(popstate.buffered_edges)
        for d, p in zip(dead, parent_list):
            # NOTE: apply "dead" flag here
            # so that we aren't giving invalid
            # indexes in the regulation step above
            popstate.parents[d].index = -1
            n0 = popstate.tables.nodes.add_row(
                time=popstate.current_generation+gen)
            n1 = popstate.tables.nodes.add_row(
                time=popstate.current_generation+gen)
            popstate.buffered_edges[p[4]][p[2]].append((0, 1, p[0], n0))
            popstate.buffered_edges[p[5]][p[3]].append((0, 1, p[1], n1))
            popstate.parents[d] = Parent(popstate.next_parent, n0, n1)
            popstate.next_parent += 1
            popstate.buffered_edges.append([[], []])

        if len(dead) > 0:
            popstate.generation_offsets.append(
                (x, len(popstate.buffered_edges)))
    popstate.current_generation = gen

    return popstate
