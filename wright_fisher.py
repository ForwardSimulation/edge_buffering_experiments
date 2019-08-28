import tskit
import numpy as np


class Parent(object):
    def __init__(self, index, n0, n1):
        self.index = index
        self.n0 = n0
        self.n1 = n1


class PopState(object):
    def __init__(self, N):
        self.edge_buffer = [[[], []] for i in range(N)]
        self.parents = [Parent(i, 2*i, 2*i+1) for i in range(N)]
        self.next_parent = N
        self.tables = tskit.TableCollection(1.0)

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
        for d,p in zip(dead,parent_list):
            n0 = popstate.tables.nodes.add_row(time=gen)
            n1 = popstate.tables.nodes.add_row(time=gen)
            popstate.tables.edges.add_row(
                left=0, right=1, parent=p[0], child=n0)
            popstate.tables.edges.add_row(
                left=0, right=1, parent=p[1], child=n1)
            popstate.parents[d] = Parent(popstate.next_parent, n0, n1)
            popstate.next_parent += 1

    return popstate
