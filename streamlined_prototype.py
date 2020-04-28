import argparse
import sys
import typing

import attr
import numpy as np
import tskit


@attr.s()
class BufferedEdge(object):
    left: float = attr.ib(converter=float)
    right: float = attr.ib(converter=float)
    child: np.int32 = attr.ib(converter=int)


@attr.s(kw_only=False)
class BufferedEdgeList(object):
    parent: np.int32 = attr.ib(converter=float)
    descendants: typing.List[BufferedEdge] = attr.ib(converter=list, factory=list)


@attr.s()
class IndexAndNodes(object):
    """
    Double-purpose class:

    1. Represent an alive individual
    2. Temporarily represent the parental nodes
       that'll replace a dead individual at ``index``
    """

    index: int = attr.ib()
    node0: int = attr.ib()
    node1: int = attr.ib()

    @property
    def nodes(self):
        return [self.node0, self.node1]


def sort_alive_at_last_simplification(alive, tables):
    alive = sorted(alive, key=lambda x: (tables.nodes.time[x], x))

    return alive


def pass_on_node(parent: IndexAndNodes):
    if np.random.uniform() < 0.5:
        return parent.node0
    return parent.node1


def get_alive_nodes(parents: typing.List[IndexAndNodes]):
    alive_nodes = []
    for p in parents:
        alive_nodes.extend(p.nodes)
    return alive_nodes


def stitch_tables(
    tables: tskit.TableCollection,
    buffered_edges: typing.List[BufferedEdgeList],
    alive_at_last_simplification: np.array,
):
    return tables


def wright_fisher(
    N: int, ngens: int, psurvival: float, simplification_period: int = 10
):
    parents = [IndexAndNodes(i, 2 * i, 2 * i + 1) for i in range(N)]
    tables = tskit.TableCollection(1.0)
    for _ in range(2 * N):
        tables.nodes.add_row(time=ngens)
    alive_at_last_simplification = np.array([], dtype=np.int32)

    buffered_edges = []
    for p in parents:
        buffered_edges.append(BufferedEdgeList(p.node0))
        buffered_edges.append(BufferedEdgeList(p.node1))

    simplified = False

    for gen in range(ngens, 0, -1):
        dead = []
        for i, _ in enumerate(parents):
            # Kill off an individual, perform a birth
            if np.random.uniform() > psurvival:
                offspring_parents = np.random.choice(len(parents), 2)
                p0node = pass_on_node(parents[offspring_parents[0]])
                p1node = pass_on_node(parents[offspring_parents[1]])
                # parent_list[i] will be replaced
                # by an individual inheriting these two nodes
                dead.append(IndexAndNodes(i, p0node, p1node))
        for d in dead:
            # Register two new nodes
            new_node_0 = tables.nodes.add_row(time=gen - 1)
            new_node_1 = tables.nodes.add_row(time=gen - 1)
            buffered_edges.append(BufferedEdgeList(new_node_0))
            buffered_edges.append(BufferedEdgeList(new_node_1))

            # Register the edges -- this is the standard method
            # tables.edges.add_row(left=0.0, right=1.0, parent=d.node0, child=new_node_0)
            # tables.edges.add_row(left=0.0, right=1.0, parent=d.node1, child=new_node_1)

            # Buffer the new edges
            buffered_edges[d.node0].descendants.append(
                BufferedEdge(left=0.0, right=1.0, child=new_node_0)
            )
            buffered_edges[d.node1].descendants.append(
                BufferedEdge(left=0.0, right=1.0, child=new_node_1)
            )

            # Make the new parent
            parents[d.index] = IndexAndNodes(d.index, new_node_0, new_node_1)

        # Simplify, if it is time to
        if gen < ngens and gen % simplification_period == 0.0:
            alive_nodes = get_alive_nodes(parents)
            tables = stitch_tables(tables, buffered_edges, alive_at_last_simplification)

    return tables, get_alive_nodes(parents)


tables, alive_nodes = wright_fisher(10, 100, 0.5)
for i in tables.edges:
    print(i, tables.nodes.time[i.parent], tables.nodes.time[i.child])
tables.sort()
tables.simplify(alive_nodes)
ts = tables.tree_sequence()
ts.first().draw(format="svg", path="foo.svg")
