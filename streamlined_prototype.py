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
    parent: np.int32 = attr.ib(converter=int)
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


def sort_alive_at_last_simplification(alive: np.ndarray, tables: tskit.TableCollection):
    """
    Sorts nodes "pastwards" in time.
    """
    alive = np.array(
        sorted(alive, key=lambda x: (tables.nodes.time[x], x)), dtype=alive.dtype
    )

    return alive


def pass_on_node(parent: IndexAndNodes):
    if np.random.uniform() < 0.5:
        return parent.node0
    return parent.node1


def get_alive_nodes(parents: typing.List[IndexAndNodes]):
    alive_nodes = []
    for p in parents:
        alive_nodes.extend(p.nodes)
    return np.array(alive_nodes, dtype=np.int32)


def add_most_ancient_edges(
    tables: tskit.TableCollection, stitched_edges: tskit.EdgeTable, time: float
):
    """
    TODO: docstring
    """
    edge_offset = 0
    # FIXME: can be more efficient via a set_columns call
    while (
        edge_offset < len(tables.edges)
        and tables.nodes.time[tables.edges.parent[edge_offset]] > time
    ):
        e = tables.edges[edge_offset]
        stitched_edges.add_row(
            left=e.left, right=e.right, parent=e.parent, child=e.child
        )
        edge_offset += 1
    return edge_offset


def handle_alive_nodes_from_last_time(
    tables: tskit.TableCollection,
    stitched_edges: tskit.EdgeTable,
    alive_at_last_simplification: np.array,
    buffered_edges: typing.List[BufferedEdgeList],
):
    """
    This is the hard part.  When generations overlap,
    samples that were alive when we last simplified
    may have birth times older than the last simplification
    and died since it happened.  Thus, they may leave new
    edges since we last simplified AND be associated
    with edges already in the edge table.

    The additional complication is that we are recording edges
    by birth time (past to present), but our pre-existing edges
    are ordered present to past by birth time.

    Almost all the code in this function is wrong/too complex.

    To make progress:

    1. process alive_at_last_simplification so that we know
       if any of these individuals have pre-existing
       edges in tables.edges
    2. if they do not, then we have a simple two-step stitch
    3. if the do, then we need to get the indexs of where those
       pre-existing edges are.
    4. once we have those indexes, we can proceed to stitch index
       group by index group, loosely following our previous,
       more complex, method in clean_implementation.py
    """
    A = 0
    # Find the first with descendants
    while A < len(alive_at_last_simplification):
        if len(buffered_edges[alive_at_last_simplification[A]].descendants) > 0:
            break
        A += 1
    # Get the range of times spanned by these samples
    most_recent = tables.nodes.time[alive_at_last_simplification].min()
    most_ancient = tables.nodes.time[alive_at_last_simplification].max()
    edge_offset = 0
    foo = None
    if edge_offset < len(tables.edges):
        foo = tables.nodes.time[tables.edges.parent[edge_offset]]
    print(most_recent, most_recent, foo)
    while (
        edge_offset < len(tables.edges)
        and tables.nodes.time[tables.edges.parent[edge_offset]] > most_recent
    ):
        edge_offset += 1
    # if edge_offset == len(tables.edges):
    #     print("early out")
    #     return edge_offset
    print("chk", edge_offset, len(tables.edges))

    while A < len(alive_at_last_simplification):
        a = alive_at_last_simplification[A]
        while edge_offset < len(tables.edges) and (
            tables.edges.parent[edge_offset] == a
            or tables.nodes.time[tables.edges.parent[edge_offset]]
            >= tables.nodes.time[a]
        ):
            e = tables.edges[edge_offset]
            stitched_edges.add_row(
                left=e.left, right=e.right, parent=e.parent, child=e.child
            )
            edge_offset -= 1
        for birth in buffered_edges[a].descendants:
            stitched_edges.add_row(
                left=birth.left, right=birth.right, parent=a, child=birth.child
            )
        A += 1

    assert A == len(alive_at_last_simplification)

    return edge_offset


def finish_initial_liftover(
    tables: tskit.TableCollection, stitched_edges: tskit.EdgeTable, edge_offset: int
):
    while edge_offset < len(tables.edges):
        stitched_edges.add_row(
            left=tables.edges.left[edge_offset],
            right=tables.edges.right[edge_offset],
            parent=tables.edges.parent[edge_offset],
            child=tables.edges.child[edge_offset],
        )
        edge_offset += 1


def add_new_edges(
    tables: tskit.TableCollection,
    stitched_edges: tskit.EdgeTable,
    buffered_edges: typing.List[BufferedEdge],
    time: float,
):
    for b in reversed(buffered_edges):
        if tables.nodes.time[b.parent] > time:
            for d in b.descendants:
                stitched_edges.add_row(
                    left=d.left, right=d.right, parent=b.parent, child=d.child
                )


def stitch_tables(
    tables: tskit.TableCollection,
    buffered_edges: typing.List[BufferedEdgeList],
    alive_at_last_simplification: np.array,
):
    """
    TODO: docstring w/details
    """
    if len(tables.edges) == 0:
        # this is our first simplification
        for b in reversed(buffered_edges):
            for d in b.descendants:
                tables.edges.add_row(
                    left=d.left, right=d.right, parent=b.parent, child=d.child
                )
        return tables

    # Get the time of the most recent node from alive_at_last_simplification
    # FIXME: this is better done by recording the last time of simplification,
    # passing that to here, and adding all elements whose parent times are
    # more recent
    time = -1
    if len(alive_at_last_simplification) > 0:
        time = tables.nodes.time[alive_at_last_simplification].min()
    stitched_edges = tskit.EdgeTable()
    for b in reversed(buffered_edges):
        if tables.nodes.time[b.parent] < time:
            for d in b.descendants:
                stitched_edges.add_row(
                    left=d.left, right=d.right, parent=b.parent, child=d.child
                )
    edge_offset = handle_alive_nodes_from_last_time(
        tables, stitched_edges, alive_at_last_simplification, buffered_edges
    )
    print(edge_offset, len(tables.edges), len(stitched_edges))
    for i in range(edge_offset):
        e = tables.edges[i]
        stitched_edges.add_row(
            left=e.left, right=e.right, parent=e.parent, child=e.child
        )
    print(edge_offset, len(tables.edges), len(stitched_edges))

    # time = -1.0
    # if len(alive_at_last_simplification) > 0:
    #     time = tables.nodes.time[alive_at_last_simplification].max()
    # edge_offset = add_most_ancient_edges(tables, stitched_edges, time)
    # edge_offset = handle_alive_nodes_from_last_time(
    #     tables,
    #     stitched_edges,
    #     alive_at_last_simplification,
    #     buffered_edges,
    #     edge_offset,
    # )
    # finish_initial_liftover(tables, stitched_edges, edge_offset)
    # add_new_edges(tables, stitched_edges, buffered_edges, time)

    tables.edges.set_columns(
        left=stitched_edges.left,
        right=stitched_edges.right,
        parent=stitched_edges.parent,
        child=stitched_edges.child,
    )

    return tables


def flush_edges_and_simplify(
    sample_nodes: np.ndarray,
    parents: typing.List[IndexAndNodes],
    alive_at_last_simplification: np.ndarray,
    buffered_edges: typing.List[BufferedEdgeList],
    tables: tskit.TableCollection,
):
    tables = stitch_tables(tables, buffered_edges, alive_at_last_simplification)
    idmap = tables.simplify(sample_nodes)
    alive_nodes = idmap[sample_nodes]
    alive_nodes = alive_nodes[np.where(alive_nodes != tskit.NULL)]
    alive_nodes = sort_alive_at_last_simplification(alive_nodes, tables)

    # One of Python's worst gotchas ;)
    buffered_edges[:] = [BufferedEdgeList(i) for i in range(len(tables.nodes))]

    for p in parents:
        p.node0 = idmap[p.node0]
        p.node1 = idmap[p.node1]
    temp = idmap[sample_nodes]
    temp = sort_alive_at_last_simplification(temp, tables)
    return temp


def simplify_classic(sample_nodes: np.ndarray, tables: tskit.TableCollection):
    tables.sort()
    tables.simplify(sample_nodes)


def pairwise_distance_branch(ts: tskit.TreeSequence, samples: np.array):
    sample_sets = []
    for i in range(len(samples)):
        for j in range(i + 1, len(samples)):
            sample_sets.append([samples[i], samples[j]])

    div = ts.diversity(sample_sets, mode="branch")
    return div


def wright_fisher(
    N: int, ngens: int, psurvival: float, simplification_period: int = 10
):
    # assert psurvival == 0.0, "Let's not get ahead of ourselves"
    parents = [IndexAndNodes(i, 2 * i, 2 * i + 1) for i in range(N)]
    tables = tskit.TableCollection(1.0)
    tables2 = tskit.TableCollection(1.0)
    for _ in range(2 * N):
        tables.nodes.add_row(time=ngens)
        tables2.nodes.add_row(time=ngens)
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
            new_node_0b = tables2.nodes.add_row(time=gen - 1)
            new_node_1b = tables2.nodes.add_row(time=gen - 1)
            assert new_node_0 == new_node_0b, f"{new_node_0} {new_node_0b}"
            assert new_node_1 == new_node_1b, f"{new_node_1} {new_node_1b}"

            buffered_edges.append(BufferedEdgeList(new_node_0))
            buffered_edges.append(BufferedEdgeList(new_node_1))

            # Register the edges -- this is the standard method
            tables2.edges.add_row(
                left=0.0, right=1.0, parent=d.node0, child=new_node_0b
            )
            tables2.edges.add_row(
                left=0.0, right=1.0, parent=d.node1, child=new_node_1b
            )

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
            sample_nodes = get_alive_nodes(parents)
            alive_at_last_simplification = flush_edges_and_simplify(
                sample_nodes,
                parents,
                alive_at_last_simplification,
                buffered_edges,
                tables,
            )
            simplify_classic(sample_nodes, tables2)
            assert len(tables.nodes) == len(
                tables2.nodes
            ), f"{len(tables.nodes)}, {len(tables2.nodes)}"
            assert len(tables.edges) == len(tables2.edges)
            assert np.array_equal(tables.nodes.time, tables2.nodes.time)
            simplified = True
            print(f"done with {gen}")
        else:
            simplified = False
    if not simplified:
        sample_nodes = get_alive_nodes(parents)
        flush_edges_and_simplify(
            sample_nodes, parents, alive_at_last_simplification, buffered_edges, tables,
        )
        simplify_classic(sample_nodes, tables2)

    return tables, tables2
