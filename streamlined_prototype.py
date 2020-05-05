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


@attr.s()
class ExistingEdges(object):
    """
    Where are the edges for this parent
    in an edge table?

    Used during stitching
    """

    parent: int = attr.ib()
    start: int = attr.ib()
    stop: int = attr.ib()


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

    alive_with_new_edges = [
        a
        for a in alive_at_last_simplification
        if len(buffered_edges[a].descendants) > 0
    ]

    # Do any in alive_with_new_edges have
    # edges already in tables.edges?
    starts = np.array([np.iinfo(np.int32).max] * len(tables.nodes), dtype=np.int32)
    stops = np.array([np.iinfo(np.int32).max] * len(tables.nodes), dtype=np.int32)
    for i, e in enumerate(tables.edges):
        if starts[e.parent] == np.iinfo(np.int32).max:
            starts[e.parent] = i
            stops[e.parent] = i
        else:
            stops[e.parent] = i

    # FIXME: a simple iteration would do better here.
    existing_edges = [
        ExistingEdges(i, starts[i], stops[i])
        for i in alive_with_new_edges
        if starts[i] != np.iinfo(np.int32).max
    ]

    num_new_births_from_old_parents = 0
    old_edges_added = 0
    if len(existing_edges) == 0:
        # Easy, so let's do less logic and get out of here
        alive_with_new_edges = sorted(
            alive_with_new_edges, key=lambda x: tables.nodes.time[x]
        )
        for a in alive_with_new_edges:
            for d in buffered_edges[a].descendants:
                num_new_births_from_old_parents += 1
                stitched_edges.add_row(
                    left=d.left, right=d.right, parent=a, child=d.child
                )
        stitched_edges.append_columns(
            tables.edges.left,
            tables.edges.right,
            tables.edges.parent,
            tables.edges.child,
        )
        return num_new_births_from_old_parents, len(tables.edges.left)

    existing_edges = [
        ExistingEdges(i, starts[i], stops[i])
        for i in alive_with_new_edges
        if len(buffered_edges[i].descendants) > 0
    ]
    # existing_edges = sorted(existing_edges, key=lambda x: x.start)
    existing_edges = sorted(
        existing_edges,
        key=lambda x: (tables.nodes.time[x.parent], x.start, x.parent)
    )
    offset = 0
    print("chk", len(stitched_edges), offset)
    for i, ex in enumerate(existing_edges):
        # Add the pre-existing edges
        # print(ex, tables.nodes.time[ex.parent])
        if ex.start != np.iinfo(np.int32).max:
            while offset < ex.start:
                e = tables.edges[offset]
                old_edges_added += 1
                assert e.parent != ex.parent
                assert tables.nodes.time[e.parent] <= tables.nodes.time[ex.parent]
                stitched_edges.add_row(
                    left=e.left, right=e.right, parent=e.parent, child=e.child
                )
                offset += 1
            for j in range(ex.start, ex.stop + 1):
                old_edges_added += 1
                e = tables.edges[j]
                assert e.parent == ex.parent
                # print(f"1: adding {j} {tables.nodes.time[e.parent]}")
                offset = j + 1
                # offset += 1
                stitched_edges.add_row(
                    left=e.left, right=e.right, parent=e.parent, child=e.child
                )
        # Any new edges are more recent than anything in the edge
        # table for this parent, and thus have larger child id values,
        # so they go in next
        for d in buffered_edges[ex.parent].descendants:
            # print(f"3: adding {ex.parent} {d.child} {tables.nodes.time[ex.parent]}")
            num_new_births_from_old_parents += 1
            stitched_edges.add_row(
                left=d.left, right=d.right, parent=ex.parent, child=d.child
            )

    print("check=", offset, len(tables.edges))
    final = offset
    for i in range(final, len(tables.edges)):
        old_edges_added += 1
        stitched_edges.add_row(
            left=tables.edges.left[i],
            right=tables.edges.right[i],
            parent=tables.edges.parent[i],
            child=tables.edges.child[i],
        )
    return num_new_births_from_old_parents, old_edges_added


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
    input_edge_table_length = len(tables.edges)
    time = -1
    if len(alive_at_last_simplification) > 0:
        time = tables.nodes.time[alive_at_last_simplification].min()
    stitched_edges = tskit.EdgeTable()
    total_births = sum([len(i.descendants) for i in buffered_edges])
    num_new_births = 0
    for b in reversed(buffered_edges):
        if tables.nodes.time[b.parent] < time:
            for d in b.descendants:
                stitched_edges.add_row(
                    left=d.left, right=d.right, parent=b.parent, child=d.child
                )
                num_new_births += 1

    (
        num_new_births_from_old_parents,
        old_edges_added,
    ) = handle_alive_nodes_from_last_time(
        tables, stitched_edges, alive_at_last_simplification, buffered_edges
    )
    print(
        f"Total = {num_new_births} + {num_new_births_from_old_parents} = {total_births}\n"
        f"\t{old_edges_added} = {input_edge_table_length}"
    )

    tables.edges.set_columns(
        left=stitched_edges.left,
        right=stitched_edges.right,
        parent=stitched_edges.parent,
        child=stitched_edges.child,
    )

    # Do some validation
    last_child = np.array([-1] * len(tables.nodes), dtype=np.int32)
    for i, e in enumerate(tables.edges):
        if last_child[e.child] != -1:
            assert last_child[e.child] == i - 1
            last_child[e.child] = i

    E = 0
    while E < len(tables.edges):
        p = tables.edges.parent[E]
        children = []
        while E < len(tables.edges) and tables.edges.parent[E] == p:
            children.append(tables.edges.child[E])
            E += 1
        assert children == sorted(
            children
        ), f"{children} {p in alive_at_last_simplification} {buffered_edges[p]}"

    return tables


def flush_edges_and_simplify(
    sample_nodes: np.ndarray,
    parents: typing.List[IndexAndNodes],
    alive_at_last_simplification: np.ndarray,
    buffered_edges: typing.List[BufferedEdgeList],
    tables: tskit.TableCollection,
):
    tables = stitch_tables(tables, buffered_edges, alive_at_last_simplification)
    print(f"stitched: {len(tables.nodes)} {len(tables.edges)}")
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
    print(f"classic: {len(tables.nodes)} {len(tables.edges)}")
    tables.simplify(sample_nodes)


def pairwise_distance_branch(ts: tskit.TreeSequence, samples: np.array):
    sample_sets = []
    indexes = []
    for i in range(len(samples)):
        sample_sets.append([i])
        for j in range(i + 1, len(samples)):
            indexes.append((i, j))

    div = ts.divergence(sample_sets, indexes=indexes, mode="branch")
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
