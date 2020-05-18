#include <vector>
#include <tuple>
#include <sstream>
#include <limits>
#include <algorithm>
#include <stdexcept>
#include "edge_buffer.hpp"

static const auto UMAX = std::numeric_limits<std::size_t>::max();

BirthData::BirthData(double l, double r, tsk_id_t c)
    : left{l}, right{r}, child{c}, next{-1}
{
    if (r <= l)
        {
            throw std::invalid_argument("BirthData: right <= left");
        }
}

EdgeBuffer::EdgeBuffer(std::size_t num_nodes) : first(num_nodes, -1), births{}
{
}

ExistingEdges::ExistingEdges(tsk_id_t p, std::size_t start_, std::size_t stop_)
    : parent{p}, start{start_}, stop{stop_}
{
}

std::int32_t
get_buffer_end(const edge_buffer_ptr& new_edges, std::size_t i)
{
    if (i >= new_edges->first.size())
        {
            throw std::runtime_error("invalid parent index");
        }
    auto f = new_edges->first[i];
    while (f != -1 && new_edges->births[f].next != -1)
        {
            f = new_edges->births[f].next;
            if (f != -1 && f >= new_edges->births.size())
                {
                    throw std::runtime_error("invalid next value");
                }
        }
    return f;
}

std::int32_t
buffer_new_edge(tsk_id_t parent, double left, double right, tsk_id_t child,
                edge_buffer_ptr& new_edges)
{
    if (parent == TSK_NULL || child == TSK_NULL)
        {
            throw std::runtime_error("bad node IDs passed to buffer_new_edge");
        }
    if (parent >= new_edges->first.size())
        {
            new_edges->first.resize(parent + 1, -1);
        }
    new_edges->births.emplace_back(left, right, child);
    if (new_edges->first[parent] == -1)
        {
            new_edges->first[parent] = new_edges->births.size() - 1;
            if (new_edges->births[new_edges->first[parent]].next != -1)
                {
                    std::ostringstream o;
                    o << "invalid next entry for first birth: " << parent << ' '
                      << new_edges->first[parent] << ' '
                      << new_edges->births[new_edges->first[parent]].next << ' '
                      << new_edges->first.size() << ' ' << new_edges->births.size();
                    throw std::runtime_error(o.str());
                }
        }
    else
        {
            auto l = get_buffer_end(new_edges, parent);
            new_edges->births[l].next = new_edges->births.size() - 1;
        }
    return new_edges->births.size() - 1;
}

std::int32_t
buffer_new_edge_at(std::int32_t loc, double left, double right, tsk_id_t child,
               edge_buffer_ptr& new_edges)
{
    if(loc>= new_edges->births.size()) { throw std::runtime_error("bad location"); }
    new_edges->births.emplace_back(left, right, child);
    new_edges->births[loc].next = new_edges->births.size() - 1;
    return new_edges->births.size() - 1;
}

std::vector<ExistingEdges>
find_pre_existing_edges(const table_collection_ptr& tables,
                        const std::vector<tsk_id_t>& alive_at_last_simplification,
                        const edge_buffer_ptr& new_edges)
// FIXME: the indexing step need go no farther than the time of the most
// recent node in alive_at_last_simplification.
{
    std::vector<tsk_id_t> alive_with_new_edges;
    for (auto a : alive_at_last_simplification)
        {
            if (new_edges->first[a] != -1)
                {
                    alive_with_new_edges.push_back(a);
                }
        }
    if (alive_with_new_edges.empty()) // get out early
        {
            return {};
        }
    // index where each node already has edges.
    std::vector<std::size_t> starts(tables->nodes.num_rows, UMAX),
        stops(tables->nodes.num_rows, UMAX);
    for (decltype(tables->edges.num_rows) i = 0; i < tables->edges.num_rows; ++i)
        {
            if (starts[tables->edges.parent[i]] == UMAX)
                {
                    starts[tables->edges.parent[i]] = i;
                    stops[tables->edges.parent[i]]
                        = i; // FIXME: idiomatically, this should be i+1
                }
            else
                {
                    stops[tables->edges.parent[i]]
                        = i; // FIXME: idiomatically, this should be i+1
                }
        }

    std::vector<ExistingEdges> existing_edges;
    for (auto a : alive_with_new_edges)
        {
            existing_edges.emplace_back(a, starts[a], stops[a]);
        }

    // Our only sort!!
    std::sort(begin(existing_edges), end(existing_edges),
              [&tables](const ExistingEdges& lhs, const ExistingEdges& rhs) {
                  // lexical comparison of tuple elements just like in Python
                  return std::tie(tables->nodes.time[lhs.parent], lhs.start, lhs.parent)
                         < std::tie(tables->nodes.time[rhs.parent], rhs.start,
                                    rhs.parent);
              });
    for (std::size_t i = 1; i < existing_edges.size(); ++i)
        {
            auto t0 = tables->nodes.time[existing_edges[i - 1].parent];
            auto t1 = tables->nodes.time[existing_edges[i].parent];
            if (t0 > t1)
                {
                    throw std::runtime_error(
                        "existing edges not properly sorted by time");
                }
        }

    return existing_edges;
}

auto
handle_pre_existing_edges(const table_collection_ptr& tables,
                          const edge_buffer_ptr& new_edges,
                          const std::vector<ExistingEdges>& existing_edges,
                          temp_edges& edge_liftover) -> decltype(tables->edges.num_rows)
{
    decltype(tables->edges.num_rows) offset = 0;
    for (const auto& ex : existing_edges)
        {
            // FIXME: this while loop is repeated 2x just w/different
            // ranges
            while (offset < tables->edges.num_rows
                   && tables->nodes.time[tables->edges.parent[offset]]
                          < tables->nodes.time[ex.parent])
                {
                    edge_liftover.add_edge(
                        tables->edges.left[offset], tables->edges.right[offset],
                        tables->edges.parent[offset], tables->edges.child[offset]);
                    ++offset;
                }
            if (ex.start != UMAX)
                {
                    while (offset < ex.start
                           && tables->nodes.time[tables->edges.parent[offset]]
                                  <= tables->nodes.time[ex.parent])
                        {
                            edge_liftover.add_edge(tables->edges.left[offset],
                                                   tables->edges.right[offset],
                                                   tables->edges.parent[offset],
                                                   tables->edges.child[offset]);
                            ++offset;
                        }
                    // FIXME: stop condition isn't idiomatic
                    for (decltype(ex.start) i = ex.start; i < ex.stop + 1; ++i)
                        {
                            edge_liftover.add_edge(
                                tables->edges.left[i], tables->edges.right[i],
                                tables->edges.parent[i], tables->edges.child[i]);
                        }
                    // NOTE: differs from python, so could be a source of error
                    offset = ex.stop + 1;
                }
            auto n = new_edges->first[ex.parent];
            while (n != -1)
                {
                    edge_liftover.add_edge(new_edges->births[n].left,
                                           new_edges->births[n].right, ex.parent,
                                           new_edges->births[n].child);
                    n = new_edges->births[n].next;
                }
        }
    return offset;
}

void
copy_births_since_last_simplification(const edge_buffer_ptr& new_edges,
                                      const table_collection_ptr& tables,
                                      double max_time, temp_edges& edge_liftover)
{

    // TODO: should validate data in new_edges
    edge_liftover.clear(); // Re-use this buffer b/c it'll get big.

    // Go backwards through new births, and add them
    // to our temporary edge table if they are newer
    // than the last simplification time

    for (auto b = new_edges->first.rbegin(); b < new_edges->first.rend(); ++b)
        {
            auto d = std::distance(new_edges->first.rbegin(), b);
            auto parent = new_edges->first.size() - d - 1;
            auto ptime = tables->nodes.time[parent];
            if (*b != -1 && ptime < max_time)
                {
                    auto n = *b;
                    while (n != -1)
                        {
                            edge_liftover.add_edge(new_edges->births[n].left,
                                                   new_edges->births[n].right, parent,
                                                   new_edges->births[n].child);
                            n = new_edges->births[n].next;
                        }
                }
            else if (*b != -1 and ptime >= max_time)
                {
                    break;
                }
        }
}

void
stitch_together_edges(const std::vector<tsk_id_t>& alive_at_last_simplification,
                      double max_time, edge_buffer_ptr& new_edges,
                      temp_edges& edge_liftover, table_collection_ptr& tables)
{
    copy_births_since_last_simplification(new_edges, tables, max_time, edge_liftover);
    auto existing_edges
        = find_pre_existing_edges(tables, alive_at_last_simplification, new_edges);
    auto offset
        = handle_pre_existing_edges(tables, new_edges, existing_edges, edge_liftover);
    for (; offset < tables->edges.num_rows; ++offset)
        {
            edge_liftover.add_edge(
                tables->edges.left[offset], tables->edges.right[offset],
                tables->edges.parent[offset], tables->edges.child[offset]);
        }
    int ret = tsk_edge_table_set_columns(
        &tables->edges, edge_liftover.size(), edge_liftover.left.data(),
        edge_liftover.right.data(), edge_liftover.parent.data(),
        edge_liftover.child.data(), nullptr, 0);
    // This resets sizes to 0, but keeps the memory allocated.
    edge_liftover.clear();
    // TODO: move this cleanup to function
    new_edges->first.resize(tables->nodes.num_rows);
    std::fill(begin(new_edges->first), end(new_edges->first), -1);
    new_edges->births.clear();
}
