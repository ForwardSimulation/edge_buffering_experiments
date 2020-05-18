#include <iostream>
#include <cstring>
#include <sstream>
#include <algorithm>
#include <tuple>
#include <limits>
#include <memory>
#include <vector>
#include <gsl/gsl_randist.h>
#include "rng.hpp"
#include "tskit_tools.hpp"
#include "edge_buffer.hpp"
#include "sort_tables.hpp"

namespace
{
    struct Parent
    {
        std::size_t index;
        tsk_id_t node0, node1;
        Parent(std::size_t i, tsk_id_t n0, tsk_id_t n1) : index(i), node0(n0), node1(n1)
        {
        }
    };

    struct Birth
    {
        std::size_t index;
        tsk_id_t p0node0, p0node1, p1node0, p1node1;
        Birth(std::size_t i, const Parent& p0, const Parent& p1)
            : index(i), p0node0(p0.node0), p0node1(p0.node1), p1node0(p1.node0),
              p1node1(p1.node1)
        {
        }
    };
}

static void
handle_tskit_return_code(int code)
{
    if (code != 0)
        {
            std::ostringstream o;
            o << tsk_strerror(code);
            throw std::runtime_error(o.str());
        }
}

static tsk_id_t
record_node(double t, table_collection_ptr& tables)
{
    return tsk_node_table_add_row(&tables->nodes,
                                  0,        // flag
                                  t,        // time
                                  TSK_NULL, // population
                                  TSK_NULL, // individual
                                  nullptr,  // metadata
                                  0         // metadata length
    );
}

static void
deaths_and_parents(const GSLrng& rng, const std::vector<Parent>& parents,
                   double psurvival, std::vector<Birth>& births)
{
    births.clear();
    for (std::size_t i = 0; i < parents.size(); ++i)
        {
            if (gsl_rng_uniform(rng.get()) > psurvival)
                {
                    std::size_t parent0 = gsl_ran_flat(rng.get(), 0, parents.size());
                    std::size_t parent1 = gsl_ran_flat(rng.get(), 0, parents.size());
                    births.emplace_back(i, parents[parent0], parents[parent1]);
                }
        }
}

void
recombination_breakpoints(const GSLrng& rng, double littler, double maxlen,
                          std::vector<double>& breakpoints)
{
    breakpoints.clear();
    auto nxovers = gsl_ran_poisson(rng.get(), littler);
    for (decltype(nxovers) i = 0; i < nxovers; ++i)
        {
            breakpoints.push_back(gsl_ran_flat(rng.get(), 0., maxlen));
        }
    std::sort(begin(breakpoints), end(breakpoints));
}

void
recombine_and_record_edges(const GSLrng& rng, double littler,
                           std::vector<double>& breakpoints, tsk_id_t parental_node0,
                           tsk_id_t parental_node1, tsk_id_t child,
                           table_collection_ptr& tables)
// NOTE: this is an improvement on what I do in fwdpp?
{
    recombination_breakpoints(rng, littler, tables->sequence_length, breakpoints);
    double left = 0.;
    std::size_t breakpoint = 1;
    auto pnode0 = parental_node0;
    auto pnode1 = parental_node1;
    for (; breakpoint < breakpoints.size(); ++breakpoint)
        {
            auto rv
                = tsk_edge_table_add_row(&tables->edges, left, breakpoints[breakpoint],
                                         pnode0, child, nullptr, 0);
            std::swap(pnode0, pnode1);
            left = breakpoints[breakpoint];
        }
    auto rv = tsk_edge_table_add_row(&tables->edges, left, tables->sequence_length,
                                     pnode0, child, nullptr, 0);
}

void
recombine_and_buffer_edges(const GSLrng& rng, double littler,
                           std::vector<double>& breakpoints, tsk_id_t parental_node0,
                           tsk_id_t parental_node1, tsk_id_t child, double maxlen,
                           edge_buffer_ptr& new_edges)
{
    recombination_breakpoints(rng, littler, maxlen, breakpoints);
    double left = 0.;
    std::size_t breakpoint = 1;
    auto pnode0 = parental_node0;
    auto pnode1 = parental_node1;
    for (; breakpoint < breakpoints.size(); ++breakpoint)
        {
            buffer_new_edge(pnode0, left, breakpoints[breakpoint], child, new_edges);
            std::swap(pnode0, pnode1);
            left = breakpoints[breakpoint];
        }
    buffer_new_edge(pnode0, left, maxlen, child, new_edges);
}

static void
generate_births(const GSLrng& rng, const std::vector<Birth>& births, double littler,
                std::vector<double>& breakpoints, double birth_time,
                bool buffer_new_edges, edge_buffer_ptr& new_edges,
                std::vector<Parent>& parents, table_collection_ptr& tables)
{
    for (auto& b : births)
        {
            auto new_node_0 = record_node(birth_time, tables);
            auto new_node_1 = record_node(birth_time, tables);
            auto p0n0 = b.p0node0;
            auto p0n1 = b.p0node1;
            if (gsl_rng_uniform(rng.get()) < 0.5)
                {
                    std::swap(p0n0, p0n1);
                }
            auto p1n0 = b.p1node0;
            auto p1n1 = b.p1node1;
            if (gsl_rng_uniform(rng.get()) < 0.5)
                {
                    std::swap(p1n0, p1n1);
                }
            if (buffer_new_edges == false)
                {
                    recombine_and_record_edges(rng, littler, breakpoints, p0n0, p0n1,
                                               new_node_0, tables);
                    recombine_and_record_edges(rng, littler, breakpoints, p1n0, p1n1,
                                               new_node_1, tables);
                }
            else
                {
                    double ptime = tables->nodes.time[p0n0];
                    double ctime = tables->nodes.time[new_node_0];
                    if (ctime >= ptime)
                        {
                            throw std::runtime_error("bad parent/child time");
                        }
                    recombine_and_buffer_edges(rng, littler, breakpoints,
                            p0n0, p0n1, new_node_0, tables->sequence_length,
                            new_edges);
                    ptime = tables->nodes.time[p1n0];
                    ctime = tables->nodes.time[new_node_1];
                    if (ctime >= ptime)
                        {
                            throw std::runtime_error("bad parent/child time");
                        }
                    recombine_and_buffer_edges(rng, littler, breakpoints,
                            p1n0, p1n1, new_node_1, tables->sequence_length,
                            new_edges);
                }
            parents[b.index] = Parent(b.index, new_node_0, new_node_1);
        }
}

// NOTE: seems like samples could/should be const?
static void
sort_n_simplify(bool cppsort, bool parallel_sort, double last_time_simplified,
                std::vector<tsk_id_t>& samples, std::vector<tsk_id_t>& node_map,
                table_collection_ptr& tables)
{
    //tsk_bookmark_t bookmark;
    //std::memset(&bookmark, 0, sizeof(bookmark));
    //tsk_id_t parent_to_sort = TSK_NULL, last_parent = TSK_NULL;
    //std::size_t last_parent_index = std::numeric_limits<std::size_t>::max();
    //for (std::size_t i = 0; i < tables->edges.num_rows && parent_to_sort == TSK_NULL;
    //     ++i)
    //    {
    //        if (tables->edges.parent[i] != last_parent)
    //            {
    //                last_parent = tables->edges.parent[i];
    //                last_parent_index = i;
    //            }
    //        if (tables->nodes.time[tables->edges.child[i]] < last_time_simplified)
    //            {
    //                parent_to_sort = tables->edges.parent[i];
    //            }
    //    }
    //if (parent_to_sort != TSK_NULL)
    //    {
    //        bookmark.edges = last_parent_index;
    //    }
    //int rv = tsk_table_collection_sort(tables.get(), &bookmark, 0);
    int rv = -1;
    if (cppsort == false)
        {
            rv = tsk_table_collection_sort(tables.get(), nullptr, 0);
            handle_tskit_return_code(rv);
        }
    else
        {
            sort_tables(tables.get(), parallel_sort);
        }
    //if (bookmark.edges > 0)
    //    {
    //        std::rotate(tables->edges.left, tables->edges.left + bookmark.edges,
    //                    tables->edges.left + tables->edges.num_rows);
    //        std::rotate(tables->edges.right, tables->edges.right + bookmark.edges,
    //                    tables->edges.right + tables->edges.num_rows);
    //        std::rotate(tables->edges.parent, tables->edges.parent + bookmark.edges,
    //                    tables->edges.parent + tables->edges.num_rows);
    //        std::rotate(tables->edges.child, tables->edges.child + bookmark.edges,
    //                    tables->edges.child + tables->edges.num_rows);
    //    }
    rv = tsk_table_collection_simplify(tables.get(), samples.data(), samples.size(), 0,
                                       node_map.data());
    handle_tskit_return_code(rv);
}

static void
flush_buffer_n_simplify(std::vector<tsk_id_t>& alive_at_last_simplification,
                        std::vector<tsk_id_t>& samples, std::vector<tsk_id_t>& node_map,
                        edge_buffer_ptr& new_edges, temp_edges& edge_liftover,
                        table_collection_ptr& tables)
{
    double max_time = std::numeric_limits<double>::max();
    for (auto a : alive_at_last_simplification)
        {
            max_time = std::min(max_time, tables->nodes.time[a]);
        }

    stitch_together_edges(alive_at_last_simplification, max_time, new_edges,
                          edge_liftover, tables);
    int rv = tsk_table_collection_simplify(tables.get(), samples.data(), samples.size(),
                                           0, node_map.data());
    handle_tskit_return_code(rv);
}

void
simulate(const GSLrng& rng, unsigned N, double psurvival, unsigned nsteps,
         unsigned simplification_interval, double rho, bool buffer_new_edges,
         bool cppsort, bool parallel_sort, table_collection_ptr& tables)
{
    std::vector<Parent> parents;
    for (unsigned i = 0; i < N; ++i)
        {
            auto id0 = record_node(nsteps, tables);
            auto id1 = record_node(nsteps, tables);
            parents.emplace_back(i, id0, id1);
        }

    // The next bits are all for buffering
    std::vector<tsk_id_t> alive_at_last_simplification;
    temp_edges edge_liftover;

    edge_buffer_ptr new_edges(nullptr);
    if (buffer_new_edges)
        {
            new_edges.reset(new EdgeBuffer(tables->nodes.num_rows));
            if (new_edges->first.size() != 2 * N)
                {
                    throw std::runtime_error("bad setup of edge_buffer_ptr");
                }
        }

    std::vector<Birth> births;
    std::vector<tsk_id_t> samples, node_map;
    bool simplified = false;
    double last_time_simplified = nsteps;
    double littler = rho / (2. * static_cast<double>(N));
    std::vector<double> breakpoints;
    for (unsigned step = 1; step <= nsteps; ++step)
        {
            deaths_and_parents(rng, parents, psurvival, births);
            generate_births(rng, births, littler, breakpoints, nsteps - step,
                            buffer_new_edges, new_edges, parents, tables);
            if (step % simplification_interval == 0.)
                {
                    samples.clear();
                    for (auto& p : parents)
                        {
                            samples.push_back(p.node0);
                            samples.push_back(p.node1);
                        }
                    node_map.resize(tables->nodes.num_rows);

                    if (buffer_new_edges == false)
                        {
                            sort_n_simplify(cppsort, parallel_sort, last_time_simplified,
                                            samples, node_map, tables);
                        }
                    else
                        {
                            flush_buffer_n_simplify(alive_at_last_simplification,
                                                    samples, node_map, new_edges,
                                                    edge_liftover, tables);
                        }
                    simplified = true;
                    last_time_simplified = nsteps - step;
                    //remap parent nodes
                    for (auto& p : parents)
                        {
                            p.node0 = node_map[p.node0];
                            p.node1 = node_map[p.node1];
                        }
                    if (buffer_new_edges == true)
                        {
                            alive_at_last_simplification.clear();
                            for (auto& p : parents)
                                {
                                    alive_at_last_simplification.push_back(p.node0);
                                    alive_at_last_simplification.push_back(p.node1);
                                }
                        }
                }
            else
                {
                    simplified = false;
                }
        }
    if (simplified == false)
        {
            samples.clear();
            for (auto& p : parents)
                {
                    samples.push_back(p.node0);
                    samples.push_back(p.node1);
                }
            node_map.resize(tables->nodes.num_rows);
            if (buffer_new_edges == false)
                {
                    sort_n_simplify(cppsort, parallel_sort, last_time_simplified,
                                    samples, node_map, tables);
                }
            else
                {
                    flush_buffer_n_simplify(alive_at_last_simplification, samples,
                                            node_map, new_edges, edge_liftover, tables);
                }
        }
}
