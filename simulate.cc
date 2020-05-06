#include <iostream>
#include <memory>
#include <vector>
#include <gsl/gsl_randist.h>
#include "rng.hpp"
#include "tskit_tools.hpp"
#include "edge_buffer.hpp"

using edge_buffer_ptr = std::unique_ptr<EdgeBuffer>;

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

// NOTE: the next couple of functions are duplicates of logic
// already implemented in fwdpp::ts::ancestry_list.  In fact,
// the entire edge buffer idea has enough overlap that we
// should create a generic template type to reuse.
std::int32_t
get_buffer_end(const edge_buffer_ptr& new_edges, std::size_t i)
{
    auto f = new_edges->first[i];
    while (f != -1 && new_edges->next[f] != -1)
        {
            f = new_edges->next[f];
        }
    return f;
}

void
buffer_new_edge(tsk_id_t parent, double left, double right, double child,
                edge_buffer_ptr& new_edges)
{
    if (parent == TSK_NULL || child == TSK_NULL)
        {
            throw std::runtime_error("bad node IDs passed to buffer_new_edge");
        }
    if (parent >= new_edges->first.size())
        {
            new_edges->first.resize(parent, -1);
        }
    if (new_edges->first[parent] == -1)
        {
            new_edges->births.emplace_back(left, right, child);
            new_edges->first[parent] = new_edges->births.size() - 1;
            if (new_edges->births.size() >= new_edges->next.size())
                {
                    new_edges->next.push_back(-1);
                }
        }
    else
        {
            new_edges->next.push_back(-1);
            auto l = get_buffer_end(new_edges, parent);
            new_edges->next[l] = -1;
        }
}

void
stitch_together_edges(const std::vector<tsk_id_t>& alive_at_last_simplification,
                      double last_simplification_time, edge_buffer_ptr& new_edges,
                      table_collection_ptr& tables)
{
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

static void
generate_births(const GSLrng& rng, const std::vector<Birth>& births, double birth_time,
                bool buffer_new_edges, edge_buffer_ptr& new_edges,
                std::vector<Parent>& parents, table_collection_ptr& tables)
{
    for (auto& b : births)
        {
            auto new_node_0 = record_node(birth_time, tables);
            auto new_node_1 = record_node(birth_time, tables);
            auto parental_node0 = b.p0node0;
            if (gsl_rng_uniform(rng.get()) < 0.5)
                {
                    parental_node0 = b.p0node1;
                }
            auto parental_node1 = b.p1node0;
            if (gsl_rng_uniform(rng.get()) < 0.5)
                {
                    parental_node1 = b.p1node1;
                }
            if (buffer_new_edges == false)
                {
                    auto rv = tsk_edge_table_add_row(
                        &tables->edges, 0., tables->sequence_length, parental_node0,
                        new_node_0, nullptr, 0);
                    rv = tsk_edge_table_add_row(&tables->edges, 0.,
                                                tables->sequence_length, parental_node1,
                                                new_node_1, nullptr, 0);
                }
            else
                {
                    buffer_new_edge(parental_node0, 0., tables->sequence_length,
                                    new_node_0, new_edges);
                    buffer_new_edge(parental_node1, 0., tables->sequence_length,
                                    new_node_1, new_edges);
                }
            parents[b.index] = Parent(b.index, new_node_0, new_node_1);
        }
}

static void
sort_n_simplify(std::vector<Parent>& parents, table_collection_ptr& tables)
{
    int rv = tsk_table_collection_sort(tables.get(), nullptr, 0);
    // FIXME: these shouldn't be allocated each time
    std::vector<tsk_id_t> samples, node_map;
    for (auto& p : parents)
        {
            samples.push_back(p.node0);
            samples.push_back(p.node1);
        }
    node_map.resize(tables->nodes.num_rows);
    rv = tsk_table_collection_simplify(tables.get(), samples.data(), samples.size(), 0,
                                       node_map.data());
    for (auto& p : parents)
        {
            p.node0 = node_map[p.node0];
            p.node1 = node_map[p.node1];
        }
}

void
simulate(const GSLrng& rng, unsigned N, double psurvival, unsigned nsteps,
         unsigned simplification_interval, double rho, bool buffer_new_edges,
         table_collection_ptr& tables)
{
    if (rho > 0.)
        {
            throw std::invalid_argument("recombination not implemented yet");
        }
    std::vector<Parent> parents;
    for (unsigned i = 0; i < N; ++i)
        {
            auto id0 = record_node(nsteps, tables);
            auto id1 = record_node(nsteps, tables);
            parents.emplace_back(i, id0, id1);
        }

    edge_buffer_ptr new_edges(nullptr);
    if (buffer_new_edges)
        {
            new_edges.reset(new EdgeBuffer(tables->nodes.num_rows));
        }

    std::vector<Birth> births;
    bool simplified = false;
    for (unsigned step = 1; step <= nsteps; ++step)
        {
            deaths_and_parents(rng, parents, psurvival, births);
            generate_births(rng, births, nsteps - step, buffer_new_edges, new_edges,
                            parents, tables);
            if (step % simplification_interval == 0.)
                {
                    sort_n_simplify(parents, tables);
                    simplified = true;
                }
            else
                {
                    simplified = false;
                }
        }
    if (simplified == false)
        {
            sort_n_simplify(parents, tables);
        }
}
