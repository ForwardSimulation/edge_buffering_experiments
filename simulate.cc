#include <iostream>
#include <vector>
#include <gsl/gsl_randist.h>
#include "rng.hpp"
#include "tskit_tools.hpp"

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
            auto rv = tsk_edge_table_add_row(&tables->edges, 0., tables->sequence_length,
                                             parental_node0, new_node_0, nullptr, 0);
            rv = tsk_edge_table_add_row(&tables->edges, 0., tables->sequence_length,
                                        parental_node1, new_node_1, nullptr, 0);
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
    std::vector<Parent> parents;
    for (unsigned i = 0; i < N; ++i)
        {
            auto id0 = record_node(nsteps, tables);
            auto id1 = record_node(nsteps, tables);
            parents.emplace_back(i, id0, id1);
        }

    std::vector<Birth> births;
    bool simplified = false;
    for (unsigned step = 1; step <= nsteps; ++step)
        {
            deaths_and_parents(rng, parents, psurvival, births);
            generate_births(rng, births, nsteps - step, parents, tables);
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
