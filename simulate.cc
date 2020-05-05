#include <vector>
#include <gsl/gsl_randist.h>
#include "rng.hpp"
#include "tskit_tools.hpp"

namespace
{
    struct IndexAndNodes
    {
        std::size_t index;
        tsk_id_t node0, node1;
        IndexAndNodes(std::size_t i, tsk_id_t n0, tsk_id_t n1)
            : index(i), node0(n0), node1(n1)
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

void
simulate(GSLrng& rng, table_collection_ptr& tables, unsigned N, double psurvival,
         unsigned nsteps, unsigned simplification_interval)
{
    std::vector<IndexAndNodes> parents;
    for (unsigned i = 0; i < 2 * N; ++i)
        {
            auto id0 = record_node(nsteps, tables);
            auto id1 = record_node(nsteps, tables);
            parents.emplace_back(i, id0, id1);
        }

    for (unsigned step = 1; step < nsteps; ++step)
        {
        }
}
