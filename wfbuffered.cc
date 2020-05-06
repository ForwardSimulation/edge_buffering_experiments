#include <cstdint>
#include <cstddef>
#include <cstdio>
#include <iostream>
#include <tskit.h>

#include "rng.hpp"
#include "tskit_tools.hpp"
#include "simulate.hpp"

int
main(int argc, char **argv)
{
    auto rng = make_rng(42);
    auto tables = make_table_collection_ptr(11.);
    simulate(rng, 10000, 0, 10000, 1000, tables);
    auto ret = tsk_table_collection_dump(tables.get(), "ts.trees", 0);
    std::cout << tables->nodes.num_rows << ' ' << tables->edges.num_rows << '\n';
}
