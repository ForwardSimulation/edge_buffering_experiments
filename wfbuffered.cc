#include <cstdint>
#include <cstddef>
#include <cstdio>
#include <tskit.h>

#include "rng.hpp"
#include "tskit_tools.hpp"

int
main(int argc, char **argv)
{
    auto rng = make_rng(42);
    auto tables = make_table_collection_ptr(11.);
}
