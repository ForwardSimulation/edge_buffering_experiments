#include <cstdint>
#include <cstddef>
#include <cstdio>
#include <tskit.h>

#include "rng.hpp"

int
main(int argc, char **argv)
{
    auto rng = make_rng(42);
}
