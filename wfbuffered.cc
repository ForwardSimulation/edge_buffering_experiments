#include <cstdint>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <tskit.h>

#include <boost/program_options.hpp>

#include "rng.hpp"
#include "tskit_tools.hpp"
#include "simulate.hpp"
#include "options.hpp"
#include "cli.hpp"

namespace po = boost::program_options;

int
main(int argc, char **argv)
{
    command_line_options options;

    auto cli = generate_main_options(options);
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, cli), vm);
    po::notify(vm);
    validate_cli(options);

    if (vm.count("help"))
        {
            std::cout << cli << '\n';
            std::exit(1);
        }
    auto rng = make_rng(options.seed);
    auto tables = make_table_collection_ptr(1.);
    simulate(rng, options.N, options.psurvival, options.nsteps,
             options.simplification_interval, options.rho, options.buffer_new_edges,
             tables);
    auto ret = tsk_table_collection_build_index(tables.get(), 0);
    ret = tsk_table_collection_dump(tables.get(), options.treefile.c_str(), 0);
}
