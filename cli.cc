#include "cli.hpp"
#include "options.hpp"

namespace po = boost::program_options;

po::options_description
generate_main_options(command_line_options &o)
{
    po::options_description options("Simulation options");
    options.add_options()("help", "Display help");
    options.add_options()("N", po::value<decltype(command_line_options::N)>(&o.N),
                          "Diploid population size. Default = 1000.");

    options.add_options()(
        "psurvival", po::value<decltype(command_line_options::psurvival)>(&o.psurvival),
        "Survival probability. Default = 0.0");
    options.add_options()("nsteps",
                          po::value<decltype(command_line_options::nsteps)>(&o.nsteps),
                          "Number of time steps to evolve. Default = 1000.");
    options.add_options()(
        "simplify",
        po::value<decltype(command_line_options::simplification_interval)>(
            &o.simplification_interval),
        "Time steps between simplifications.  Default = 100.");
    options.add_options()("rho", po::value<decltype(command_line_options::rho)>(&o.rho),
                          "Scaled recombination rate, 4Nr.  Default=0.");
    options.add_options()(
        "treefile", po::value<decltype(command_line_options::treefile)>(&o.treefile),
        "Ouput file name.  Default = treefile.trees");
    options.add_options()(
        "buffer",
        po::value<decltype(command_line_options::buffer_new_edges)>(&o.buffer_new_edges),
        "If true, use edge buffering algorithm. If not, sort and simplify. Default = "
        "false");
    options.add_options()(
        "seed",
        po::value<decltype(command_line_options::seed)>(&o.seed),
        "Random number seed.  Default = 42.");


    return options;
}

