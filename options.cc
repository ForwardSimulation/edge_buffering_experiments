#include <cmath>
#include <stdexcept>
#include "options.hpp"

command_line_options::command_line_options()
    : N{1000}, psurvival{0.}, nsteps{1000},
      simplification_interval{100}, rho{0.}, treefile{"treefile.trees"},
      buffer_new_edges{false}, cppsort{false}, parallel_sort{false}, seed{42}
{
}

void
validate_cli(const command_line_options& options)
{
    if (options.N == 0)
        {
            throw std::invalid_argument("Population size must be > 0");
        }

    if (options.psurvival < 0. || options.psurvival >= 1.
        || std::isfinite(options.psurvival) == false)
        {
            throw std::invalid_argument("psurvival must be 0.0 <= p < 1.0");
        }

    if (options.rho < 0.0 || std::isfinite(options.rho) == false)
        {
            throw std::invalid_argument("rho must be >= 0.0");
        }

    if (options.treefile.empty())
        {
            throw std::invalid_argument("treefile must not be an empty string");
        }
}
