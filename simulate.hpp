#pragma once

#include "rng.hpp"
#include "tskit_tools.hpp"

void simulate(const GSLrng& rng, unsigned N, double psurvival, unsigned nsteps,
              unsigned simplification_interval, double rho, bool buffer_new_edges,
              table_collection_ptr& tables);
