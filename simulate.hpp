#pragma once

class GSLrng;
class table_collection_ptr;

void simulate(const GSLrng& rng, unsigned N, double psurvival, unsigned nsteps,
              unsigned simplification_interval, table_collection_ptr& tables);
