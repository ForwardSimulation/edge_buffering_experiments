#pragma once

class GSLrng;
class table_collection_ptr;

void simulate(GSLrng& rng, table_collection_ptr& tables, unsigned N, double psurvival,
              unsigned nsteps, unsigned simplification_interval);
