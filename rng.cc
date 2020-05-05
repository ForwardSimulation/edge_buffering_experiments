#include "rng.hpp"

GSLrng
make_rng(unsigned seed)
{
    GSLrng rng(gsl_rng_alloc(gsl_rng_mt19937), [](gsl_rng* rng) { gsl_rng_free(rng); });
    gsl_rng_set(rng.get(), seed);
    return rng;
}
