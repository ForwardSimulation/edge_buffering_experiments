#pragma once

#include <memory>
#include <functional>
#include <gsl/gsl_rng.h>

using GSLrng = std::unique_ptr<gsl_rng, std::function<void(gsl_rng *)>>;


GSLrng make_rng(unsigned seed);
