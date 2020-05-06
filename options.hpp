#pragma once

#include <string>

struct command_line_options
{
    unsigned N;
    double psurvival;
    unsigned nsteps;
    unsigned simplification_interval;
    double rho;
    std::string treefile;
    bool buffer_new_edges;
    unsigned seed;

    command_line_options();
};

void validate_cli(const command_line_options &);
