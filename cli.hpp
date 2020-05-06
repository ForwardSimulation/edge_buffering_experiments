#pragma once

#include <boost/program_options.hpp>

class command_line_options;

boost::program_options::options_description
generate_main_options(command_line_options &o);

