#ifndef CGAL_MWT_SRC_OPTIONS_H_INCLUDED_
#define CGAL_MWT_SRC_OPTIONS_H_INCLUDED_

#include <algorithm>
#include <boost/program_options.hpp>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <random>
#include <string>

namespace mwt {

using OptDescription = boost::program_options::options_description;

template<typename T> void option_with_default(OptDescription &desc, const char *name, T &output, const char *help) {
    desc.add_options()(name, boost::program_options::value<T>(&output)->default_value(output), help);
}

inline void bool_switch(OptDescription &desc, bool &b, const char *name, const char *help) {
    desc.add_options()(name, boost::program_options::bool_switch(&b), help);
}

inline void random_seed_option(OptDescription &desc, std::size_t &seed) {
    seed = std::numeric_limits<std::size_t>::max();
    option_with_default(desc, "seed", seed, "seed for random number generator (random if not specified)");
}

inline void num_points_option(OptDescription &desc, std::size_t &num_points) {
    option_with_default(desc, "num-points,n", num_points, "number of points to generate");
}

inline void grid_size_option(OptDescription &desc, std::uint32_t &grid_size) {
    option_with_default(desc, "grid-size,s", grid_size, "grid size to generate points in");
}

inline void print_step_times_option(OptDescription &desc, bool &print_times) {
    bool_switch(desc, print_times, "print-step-times", "print times for each main step of the algorithm");
}

inline void instance_file_option(OptDescription &desc, std::string &filename) {
    filename = "";
    option_with_default(desc, "instance-file,i", filename,
                        "file to read instance from, instead of generating a random one");
}

inline void output_file_option(OptDescription &desc, std::string &filename) {
    filename = "";
    option_with_default(desc, "output-file,o", filename, "file to write output to, instead of stdout");
}

inline void step_file_option(OptDescription &desc, const char *optname, std::string &filename, const char *help) {
    filename = "";
    option_with_default(desc, optname, filename, help);
}

inline void verification_option(OptDescription &desc, std::string &verification_mode) {
    verification_mode = "quick";
    option_with_default(desc, "verification", verification_mode, "verification mode ('none', 'quick', or 'full')");
}

inline void nonsimple_face_option(OptDescription &desc, std::string &nonsimple_face_mode) {
    nonsimple_face_mode = "inexact_with_gap";
    option_with_default(desc, "nonsimple-face-mode", nonsimple_face_mode,
                        "nonsimple face triangulation mode ('inexact_with_gap', 'exact')");
}

inline void only_dump_instance_option(OptDescription &desc, bool &only_dump_instance) {
    bool_switch(desc, only_dump_instance, "only-dump-instance",
                "only dump the instance in JSON format and exit without solving");
}

inline void lp_verbose_option(OptDescription &desc, bool &lp_verbose) {
    bool_switch(desc, lp_verbose, "lp-verbose", "print output from LP/MIP solver backend");
}

inline void resource_monitor_option(OptDescription &desc, bool &monitor_resources) {
    bool_switch(desc, monitor_resources, "monitor-resources",
                "monitor and record maximum resource usage during execution; support depends on platform. "
                "Unsupported values are set to -1 in the output.");
}

inline void use_cuts_option(OptDescription &desc, bool &use_cuts) {
    bool_switch(desc, use_cuts, "use-cuts", "use cuts in the LP before resorting to MIP");
}

inline void serial_lmt_option(OptDescription &desc, bool &serial_lmt) {
    bool_switch(desc, serial_lmt, "serial-lmt", "use serial LMT instead of parallel LMT");
}

inline std::mt19937_64 initialize_rng(std::size_t seed) {
    if(seed == std::numeric_limits<std::size_t>::max()) {
        std::random_device rd;
        seed = rd();
    }
    return std::mt19937_64(seed);
}

inline bool have_option(const OptDescription &desc, const char *name) {
    return std::any_of(desc.options().begin(), desc.options().end(),
                       [name](const auto &opt) { return opt->long_name() == name; });
}

inline void parse_options(int argc, char **argv, OptDescription &desc) {
    const auto &opts = desc.options();
    bool help_flag = false;
    if(!have_option(desc, "help")) {
        bool_switch(desc, help_flag, "help", "produce this help message");
    }
    bool have_verification_mode = have_option(desc, "verification");
    bool have_nonsimple_face_mode = have_option(desc, "nonsimple-face-mode");
    try {
        boost::program_options::variables_map vm;
        boost::program_options::store(boost::program_options::parse_command_line(argc, argv, desc), vm);
        boost::program_options::notify(vm);
        if(have_verification_mode && vm.count("verification")) {
            std::string mode = vm["verification"].as<std::string>();
            if(mode != "none" && mode != "quick" && mode != "full") {
                std::cerr << "Unknown verification mode '" << mode << "'!" << std::endl;
                std::cerr << desc << std::endl;
                std::exit(EXIT_FAILURE);
            }
        }
        if(have_nonsimple_face_mode && vm.count("nonsimple-face-mode")) {
            std::string mode = vm["nonsimple-face-mode"].as<std::string>();
            if(mode != "inexact_with_gap" && mode != "exact") {
                std::cerr << "Unknown nonsimple face triangulation mode '" << mode << "'!" << std::endl;
                std::cerr << desc << std::endl;
                std::exit(EXIT_FAILURE);
            }
        }
        if(vm.count("help") && vm["help"].as<bool>()) {
            std::cerr << desc << std::endl;
            std::exit(EXIT_SUCCESS);
        }
    } catch(std::exception &ex) {
        std::cerr << "Could not parse commandline options: " << ex.what() << std::endl;
        std::cerr << desc << std::endl;
        std::exit(EXIT_FAILURE);
    }
}

} // namespace mwt

#endif
