#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL_MWT/Read_instance.h>
#include <CGAL_MWT/Validate.h>
#include <CGAL_MWT/time_util.h>
#include <boost/program_options.hpp>
#include <cstdlib>
#include <functional>

template<typename Points, typename Edges>
void dump_instance(const Points &points, const Edges &edges, const std::string &filename) {
    std::ofstream output(filename.c_str(), std::ios::out | std::ios::trunc);
    output << "{\n";
    output << "  \"points\": [\n";
    output << std::setprecision(19);
    bool first = true;
    for(const auto &p : points) {
        if(!first) {
            output << ",\n";
        }
        first = false;
        output << "    [" << p.x() << ", " << p.y() << "]";
    }
    output << "\n  ],\n  \"edges\": [\n";
    first = true;
    for(const auto &e : edges) {
        if(!first) {
            output << ",\n";
        }
        first = false;
        output << "    [" << e[0] << ", " << e[1] << "]";
    }
    output << "\n  ]\n";
    output << "}\n";
}

int main(int argc, char **argv) {
    using namespace mwt;
    namespace po = boost::program_options;
    po::options_description desc("Allowed options");
    std::string xfree_option, edge_count_mode;
    std::string input_filename, dump_json_path;
    bool print_steps;
    bool expect_line_format;
    desc.add_options()("help", "produce this help message")("triangulation-file,i",
                                                            po::value<std::string>(&input_filename),
                                                            "specify the file containing a triangulation to load")(
        "xfree-mode", po::value<std::string>(&xfree_option)->default_value("auto"),
        "specify how crossing-freeness shall be detected ('auto', 'sweep', 'cgal', 'both', 'none')")(
        "print-steps", "print times for each main step of the algorithm")(
        "edge-count-mode", po::value<std::string>(&edge_count_mode)->default_value("cgal"),
        "specify how to check the number of edges in the triangulation ('cgal', 'none')")(
        "expect-line-format", "expect a line-based format instead of JSON")(
        "dump-json-to", po::value<std::string>(&dump_json_path)->default_value(""),
        "output the triangulation in JSON format to the indicated file (default: do not output)");

    try {
        po::variables_map vm;
        po::store(po::parse_command_line(argc, const_cast<const char *const *>(argv), desc), vm);
        po::notify(vm);

        if(vm.count("help")) {
            std::cout << desc << std::endl;
            std::exit(EXIT_SUCCESS);
        }

        expect_line_format = vm.count("expect-line-format");

        print_steps = vm.count("print-steps");

        std::vector<std::string> xfree_options{"auto", "sweep", "cgal", "both", "none"};
        if(std::find(xfree_options.begin(), xfree_options.end(), xfree_option) == xfree_options.end()) {
            throw po::validation_error(po::validation_error::invalid_option_value, "xfree-mode");
        }
        std::vector<std::string> edge_count_modes{"cgal", "none"};
        if(std::find(edge_count_modes.begin(), edge_count_modes.end(), edge_count_mode) == edge_count_modes.end()) {
            throw po::validation_error(po::validation_error::invalid_option_value, "edge-count-mode");
        }
    } catch(std::exception &ex) {
        std::cerr << "Error: " << ex.what() << std::endl;
        std::cerr << desc << std::endl;
        std::exit(EXIT_FAILURE);
    }

    auto solution = !expect_line_format ? read_solution<CGAL::Epick>(input_filename)
                                        : read_text_solution<CGAL::Epick>(input_filename);
    if(print_steps) {
        std::cout << "Read solution with " << solution.first.size() << " points and " << solution.second.size()
                  << " edges" << std::endl;
    }

    LoadedMWTVerifier<CGAL::Epick> verifier(std::move(solution.first), std::move(solution.second));
    if(print_steps) {
        std::cout << "Created verifier..." << std::endl;
    }

    std::map<std::string, std::function<void()>> xfree_verifiers{
        {"none", []() {}},
        {"sweep",
         [&]() {
             double t = measure_time([&]() { verifier.check_no_intersections(); });
             if(print_steps) {
                 std::cout << "No intersections found using sweep-line algorithm (" << t << " seconds)" << std::endl;
             }
         }},
        {"cgal",
         [&]() {
             double t = measure_time([&]() { verifier.check_no_intersections_cgal_exact(); });
             if(print_steps) {
                 std::cout << "No intersections found using CGAL's exact intersection check (" << t << " seconds)"
                           << std::endl;
             }
         }},
        {"both",
         [&]() {
             double t1 = measure_time([&]() { verifier.check_no_intersections(); });
             if(print_steps) {
                 std::cout << "No intersections found using sweep-line algorithm (" << t1 << " seconds)" << std::endl;
             }
             double t2 = measure_time([&]() { verifier.check_no_intersections_cgal_exact(); });
             if(print_steps) {
                 std::cout << "No intersections found using CGAL's exact intersection check (" << t2 << " seconds)"
                           << std::endl;
             }
         }},
        {"auto", [&]() {
             double t = measure_time([&verifier] { verifier.check_no_intersections(); });
             if(print_steps) {
                 std::cout << "No intersections found using sweep-line algorithm (" << t << " seconds)" << std::endl;
             }
             if(t < 0.5) {
                 if(print_steps) {
                     std::cout << "Double-checking using CGAL..." << std::endl;
                 }
                 t = measure_time([&verifier] { verifier.check_no_intersections_cgal_exact(); });
                 if(print_steps) {
                     std::cout << "No intersections found using CGAL's exact intersection check (" << t << " seconds)"
                               << std::endl;
                 }
             }
         }}};

    try {
        if(edge_count_mode == "cgal") {
            double tdel = measure_time([&]() { verifier.check_edge_count(); });
            if(print_steps) {
                std::cout << "Computed Delaunay triangulation; edge count matches expected value (" << tdel
                          << " seconds)" << std::endl;
            }
        }
        xfree_verifiers.at(xfree_option)();
        if(!dump_json_path.empty()) {
            dump_instance(verifier.get_points(), verifier.get_edges(), dump_json_path);
        }
        std::exit(EXIT_SUCCESS);
    } catch(const std::runtime_error &err) {
        std::cerr << err.what() << std::endl;
    }
    std::exit(EXIT_FAILURE);
}
