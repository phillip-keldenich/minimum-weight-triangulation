#include "options.h"
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/FPU.h>
#include <CGAL_MWT/Are_all_collinear.h>
#include <CGAL_MWT/Diamond_filter.h>
#include <CGAL_MWT/Exact_simple_face_triangulator.h>
#include <CGAL_MWT/Face_collector.h>
#include <CGAL_MWT/Generate_random_integral.h>
#include <CGAL_MWT/LMT_skeleton.h>
#include <CGAL_MWT/Mwt_traits.h>
#include <CGAL_MWT/Read_instance.h>
#include <CGAL_MWT/Resource_monitor.h>
#include <CGAL_MWT/Select_LP_backend.h>
#include <CGAL_MWT/Static_quadtree.h>
#include <CGAL_MWT/Validate.h>
#include <CGAL_MWT/output.h>
#include <CGAL_MWT/time_util.h>
#include <chrono>
#include <iostream>
#include <string>

using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using Traits = mwt::Mwt_traits_2<Kernel>;
using Point = Traits::Point_2;
using Rect = Kernel::Iso_rectangle_2;
using Points = std::vector<Point>;
using PointIterator = Points::iterator;
using Tree = mwt::Point_quadtree<Traits, PointIterator>;
using Filter = mwt::Diamond_edge_filter<Tree, true>;

void print_time(const std::string &name, double time, bool print_times) {
    if(print_times) {
        std::cout << name << " took " << time << " seconds" << std::endl;
    }
}

template<typename Skeleton>
void skeleton_to_file(const std::string &name, const Skeleton &skeleton, bool with_possible,
                      const mwt::InfoMap &info_map = {}) {
    if(!name.empty()) {
        std::ofstream output;
        output.exceptions(std::ofstream::failbit | std::ofstream::badbit);
        output.open(name.c_str(), std::ios::out | std::ios::trunc);
        mwt::output_triangulation(output, skeleton.get_tree(), skeleton, with_possible, info_map);
    }
}

// make sure we reset the rounding mode at the end
// to avoid bad CGAL assertions.
CGAL::Protect_FPU_rounding round_protect;

int main(int argc, char **argv) {
    using namespace mwt;

    std::size_t seed;
    std::size_t num_points = 0;
    std::uint32_t grid_size = 1 << 27;
    bool print_times = false;
    bool only_dump_instance = false;
    bool lp_verbose = false;
    bool use_cuts = false;
    bool serial_lmt = false;
    bool monitor_resources = false;
    std::string filename, out_filename;
    std::string diamond_filename, lmt_filename, advanced_lmt_filename;
    std::string nonsimple_face_mode;
    std::string verification_mode;
    OptDescription desc("Allowed options");

    random_seed_option(desc, seed);
    num_points_option(desc, num_points);
    grid_size_option(desc, grid_size);
    print_step_times_option(desc, print_times);
    instance_file_option(desc, filename);
    output_file_option(desc, out_filename);
    serial_lmt_option(desc, serial_lmt);
    if(mwt::CGAL_MWT_RESOURCE_MONITOR_SUPPORTED) {
        resource_monitor_option(desc, monitor_resources);
    }
    step_file_option(desc, "diamond", diamond_filename, "Dump the diamond filtered edges to the indicated file");
    step_file_option(desc, "lmt", lmt_filename,
                     "Dump the LMT skeleton edges (prior to advanced LMT) to the indicated file");
    step_file_option(desc, "advanced-lmt", advanced_lmt_filename,
                     "Dump the LMT skeleton edges (after advanced LMT) to the indicated file");
    nonsimple_face_option(desc, nonsimple_face_mode);
    verification_option(desc, verification_mode);
    only_dump_instance_option(desc, only_dump_instance);
    lp_verbose_option(desc, lp_verbose);
    use_cuts_option(desc, use_cuts);
    parse_options(argc, argv, desc);
    mwt::ResourceMonitor resource_monitor;
    if(monitor_resources) {
        resource_monitor.start();
    }
    std::mt19937_64 rng = initialize_rng(seed);

    // point generation or instance reading
    Points points;
    double gen_time = -1.0;
    if(filename.empty()) {
        if(num_points == 0) {
            std::cerr << "You must specify an instance file name or a number of points to generate!\n";
            std::cerr << desc << std::endl;
            std::exit(EXIT_FAILURE);
        }

        // generate points
        gen_time = measure_time([&]() {
            points = generate_random_integral<Kernel>(Rect{Point{0, 0}, Point{grid_size, grid_size}}, num_points, rng);
        });
        print_time("Generating points", gen_time, print_times);
    } else {
        try {
            // read points from file
            points = read_instance<Kernel>(filename);
        } catch(std::exception &error) {
            std::cerr << "Failed to read instance file '" << filename << "': " << error.what() << std::endl;
            return EXIT_FAILURE;
        }
    }

    // option to only dump the instance; useful to convert instances to JSON format
    if(only_dump_instance) {
        std::ostream *output = &std::cout;
        std::ofstream fout;
        fout.exceptions(std::ofstream::failbit | std::ofstream::badbit);
        if(!out_filename.empty()) {
            fout.open(out_filename.c_str(), std::ios::out | std::ios::trunc);
            output = &fout;
        }
        output_only_points(*output, points.begin(), points.end());
        return EXIT_SUCCESS;
    }

    // does the instance have a triangulation?
    if(mwt::are_all_points_collinear(points.begin(), points.end())) {
        std::cerr << "You asked for an MWT for a point set of only collinear points!\n";
        std::exit(EXIT_FAILURE);
    }

    // build tree
    std::optional<Tree> tree_;
    double tree_time = measure_time([&]() { tree_.emplace(points.begin(), points.end()); });
    Tree &tree = *tree_;
    print_time("Building quadtree", tree_time, print_times);

    // run diamond filter
    std::optional<Filter> filter_;
    double filter_time = measure_time([&]() {
        filter_.emplace(&tree);
        filter_->compute_remaining_edges();
    });
    Filter &filter = *filter_;
    std::size_t num_diamond_edges = filter.get_edges().size();
    print_time("Diamond filter (" + std::to_string(num_diamond_edges) + " edges resulting)", filter_time, print_times);

    // run normal LMT
    using Skeleton = LMTSkeleton<Traits, Tree>;
    std::optional<Skeleton> skeleton_builder_;
    double lmt_time = measure_time([&]() {
        skeleton_builder_.emplace(&tree, std::move(filter.get_edges()));
        skeleton_to_file(diamond_filename, *skeleton_builder_, true);
        skeleton_builder_->lmt_main_loop(!serial_lmt);
    });
    Skeleton &skeleton_builder = *skeleton_builder_;
    std::size_t lmt_fixed_edges = skeleton_builder.skeleton().size();
    std::size_t lmt_candidate_edges = skeleton_builder.candidates().size();
    print_time("Normal LMT loop (" + std::to_string(lmt_fixed_edges) + " edges fixed, " +
                   std::to_string(lmt_candidate_edges) + " remaining candidates; full triangulation must have " +
                   std::to_string(skeleton_builder.expected_total_triangulation_edges()) + " edges)",
               lmt_time, print_times);
    skeleton_to_file(lmt_filename, skeleton_builder, true);

    // advanced LMT loop
    double adv_lmt_time = measure_time([&]() { skeleton_builder.advanced_lmt_loop(); });
    std::size_t adv_lmt_fixed_edges = skeleton_builder.skeleton().size();
    std::size_t adv_lmt_candidate_edges = skeleton_builder.candidates().size();
    print_time("Advanced LMT loop (" + std::to_string(adv_lmt_fixed_edges) + " edges fixed, " +
                   std::to_string(adv_lmt_candidate_edges) + " remaining candidates)",
               adv_lmt_time, print_times);
    skeleton_builder.advanced_lmt_loop();
    skeleton_to_file(advanced_lmt_filename, skeleton_builder, true);

    // triangulate simple faces
    double max_gap = 0.0;
    nlohmann::json face_triangulation_stats;
    mwt::FaceTriangulatorOptions face_triangulator_options;
    face_triangulator_options.use_cuts = use_cuts;
    face_triangulator_options.lp_verbose = lp_verbose;
    double face_triangulation_time = measure_time([&]() {
        if(nonsimple_face_mode == "exact") {
            face_triangulation_stats = triangulate_nonsimple_faces_with_lp<Face_collector<Skeleton>>(
                skeleton_builder, face_triangulator_options);
        } else {
            face_triangulation_stats = triangulate_nonsimple_faces_with_gap<Face_collector<Skeleton>>(
                skeleton_builder, &max_gap, face_triangulator_options);
        }
    });
    print_time("Triangulating faces", face_triangulation_time, print_times);
    mwt::InfoMap info_map;
    info_map["num_points"] = std::int64_t(tree.size());
    info_map["max_gap"] = max_gap;
    info_map["gen_time"] = gen_time;
    info_map["diamond_time"] = filter_time;
    info_map["edges_passing_diamond_filter"] = std::int64_t(num_diamond_edges);
    info_map["lmt_time"] = lmt_time;
    info_map["lmt_fixed_edges"] = std::int64_t(lmt_fixed_edges);
    info_map["lmt_candidate_edges"] = std::int64_t(lmt_candidate_edges);
    info_map["advanced_lmt_time"] = adv_lmt_time;
    info_map["advanced_lmt_fixed_edges"] = std::int64_t(adv_lmt_fixed_edges);
    info_map["advanced_lmt_candidate_edges"] = std::int64_t(adv_lmt_candidate_edges);
    info_map["face_triangulation_time"] = face_triangulation_time;
    info_map["face_triangulation_stats"] = face_triangulation_stats;
    info_map["verification_mode"] = verification_mode;
    info_map["instance_file"] = filename;
    CGAL::Interval_nt_advanced weight = mwt::total_weight<CGAL::Interval_nt_advanced>(skeleton_builder);
    info_map["total_weight_lb"] = weight.inf();
    info_map["total_weight_ub"] = weight.sup();

    if(verification_mode != "none") {
        mwt::MWTVerifier<Skeleton> verifier{&skeleton_builder};
        if(verification_mode == "full") {
            verifier.strong_check();
            print_time("Verifier: Delaunay triangulation", verifier.get_delaunay_time(), print_times);
            info_map["verification_delaunay_time"] = verifier.get_delaunay_time();
            print_time("Verifier: Intersection check", verifier.get_intersection_time(), print_times);
            info_map["verification_intersection_time"] = verifier.get_intersection_time();
        } else {
            verifier.quick_check();
            print_time("Verifier: Delaunay triangulation", verifier.get_delaunay_time(), print_times);
            info_map["verification_delaunay_time"] = verifier.get_delaunay_time();
        }
    }
    if(monitor_resources) {
        resource_monitor.stop();
        info_map["max_resource_usage"] = resource_monitor.get_max_usage().to_json();
    }
    skeleton_to_file(out_filename, skeleton_builder, false, info_map);
    return 0;
}
