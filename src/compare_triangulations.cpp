#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/FPU.h>
#include <CGAL/Interval_nt.h>
#include <CGAL_MWT/Compare_weights.h>
#include <CGAL_MWT/Read_instance.h>
#include <CGAL_MWT/time_util.h>
#include <algorithm>
#include <boost/program_options.hpp>
#include <boost/unordered/unordered_flat_set.hpp>
#include <iomanip>
#include <iostream>
#include <vector>

/**
 * Exit codes. 0 means that better <= worse (weight-wise).
 */
constexpr int COMPARE_SUCCESS = EXIT_SUCCESS;

/**
 * Exit code for invalid arguments.
 */
constexpr int INVALID_ARGUMENTS = 1;

/**
 * Exit code for unequal point sets.
 */
constexpr int UNEQUAL_POINT_SETS = 2;

/**
 * Exit code for file read errors.
 */
constexpr int FILE_READ_ERROR = 3;

/**
 * Exit code for unequal edge counts on
 * the same point set, indicating one of the
 * triangulations must be invalid.
 */
constexpr int UNEQUAL_EDGE_COUNT = 4;

/**
 * The triangulation comparator has encountered an unexpected internal error.
 */
constexpr int INTERNAL_ERROR = 5;

/**
 * The supposedly worse triangulation is actually better.
 */
constexpr int WORSE_IS_BETTER = 6;

/**
 * Epick should cover all our needs.
 */
using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point_2 = Kernel::Point_2;
using Edge = std::array<std::uint64_t, 2>;
using EdgeList = std::vector<Edge>;

/**
 * Normalize the points and edges of the triangulations
 * to have the same order.
 */
template<typename PointType> void normalize_points(std::vector<PointType> &points, EdgeList &edges) {
    std::vector<std::size_t> ordered_index;
    ordered_index.resize(points.size());
    std::iota(ordered_index.begin(), ordered_index.end(), std::size_t(0));
    std::sort(ordered_index.begin(), ordered_index.end(),
              [&points](std::size_t i, std::size_t j) { return points[i] < points[j]; });
    std::vector<std::size_t> to_ordered_index;
    to_ordered_index.resize(points.size(), 0);
    for(std::size_t i = 0, n = ordered_index.size(); i < n; ++i) {
        to_ordered_index[ordered_index[i]] = i;
    }
    std::vector<PointType> tmp_points;
    tmp_points.reserve(points.size());
    for(std::size_t oi : ordered_index) {
        tmp_points.push_back(points[oi]);
    }
    points = std::move(tmp_points);
    for(auto &edge : edges) {
        edge[0] = to_ordered_index[edge[0]];
        edge[1] = to_ordered_index[edge[1]];
    }
}

/**
 * Compute set of edges only in one of the triangulations.
 * Returns two sets (only in better, only in worse).
 * The edges are normalized (source index < target index).
 */
std::pair<EdgeList, EdgeList> compute_symmetric_difference(const EdgeList &better_edges, const EdgeList &worse_edges) {
    EdgeList only_in_better;
    EdgeList only_in_worse;
    boost::unordered_flat_set<Edge> better_set;
    better_set.reserve(better_edges.size());
    auto normalize_edge = [](Edge edge) {
        if(edge[0] > edge[1]) {
            std::swap(edge[0], edge[1]);
        }
        return edge;
    };
    for(const Edge &edge_ : better_edges) {
        Edge edge = normalize_edge(edge_);
        better_set.insert(edge);
    }
    for(const Edge &edge_ : worse_edges) {
        Edge edge = normalize_edge(edge_);
        auto it = better_set.find(edge);
        if(it == better_set.end()) {
            only_in_worse.push_back(edge);
        } else {
            better_set.erase(it);
        }
    }
    only_in_better.assign(better_set.begin(), better_set.end());
    return {std::move(only_in_better), std::move(only_in_worse)};
}

/**
 * Compute an interval containing the difference of the two edge weights
 * of better and worse on the given point set.
 */
CGAL::Interval_nt_advanced compute_interval_difference(const EdgeList &better, const EdgeList &worse,
                                                       const std::vector<Point_2> &points) {
    CGAL::Protect_FPU_rounding pfr;
    std::vector<CGAL::Interval_nt_advanced> values;
    values.reserve(better.size() + worse.size());
    auto iv_dist = [](const Point_2 &p1, const Point_2 &p2) -> CGAL::Interval_nt_advanced {
        CGAL::Interval_nt_advanced p1x(p1.x()), p1y(p1.y());
        p1x -= p2.x();
        p1y -= p2.y();
        p1x *= p1x;
        p1y *= p1y;
        p1x += p1y;
        p1x = CGAL::sqrt(p1x);
        return p1x;
    };
    for(const Edge &e : better) {
        values.push_back(iv_dist(points[e[0]], points[e[1]]));
    }
    for(const Edge &e : worse) {
        values.push_back(-iv_dist(points[e[0]], points[e[1]]));
    }
    std::sort(values.begin(), values.end(), [](CGAL::Interval_nt_advanced i1, CGAL::Interval_nt_advanced i2) {
        return CGAL::abs(i1).inf() < CGAL::abs(i2).inf();
    });
    CGAL::Interval_nt_advanced iv_diff(0);
    for(const auto &iv : values) {
        iv_diff += iv;
    }
    return iv_diff;
}

int main(int argc, char **argv) {
    std::cerr << std::setprecision(19);
    std::cout << std::setprecision(19);

    /**
     * Command-line parsing.
     */
    namespace po = boost::program_options;
    po::options_description desc("Allowed options");
    std::string better_input_filename;
    bool better_is_line_format = false;
    std::string worse_input_filename;
    bool worse_is_line_format = false;
    bool precise_interval_difference = false;
    desc.add_options()("help", "produce this help message")(
        "better-triangulation-file,i", po::value<std::string>(&better_input_filename),
        "specify the file containing a better triangulation to load")(
        "better-expect-line-format", "expect a line-based format instead of JSON for the better triangulation")(
        "worse-triangulation-file,j", po::value<std::string>(&worse_input_filename),
        "specify the file containing a worse triangulation to load")(
        "worse-expect-line-format", "expect a line-based format instead of JSON for the worse triangulation")(
        "precise-interval-difference",
        "compute a precise interval difference between the two triangulations if it is not 0");
    po::variables_map vmap;
    try {
        po::store(po::parse_command_line(argc, argv, desc), vmap);
        po::notify(vmap);
    } catch(std::exception &ex) {
        std::cerr << "Could not parse command line: " << ex.what() << std::endl;
        return INVALID_ARGUMENTS;
    }
    if(vmap.count("help")) {
        std::cout << desc << std::endl;
        return COMPARE_SUCCESS;
    }
    if(!vmap.count("better-triangulation-file") || !vmap.count("worse-triangulation-file")) {
        std::cerr << "Both better and worse triangulation files must be specified." << std::endl;
        return INVALID_ARGUMENTS;
    }
    better_is_line_format = vmap.count("better-expect-line-format");
    worse_is_line_format = vmap.count("worse-expect-line-format");
    precise_interval_difference = vmap.count("precise-interval-difference");

    /**
     * Reading the triangulations.
     */
    std::vector<Point_2> better_points, worse_points;
    std::vector<std::array<std::uint64_t, 2>> better_edges, worse_edges;
    try {
        if(better_is_line_format) {
            std::tie(better_points, better_edges) = mwt::read_text_solution<Kernel>(better_input_filename);
        } else {
            std::tie(better_points, better_edges) = mwt::read_solution<Kernel>(better_input_filename);
        }
        if(worse_is_line_format) {
            std::tie(worse_points, worse_edges) = mwt::read_text_solution<Kernel>(worse_input_filename);
        } else {
            std::tie(worse_points, worse_edges) = mwt::read_solution<Kernel>(worse_input_filename);
        }
    } catch(std::exception &ex) {
        std::cerr << "Failed to read input files: " << ex.what() << std::endl;
        return FILE_READ_ERROR;
    }

    /**
     * Make sure we have a common point set with common order.
     */
    if(worse_points.size() != better_points.size()) {
        std::cerr << "The number of points in the better and worse triangulations must be the same." << std::endl;
        return UNEQUAL_POINT_SETS;
    }
    normalize_points(better_points, better_edges);
    normalize_points(worse_points, worse_edges);
    auto [iter_b, iter_w] = std::mismatch(better_points.begin(), better_points.end(), worse_points.begin());
    if(iter_b != better_points.end()) {
        std::cerr << "The point sets are not the same." << std::endl;
        std::cerr << "The better point set contains the point (" << *iter_b << ") compared to (" << *iter_w
                  << ") in the worse point set!" << std::endl;
        return UNEQUAL_POINT_SETS;
    }
    if(worse_edges.size() != better_edges.size()) {
        std::cerr << "The number of edges in the better and worse triangulations must be the same." << std::endl;
        return UNEQUAL_EDGE_COUNT;
    }

    /**
     * Compute the symmetric difference and check for equality.
     */
    auto [only_in_better, only_in_worse] = compute_symmetric_difference(better_edges, worse_edges);
    if(only_in_better.empty() != only_in_worse.empty()) {
        std::cerr << "Internal error in symmetric difference computation!" << std::endl;
        std::cerr << "Better edge count: " << better_edges.size() << ", worse edge count: " << worse_edges.size()
                  << std::endl;
        std::cerr << "Only in better: " << only_in_better.size() << ", only in worse: " << only_in_worse.size()
                  << std::endl;
        return INTERNAL_ERROR;
    }
    if(only_in_better.empty()) {
        std::cout << "The triangulations are exactly the same." << std::endl;
        std::cout << "Interval difference: [0.0, 0.0]" << std::endl;
        return COMPARE_SUCCESS;
    }

    /**
     * Do the actual weight comparison.
     */
    mwt::CompareWeights<Kernel> compare_weights;
    for(const Edge &e : only_in_better) {
        auto i1 = e[0], i2 = e[1];
        compare_weights.add_lhs(better_points[i1], better_points[i2]);
    }
    for(const Edge &e : only_in_worse) {
        auto i1 = e[0], i2 = e[1];
        compare_weights.add_rhs(better_points[i1], better_points[i2]);
    }
    auto sign = compare_weights.sign();
    CGAL::Interval_nt_advanced iv_diff;
    if(precise_interval_difference) {
        iv_diff = compare_weights.compute_interval();
    } else {
        iv_diff = compute_interval_difference(only_in_better, only_in_worse, better_points);
    }
    if(sign == CGAL::NEGATIVE) {
        std::cout << "The better triangulation is strictly better than the worse triangulation." << std::endl;
        std::cout << "Interval difference: [" << iv_diff.inf() << ", " << iv_diff.sup() << "]" << std::endl;
        return COMPARE_SUCCESS;
    }
    if(sign == CGAL::ZERO) {
        std::cout << "The triangulations differ, but have exactly the same weight." << std::endl;
        std::cout << "Interval difference: [0.0, 0.0]" << std::endl;
        return COMPARE_SUCCESS;
    }
    std::cout << "The triangulation given as 'worse' is better than the triangulation given as 'better'." << std::endl;
    std::cout << "Interval difference: [" << iv_diff.inf() << ", " << iv_diff.sup() << "]" << std::endl;
    return WORSE_IS_BETTER;
}
