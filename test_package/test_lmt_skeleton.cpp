#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/FPU.h>
#include <CGAL/Interval_arithmetic.h>
#include <CGAL_MWT/Diamond_filter.h>
#include <CGAL_MWT/Directional_filter.h>
#include <CGAL_MWT/Exact_diamond_filtered_search_driver.h>
#include <CGAL_MWT/Generate_random_integral.h>
#include <CGAL_MWT/LMT_skeleton.h>
#include <CGAL_MWT/Mwt_traits.h>
#include <CGAL_MWT/Static_quadtree.h>
#include <chrono>
#include <doctest/doctest.h>
#include <iomanip>
#include <random>

namespace mwt {

namespace test {

class LMTSkeletonTester {
  public:
    template<typename Kernel>
    static void random_initialization_test(std::size_t grid_dim, std::size_t num_points, std::size_t seed) {
        using Traits = mwt::Mwt_traits_2<Kernel>;
        using Point = typename Traits::Point_2;
        using Points = std::vector<Point>;
        using PointIterator = typename Points::iterator;
        using Tree = mwt::Point_quadtree<Kernel, PointIterator>;
        using LMTSkeleton = mwt::LMTSkeleton<Traits, Tree>;

        CGAL::FPU_set_cw(CGAL_FE_UPWARD);
        Points points =
            mwt::generate_random_integral<Kernel>({{-grid_dim, -grid_dim}, {grid_dim, grid_dim}}, num_points, seed);
        Tree tree(points.begin(), points.end());
        mwt::Diamond_edge_filter<Tree, true> diamond_filter(&tree);
        diamond_filter.compute_remaining_edges();
        LMTSkeleton skeleton(&tree, std::move(diamond_filter.get_edges()));
        CHECK(diamond_filter.get_edges().size() == 0);
        skeleton.p_validate_initial_invariants();
    }

    template<typename Kernel>
    static void random_lmt_test(std::size_t grid_dim, std::size_t num_points, std::size_t seed) {
        using Traits = mwt::Mwt_traits_2<Kernel>;
        using Point = typename Traits::Point_2;
        using Points = std::vector<Point>;
        using PointIterator = typename Points::iterator;
        using Tree = mwt::Point_quadtree<Kernel, PointIterator>;
        using LMTSkeleton = mwt::LMTSkeleton<Traits, Tree>;

        CGAL::FPU_set_cw(CGAL_FE_UPWARD);
        Points points =
            mwt::generate_random_integral<Kernel>({{-grid_dim, -grid_dim}, {grid_dim, grid_dim}}, num_points, seed);
        Tree tree(points.begin(), points.end());
        mwt::Diamond_edge_filter<Tree, true> diamond_filter(&tree);
        diamond_filter.compute_remaining_edges();
        LMTSkeleton skeleton(&tree, std::move(diamond_filter.get_edges()));
        CHECK(diamond_filter.get_edges().size() == 0);
        skeleton.p_validate_initial_invariants();
        skeleton.lmt_main_loop();
        skeleton.p_validate_twin_relationship();
        skeleton.p_validate_candidates_and_skeleton();
        skeleton.advanced_lmt_loop();
        skeleton.p_validate_twin_relationship();

        /**
         * Check that there still is a triangulation
         * with the forced edges as constraints.
         */
        CGAL::Constrained_Delaunay_triangulation_2<Kernel, CGAL::Default, CGAL::Exact_predicates_tag>
            constrained_triangulation;
        constrained_triangulation.insert(points.begin(), points.end());
        const auto &forced_edges = skeleton.skeleton();
        std::vector<std::pair<Point, Point>> forced_pairs;
        for(const auto *edge : forced_edges) {
            CHECK((edge->status == LMTStatus::CH || edge->status == LMTStatus::Certain));
            forced_pairs.emplace_back(edge->source(), edge->target());
        }
        constrained_triangulation.insert_constraints(forced_pairs.begin(), forced_pairs.end());
        CHECK(constrained_triangulation.is_valid());
    }

    template<typename Kernel>
    static void random_lmt_parallel_vs_sequential(std::size_t grid_dim, std::size_t num_points, std::size_t seed) {
        using Traits = mwt::Mwt_traits_2<Kernel>;
        using Point = typename Traits::Point_2;
        using Points = std::vector<Point>;
        using PointIterator = typename Points::iterator;
        using Tree = mwt::Point_quadtree<Kernel, PointIterator>;
        using LMTSkeleton = mwt::LMTSkeleton<Traits, Tree>;
        using Halfedge = typename LMTSkeleton::Halfedge;

        CGAL::FPU_set_cw(CGAL_FE_UPWARD);
        Points points =
            mwt::generate_random_integral<Kernel>({{-grid_dim, -grid_dim}, {grid_dim, grid_dim}}, num_points, seed);
        Tree tree(points.begin(), points.end());

        std::vector<Halfedge *> parallel_result;
        std::vector<Halfedge *> sequential_result;
        std::unique_ptr<LMTSkeleton> parallel_skeleton;
        std::unique_ptr<LMTSkeleton> sequential_skeleton;

        {
            /* parallel */
            mwt::Diamond_edge_filter<Tree, true> diamond_filter(&tree);
            diamond_filter.compute_remaining_edges();
            parallel_skeleton = std::make_unique<LMTSkeleton>(&tree, std::move(diamond_filter.get_edges()));
            parallel_skeleton->p_validate_initial_invariants();
            parallel_skeleton->lmt_main_loop(/*parallel=*/true);
            parallel_skeleton->p_validate_twin_relationship();
            parallel_skeleton->p_validate_candidates_and_skeleton();
            parallel_skeleton->advanced_lmt_loop();
            parallel_skeleton->p_validate_twin_relationship();
            parallel_skeleton->p_validate_candidates_and_skeleton();
        }
        {
            /* sequential */
            mwt::Diamond_edge_filter<Tree, false> diamond_filter(&tree);
            diamond_filter.compute_remaining_edges();
            sequential_skeleton = std::make_unique<LMTSkeleton>(&tree, std::move(diamond_filter.get_edges()));
            sequential_skeleton->p_validate_initial_invariants();
            sequential_skeleton->lmt_main_loop(/*parallel=*/false);
            sequential_skeleton->p_validate_twin_relationship();
            sequential_skeleton->p_validate_candidates_and_skeleton();
            sequential_skeleton->advanced_lmt_loop();
            sequential_skeleton->p_validate_twin_relationship();
            sequential_skeleton->p_validate_candidates_and_skeleton();
        }
        parallel_result = parallel_skeleton->skeleton();
        sequential_result = sequential_skeleton->skeleton();
        CHECK(parallel_result.size() == sequential_result.size());
        std::sort(parallel_result.begin(), parallel_result.end());
        std::sort(sequential_result.begin(), sequential_result.end());
        for(std::size_t i = 0, nres = parallel_result.size(); i < nres; ++i) {
            CHECK(parallel_result[i]->status == sequential_result[i]->status);
            CHECK(parallel_result[i]->source() == sequential_result[i]->source());
            CHECK(parallel_result[i]->target() == sequential_result[i]->target());
        }
    }

    template<typename Skeleton> static void validate_initial_invariants(Skeleton &skeleton) {
        skeleton.p_validate_initial_invariants();
    }
};

} // namespace test

} // namespace mwt

TEST_CASE_TEMPLATE("LMT Skeleton Initialization", Kernel, CGAL::Exact_predicates_inexact_constructions_kernel,
                   CGAL::Exact_predicates_exact_constructions_kernel) {
    CGAL::Protect_FPU_rounding rounder;
    mwt::test::LMTSkeletonTester::random_initialization_test<Kernel>(10000, 40000, 1338);
    mwt::test::LMTSkeletonTester::random_initialization_test<Kernel>(10000, 40000, 1339);
}

TEST_CASE_TEMPLATE("Random LMT Skeleton", Kernel, CGAL::Exact_predicates_inexact_constructions_kernel,
                   CGAL::Exact_predicates_exact_constructions_kernel) {
    CGAL::Protect_FPU_rounding rounder;
    mwt::test::LMTSkeletonTester::random_lmt_test<Kernel>(10000, 40000, 13337);
}

TEST_CASE_TEMPLATE("Random LMT Skeleton - parallel vs. sequential", Kernel,
                   CGAL::Exact_predicates_inexact_constructions_kernel,
                   CGAL::Exact_predicates_exact_constructions_kernel) {
    CGAL::Protect_FPU_rounding rounder;
    mwt::test::LMTSkeletonTester::random_lmt_parallel_vs_sequential<Kernel>(10000, 10000, 155);
}

TEST_CASE("LMT skeleton regression - 16 point reduced instance") {
    using Kernel = CGAL::Epick;
    using Point = Kernel::Point_2;
    using Traits = mwt::Mwt_traits_2<Kernel>;
    using Points = std::vector<Point>;
    using PointIterator = Points::iterator;
    using Tree = mwt::Point_quadtree<Traits, PointIterator>;
    using Filter = mwt::Diamond_edge_filter<Tree, true>;
    using Skeleton = mwt::LMTSkeleton<Traits, Tree>;
    using Halfedge = Skeleton::Halfedge;

    Points points{{1093684.7861247056, 141422460.98395497}, {231745774.85289013, 55694386.866407484},
                  {231745774.85288987, 55694386.86640749},  {101736053.88626969, 385854957.7550138},
                  {226847995.81920564, 19180617.15081631},  {226847995.81920564, 19180617.150816284},
                  {224574292.12692818, 887386249.6478022},  {79279504.69066752, 537771864.0459071},
                  {59613341.37232336, 508445775.9244265},   {238677687.76589078, 519486350.9613699},
                  {413136606.88896614, 878380637.4216712},  {775435646.4357706, 873741531.0111102},
                  {424559431.10804564, 771978212.8797358},  {529402951.87850744, 584151948.2797536},
                  {533427355.5398562, 201913953.64625964},  {421192840.84122115, 324686862.71105677}};

    CGAL::Protect_FPU_rounding rounder;
    CHECK(points.size() == 16);
    std::sort(points.begin(), points.end());
    points.erase(std::unique(points.begin(), points.end()), points.end());
    CHECK(points.size() == 16);
    Tree tree(points.begin(), points.end());
    Filter diamond_filter(&tree);
    diamond_filter.compute_remaining_edges();
    CHECK(diamond_filter.get_edges().size() == 83);
    Skeleton skeleton(&tree, std::move(diamond_filter.get_edges()));
    mwt::test::LMTSkeletonTester::validate_initial_invariants(skeleton);
    CHECK(skeleton.num_convex_hull_edges() == 6);
    auto find_edge = [&](const Point &src, const Point &tgt) -> Halfedge * {
        for(auto &edge : skeleton.get_all_halfedges()) {
            if(edge.source() == src && edge.target() == tgt) {
                return &edge;
            }
        }
        REQUIRE(false);
        return nullptr;
    };

    Halfedge *last_edge =
        find_edge({231745774.8528898656, 55694386.86640749127}, {231745774.8528901339, 55694386.86640748382});
    CHECK(last_edge != nullptr);
    CHECK(last_edge->status == mwt::LMTStatus::Possible);
    std::vector<Point> targets_around_src;
    auto e = last_edge;
    do {
        targets_around_src.push_back(e->target());
        CHECK(e->status == mwt::LMTStatus::Possible);
        CHECK(e->source() == last_edge->source());
        e = e->next;
    } while(e != last_edge);
    std::vector<Point> expected_targets = {
        {231745774.85289013, 55694386.866407484}, {533427355.5398562, 201913953.64625964},
        {421192840.84122115, 324686862.71105677}, {529402951.87850744, 584151948.27975357},
        {238677687.76589078, 519486350.96136987}, {79279504.690667525, 537771864.04590714},
        {59613341.372323357, 508445775.9244265},  {101736053.88626969, 385854957.75501382},
        {1093684.7861247056, 141422460.98395497}, {226847995.81920564, 19180617.15081631},
        {226847995.81920564, 19180617.150816284}};
    CHECK(targets_around_src == expected_targets);

    Halfedge *first_deleted =
        find_edge({1093684.7861247056, 141422460.98395497}, {533427355.5398562, 201913953.64625964});
    Halfedge *first_deleted_twn = first_deleted->twin();
    CHECK(first_deleted->status == mwt::LMTStatus::Possible);
    CHECK(first_deleted_twn->status == mwt::LMTStatus::Possible);
    first_deleted_twn->reset();
    CHECK(first_deleted_twn->next_triangle());
    CHECK(first_deleted_twn->i->target() == first_deleted_twn->j->target());
    CHECK(first_deleted_twn->i->target() == Point(231745774.85288987, 55694386.866407491));
    CHECK(first_deleted_twn->next_triangle());
    CHECK(first_deleted_twn->i->target() == Point(231745774.85289013, 55694386.866407484));
    CHECK(first_deleted_twn->i->target() == first_deleted_twn->j->target());
    CHECK(!first_deleted_twn->next_triangle());
    first_deleted_twn->reset();
    first_deleted->reset();
    CHECK(first_deleted->next_triangle());
    CHECK(first_deleted->i->target() == first_deleted->j->target());
    CHECK(first_deleted->i->target() == Point(421192840.84122115, 324686862.71105677));
    CHECK(first_deleted->next_triangle());
    CHECK(first_deleted->i->target() == first_deleted->j->target());
    CHECK(first_deleted->i->target() == Point(238677687.76589078, 519486350.96136987));
    first_deleted->reset();
    skeleton.lmt_main_loop(/*parallel=*/false);
}
