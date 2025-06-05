#include <CGAL/FPU.h>
#include <CGAL_MWT/Diamond_filter.h>
#include <CGAL_MWT/Exact_simple_face_triangulator.h>
#include <CGAL_MWT/Face_collector.h>
#include <CGAL_MWT/LMT_skeleton.h>
#include <CGAL_MWT/Mwt_traits.h>
#include <CGAL_MWT/Validate.h>
#include <CGAL_MWT/output.h>
#include <doctest/doctest.h>

using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using Traits = mwt::Mwt_traits_2<Kernel>;
using Point = Traits::Point_2;
using Points = std::vector<Point>;
using PointIterator = Points::iterator;
using Tree = mwt::Point_quadtree<Traits, PointIterator>;
using Filter = mwt::Diamond_edge_filter<Tree, true>;
using Skeleton = mwt::LMTSkeleton<Traits, Tree>;
using Halfedge = Skeleton::Halfedge;
using FaceCollector = mwt::Face_collector<Skeleton>;
using Verifier = mwt::MWTVerifier<Skeleton>;


extern std::vector<Halfedge> halfedges_from_graph(const std::vector<Point> &points,
                                                  const std::vector<std::pair<std::size_t, std::size_t>> &skeleton_edges,
                                                  const std::vector<std::pair<std::size_t, std::size_t>> &possible_edges);


TEST_CASE("Test triangulation of a rectangle via DP") {
    /**
     * Due to the explicit check of length-3 edges,
     * DP with intervals can actually fully triangulate
     * empty rectangles.
     */
    CGAL::Protect_FPU_rounding round_protect;
    Points points{{0, 0}, {3, 0}, {3, 1}, {0, 1}};
    Tree tree(points.begin(), points.end());
    Filter filter(&tree);
    filter.compute_remaining_edges();
    Skeleton skeleton(&tree, std::move(filter.get_edges()));
    skeleton.lmt_main_loop();
    skeleton.advanced_lmt_loop();
    FaceCollector collector(skeleton);
    CHECK(collector.next());
    auto face = collector.get_current_face_copy();
    CHECK(!collector.next());
    mwt::Interval_simple_face_triangulator<Skeleton> triangulator;
    triangulator(std::move(face));
    CHECK(triangulator.is_unique());
    Verifier verifier{&skeleton};
    verifier.strong_check();
}

Points make_regular_ngon(std::size_t n) {
    Points points;
    double recn = (1.0 / n) * 2 * CGAL_PI;
    for(std::size_t i = 0; i < n; ++i) {
        points.push_back({std::cos(recn * i), std::sin(recn * i)});
    }
    return points;
}

TEST_CASE("Test triangulation of an (almost) regular hexagon via DP") {
    CGAL::Protect_FPU_rounding round_protect;
    Points points = make_regular_ngon(6);
    Tree tree(points.begin(), points.end());
    Filter filter(&tree);
    filter.compute_remaining_edges();
    Skeleton skeleton(&tree, std::move(filter.get_edges()));
    skeleton.lmt_main_loop();
    skeleton.advanced_lmt_loop();
    FaceCollector collector(skeleton);
    CHECK(collector.next());
    auto face = collector.get_current_face_copy();
    CHECK(!collector.next());
    mwt::Interval_simple_face_triangulator<Skeleton> triangulator;
    triangulator(std::move(face));
    CHECK(triangulator.is_unique());
    Verifier verifier{&skeleton};
    verifier.strong_check();
}

TEST_CASE("Test triangulation of an (almost) regular 33-gon via DP") {
    CGAL::Protect_FPU_rounding round_protect;
    Points points = make_regular_ngon(33);
    Tree tree(points.begin(), points.end());
    Filter filter(&tree);
    filter.compute_remaining_edges();
    Skeleton skeleton(&tree, std::move(filter.get_edges()));
    skeleton.lmt_main_loop();
    skeleton.advanced_lmt_loop();
    FaceCollector collector(skeleton);
    CHECK(collector.next());
    auto face = collector.get_current_face_copy();
    CHECK(!collector.next());
    mwt::Interval_simple_face_triangulator<Skeleton> triangulator;
    triangulator(std::move(face));
    CHECK(!triangulator.is_unique());
    fix_nonunique_interval_dp(triangulator);
    CHECK(triangulator.is_unique());
    Verifier verifier{&skeleton};
    verifier.strong_check();
}

TEST_CASE("Test triangulation of a regular grid") {
    CGAL::Protect_FPU_rounding round_protect;
    Points points;
    for(int i = 0; i < 20; ++i) {
        for(int j = 0; j < 20; ++j) {
            points.emplace_back(i, j);
        }
    }
    Tree tree(points.begin(), points.end());
    Filter filter(&tree);
    filter.compute_remaining_edges();
    Skeleton skeleton(&tree, std::move(filter.get_edges()));
    skeleton.lmt_main_loop();
    skeleton.advanced_lmt_loop();
    FaceCollector collector(skeleton);
    while(collector.next()) {
        auto face = collector.get_current_face_copy();
        mwt::Interval_simple_face_triangulator<Skeleton> triangulator;
        triangulator(std::move(face));
        CHECK(triangulator.is_unique());
    }
    Verifier verifier{&skeleton};
    verifier.strong_check();
}

TEST_CASE("Test triangulation of a slightly irregular grid") {
    /**
     * The solution here is unique (and already done by the LMT skeleton).
     */
    CGAL::Protect_FPU_rounding round_protect;
    Points points;
    for(int y = 0; y < 20; ++y) {
        for(int x = 0; x < 20; ++x) {
            double rx = x;
            if(y % 2) {
                rx = std::nextafter(rx, rx - 1);
            } else {
                rx = std::nextafter(rx, rx + 1);
            }
            points.emplace_back(rx, double(y));
        }
    }
    Tree tree(points.begin(), points.end());
    Filter filter(&tree);
    filter.compute_remaining_edges();
    Skeleton skeleton(&tree, std::move(filter.get_edges()));
    skeleton.lmt_main_loop();
    skeleton.advanced_lmt_loop();
    FaceCollector collector(skeleton);
    while(collector.next()) {
        CHECK(false);
        auto face = collector.get_current_face_copy();
        mwt::Interval_simple_face_triangulator<Skeleton> triangulator;
        triangulator(std::move(face));
        CHECK(triangulator.is_unique());
    }
    Verifier verifier{&skeleton};
    verifier.strong_check();
}

/**
 * A problematic case of a face from an instance
 * which, depending on which boundary edge is chosen in the DP,
 * uses a halfedge with a non-unique left-hand side triangulation.
 */
TEST_CASE("Test problematic face #1") {
    CGAL::Protect_FPU_rounding round_protect;
    Points points{
        {3497.6, 1028.2}, // 0
        {3516.6, 1072.6}, // 1  --  forbidden (close to optimal) edge 1 -- 3
        {3567.4, 1104.4}, // 2
        {3554.7, 1123.4}, // 3
        {3548.4, 1142.5}, // 4
        {3529.3, 1155.2}, // 5  --  correct edge 5 -- 2
        {3472.2, 1174.2}, // 6
        {3408.7, 1161.5}, // 7
        {3396.0, 1123.4}, // 8
        {3396.0, 1053.6}  // 9
    };
    Tree tree(points.begin(), points.end());
    Filter filter(&tree);
    filter.compute_remaining_edges();
    Skeleton skeleton(&tree, std::move(filter.get_edges()));
    skeleton.lmt_main_loop();
    skeleton.advanced_lmt_loop();
    FaceCollector collector(skeleton);
    CHECK(collector.next());
    auto face = collector.get_current_face_copy();
    CHECK(face.is_simple());
    CHECK(face.boundary.size() == 10);
    CHECK(face.inner_edges.size() == 30);
    CHECK(!collector.next());
    mwt::Interval_simple_face_triangulator<Skeleton> triangulator;
    triangulator(std::move(face));
    std::vector<Halfedge*> certain;
    bool found = false;
    for(Halfedge &h : skeleton.get_all_halfedges()) {
        if(h.status == mwt::LMTStatus::Certain || h.status == mwt::LMTStatus::CH) {
            certain.push_back(&h);
            if(h.source() == points[1] && h.target() == points[3] ||
               h.source() == points[3] && h.target() == points[1]) {
                CHECK(h.status != mwt::LMTStatus::Certain);
                CHECK(h.status != mwt::LMTStatus::CH);
            }
            if(h.source() == points[5] && h.target() == points[2] ||
               h.source() == points[2] && h.target() == points[5]) {
                found = true;
            }
        }
    }
    CHECK(triangulator.is_unique());
    CHECK(certain.size() == 38);
    CHECK(found);
}

/**
 * The same case as above, but with the points permuted
 * such as to enforce the usage of the non-unique edge.
 */
TEST_CASE("Test problematic face #1 - permuted") {
    CGAL::Protect_FPU_rounding round_protect;
    Points points{
        {3497.6, 1028.2},
        {3516.6, 1072.6},
        {3567.4, 1104.4},
        {3554.7, 1123.4},
        {3548.4, 1142.5},
        {3529.3, 1155.2},
        {3472.2, 1174.2},
        {3408.7, 1161.5},
        {3396.0, 1123.4},
        {3396.0, 1053.6}
    };
    std::vector<std::pair<std::size_t, std::size_t>> certain_edges{
        {0, 1}, {1, 2}, {2, 3}, {3, 4}
    };
    std::vector<std::pair<std::size_t, std::size_t>> possible_edges{
        {0, 5}, {0, 7}, {1, 3}, {1, 4}, {1, 5}, {1, 6}, {1, 7}, {1, 8}, {1, 9},
        {2, 5}, {3, 5}, {4, 6}, {5, 7}, {6, 8}, {7, 9}
    };
    std::vector<std::size_t> permutation{
        3, 9, 7, 8, 6, 1, 2, 5, 4, 0
    };
    Points permuted_points;
    std::vector<std::size_t> new_index_of(10, 0);
    for(std::size_t i : permutation) {
        permuted_points.push_back(points[i]);
        new_index_of[i] = permuted_points.size() - 1;
    }
    for(auto& [s, t] : certain_edges) {
        s = new_index_of[s];
        t = new_index_of[t];
    }
    for(auto& [s, t] : possible_edges) {
        s = new_index_of[s];
        t = new_index_of[t];
    }
    std::vector<Halfedge> halfedges = halfedges_from_graph(permuted_points, certain_edges, possible_edges);
    std::vector<Halfedge *> lmt_skeleton;
    for(std::size_t i = 0; i < halfedges.size(); ++i) {
        if(halfedges[i].status == mwt::LMTStatus::Certain || halfedges[i].status == mwt::LMTStatus::CH) {
            lmt_skeleton.push_back(&halfedges[i]);
        }
    }
    mwt::Face_collector<Skeleton> face_collector(lmt_skeleton.begin(), lmt_skeleton.end());
    CHECK(face_collector.next());
    auto f1 = face_collector.get_current_face_copy();
    CHECK(f1.is_simple());
    CHECK(f1.has_inner_edges());
    CHECK(!face_collector.next());
    mwt::Interval_simple_face_triangulator<Skeleton> triangulator;
    triangulator(std::move(f1));
    std::vector<Halfedge*> certain;
    if(!triangulator.is_unique()) {
        fix_nonunique_interval_dp(triangulator);
        CHECK(triangulator.is_unique());
    }
    bool found = false;
    for(Halfedge &h : halfedges) {
        if(h.status == mwt::LMTStatus::Certain || h.status == mwt::LMTStatus::CH) {
            certain.push_back(&h);
            if(h.source() == points[1] && h.target() == points[3] ||
               h.source() == points[3] && h.target() == points[1]) {
                CHECK(h.status != mwt::LMTStatus::Certain);
                CHECK(h.status != mwt::LMTStatus::CH);
            }
            if(h.source() == points[5] && h.target() == points[2] ||
               h.source() == points[2] && h.target() == points[5]) {
                found = true;
            }
        }
    }
    CHECK(found);
    CHECK(certain.size() == 38);
}
