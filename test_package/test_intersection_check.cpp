#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Kernel/global_functions.h>
#include <CGAL/Point_2.h>
#include <CGAL/Segment_2.h>
#include <CGAL/Sweep_line_2_algorithms.h>
#include <CGAL_MWT/Validate.h>
#include <doctest/doctest.h>

TEST_CASE_TEMPLATE("[IntersectionChecker] test CGAL behavior", Kernel, CGAL::Epeck, CGAL::Epick) {
    using Point = typename Kernel::Point_2;
    using Segment = typename Kernel::Segment_2;
    using FT = typename Kernel::FT;

    Point below_vert(0, -1);
    Point at_lower_vert(0, 0);
    Point in_vert(0, 0.5);
    Point at_upper_vert(0, 1);
    Point above_vert(0, 2);
    Segment vert(at_lower_vert, at_upper_vert);
    Segment vert2(at_upper_vert, at_lower_vert);
    Segment s3{Point{1, 1}, Point{2, 2}};
    Segment s4{Point{0, 0}, Point{1, 1}};

    auto comp = [](const Point &p, const Segment &s) { return CGAL::compare_y_at_x(p, s); };
    CHECK(comp(below_vert, vert) == CGAL::SMALLER);
    CHECK(comp(at_lower_vert, vert) == CGAL::EQUAL);
    CHECK(comp(in_vert, vert) == CGAL::EQUAL);
    CHECK(comp(at_upper_vert, vert) == CGAL::EQUAL);
    CHECK(comp(above_vert, vert) == CGAL::LARGER);
    CHECK(comp(below_vert, vert2) == CGAL::SMALLER);
    CHECK(comp(at_lower_vert, vert2) == CGAL::EQUAL);
    CHECK(comp(in_vert, vert2) == CGAL::EQUAL);
    CHECK(comp(at_upper_vert, vert2) == CGAL::EQUAL);
    CHECK(comp(above_vert, vert2) == CGAL::LARGER);
    CHECK(CGAL::do_intersect(vert, vert2));
    CHECK(CGAL::do_intersect(vert2, vert));
    CHECK(CGAL::do_intersect(s3, s4));
}

TEST_CASE_TEMPLATE("[IntersectionChecker] test intersection checker", Kernel, CGAL::Epeck, CGAL::Epick) {
    using Point = typename Kernel::Point_2;
    using Segment = typename Kernel::Segment_2;
    using EPoint = CGAL::Epeck::Point_2;
    using ESegment = CGAL::Epeck::Segment_2;
    using FT = typename Kernel::FT;

    std::vector<Point> points = {
        {112, 64},  // 0
        {144, 48},  // 1
        {144, 64},  // 2
        {144, 80},  // 3
        {144, 96},  // 4
        {144, 112}, // 5
        {176, 64},  // 6
        {176, 96},  // 7
        {176, 128}  // 8
    };
    std::vector<std::pair<std::size_t, std::size_t>> edges = {{0, 1}, {0, 2}, {0, 3}, {0, 4}, {0, 5}, {1, 2},
                                                              {2, 3}, {3, 4}, {4, 5}, {6, 7}, {7, 8}, {1, 6},
                                                              {2, 6}, {3, 6}, {3, 7}, {4, 7}, {4, 8}, {5, 8}};
    std::vector<Segment> segments;
    std::vector<ESegment> exact_segments;
    for(auto [i1, i2] : edges) {
        segments.push_back({points[i1], points[i2]});
        exact_segments.push_back({EPoint(points[i1].x(), points[i1].y()), EPoint(points[i2].x(), points[i2].y())});
    }
    CHECK(!CGAL::do_curves_intersect(exact_segments.begin(), exact_segments.end()));
    mwt::CrossingFreenessVerifier<Kernel> verifier(std::move(segments));
    CHECK(!verifier());
    for(std::size_t v1 = 0; v1 < 8; ++v1) {
        for(std::size_t v2 = v1 + 1; v2 < 9; ++v2) {
            std::vector<std::pair<std::size_t, std::size_t>> edge_cpy(edges);
            std::pair<std::size_t, std::size_t> edge(v1, v2);
            if(std::find(edge_cpy.begin(), edge_cpy.end(), edge) != edge_cpy.end())
                continue;
            edge_cpy.push_back(edge);
            segments.clear();
            exact_segments.clear();
            for(auto [i1, i2] : edge_cpy) {
                segments.push_back({points[i1], points[i2]});
                exact_segments.push_back(
                    {EPoint(points[i1].x(), points[i1].y()), EPoint(points[i2].x(), points[i2].y())});
            }
            CHECK(CGAL::do_curves_intersect(exact_segments.begin(), exact_segments.end()));
            mwt::CrossingFreenessVerifier<Kernel> verifier(std::move(segments));
            CHECK(verifier());
            auto segs = verifier.get_intersecting_segments();
            CHECK(CGAL::do_intersect(segs.first, segs.second));
        }
    }
}
