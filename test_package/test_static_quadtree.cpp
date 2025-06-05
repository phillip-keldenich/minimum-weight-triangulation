#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/FPU.h>
#include <CGAL_MWT/Mwt_traits.h>
#include <CGAL_MWT/Static_quadtree.h>
#include <doctest/doctest.h>
#include <vector>

using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point = Kernel::Point_2;
using Traits = mwt::Mwt_traits_2<Kernel>;

TEST_CASE("Static_quadtree") {
    using Vec = std::vector<Point>;
    using Iter = std::vector<Point>::iterator;

    CGAL::Protect_FPU_rounding rounder;

    SUBCASE("Small") {
        std::vector<Point> points{{0, 0},  {1, 1},  {1, 2},  {1, 3},  {2, 2}, {2, 3},  {5, 9},   {7, 17}, {13, 13},
                                  {11, 1}, {7, 13}, {2, 22}, {8, 19}, {5, 8}, {9, 10}, {10, 10}, {3, 7},  {900, 100}};
        mwt::Static_quadtree<Traits, Iter> quadtree(points.begin(), points.end());
        quadtree.verify_tree();
    }
}
