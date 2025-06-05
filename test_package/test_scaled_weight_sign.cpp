#include <CGAL_MWT/Scaled_weight_sign.h>
#include <CGAL_MWT/Mwt_traits.h>
#include <CGAL_MWT/LMT_skeleton.h>
#include <CGAL_MWT/Static_quadtree.h>
#include <CGAL_MWT/Rational_or_int.h>
#include <CGAL/Exact_rational.h>
#include <CGAL/FPU.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <doctest/doctest.h>


TEST_CASE("problematic weight sign sequence #1") {
    using namespace mwt;
    using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
    using Traits = Mwt_traits_2<Kernel>;
    using Point = Traits::Point_2;
    using Points = std::vector<Point>;
    using PointIterator = Points::iterator;
    using Tree = Point_quadtree<Traits, PointIterator>;
    using Skeleton = LMTSkeleton<Traits, Tree>;
    using Halfedge = Skeleton::Halfedge;

    std::vector<Point> points{
        {212500773.17867965, 572648488.0832958},
        {249743464.13797587, 581217155.4036571},
        {212500773.17867976, 572648488.0832964},
        {249743464.1379756, 581217155.4036584},
        {212500773.1786797, 572648488.0832971}
    };
    auto p = [&] (int i) -> Point* { return &points[i]; };
    // setup the following edges:
    // 0 -> 1, 3
    // 1 -> 2, 4
    // 2 -> 3
    std::vector<Halfedge> halfedges{10, Halfedge{}};
    halfedges[0].tar = p(1);
    halfedges[0].twn = &halfedges[1];
    halfedges[1].tar = p(0);
    halfedges[1].twn = &halfedges[0];
    halfedges[2].tar = p(3);
    halfedges[2].twn = &halfedges[3];
    halfedges[3].tar = p(0);
    halfedges[3].twn = &halfedges[2];
    halfedges[4].tar = p(2);
    halfedges[4].twn = &halfedges[5];
    halfedges[5].tar = p(1);
    halfedges[5].twn = &halfedges[4];
    halfedges[6].tar = p(4);
    halfedges[6].twn = &halfedges[7];
    halfedges[7].tar = p(1);
    halfedges[7].twn = &halfedges[6];
    halfedges[8].tar = p(3);
    halfedges[8].twn = &halfedges[9];
    halfedges[9].tar = p(2);
    halfedges[9].twn = &halfedges[8];
        /*
        he[0] 1 
        he[4] 1
        he[2] -1
        he[9] -1
        he[0] -1
        he[3] -1
        he[7] 1
        he[8] -1
        he[5] 1
        he[6] 1*/

    using Scaled_weight_sign = mwt::Scaled_weight_sign<Halfedge, CGAL::Exact_rational, CGAL::Lazy_exact_nt<CGAL::Exact_rational>>;
    Scaled_weight_sign sign;
    sign.clear();

    // these should cancel
    sign.add_halfedge_with_coefficient(1, &halfedges[0]);
    sign.add_halfedge_with_coefficient(-1, &halfedges[0]);

    // add up to -2
    sign.add_halfedge_with_coefficient(-1, &halfedges[2]);
    sign.add_halfedge_with_coefficient(-1, &halfedges[3]);

    // add up to 2
    sign.add_halfedge_with_coefficient(1, &halfedges[4]);
    sign.add_halfedge_with_coefficient(1, &halfedges[5]);
    
    // add up to -2
    sign.add_halfedge_with_coefficient(-1, &halfedges[9]);
    sign.add_halfedge_with_coefficient(-1, &halfedges[8]);
    
    // add up to 2
    sign.add_halfedge_with_coefficient(1, &halfedges[7]);
    sign.add_halfedge_with_coefficient(1, &halfedges[6]);
    CHECK(sign.compute_sign() == CGAL::NEGATIVE);
}
