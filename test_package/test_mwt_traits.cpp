#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/FPU.h>
#include <CGAL/Interval_arithmetic.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL_MWT/Mwt_traits.h>
#include <doctest/doctest.h>

TEST_CASE("MWT Traits with Epeck") {
    using Kernel = CGAL::Exact_predicates_exact_constructions_kernel;
    using Traits = mwt::Mwt_traits_2<Kernel>;
    using Point = Kernel::Point_2;
    using Segment = Kernel::Segment_2;
    using FT = Kernel::FT;

    static_assert(std::is_same_v<Point, Traits::Point_2>, "Traits do not define correct point type!");

    CGAL::Protect_FPU_rounding rounder;
    Traits trait_object;

    Point very_close{1.0000000000000002220446049250313080847263336181640625,
                     1.0000000000000002220446049250313080847263336181640625};
    Point very_far_coll{1.0e100, 1.0e100};
    Point very_far_below{1.0e100, 9.999999999999998e+99};
    Point very_far_above{1.0e100, 1.0000000000000002e+100};

    SUBCASE("MWT Traits: Less_angle_2 (around origin = (1,1))") {
        auto less_angle_2 = trait_object.angle_less_2_object(Point{1, 1});
        CHECK(less_angle_2(Point{2, 2}, Point{-2, 2}));
        CHECK(!less_angle_2(Point{2, 2}, Point{2, 1}));
        CHECK(less_angle_2(Point{2, 1}, Point{2, 2}));
        CHECK(!less_angle_2(very_close, very_far_coll));
        CHECK(!less_angle_2(very_far_coll, very_close));
        CHECK(less_angle_2(very_far_below, very_close));
        CHECK(less_angle_2(very_close, very_far_above));
        CHECK(less_angle_2(very_far_below, very_far_above));
    }

    SUBCASE("MWT Traits: Compute_angle (around origin = (1,1))") {
        auto compute_angle_2 = trait_object.compute_polar_angle_2_object(Point{1, 1});
        auto compare_angle_2 = trait_object.less_polar_angle_2_object(Point{1, 1});
        auto less_angle_2 = trait_object.angle_less_2_object(Point{1, 1});

        std::vector<Point> points{{2, 2}, {-2, 2}, {2, 1}, very_close, very_far_coll, very_far_below, very_far_above};
        std::vector<FT> angles;
        for(Point p : points) {
            angles.push_back(compute_angle_2(p));
        }
        for(std::size_t i = 0, n = points.size(); i < n; ++i) {
            for(std::size_t j = 0; j < i; ++j) {
                const FT &angle_i = angles[i];
                const FT &angle_j = angles[j];
                if(angle_i < angle_j) {
                    CHECK(less_angle_2(points[i], points[j]));
                    CHECK(compare_angle_2(points[i], points[j]));
                } else if(angle_j < angle_i) {
                    CHECK(less_angle_2(points[j], points[i]));
                    CHECK(compare_angle_2(points[j], points[i]));
                } else {
                    CHECK(!less_angle_2(points[j], points[i]));
                    CHECK(!compare_angle_2(points[j], points[i]));
                    CHECK(!less_angle_2(points[i], points[j]));
                    CHECK(!compare_angle_2(points[i], points[j]));
                }
            }
        }
    }

    SUBCASE("MWT Traits: Diamond test") {
        auto diamond_test = trait_object.diamond_test_2_object();
        Segment segment{Point{192, 576}, Point{288, 576}};
        auto mirror_vertical = [](const Point &p) {
            auto dy = p.y() - 576;
            return Point{p.x(), 576 - dy};
        };
        std::vector<Point> points_true{Point{224, 592}, Point{240, 592}, Point{256, 592},
                                       Point{240, 608}, Point{240, 615}, Point{207, 588}};
        std::vector<Point> points_false{Point{208, 624}, Point{192, 608}, Point{224, 672}, Point{240, 720},
                                        Point{272, 592}, Point{208, 608}, Point{224, 608}, Point{256, 608},
                                        Point{272, 608}, Point{224, 624}, Point{240, 624}, Point{256, 624}};
        for(Point l : points_true) {
            Point r = mirror_vertical(l);
            CHECK(diamond_test(segment, l, r));
        }
        for(Point l : points_false) {
            Point r = mirror_vertical(l);
            CHECK(!diamond_test(segment, l, r));
        }
    }

    auto is_in_sector = [](const auto &s, const auto &t, const auto &l, const auto &r, const auto &sector) -> bool {
        if(!sector)
            return false;
        const auto &sec = *sector;
        auto distance = CGAL::max(CGAL::squared_distance(s, l), CGAL::squared_distance(s, r));
        auto activation_distance = distance * Traits::Construct_dead_sector_2::squared_distance_coefficient();
        Traits trait_object;
        auto polar = trait_object.compute_polar_angle_2_object(s);
        auto pol_t = polar(t);
        if(pol_t < sec.first)
            return false;
        if(pol_t > sec.second)
            return false;
        if(CGAL::squared_distance(s, t) < activation_distance)
            return false;
        auto diamond_test = trait_object.diamond_test_2_object();
        CHECK(diamond_test(Segment{s, t}, l, r));
        return true;
    };

    SUBCASE("MWT Traits: Dead sector computation (case 1: no sector)") {
        auto dead_sector = trait_object.construct_dead_sector_2_object();
        Point s{384, 400}, l{385.321, 524.957}, r{469.869, 416.683};
        auto sector = dead_sector(s, l, r);
        CHECK(!sector);
    }

    SUBCASE("MWT Traits: Dead sector computation (case 2: delimited by points)") {
        auto dead_sector = trait_object.construct_dead_sector_2_object();
        Point s{384, 448}, l{428.8, 504}, r{484.083, 469.895};
        auto sector = dead_sector(s, l, r);
        CHECK(!!sector);
        std::vector<Point> positive_points{{512, 576}, {528, 544}, {544, 528}, {576, 496}, {528, 624}};
        std::vector<Point> negative_points{{480, 544}, {512, 512}, {480, 592}, {560, 464}, {480, 592}};
        for(const auto &t : positive_points) {
            CHECK(is_in_sector(s, t, l, r, sector));
        }
        for(const auto &t : negative_points) {
            CHECK(!is_in_sector(s, t, l, r, sector));
        }
    }

    SUBCASE("MWT Traits: Dead sector computation (case 3: not delimited by points)") {
        auto dead_sector = trait_object.construct_dead_sector_2_object();
        Point s{176, 560}, l{192, 640}, r{256, 576};
        auto sector = dead_sector(s, l, r);
        CHECK(!!sector);
        std::vector<Point> positive_points{{320, 704}, {336, 704}, {328, 689}};
        std::vector<Point> negative_points{{336, 688}, {320, 672}, {257, 653}, {288, 720}, {288, 730}};
        for(const auto &t : positive_points) {
            CHECK(is_in_sector(s, t, l, r, sector));
        }
        for(const auto &t : negative_points) {
            CHECK(!is_in_sector(s, t, l, r, sector));
        }
    }
}
