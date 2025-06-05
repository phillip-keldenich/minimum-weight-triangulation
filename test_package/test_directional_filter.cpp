#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/FPU.h>
#include <CGAL_MWT/Directional_filter.h>
#include <CGAL_MWT/Mwt_traits.h>
#include <doctest/doctest.h>

using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using Traits = mwt::Mwt_traits_2<Kernel>;
using Point = Traits::Point_2;
using Rect = Traits::Iso_rectangle_2;

namespace mwt {

namespace test {

template<typename Kernel> class Test_exact_diamond_filter {
  public:
    void test1() {
        // for point coordinates, see test_package/test_directional_filter_test1.ipe
        using Point = typename Kernel::Point_2;
        Exact_diamond_filter<Kernel> filter;
        filter.reset(Point{256, 704});

        std::vector<Point> points{{272, 768}, {308.711, 743.668},       {336, 688}, {352, 736}, {208, 816}, {208, 576},
                                  {240, 576}, {336, 687.9999999999999}, {352, 624}, {288, 576}};
        auto pi = points.begin();

        CHECK(CGAL::right_turn(Point{256, 704}, pi[0], pi[1]));
        filter.insert_sector(pi + 0, pi + 1);
        CHECK(filter.m_dead_sectors.size() == 3);
        CHECK(filter.m_dead_sectors[1].low_point == &points[1]);
        CHECK(filter.m_dead_sectors[1].high_point == &points[0]);

        CHECK(CGAL::right_turn(Point{256, 704}, pi[1], pi[2]));
        filter.insert_sector(pi + 1, pi + 2);
        CHECK(filter.m_dead_sectors.size() == 4);
        CHECK(filter.m_dead_sectors[1].low_point == nullptr);
        CHECK(filter.m_dead_sectors[1].high_point == nullptr);
        CHECK(filter.m_dead_sectors[2].low_point == &points[1]);
        CHECK(filter.m_dead_sectors[2].high_point == &points[0]);

        CHECK(CGAL::right_turn(Point{256, 704}, pi[1], pi[3]));
        filter.insert_sector(pi + 1, pi + 3);
        CHECK(filter.m_dead_sectors.size() == 3);
        CHECK(filter.m_dead_sectors[1].low_point == &points[3]);
        CHECK(filter.m_dead_sectors[1].high_point == &points[0]);
        CHECK(!filter(Point{368, 688}));

        CHECK(CGAL::right_turn(Point{256, 704}, pi[3], pi[2]));
        filter.insert_sector(pi + 3, pi + 2);
        CHECK(filter.m_dead_sectors.size() == 3);
        CHECK(filter.m_dead_sectors[1].low_point == &points[2]);
        CHECK(filter.m_dead_sectors[1].high_point == &points[0]);
        CHECK(filter(Point{368, 688}));
        CHECK(!filter(Point{336, 687.9999999999999}));

        CHECK(CGAL::right_turn(Point{256, 704}, pi[4], pi[0]));
        filter.insert_sector(pi + 4, pi + 0);
        CHECK(filter.m_dead_sectors.size() == 2);
        CHECK(filter.m_dead_sectors[1].low_point == &points[2]);
        CHECK(filter.m_dead_sectors[1].high_point == &points[4]);
        CHECK(filter.m_dead_sectors[1].pseudoangle_ub > 5000.0);
        CHECK(filter.m_dead_sectors[1].pseudoangle_lb < 0);
        CHECK(!filter(Point{224, 560}));

        CHECK(CGAL::right_turn(Point{256, 704}, pi[6], pi[5]));
        filter.insert_sector(pi + 6, pi + 5);
        CHECK(filter.m_dead_sectors.size() == 2);
        CHECK(filter.m_dead_sectors[0].low_point == &points[5]);
        CHECK(filter.m_dead_sectors[0].high_point == &points[6]);
        CHECK(filter(Point{224, 560}));

        CHECK(CGAL::right_turn(Point{256, 704}, pi[2], pi[7]));
        filter.insert_sector(pi + 2, pi + 7);
        CHECK(filter.m_dead_sectors.size() == 2);
        CHECK(filter.m_dead_sectors[1].low_point == &points[7]);
        CHECK(filter.m_dead_sectors[1].high_point == &points[4]);

        filter.insert_sector(pi + 4, pi + 8);
        CHECK(filter.m_dead_sectors.size() == 2);
        CHECK(filter.m_dead_sectors[1].low_point == &points[7]);
        CHECK(filter.m_dead_sectors[1].high_point == &points[4]);
        CHECK(filter.m_dead_sectors[0].low_point == &points[5]);
        CHECK(filter.m_dead_sectors[0].high_point == &points[6]);

        filter.insert_sector(pi + 9, pi + 6);
        CHECK(filter.m_dead_sectors.size() == 2);
        CHECK(filter.m_dead_sectors[0].low_point == &points[5]);
        CHECK(filter.m_dead_sectors[0].high_point == &points[9]);
        CHECK(filter.m_dead_sectors[0].pseudoangle_ub < 0);
        CHECK(filter.m_dead_sectors[0].pseudoangle_lb < -5000.0);
        CHECK(CGAL::right_turn(Point{256, 704}, pi[8], pi[9]));
        filter.insert_sector(pi + 8, pi + 9);
        CHECK(filter.m_dead_sectors.size() == 2);
        CHECK(filter.m_dead_sectors[0].low_point == &points[5]);
        CHECK(filter.m_dead_sectors[0].high_point == &points[8]);
        CHECK(filter.m_dead_sectors[0].pseudoangle_ub < 0);
        CHECK(filter.m_dead_sectors[0].pseudoangle_lb < -5000.0);
        CHECK(!filter());
        filter.insert_sector(pi + 7, pi + 8);
        CHECK(filter.m_dead_sectors.size() == 1);
        CHECK(filter());
        CHECK(filter.m_dead_sectors[0].pseudoangle_lb < -5000.0);
        CHECK(filter.m_dead_sectors[0].pseudoangle_ub > 5000.0);
        CHECK(filter.m_dead_sectors[0].low_point == &points[5]);
        CHECK(filter.m_dead_sectors[0].high_point == &points[4]);
    }

    void test_issue_priv3() {
        using Point = typename Kernel::Point_2;
        Point query{-25441084, -31235852};
        Point points[4] = {
            {-25464824, -31397672}, {-25513251, -31035362}, {-25485864, -31002255}, {-25449630, -31100149}};
        Exact_diamond_filter<Kernel> filter;
        filter.reset(query);
        filter.insert_sector(-2, -1.25, nullptr, nullptr);
        filter.insert_sector(-1.127937055399870658, 2, &points[0], &points[1]);
        CHECK(filter.size() == 2);
        CHECK(filter.contains_pseudo_angle_interval(-1.127937055399870658, 2));

        filter.insert_sector(1.160860990670924631, 1.264680532683921665, &points[2], &points[1]);
        CHECK(filter.size() == 2);
        CHECK(filter.contains_pseudo_angle_interval(-1.127937055399870658, 2));

        filter.insert_sector(1.059244778126711362, 1.160860990670924631, &points[3], &points[2]);
        CHECK(filter.size() == 2);
        CHECK(filter.contains_pseudo_angle_interval(-1.127937055399870658, 2));
    }

    void test_issue_priv4() {
        using Point = typename Kernel::Point_2;
        using Rect = typename Kernel::Iso_rectangle_2;

        std::vector<std::pair<double, double>> sectors[] = {
            {{-2.0, -1.25}, {1.25, 2.0}}, {{-2.0, -1.0}, {1.0, 2.0}}, {{1.25, 1.75}}, {{-2.0, 0.0}, {1.9999, 2.0}}};

        CGAL::FPU_set_cw(CGAL_FE_UPWARD);
        auto angle_interval = mwt::Compute_pseudo_angle<CGAL::Interval_nt_advanced>();

        for(const auto &secs : sectors) {
            Exact_diamond_filter<Kernel> filter;
            filter.reset(Point{0, 0});
            for(const auto &sec : secs) {
                filter.insert_sector(sec.first, sec.second, nullptr, nullptr);
            }
            for(int i = 0; i < 3600; ++i) {
                double angle = (2 * 3.141592653589793) * i / 3600.0;
                double cx = std::cos(angle);
                double cy = std::sin(angle);
                double xmin = cx - 0.375, xmax = cx + 0.375;
                double ymin = cy - 0.375, ymax = cy + 0.375;
                Rect rect{xmin, ymin, xmax, ymax};
                if(filter(rect)) {
                    Point points[4] = {{xmin, ymin}, {xmax, ymin}, {xmax, ymax}, {xmin, ymax}};
                    double min_a = 100.0, max_a = -100.0;
                    bool below_break = false, above_break = false;
                    for(const auto &p : points) {
                        auto a = angle_interval(p.x(), p.y());
                        if(a.inf() < min_a)
                            min_a = a.inf();
                        if(a.sup() > max_a)
                            max_a = a.sup();

                        if(a.inf() < -1.25)
                            below_break = true;
                        if(a.sup() > 1.25)
                            above_break = true;
                        CHECK(filter.contains_pseudo_angle_interval(a.inf(), a.sup()));
                        if(!below_break || !above_break) {
                            CHECK(filter.contains_pseudo_angle_interval(min_a, max_a));
                        }
                    }
                }
            }
        }
    }

    void test_issue_priv5() {
        using Point = typename Kernel::Point_2;
        using Rect = typename Kernel::Iso_rectangle_2;

        Exact_diamond_filter<Kernel> filter;
        filter.reset(Point{448, 320});

        CGAL::FPU_set_cw(CGAL_FE_UPWARD);
        auto angle_interval = mwt::Compute_pseudo_angle<CGAL::Interval_nt_advanced>();

        {
            Point p1{384, 400}, p2{416, 432};
            Rect r1{p1, p2};
            CHECK(!filter(r1));

            Point p3{368, 383}, p4{400, 415};
            Rect r2{p3, p4};
            CHECK(filter(r2));

            Point p5{397, 368}, p6{429, 400};
            Rect r3{p5, p6};
            CHECK(!filter(r3));
        }
    }
};

} // namespace test

} // namespace mwt

TEST_CASE_TEMPLATE("Exact_diamond_filter", Kernel, CGAL::Exact_predicates_exact_constructions_kernel,
                   CGAL::Exact_predicates_inexact_constructions_kernel) {
    CGAL::Protect_FPU_rounding rounder;
    mwt::test::Test_exact_diamond_filter<Kernel> tester;
    tester.test1();
}

TEST_CASE_TEMPLATE("Exact_diamond_filter_issue_priv3", Kernel, CGAL::Exact_predicates_exact_constructions_kernel,
                   CGAL::Exact_predicates_inexact_constructions_kernel) {
    CGAL::Protect_FPU_rounding rounder;
    mwt::test::Test_exact_diamond_filter<Kernel> tester;
    tester.test_issue_priv3();
}

TEST_CASE("Exact_diamond_filter_issue_priv4") {
    CGAL::Protect_FPU_rounding rounder;
    mwt::test::Test_exact_diamond_filter<CGAL::Exact_predicates_inexact_constructions_kernel> tester;
    tester.test_issue_priv4();
}

TEST_CASE("Exact_diamond_filter_issue_priv5") {
    CGAL::Protect_FPU_rounding rounder;
    mwt::test::Test_exact_diamond_filter<CGAL::Exact_predicates_inexact_constructions_kernel> tester;
    tester.test_issue_priv5();
}

TEST_CASE("Exact_diamond_filter_issue_priv6") {
    CGAL::Protect_FPU_rounding rounder;
    auto angle_interval = mwt::Compute_pseudo_angle<CGAL::Interval_nt_advanced>();
    CGAL::Interval_nt_advanced x(1.0, 1.0);
    CGAL::Interval_nt_advanced y(-0.25, 0.25);
    auto aiv = angle_interval(x, y);
    CHECK(aiv.inf() < -0.19);
    CHECK(aiv.sup() > 0.19);
}
