#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL_MWT/Read_instance.h>
#include <doctest/doctest.h>

TEST_CASE_TEMPLATE("[ReadInstance] Integer pair with comments", Kernel, CGAL::Epick, CGAL::Epeck) {
    using Point_2 = typename Kernel::Point_2;
    std::string instance = R"(
        # this is a comment
        0 0
        # this is another comment
        1 1
        1 2
        2 1
        3 3
    )";
    std::istringstream buffer(instance);
    auto points = mwt::read_instance<Kernel>(buffer);
    CHECK(points == std::vector<Point_2>{{0, 0}, {1, 1}, {1, 2}, {2, 1}, {3, 3}});
}

TEST_CASE_TEMPLATE("[ReadInstance] Double pair with comments", Kernel, CGAL::Epick, CGAL::Epeck) {
    std::string instance = R"(
        # this is a comment
        # multiple instances of the same points are removed, points are sorted.
        0.0 -0.0
        -1.0 2.0
        -3.0 3.0
        2 1
        -3.0 3.0
        -3.0 3.0
        EOF
    )";
    std::istringstream buffer(instance);
    auto points = mwt::read_instance<Kernel>(buffer);
    CHECK(points == std::vector<typename Kernel::Point_2>{{-3, 3}, {-1, 2}, {0, 0}, {2, 1}});
}

TEST_CASE_TEMPLATE("[ReadInstance] Double triple", Kernel, CGAL::Epick, CGAL::Epeck) {
    std::string instance = R"(
        # Epeck reads decimal values differently from Epick!
        0 -0.0 1.0
        1 -2.0 1.0
        2 -3.0 3.0
        3 2 4.940656458412465442e-324
        4 4.940656458412465442e-324 4.940656458412465442e-324
    )";
    std::istringstream buffer(instance);
    auto points = mwt::read_instance<Kernel>(buffer);
    if constexpr(std::is_same_v<Kernel, CGAL::Epeck>) {
        typename Kernel::FT val;
        std::istringstream val_i("4.940656458412465442e-324");
        val_i >> val;
        CHECK(points == std::vector<typename Kernel::Point_2>{{-3, 3}, {-2, 1}, {0, 1}, {val, val}, {2, val}});
    } else {
        CHECK(points == std::vector<typename Kernel::Point_2>{{-3, 3},
                                                              {-2, 1},
                                                              {0, 1},
                                                              {4.940656458412465442e-324, 4.940656458412465442e-324},
                                                              {2, 4.940656458412465442e-324}});
    }
}

TEST_CASE_TEMPLATE("[ReadInstance] JSON instance", Kernel, CGAL::Epick, CGAL::Epeck) {
    std::string instance = R"(
        {
            "points": [
                [0, 0],
                [1, 1],
                [1, 2],
                [2, 1],
                [3, 3.9]
            ],
            "remark": "epeck reads decimal values as floats like epick!"
        }
    )";
    std::istringstream buffer(instance);
    auto points = mwt::read_instance<Kernel>(buffer);
    CHECK(points == std::vector<typename Kernel::Point_2>{{0, 0}, {1, 1}, {1, 2}, {2, 1}, {3, 3.9}});
}
