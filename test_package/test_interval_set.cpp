#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/FPU.h>
#include <CGAL_MWT/Interval_set.h>
#include <doctest/doctest.h>

using Epeck = CGAL::Exact_predicates_exact_constructions_kernel;
using EPECKFT = Epeck::FT;
using Epick = CGAL::Exact_predicates_inexact_constructions_kernel;
using EPICKFT = Epeck::FT;
using namespace mwt;

TEST_CASE_TEMPLATE("Interval Set", T, int, float, double, EPECKFT, EPICKFT) {
    CGAL::Protect_FPU_rounding rounder;
    Interval_set<T> s1(T(-3), T(10));
    CHECK(s1.begin() == s1.end());
    CHECK(s1.get_range_min() == T(-3));
    CHECK(s1.get_range_max() == T(10));
    CHECK(!s1.contains(T(1)));
    CHECK(!s1.contains(T(2)));
    CHECK(!s1.contains(T(3)));
    CHECK(s1.size() == 0);

    s1.insert(T(1), T(3));
    CHECK(!s1.completely_covered());
    CHECK(s1.contains(T(1)));
    CHECK(s1.contains(T(2)));
    CHECK(s1.contains(T(3)));
    CHECK(!s1.contains(T(4)));
    CHECK(!s1.contains(T(0)));
    CHECK(s1.contains(T(2), T(3)));
    CHECK(s1.size() == 1);

    s1.insert(T(0), T(2));
    CHECK(s1.contains(T(0), T(3)));
    CHECK(!s1.completely_covered());
    CHECK(!s1.contains(T(-2), T(2)));

    s1.insert(T(8), T(10));
    s1.insert(T(-3), T(2));
    CHECK(!s1.completely_covered());
    CHECK(s1.contains(T(-3), T(3)));
    CHECK(!s1.contains(T(-3), T(4)));
    CHECK(s1.contains(T(-2), T(3)));
    CHECK(!s1.contains(T(1), T(5)));
    CHECK(!s1.contains(T(3), T(10)));
    CHECK(s1.contains(T(9)));

    s1.insert(T(4), T(7));
    CHECK(!s1.completely_covered());
    s1.insert(T(3), T(8));
    CHECK(s1.completely_covered());
}

TEST_CASE_TEMPLATE("Circular Interval Set", T, int, float, double, EPECKFT, EPICKFT) {
    CGAL::Protect_FPU_rounding rounder;
    Circular_interval_set<T> s1(T(-3), T(10));
    CHECK(!s1.contains(T(1)));
    CHECK(!s1.contains(T(2)));
    CHECK(!s1.contains(T(3)));
    CHECK(s1.size() == 0);

    s1.insert(T(1), T(3));
    CHECK(!s1.completely_covered());
    CHECK(s1.contains(T(1)));
    CHECK(s1.contains(T(2)));
    CHECK(s1.contains(T(3)));
    CHECK(!s1.contains(T(4)));
    CHECK(!s1.contains(T(0)));
    CHECK(s1.contains(T(2), T(3)));
    CHECK(s1.size() == 1);

    s1.insert(T(0), T(2));
    CHECK(s1.contains(T(0), T(3)));
    CHECK(!s1.completely_covered());
    CHECK(!s1.contains(T(-2), T(2)));

    s1.insert(T(8), T(10));
    s1.insert(T(-3), T(2));
    CHECK(!s1.completely_covered());
    CHECK(s1.contains(T(-3), T(3)));
    CHECK(!s1.contains(T(-3), T(4)));
    CHECK(s1.contains(T(-2), T(3)));
    CHECK(!s1.contains(T(1), T(5)));
    CHECK(!s1.contains(T(3), T(10)));
    CHECK(s1.contains(T(9)));

    s1.insert(T(4), T(7));
    CHECK(!s1.completely_covered());
    s1.insert(T(10), T(8));
    CHECK(s1.completely_covered());
}
