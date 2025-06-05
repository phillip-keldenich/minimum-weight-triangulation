#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/FPU.h>
#include <CGAL_MWT/Directional_filter.h>
#include <CGAL_MWT/Filtered_incremental_search.h>
#include <CGAL_MWT/Mwt_traits.h>
#include <doctest/doctest.h>
#include <random>

inline std::mt19937_64 &default_rng() {
    thread_local std::mt19937_64 rng{42};
    return rng;
}

template<typename Point, typename RNG = std::mt19937_64>
void generate_exactness_test_points(std::vector<Point> &points, const Point &query, int group_size,
                                    RNG &rng = default_rng()) {
    using FT = std::remove_cv_t<std::remove_reference_t<decltype(query.x())>>;
    FT large = 3.094850098213451e+16;
    FT tiny = 2.938735877055719e-39;

    std::uniform_int_distribution<int> dist{0, 1};
    auto randomize_sign = [&](auto d) {
        int value = 1 - 2 * dist(rng);
        return d * value;
    };

    points.push_back(query);
    for(int i = 0; i < group_size; ++i) {
        points.emplace_back(large, randomize_sign(i * tiny));
        points.emplace_back(-large, randomize_sign(i * tiny));
        points.emplace_back(randomize_sign(i * tiny), large);
        points.emplace_back(randomize_sign(i * tiny), -large);
    }

    std::sort(points.begin(), points.end(), [](const auto &lhs, const auto &rhs) {
        return lhs.x() < rhs.x() || (lhs.x() == rhs.x() && lhs.y() < rhs.y());
    });
    points.erase(std::unique(points.begin(), points.end()), points.end());
    std::shuffle(points.begin(), points.end(), rng);
}

TEST_CASE_TEMPLATE("Filtered incremental search (new exact version), counting filter, exactness test", Kernel,
                   CGAL::Exact_predicates_exact_constructions_kernel,
                   CGAL::Exact_predicates_inexact_constructions_kernel) {
    using Point = typename Kernel::Point_2;
    using Rect = typename Kernel::Iso_rectangle_2;

    CGAL::Protect_FPU_rounding rounder;
    constexpr std::size_t NEIGHBORS_PER_POINT = 1000;

    struct CountingFilter {
        std::size_t point_count{0};

        void reset(const Point & /*query*/) { point_count = 0; }

        bool operator()() const { return point_count >= NEIGHBORS_PER_POINT; }
        bool operator()(const Point & /*point*/) const { return point_count >= NEIGHBORS_PER_POINT; }
        bool operator()(const Rect & /*rect*/) const { return point_count >= NEIGHBORS_PER_POINT; }
    };

    CGAL::FPU_set_cw(CGAL_FE_UPWARD);
    std::vector<Point> points;
    std::size_t group_size = 400;
    generate_exactness_test_points(points, Point{0, 0}, group_size);

    CHECK(points.size() == 1 + 4 * 400);

    auto distance_rank = [&](const Point &point) {
        const double tiny = 2.938735877055719e-39;
        const double large = 3.094850098213451e+16;
        auto ax = CGAL::abs(point.x());
        auto ay = CGAL::abs(point.y());
        if(ax < ay) {
            double ddist = CGAL::to_double(ax / tiny);
            int idist = static_cast<int>(ddist);
            CHECK(ddist == idist);
            return idist;
        } else {
            double ddist = CGAL::to_double(ay / tiny);
            int idist = static_cast<int>(ddist);
            CHECK(ddist == idist);
            return idist;
        }
    };

    Point query{0, 0};
    for(std::size_t i = 0; i + 1 < points.size(); ++i) {
        Point cur = points[i];
        Point next = points[i + 1];
        if(cur == query || next == query)
            continue;
        int dc = distance_rank(cur);
        int dn = distance_rank(next);
        if(dc < dn) {
            CHECK(CGAL::has_larger_distance_to_point(query, next, cur));
        } else if(dn < dc) {
            CHECK(CGAL::has_larger_distance_to_point(query, cur, next));
        } else {
            CHECK(!CGAL::has_larger_distance_to_point(query, cur, next));
            CHECK(!CGAL::has_larger_distance_to_point(query, next, cur));
        }
    }

    using PointIterator = typename decltype(points)::iterator;
    using Tree = mwt::Point_quadtree<Kernel, PointIterator>;
    using Search = typename Tree::template Incremental_search<CountingFilter>;
    Tree tree(points.begin(), points.end());
    tree.verify();
    Search search(&tree, Point{0, 0}, CountingFilter{});
    CountingFilter &filter = search.filter();
    CHECK(search.start_search(Point(0, -3.094850098213451e+16)));
    CHECK(*search.current().handle() == Point(0, -3.094850098213451e+16));
    CHECK(search.start_search(Point{0, 0}));
    CHECK(*search.current().handle() == Point{0, 0});

    int prev_distance_rank = -1;
    int at_that_rank = 4;
    for(;;) {
        if(!filter()) {
            CHECK(search.next());
        } else {
            CHECK(!search.next());
            break;
        }
        Point current = *search.current().handle();
        int drank = distance_rank(current);
        if(drank == prev_distance_rank) {
            CHECK(at_that_rank < 4);
            CHECK(prev_distance_rank >= 0);
            ++at_that_rank;
        } else {
            CHECK(drank == prev_distance_rank + 1);
            CHECK(at_that_rank == 4);
            prev_distance_rank = drank;
            at_that_rank = 1;
        }
        filter.point_count += 1;
    }

    for(const Point &p : points) {
        CHECK(search.start_search(p));
        CHECK(*search.current().handle() == p);
    }
}

struct TrivialFilter {
    template<typename T> void reset(const T & /*query*/) {}
    bool operator()() const { return false; }
    template<typename T> bool operator()(const T &) { return false; }
};

TEST_CASE_TEMPLATE("Exact incremental search - empty filter", Kernel, CGAL::Exact_predicates_exact_constructions_kernel,
                   CGAL::Exact_predicates_inexact_constructions_kernel) {
    using Point = typename Kernel::Point_2;
    using Points = std::vector<Point>;
    using PointIterator = typename Points::iterator;
    using Tree = mwt::Point_quadtree<Kernel, PointIterator>;
    using Search = typename Tree::template Incremental_search<TrivialFilter>;
    using Traits = mwt::Mwt_traits_2<Kernel>;

    CGAL::Protect_FPU_rounding rounder;

    Points points{{76, 684},  {292, 616}, {132, 548}, {224, 748}, {228, 520}, {176, 700}, {552, 644}, {440, 524},
                  {416, 580}, {404, 712}, {308, 700}, {324, 492}, {92, 448},  {76, 624},  {196, 600}, {188, 444},
                  {356, 424}, {356, 552}, {264, 552}, {432, 632}, {432, 632}, {332, 756}, {456, 772}, {460, 640},
                  {344, 644}, {332, 408}, {196, 488}, {248, 464}, {128, 612}, {240, 612}, {84, 512},  {72, 568},
                  {28, 580},  {28, 696},  {80, 776},  {156, 772}, {112, 700}, {256, 684}, {224, 664}, {272, 652},
                  {248, 664}, {296, 648}, {268, 592}};
    Tree tree(points.begin(), points.end());
    Point query{272, 652};
    Search search(&tree, query, TrivialFilter{});

    CHECK(search.start_search(query));
    Point last_point{query};
    Traits trait_object;
    auto less_distance = trait_object.less_distance_2_object();
    std::size_t counter = 0;
    while(search.next()) {
        auto handle = search.current().handle();
        Point cur = *handle;
        CHECK(!CGAL::has_larger_distance_to_point(query, last_point, cur));
        last_point = cur;
        ++counter;
    }
    CHECK(counter == points.size() - 1);
}
