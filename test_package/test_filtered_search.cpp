#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/FPU.h>
#include <CGAL_MWT/Directional_filter.h>
#include <CGAL_MWT/Exact_diamond_filtered_search_driver.h>
#include <CGAL_MWT/Generate_random_integral.h>
#include <CGAL_MWT/Mwt_traits.h>
#include <CGAL_MWT/Static_quadtree.h>
#include <chrono>
#include <doctest/doctest.h>
#include <iomanip>
#include <random>

/**
 * Test the diamond filtered search with a simple
 * example.
 */
TEST_CASE_TEMPLATE("Simple diamond filter test", Kernel, CGAL::Exact_predicates_exact_constructions_kernel,
                   CGAL::Exact_predicates_inexact_constructions_kernel) {
    using Traits = mwt::Mwt_traits_2<Kernel>;
    using Point = typename Traits::Point_2;
    using Points = std::vector<Point>;
    using PointIterator = typename Points::iterator;
    using Tree = mwt::Point_quadtree<Kernel, PointIterator>;

    CGAL::Protect_FPU_rounding rounder;
    Points points{{256, 512}, {304, 544}, {304, 592}, {272, 624}, {368, 624}, {368, 544}, {368, 464},
                  {336, 480}, {432, 496}, {432, 560}, {320, 640}, {416, 608}, {352, 576}};
    Points expected{{304, 544}, {336, 480}, {304, 592}, {272, 624}, {368, 544}, {368, 464}, {320, 640}};
    Tree tree(points.begin(), points.end());
    mwt::Exact_diamond_filtered_search_driver<Traits, PointIterator> search_driver(&tree);
    std::size_t count = 0;
    search_driver.enumerate_filtered_neighbors_of(Point{256, 512}, [&](PointIterator it) {
        CHECK(count < expected.size());
        CHECK(*it == expected[count]);
        ++count;
    });
    CHECK(count == expected.size());
}

/**
 * Simple check for the diamond predicate object.
 */
TEST_CASE_TEMPLATE("Diamond filter debug check", Kernel, CGAL::Exact_predicates_exact_constructions_kernel,
                   CGAL::Exact_predicates_inexact_constructions_kernel) {
    using Traits = mwt::Mwt_traits_2<Kernel>;
    using Point = typename Traits::Point_2;
    CGAL::Protect_FPU_rounding rounder;
    auto dtest = Traits{}.diamond_test_2_object();
    Point r{105, 265};
    Point q{104, 332};
    Point u{125, 56};
    Point l{128, 264};
    CHECK(dtest(q, u, l));
    CHECK(dtest(u, q, r));
}

/**
 * Test the diamond filtered search with a
 * randomly chosen example.
 */
TEST_CASE_TEMPLATE("Random diamond filter test vs. diamond predicate", Kernel,
                   CGAL::Exact_predicates_exact_constructions_kernel,
                   CGAL::Exact_predicates_inexact_constructions_kernel) {
    using Traits = mwt::Mwt_traits_2<Kernel>;
    using Point = typename Traits::Point_2;
    using Points = std::vector<Point>;
    using PointIterator = typename Points::iterator;
    using Tree = mwt::Point_quadtree<Kernel, PointIterator>;

    CGAL::Protect_FPU_rounding rounder;
    constexpr int NUM_POINTS_GENERATED = 1'000;
    constexpr int NUM_UNIQUE_EXPECTED = 950;
    constexpr int NUM_RANDOM_SEARCHES = 20;
    CHECK(NUM_POINTS_GENERATED >= NUM_UNIQUE_EXPECTED);
    std::mt19937_64 rng(42);
    auto points = mwt::generate_random_integral<Kernel>({0, 0, 1024, 1024}, NUM_POINTS_GENERATED, rng);
    CHECK(points.size() > NUM_UNIQUE_EXPECTED);

    auto before_tree = std::chrono::steady_clock::now();
    Tree tree(points.begin(), points.end());
    auto after_tree = std::chrono::steady_clock::now();
    double seconds = std::chrono::duration_cast<std::chrono::duration<double>>(after_tree - before_tree).count();
    CHECK(seconds < 1.0);

    mwt::Exact_diamond_filtered_search_driver<Traits, PointIterator> search_driver(&tree);
    std::uniform_int_distribution<std::size_t> qpdist(0, points.size() - 1);
    struct SingleSearch {
        PointIterator query;
        std::vector<PointIterator> result;
    };

    std::vector<SingleSearch> searches;
    auto before_searches = std::chrono::steady_clock::now();
    for(int i = 0; i < NUM_RANDOM_SEARCHES; ++i) {
        PointIterator query = points.begin() + qpdist(rng);
        std::vector<PointIterator> result;
        search_driver.enumerate_filtered_neighbors_of(*query, [&](PointIterator it) { result.push_back(it); });
        searches.emplace_back(SingleSearch{query, std::move(result)});
    }
    auto after_searches = std::chrono::steady_clock::now();
    double search_seconds =
        std::chrono::duration_cast<std::chrono::duration<double>>(after_searches - before_searches).count();
    CHECK(search_seconds < 1.0);

    auto fails_diamond_test = [&](PointIterator unreported, const SingleSearch &search) {
        const Point q = *search.query;
        const Point u = *unreported;
        if(!typename Kernel::Less_xy_2{}(q, u)) {
            return;
        }
        std::vector<Point> reported_left, reported_right;
        for(const auto piter : search.result) {
            auto o = CGAL::orientation(q, u, *piter);
            if(o == CGAL::COLLINEAR) {
                CHECK(CGAL::collinear_are_ordered_along_line(q, *piter, u));
                return;
            }
            if(o == CGAL::LEFT_TURN) {
                reported_left.push_back(*piter);
            } else {
                reported_right.push_back(*piter);
            }
        }

        auto dtest = Traits{}.diamond_test_2_object();
        bool in_left = std::any_of(reported_left.begin(), reported_left.end(), [&](Point p) { return dtest(q, u, p); });
        bool in_right =
            std::any_of(reported_right.begin(), reported_right.end(), [&](Point p) { return dtest(u, q, p); });
        std::vector<Point> all_left, all_right;
        for(const auto p : points) {
            if(p == q || p == u) {
                continue;
            }
            auto o = CGAL::orientation(q, u, p);
            if(o == CGAL::COLLINEAR) {
                if(CGAL::collinear_are_ordered_along_line(q, p, u)) {
                    return;
                }
            } else {
                if(o == CGAL::LEFT_TURN) {
                    all_left.push_back(p);
                } else {
                    all_right.push_back(p);
                }
            }
        }
        if(!in_left) {
            CHECK(std::any_of(all_left.begin(), all_left.end(), [&](Point p) { return dtest(q, u, p); }));
        }
        if(!in_right) {
            CHECK(std::any_of(all_right.begin(), all_right.end(), [&](Point p) { return dtest(u, q, p); }));
        }
    };

    auto verify_diamonds = [&](const SingleSearch &search) {
        // we do not want the query to be reported
        CHECK(std::find(search.result.begin(), search.result.end(), search.query) == search.result.end());

        // points must be in order of non-decreasing distance
        PointIterator previous = search.query;
        for(const auto point : search.result) {
            CHECK(point != previous);
            CHECK(!CGAL::has_larger_distance_to_point(*search.query, *previous, *point));
            previous = point;
        }

        std::vector<bool> indices_reported(points.size(), false);
        for(const auto point : search.result) {
            indices_reported[point - points.begin()] = true;
        }
        std::vector<PointIterator> unreported;
        for(std::size_t i = 0; i < points.size(); ++i) {
            if(!indices_reported[i]) {
                unreported.push_back(points.begin() + i);
            }
        }

        for(PointIterator unrep : unreported) {
            fails_diamond_test(unrep, search);
        }
    };

    for(const auto &search : searches) {
        verify_diamonds(search);
    }
}

namespace mwt {

namespace test {

class Test_exact_diamond_filtered_search_driver {
  public:
    template<typename Kernel> static void test_point_set() {
        using Traits = Mwt_traits_2<Kernel>;
        using Point = typename Traits::Point_2;
        using Points = std::vector<Point>;
        using PointIterator = typename Points::iterator;
        using PointSet = typename mwt::Exact_diamond_filtered_search_driver<Traits, PointIterator>::PointSet;
        using PointEntry = typename PointSet::value_type;

        auto pseudo_angle_ = mwt::Compute_pseudo_angle<CGAL::Interval_nt_advanced>{};
        auto pseudo_angle = [&](const Point &query, const Point &point) {
            CGAL::Interval_nt_advanced px{CGAL::to_double(point.x())};
            CGAL::Interval_nt_advanced py{CGAL::to_double(point.y())};
            px -= CGAL::to_double(query.x());
            py -= CGAL::to_double(query.y());
            return pseudo_angle_(px, py);
        };
        constexpr int NUM_POINTS_GENERATED = 1'000;

        std::mt19937_64 rng(4211);
        std::uniform_int_distribution<int> dist(-5000, 5000);
        Point query{0, 0};
        Points points;
        for(int i = 0; i < NUM_POINTS_GENERATED; ++i) {
            points.emplace_back(dist(rng), dist(rng));
        }
        std::sort(points.begin(), points.end(), typename Kernel::Less_xy_2{});
        points.erase(std::unique(points.begin(), points.end()), points.end());

        PointSet point_set;
        for(auto it = points.begin(); it != points.end(); ++it) {
            auto res = point_set.emplace(it, query);
            if(res.second) {
                if(res.first != point_set.begin()) {
                    auto prev = std::prev(res.first);
                    CHECK(CGAL::possibly(pseudo_angle(query, *prev->point) < pseudo_angle(query, *res.first->point)));
                }
                if(std::next(res.first) != point_set.end()) {
                    auto next = std::next(res.first);
                    CHECK(CGAL::possibly(pseudo_angle(query, *res.first->point) < pseudo_angle(query, *next->point)));
                }
            } else {
                CHECK(CGAL::collinear(query, *res.first->point, *it));
            }
        }
    }
};

} // namespace test

} // namespace mwt

TEST_CASE("Diamond filter point set") {
    CGAL::Protect_FPU_rounding rounder;
    mwt::test::Test_exact_diamond_filtered_search_driver::test_point_set<
        CGAL::Exact_predicates_inexact_constructions_kernel>();
}

TEST_CASE_TEMPLATE("Diamond filter/Dead sector test (issue #x1)", Kernel,
                   CGAL::Exact_predicates_exact_constructions_kernel,
                   CGAL::Exact_predicates_inexact_constructions_kernel) {
    using Traits = mwt::Mwt_traits_2<Kernel>;
    using Point = typename Traits::Point_2;
    using Points = std::vector<Point>;
    using PointIterator = typename Points::iterator;

    CGAL::Protect_FPU_rounding rounder;
    mwt::Exact_diamond_filter<Traits> filter;
    Point query{-25441084, -31235852};
    filter.reset(query);
    Points early_points{{-25337753, -31444845}, {-25341638, -30272568}, {-25488203, -30864405}, {-25292734, -31241373},
                        {-25485864, -31002255}, {-25452033, -32451733}, {-25557599, -32003402}, {-25271848, -31442036},
                        {-25278292, -31233904}, {-25639231, -30714980}, {-25264940, -30164614}, {-25318547, -30307017},
                        {-25687290, -30873973}, {-25388460, -31083377}, {-25522465, -30496683}, {-25606301, -31617244},
                        {-25264012, -31768132}, {-25516774, -30661322}, {-25206576, -31176077}, {-25340284, -31157452},
                        {-25205267, -29993821}, {-25501794, -30426132}, {-25542914, -32337100}, {-25267836, -31946682},
                        {-25444756, -31393350}, {-25671310, -31823070}, {-25187991, -31876165}, {-25670132, -30861392},
                        {-25398134, -30973998}, {-25330167, -31322046}, {-25477938, -30426057}, {-25449960, -30660692},
                        {-25182562, -31614425}, {-25589750, -32188997}, {-25471441, -30173824}, {-25671400, -32399976},
                        {-25437936, -32154154}, {-25346591, -30324113}, {-25394944, -31043446}, {-25441178, -30692257},
                        {-25251821, -31522451}, {-25278554, -31208819}, {-25360253, -30176939}, {-25533988, -31592127},
                        {-25533698, -30879046}, {-25346441, -31289979}, {-25319108, -31037453}, {-25377202, -30181220},
                        {-25302661, -31171434}, {-25468577, -30553926}, {-25452994, -32028519}, {-25470457, -32407370},
                        {-25481124, -31894144}, {-25372253, -31575531}, {-25175654, -31682966}, {-25299927, -30189861},
                        {-25464824, -31397672}, {-25608783, -31761705}, {-25484094, -30232756}, {-25506282, -31816541},
                        {-25682514, -30799307}, {-25513251, -31035362}, {-25283124, -31633767}, {-25350477, -31460579},
                        {-25260266, -31129435}, {-25192365, -30946942}, {-25183662, -29968829}, {-25501057, -31975079},
                        {-25665493, -32240982}, {-25263151, -31711225}, {-25337078, -31061371}, {-25277556, -31439656},
                        {-25421822, -31678084}, {-25671062, -32324993}, {-25214092, -30312435}, {-25242349, -30340303},
                        {-25372436, -31153449}, {-25214852, -29716788}, {-25291454, -29691837}, {-25554366, -31767627},
                        {-25626167, -31005309}, {-25382431, -31496728}, {-25467804, -30282363}, {-25294687, -30393257},
                        {-25189259, -31382887}, {-25265333, -31638053}, {-25578181, -30777803}, {-25632658, -31580890},
                        {-25185153, -30246881}, {-25483575, -31897836}, {-25334075, -31075203}, {-25223192, -29841367},
                        {-25370863, -31595796}, {-25617440, -32063503}, {-25447825, -30456764}, {-25233085, -29847966},
                        {-25198907, -29804586}, {-25439919, -30579809}, {-25640419, -32097839}, {-25438224, -32176058},
                        {-25505647, -30617729}, {-25386145, -31045820}, {-25269176, -30329709}, {-25465297, -30526491},
                        {-25463900, -32216924}, {-25178279, -29898432}, {-25631696, -30749905}, {-25437870, -32096488},
                        {-25343704, -31442048}, {-25551637, -31804291}, {-25237204, -29882907}, {-25503656, -32396869},
                        {-25277819, -31790359}, {-25498140, -30606835}, {-25240077, -31136790}, {-25285766, -30946603},
                        {-25486853, -30563001}, {-25637980, -32410557}, {-25393327, -31274766}, {-25525157, -31904107},
                        {-25360294, -30300505}, {-25649999, -32243166}, {-25210474, -30207635}, {-25244258, -29727644},
                        {-25306165, -31079955}, {-25569738, -32283529}, {-25449630, -31100149}, {-25393021, -31486707},
                        {-25634397, -30916434}, {-25427159, -31604187}, {-25185816, -31958026}, {-25288271, -31665605},
                        {-25687149, -31556153}, {-25540850, -32100687}, {-25359780, -30339311}, {-25274348, -31565406},
                        {-25184087, -29645082}, {-25387377, -30214797}, {-25450917, -30398112}, {-25640933, -30695172},
                        {-25454322, -32254793}, {-25479484, -30254958}, {-25471613, -30807768}, {-25657009, -31920112},
                        {-25579559, -32323596}, {-25274546, -31557160}, {-25175804, -30009219}, {-25350613, -31530763},
                        {-25500947, -30776903}, {-25208750, -31188964}, {-25449597, -31339161}, {-25486020, -31818902},
                        {-25269666, -31837398}, {-25263642, -31133239}, {-25231353, -29728477}, {-25542147, -30323960},
                        {-25627616, -30729635}, {-25457189, -30159946}, {-25375039, -31477972}, {-25169476, -31211139},
                        {-25191165, -31720369}, {-25285268, -29891648}, {-25492608, -32004788}, {-25661862, -31896372},
                        {-25298866, -31092513}, {-25232336, -31955385}, {-25563077, -31622943}, {-25440899, -32434984},
                        {-25411759, -31508173}, {-25284462, -29769014}, {-25651466, -31951313}, {-25548040, -32006385},
                        {-25523685, -32346085}, {-25505018, -30652019}, {-25318474, -31331988}, {-25221252, -31762718},
                        {-25552727, -30971600}, {-25643370, -30920742}, {-25605201, -32476268}, {-25384188, -30174175},
                        {-25304411, -31407289}, {-25208529, -30077864}, {-25558325, -30782269}, {-25501329, -31976769},
                        {-25167472, -29655181}, {-25602861, -31815972}, {-25447265, -30309184}, {-25550189, -30421721}};
    // check point uniqueness
    struct PHash {
        std::size_t operator()(const Point &p) const {
            auto hasher = std::hash<double>{};
            return hasher(CGAL::to_double(p.x())) ^ hasher(CGAL::to_double(p.y()));
        }
    };
    std::unordered_set<Point, PHash> point_set;
    for(const Point &p : early_points) {
        if(!point_set.insert(p).second) {
            CHECK(p == Point{0, 0});
        }
    }

    // sorting
    auto ordering = [&](const Point &a, const Point &b) -> bool {
        bool abelow = (a.y() <= query.y());
        bool bbelow = (b.y() <= query.y());
        if(abelow != bbelow) {
            return !abelow;
        }
        return a != b && CGAL::right_turn(query, a, b);
    };
    std::sort(early_points.begin(), early_points.end(), ordering);
    for(auto cur = early_points.begin() + 1, prev = early_points.begin(); cur != early_points.end(); ++cur, ++prev) {
        CHECK(CGAL::right_turn(query, *prev, *cur));
        CHECK(ordering(*prev, *cur));
        filter.insert_sector(prev, cur);
    }
    CGAL::Iso_rectangle_2<Kernel> box1{-16777220, -31457262, -16252933, -30932974};
    CHECK(filter(box1));
    CHECK(filter());
}
