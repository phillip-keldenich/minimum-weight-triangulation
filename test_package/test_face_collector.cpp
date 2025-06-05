#include <CGAL/Convex_hull_traits_adapter_2.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/FPU.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/property_map.h>
#include <CGAL_MWT/Diamond_filter.h>
#include <CGAL_MWT/Face_analyzer.h>
#include <CGAL_MWT/Face_collector.h>
#include <CGAL_MWT/LMT_skeleton.h>
#include <CGAL_MWT/Mwt_traits.h>
#include <CGAL_MWT/Static_quadtree.h>
#include <boost/iterator/counting_iterator.hpp>
#include <doctest/doctest.h>
#include <random>

namespace mwt {

namespace test {

class Face_collector_tester {
  public:
    template<typename S> static void verify_initial(Face_collector<S> &c) { c.p_verify_initial(); }
};

} // namespace test

} // namespace mwt

struct HalfedgeInfo {
    std::size_t source_id, target_id;
    std::size_t next_id, twin_id;
    mwt::LMTStatus status;

    HalfedgeInfo(std::size_t source_id, std::size_t target_id, std::size_t next_id, std::size_t twin_id,
                 mwt::LMTStatus status)
        : source_id{source_id}, target_id{target_id}, next_id{next_id}, twin_id{twin_id}, status{status} {}
};

template<typename Halfedge, typename Point>
std::vector<Halfedge> halfedges_from_info(const std::vector<Point> &points, const std::vector<HalfedgeInfo> &info) {
    std::vector<Halfedge> halfedges;
    halfedges.resize(info.size());
    for(std::size_t i = 0; i < info.size(); ++i) {
        Halfedge &e = halfedges[i];
        const auto &ei = info[i];
        e.tar = &points[ei.target_id];
        e.twn = &halfedges[ei.twin_id];
        e.next = &halfedges[ei.next_id];
        e.flagged = false;
        e.status = ei.status;
    }
    return halfedges;
}

using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point = Kernel::Point_2;
using Points = std::vector<Point>;
using Traits = mwt::Mwt_traits_2<Kernel>;
using Tree = mwt::Point_quadtree<Traits, Points::iterator>;
using Skeleton = mwt::LMTSkeleton<Traits, Tree>;
using Halfedge = Skeleton::Halfedge;
using Filter = mwt::Diamond_edge_filter<Tree, false>;

static Points fake_lmt_test1_points() {
    return Points{{224, 704}, {416, 704}, {192, 672}, {448, 672}, {256, 640}, {320, 640},
                  {384, 640}, {256, 576}, {384, 576}, {256, 512}, {320, 512}, {384, 512},
                  {192, 480}, {448, 480}, {224, 448}, {416, 448}, {320, 416}, {208, 576}};
}

template<typename Range1Iter, typename Range2Iter>
static bool same_cyclic_order(Range1Iter begin1, Range1Iter end1, Range2Iter begin2, Range2Iter end2) {
    if(std::distance(begin1, end1) != std::distance(begin2, end2)) {
        return false; // different lengths
    }
    if(begin1 == end1) {
        return true; // both empty
    }

    // go (cyclically) from one element in sequence 2 to the next
    auto next2 = [&](auto it) {
        if(++it == end2) {
            it = begin2;
        }
        return it;
    };

    // for each occurrence of *begin1 in sequence 2, check if the rest of the elements match
    auto candidate_begin = begin2;
    while((candidate_begin = std::find(candidate_begin, end2, *begin1)) != end2) {
        auto cur2 = candidate_begin;
        bool are_same = true;
        for(auto it = std::next(begin1); it != end1; ++it) {
            cur2 = next2(cur2);
            if(*it != *cur2) {
                are_same = false;
                break;
            }
        }
        if(are_same)
            return true;
        ++candidate_begin;
    }
    return false;
}

TEST_CASE("Check same_cyclic_order") {
    std::vector<int> v1{1, 2, 3, 4, 5};
    std::vector<int> v2{3, 4, 5, 1, 2};
    std::vector<int> v3{3, 4, 5, 2, 1};
    CHECK(same_cyclic_order(v1.begin(), v1.end(), v2.begin(), v2.end()));
    CHECK(!same_cyclic_order(v1.begin(), v1.end(), v3.begin(), v3.end()));
}

std::vector<std::size_t> convex_hull_indices(const std::vector<Point> &points) {
    using Adapter = CGAL::Convex_hull_traits_adapter_2<Kernel, CGAL::Pointer_property_map<Point>::const_type>;
    using IndexIterator = boost::counting_iterator<std::size_t>;
    std::vector<std::size_t> result;
    CGAL::convex_hull_2(IndexIterator(0), IndexIterator(points.size()), std::back_inserter(result),
                        Adapter(CGAL::make_property_map(points.data())));
    return result;
}

TEST_CASE("Check convex_hull_indices") {
    CGAL::Protect_FPU_rounding rounder;
    std::vector<std::size_t> indices = convex_hull_indices(fake_lmt_test1_points());
    std::vector<std::size_t> expected = {0, 2, 12, 14, 16, 15, 13, 3, 1};
    CHECK(indices.size() == expected.size());
    CHECK(same_cyclic_order(indices.begin(), indices.end(), expected.begin(), expected.end()));
}

/**
 * Construct halfedges (with correctly set next and twin pointers)
 * from a list of points and lists of certain and possible edges.
 * The convex hull is added automatically and must not be in the skeleton_edges list.
 */
std::vector<Halfedge> halfedges_from_graph(const std::vector<Point> &points,
                                           const std::vector<std::pair<std::size_t, std::size_t>> &skeleton_edges,
                                           const std::vector<std::pair<std::size_t, std::size_t>> &possible_edges) {
    using LMTStatus = mwt::LMTStatus;
    std::vector<std::pair<std::size_t, std::size_t>> certain_edges = skeleton_edges;
    std::vector<std::size_t> hull_indices = convex_hull_indices(points);
    std::size_t prev = hull_indices.back();
    std::vector<std::pair<std::size_t, std::size_t>> convex_hull_edges;
    for(std::size_t i : hull_indices) {
        convex_hull_edges.emplace_back(prev, i);
        prev = i;
    }
    std::vector<std::size_t> out_edge_begin;
    std::vector<std::size_t> edge_counts(points.size(), 0);
    auto count_edge = [&](std::pair<std::size_t, std::size_t> edge) {
        ++edge_counts[edge.first];
        ++edge_counts[edge.second];
    };
    for(const auto &e : certain_edges) {
        count_edge(e);
    }
    for(const auto &e : convex_hull_edges) {
        count_edge(e);
    }
    for(const auto &e : possible_edges) {
        count_edge(e);
    }
    std::vector<std::size_t> edge_offsets;
    edge_offsets.push_back(0);
    for(std::size_t i = 0; i < edge_counts.size(); ++i) {
        edge_offsets.push_back(edge_offsets.back() + edge_counts[i]);
    }
    std::vector<Halfedge> result{2 * certain_edges.size() + 2 * possible_edges.size() + 2 * convex_hull_edges.size()};
    edge_counts.assign(points.size(), 0);
    auto append_edge = [&](std::pair<std::size_t, std::size_t> edge, LMTStatus status) {
        std::size_t o1 = edge_offsets[edge.first] + edge_counts[edge.first]++;
        std::size_t o2 = edge_offsets[edge.second] + edge_counts[edge.second]++;
        result[o1] = Halfedge{&points[edge.second], &result[o2]};
        result[o1].flagged = false;
        result[o1].status = status;
        result[o2] = Halfedge{&points[edge.first], &result[o1]};
        result[o2].flagged = false;
        result[o2].status = status;
    };
    for(const auto &e : certain_edges) {
        append_edge(e, LMTStatus::Certain);
    }
    for(const auto &e : possible_edges) {
        append_edge(e, LMTStatus::Possible);
    }
    for(const auto &e : convex_hull_edges) {
        append_edge(e, LMTStatus::CH);
    }
    for(std::size_t i = 0; i < points.size(); ++i) {
        std::size_t begin_offset = edge_offsets[i];
        std::size_t end_offset = edge_offsets[i + 1];
        mwt::detail::sort_halfedges_around_point(points[i], &result[begin_offset], &result[end_offset]);
        mwt::detail::init_start_pointers(&result[begin_offset], &result[end_offset]);
    }
    return result;
}

std::pair<Points, std::vector<Halfedge>> test_case1() {
    Points ordered_points{
        Point{-364320, -2427080},  Point{-3885180, -2102440}, Point{-3834044, -3092987}, Point{-648985, -2381971},
        Point{-974585, -1989286},  Point{-1114909, -3820815}, Point{-3666032, -3669508}, Point{-1903415, -269212},
        Point{-2041361, -2425846}, Point{-3271404, -744881},  Point{-3154414, -2437543}, Point{-2280531, -4094692},
        Point{-2594934, -1547983}, Point{-2494082, 581525},   Point{-2409929, 1450935},  Point{-2874937, 1748843},
        Point{-2154519, 1132808},  Point{-3617374, 1923819},  Point{-3683765, 857970},   Point{-4097386, 1944838},
        Point{-3797937, 3422891},  Point{-3701726, 3655770},  Point{-3641864, 3038848},  Point{-217937, 3797519},
        Point{-1480722, 3559245},  Point{-1550583, 2918217},  Point{-581460, 2108960},   Point{-1509898, 2537066},
        Point{-809745, 1307314},   Point{-969261, 1317646},   Point{-1200241, 608172},   Point{4133122, 4058162},
        Point{4062611, 760321},    Point{3425885, 3419553},   Point{3086876, 3190494},   Point{875484, 2519280},
        Point{785569, -403682},    Point{896404, -2292982},   Point{903090, -967665},    Point{1024596, -3106057},
        Point{1295133, -57481},    Point{1350785, -1415805},  Point{2525839, -1112449},  Point{2682794, -576839},
        Point{477766, -2778850},   Point{3124081, -3944745},  Point{476973, -61158},     Point{4021511, -2666645},
        Point{451064, -3053294},   Point{300657, -697448},
    };

    using LMTStatus = mwt::LMTStatus;
    std::vector<HalfedgeInfo> info{{0, 37, 4, 575, LMTStatus::Possible},
                                   {0, 41, 4, 636, LMTStatus::Impossible},
                                   {0, 38, 4, 593, LMTStatus::Impossible},
                                   {0, 36, 4, 556, LMTStatus::Impossible},
                                   {0, 49, 8, 760, LMTStatus::Possible},
                                   {0, 46, 8, 717, LMTStatus::Impossible},
                                   {0, 13, 8, 217, LMTStatus::Impossible},
                                   {0, 7, 8, 119, LMTStatus::Impossible},
                                   {0, 4, 9, 76, LMTStatus::Certain},
                                   {0, 3, 15, 56, LMTStatus::Certain},
                                   {0, 8, 15, 144, LMTStatus::Impossible},
                                   {0, 10, 15, 159, LMTStatus::Impossible},
                                   {0, 2, 15, 29, LMTStatus::Impossible},
                                   {0, 11, 15, 180, LMTStatus::Impossible},
                                   {0, 5, 15, 83, LMTStatus::Impossible},
                                   {0, 48, 16, 739, LMTStatus::Certain},
                                   {0, 44, 0, 681, LMTStatus::Certain},
                                   {1, 12, 19, 196, LMTStatus::Certain},
                                   {1, 13, 19, 213, LMTStatus::Impossible},
                                   {1, 9, 23, 153, LMTStatus::Certain},
                                   {1, 14, 23, 230, LMTStatus::Impossible},
                                   {1, 18, 23, 292, LMTStatus::Impossible},
                                   {1, 17, 23, 278, LMTStatus::Impossible},
                                   {1, 19, 24, 303, LMTStatus::Certain},
                                   {1, 2, 27, 34, LMTStatus::Certain},
                                   {1, 6, 27, 99, LMTStatus::Impossible},
                                   {1, 5, 27, 89, LMTStatus::Impossible},
                                   {1, 10, 17, 170, LMTStatus::Certain},
                                   {1, 8, 17, 135, LMTStatus::Impossible},
                                   {2, 0, 32, 12, LMTStatus::Impossible},
                                   {2, 3, 32, 50, LMTStatus::Impossible},
                                   {2, 8, 32, 137, LMTStatus::Impossible},
                                   {2, 10, 34, 171, LMTStatus::Certain},
                                   {2, 12, 34, 197, LMTStatus::Impossible},
                                   {2, 1, 35, 24, LMTStatus::Certain},
                                   {2, 19, 36, 302, LMTStatus::CH},
                                   {2, 6, 32, 100, LMTStatus::CH},
                                   {2, 11, 32, 186, LMTStatus::Impossible},
                                   {2, 5, 32, 90, LMTStatus::Impossible},
                                   {3, 37, 47, 574, LMTStatus::Impossible},
                                   {3, 41, 47, 635, LMTStatus::Impossible},
                                   {3, 38, 47, 592, LMTStatus::Impossible},
                                   {3, 36, 47, 555, LMTStatus::Impossible},
                                   {3, 49, 47, 759, LMTStatus::Impossible},
                                   {3, 46, 47, 716, LMTStatus::Impossible},
                                   {3, 30, 47, 472, LMTStatus::Impossible},
                                   {3, 7, 47, 118, LMTStatus::Impossible},
                                   {3, 4, 49, 75, LMTStatus::Certain},
                                   {3, 10, 49, 161, LMTStatus::Impossible},
                                   {3, 8, 53, 124, LMTStatus::Certain},
                                   {3, 2, 53, 30, LMTStatus::Impossible},
                                   {3, 6, 53, 94, LMTStatus::Impossible},
                                   {3, 11, 53, 181, LMTStatus::Impossible},
                                   {3, 5, 54, 84, LMTStatus::Certain},
                                   {3, 48, 56, 740, LMTStatus::Certain},
                                   {3, 44, 56, 682, LMTStatus::Impossible},
                                   {3, 0, 47, 9, LMTStatus::Certain},
                                   {4, 41, 60, 634, LMTStatus::Impossible},
                                   {4, 38, 60, 591, LMTStatus::Impossible},
                                   {4, 36, 60, 554, LMTStatus::Impossible},
                                   {4, 49, 66, 758, LMTStatus::Certain},
                                   {4, 46, 66, 715, LMTStatus::Impossible},
                                   {4, 28, 66, 431, LMTStatus::Impossible},
                                   {4, 29, 66, 452, LMTStatus::Impossible},
                                   {4, 30, 66, 471, LMTStatus::Impossible},
                                   {4, 16, 66, 269, LMTStatus::Impossible},
                                   {4, 7, 69, 117, LMTStatus::Certain},
                                   {4, 13, 69, 215, LMTStatus::Impossible},
                                   {4, 9, 69, 157, LMTStatus::Impossible},
                                   {4, 12, 71, 203, LMTStatus::Certain},
                                   {4, 10, 71, 162, LMTStatus::Impossible},
                                   {4, 8, 75, 125, LMTStatus::Certain},
                                   {4, 6, 75, 95, LMTStatus::Impossible},
                                   {4, 11, 75, 182, LMTStatus::Impossible},
                                   {4, 5, 75, 85, LMTStatus::Impossible},
                                   {4, 3, 76, 47, LMTStatus::Certain},
                                   {4, 0, 77, 8, LMTStatus::Certain},
                                   {4, 44, 60, 680, LMTStatus::Possible},
                                   {4, 37, 60, 573, LMTStatus::Impossible},
                                   {5, 47, 80, 731, LMTStatus::Impossible},
                                   {5, 39, 81, 608, LMTStatus::Certain},
                                   {5, 48, 84, 743, LMTStatus::Certain},
                                   {5, 44, 84, 685, LMTStatus::Impossible},
                                   {5, 0, 84, 14, LMTStatus::Impossible},
                                   {5, 3, 87, 53, LMTStatus::Certain},
                                   {5, 4, 87, 74, LMTStatus::Impossible},
                                   {5, 12, 87, 202, LMTStatus::Impossible},
                                   {5, 8, 92, 140, LMTStatus::Certain},
                                   {5, 10, 92, 174, LMTStatus::Impossible},
                                   {5, 1, 92, 26, LMTStatus::Impossible},
                                   {5, 2, 92, 38, LMTStatus::Impossible},
                                   {5, 6, 92, 102, LMTStatus::Impossible},
                                   {5, 11, 93, 177, LMTStatus::Certain},
                                   {5, 45, 80, 698, LMTStatus::Certain},
                                   {6, 3, 98, 51, LMTStatus::Impossible},
                                   {6, 4, 98, 72, LMTStatus::Impossible},
                                   {6, 8, 98, 138, LMTStatus::Impossible},
                                   {6, 12, 98, 199, LMTStatus::Impossible},
                                   {6, 10, 100, 172, LMTStatus::Certain},
                                   {6, 1, 100, 25, LMTStatus::Impossible},
                                   {6, 2, 101, 36, LMTStatus::CH},
                                   {6, 11, 98, 187, LMTStatus::CH},
                                   {6, 5, 98, 91, LMTStatus::Impossible},
                                   {7, 46, 105, 711, LMTStatus::Impossible},
                                   {7, 35, 105, 535, LMTStatus::Impossible},
                                   {7, 30, 110, 467, LMTStatus::Certain},
                                   {7, 28, 110, 428, LMTStatus::Impossible},
                                   {7, 29, 110, 449, LMTStatus::Impossible},
                                   {7, 26, 110, 392, LMTStatus::Impossible},
                                   {7, 16, 110, 268, LMTStatus::Impossible},
                                   {7, 13, 113, 216, LMTStatus::Certain},
                                   {7, 17, 113, 280, LMTStatus::Impossible},
                                   {7, 18, 113, 296, LMTStatus::Impossible},
                                   {7, 9, 115, 147, LMTStatus::Certain},
                                   {7, 10, 115, 166, LMTStatus::Impossible},
                                   {7, 12, 117, 192, LMTStatus::Certain},
                                   {7, 8, 117, 132, LMTStatus::Impossible},
                                   {7, 4, 122, 66, LMTStatus::Certain},
                                   {7, 3, 122, 46, LMTStatus::Impossible},
                                   {7, 0, 122, 7, LMTStatus::Impossible},
                                   {7, 41, 122, 632, LMTStatus::Impossible},
                                   {7, 38, 122, 589, LMTStatus::Impossible},
                                   {7, 49, 105, 754, LMTStatus::Certain},
                                   {7, 36, 105, 551, LMTStatus::Impossible},
                                   {8, 3, 125, 49, LMTStatus::Certain},
                                   {8, 4, 133, 71, LMTStatus::Certain},
                                   {8, 38, 133, 590, LMTStatus::Impossible},
                                   {8, 49, 133, 757, LMTStatus::Impossible},
                                   {8, 46, 133, 714, LMTStatus::Impossible},
                                   {8, 28, 133, 430, LMTStatus::Impossible},
                                   {8, 29, 133, 451, LMTStatus::Impossible},
                                   {8, 30, 133, 470, LMTStatus::Impossible},
                                   {8, 7, 133, 116, LMTStatus::Impossible},
                                   {8, 12, 136, 201, LMTStatus::Certain},
                                   {8, 9, 136, 155, LMTStatus::Impossible},
                                   {8, 1, 136, 28, LMTStatus::Impossible},
                                   {8, 10, 139, 160, LMTStatus::Certain},
                                   {8, 2, 139, 31, LMTStatus::Impossible},
                                   {8, 6, 139, 96, LMTStatus::Impossible},
                                   {8, 11, 140, 183, LMTStatus::Certain},
                                   {8, 5, 124, 87, LMTStatus::Certain},
                                   {8, 45, 124, 697, LMTStatus::Impossible},
                                   {8, 48, 124, 741, LMTStatus::Impossible},
                                   {8, 44, 124, 683, LMTStatus::Impossible},
                                   {8, 0, 124, 10, LMTStatus::Impossible},
                                   {9, 49, 147, 755, LMTStatus::Impossible},
                                   {9, 46, 147, 712, LMTStatus::Impossible},
                                   {9, 7, 149, 113, LMTStatus::Certain},
                                   {9, 16, 149, 267, LMTStatus::Impossible},
                                   {9, 13, 151, 212, LMTStatus::Certain},
                                   {9, 14, 151, 231, LMTStatus::Impossible},
                                   {9, 18, 152, 294, LMTStatus::Certain},
                                   {9, 19, 153, 305, LMTStatus::Certain},
                                   {9, 1, 156, 19, LMTStatus::Certain},
                                   {9, 10, 156, 167, LMTStatus::Impossible},
                                   {9, 8, 156, 134, LMTStatus::Impossible},
                                   {9, 12, 147, 195, LMTStatus::Certain},
                                   {9, 4, 147, 68, LMTStatus::Impossible},
                                   {9, 37, 147, 571, LMTStatus::Impossible},
                                   {10, 0, 160, 11, LMTStatus::Impossible},
                                   {10, 8, 165, 136, LMTStatus::Certain},
                                   {10, 3, 165, 48, LMTStatus::Impossible},
                                   {10, 4, 165, 70, LMTStatus::Impossible},
                                   {10, 41, 165, 633, LMTStatus::Impossible},
                                   {10, 30, 165, 469, LMTStatus::Impossible},
                                   {10, 12, 170, 198, LMTStatus::Certain},
                                   {10, 7, 170, 114, LMTStatus::Impossible},
                                   {10, 9, 170, 154, LMTStatus::Impossible},
                                   {10, 18, 170, 293, LMTStatus::Impossible},
                                   {10, 19, 170, 304, LMTStatus::Impossible},
                                   {10, 1, 171, 27, LMTStatus::Certain},
                                   {10, 2, 172, 32, LMTStatus::Certain},
                                   {10, 6, 173, 98, LMTStatus::Certain},
                                   {10, 11, 160, 185, LMTStatus::Certain},
                                   {10, 5, 160, 88, LMTStatus::Impossible},
                                   {11, 45, 177, 699, LMTStatus::CH},
                                   {11, 47, 177, 732, LMTStatus::Impossible},
                                   {11, 5, 183, 92, LMTStatus::Certain},
                                   {11, 48, 183, 742, LMTStatus::Impossible},
                                   {11, 44, 183, 684, LMTStatus::Impossible},
                                   {11, 0, 183, 13, LMTStatus::Impossible},
                                   {11, 3, 183, 52, LMTStatus::Impossible},
                                   {11, 4, 183, 73, LMTStatus::Impossible},
                                   {11, 8, 185, 139, LMTStatus::Certain},
                                   {11, 12, 185, 200, LMTStatus::Impossible},
                                   {11, 10, 187, 173, LMTStatus::Certain},
                                   {11, 2, 187, 37, LMTStatus::Impossible},
                                   {11, 6, 175, 101, LMTStatus::CH},
                                   {12, 49, 192, 756, LMTStatus::Impossible},
                                   {12, 36, 192, 552, LMTStatus::Impossible},
                                   {12, 46, 192, 713, LMTStatus::Impossible},
                                   {12, 30, 192, 468, LMTStatus::Impossible},
                                   {12, 7, 195, 115, LMTStatus::Certain},
                                   {12, 13, 195, 214, LMTStatus::Impossible},
                                   {12, 18, 195, 295, LMTStatus::Impossible},
                                   {12, 9, 196, 156, LMTStatus::Certain},
                                   {12, 1, 198, 17, LMTStatus::Certain},
                                   {12, 2, 198, 33, LMTStatus::Impossible},
                                   {12, 10, 201, 165, LMTStatus::Certain},
                                   {12, 6, 201, 97, LMTStatus::Impossible},
                                   {12, 11, 201, 184, LMTStatus::Impossible},
                                   {12, 8, 203, 133, LMTStatus::Certain},
                                   {12, 5, 203, 86, LMTStatus::Impossible},
                                   {12, 4, 192, 69, LMTStatus::Certain},
                                   {12, 37, 192, 572, LMTStatus::Impossible},
                                   {13, 30, 207, 466, LMTStatus::Certain},
                                   {13, 29, 207, 448, LMTStatus::Impossible},
                                   {13, 16, 208, 266, LMTStatus::Certain},
                                   {13, 14, 209, 232, LMTStatus::Certain},
                                   {13, 15, 211, 247, LMTStatus::Certain},
                                   {13, 17, 211, 281, LMTStatus::Impossible},
                                   {13, 18, 212, 297, LMTStatus::Certain},
                                   {13, 9, 216, 149, LMTStatus::Certain},
                                   {13, 1, 216, 18, LMTStatus::Impossible},
                                   {13, 12, 216, 193, LMTStatus::Impossible},
                                   {13, 4, 216, 67, LMTStatus::Impossible},
                                   {13, 7, 205, 110, LMTStatus::Certain},
                                   {13, 0, 205, 6, LMTStatus::Impossible},
                                   {13, 46, 205, 710, LMTStatus::Impossible},
                                   {14, 26, 220, 390, LMTStatus::Impossible},
                                   {14, 27, 226, 410, LMTStatus::Certain},
                                   {14, 25, 226, 374, LMTStatus::Impossible},
                                   {14, 24, 226, 360, LMTStatus::Impossible},
                                   {14, 21, 226, 327, LMTStatus::Impossible},
                                   {14, 20, 226, 318, LMTStatus::Impossible},
                                   {14, 22, 226, 338, LMTStatus::Impossible},
                                   {14, 15, 232, 249, LMTStatus::Certain},
                                   {14, 17, 232, 283, LMTStatus::Impossible},
                                   {14, 19, 232, 308, LMTStatus::Impossible},
                                   {14, 18, 232, 288, LMTStatus::Impossible},
                                   {14, 1, 232, 20, LMTStatus::Impossible},
                                   {14, 9, 232, 150, LMTStatus::Impossible},
                                   {14, 13, 233, 208, LMTStatus::Certain},
                                   {14, 16, 235, 261, LMTStatus::Certain},
                                   {14, 30, 235, 464, LMTStatus::Impossible},
                                   {14, 29, 220, 446, LMTStatus::Certain},
                                   {14, 28, 220, 425, LMTStatus::Impossible},
                                   {15, 26, 238, 389, LMTStatus::Impossible},
                                   {15, 27, 239, 409, LMTStatus::Certain},
                                   {15, 25, 241, 373, LMTStatus::Certain},
                                   {15, 24, 241, 359, LMTStatus::Impossible},
                                   {15, 21, 243, 326, LMTStatus::Possible},
                                   {15, 20, 243, 317, LMTStatus::Impossible},
                                   {15, 22, 244, 337, LMTStatus::Certain},
                                   {15, 17, 246, 284, LMTStatus::Certain},
                                   {15, 19, 246, 309, LMTStatus::Impossible},
                                   {15, 18, 247, 289, LMTStatus::Certain},
                                   {15, 13, 249, 209, LMTStatus::Certain},
                                   {15, 16, 249, 262, LMTStatus::Impossible},
                                   {15, 14, 238, 226, LMTStatus::Certain},
                                   {15, 29, 238, 445, LMTStatus::Impossible},
                                   {15, 28, 238, 424, LMTStatus::Impossible},
                                   {16, 28, 253, 427, LMTStatus::Impossible},
                                   {16, 29, 261, 447, LMTStatus::Certain},
                                   {16, 26, 261, 391, LMTStatus::Impossible},
                                   {16, 27, 261, 411, LMTStatus::Impossible},
                                   {16, 25, 261, 375, LMTStatus::Impossible},
                                   {16, 24, 261, 361, LMTStatus::Impossible},
                                   {16, 21, 261, 328, LMTStatus::Impossible},
                                   {16, 20, 261, 319, LMTStatus::Impossible},
                                   {16, 22, 261, 339, LMTStatus::Impossible},
                                   {16, 14, 266, 233, LMTStatus::Certain},
                                   {16, 15, 266, 248, LMTStatus::Impossible},
                                   {16, 17, 266, 282, LMTStatus::Impossible},
                                   {16, 19, 266, 307, LMTStatus::Impossible},
                                   {16, 18, 266, 287, LMTStatus::Impossible},
                                   {16, 13, 271, 207, LMTStatus::Certain},
                                   {16, 9, 271, 148, LMTStatus::Impossible},
                                   {16, 7, 271, 109, LMTStatus::Impossible},
                                   {16, 4, 271, 65, LMTStatus::Impossible},
                                   {16, 49, 271, 753, LMTStatus::Impossible},
                                   {16, 30, 253, 465, LMTStatus::Certain},
                                   {17, 27, 275, 408, LMTStatus::Impossible},
                                   {17, 25, 275, 372, LMTStatus::Impossible},
                                   {17, 24, 275, 358, LMTStatus::Impossible},
                                   {17, 22, 277, 336, LMTStatus::Certain},
                                   {17, 20, 277, 315, LMTStatus::Impossible},
                                   {17, 19, 279, 310, LMTStatus::Certain},
                                   {17, 1, 279, 22, LMTStatus::Impossible},
                                   {17, 18, 284, 290, LMTStatus::Certain},
                                   {17, 7, 284, 111, LMTStatus::Impossible},
                                   {17, 13, 284, 210, LMTStatus::Impossible},
                                   {17, 16, 284, 263, LMTStatus::Impossible},
                                   {17, 14, 284, 227, LMTStatus::Impossible},
                                   {17, 15, 275, 244, LMTStatus::Certain},
                                   {17, 29, 275, 444, LMTStatus::Impossible},
                                   {17, 28, 275, 423, LMTStatus::Impossible},
                                   {18, 16, 289, 265, LMTStatus::Impossible},
                                   {18, 14, 289, 229, LMTStatus::Impossible},
                                   {18, 15, 290, 246, LMTStatus::Certain},
                                   {18, 17, 291, 279, LMTStatus::Certain},
                                   {18, 19, 294, 306, LMTStatus::Certain},
                                   {18, 1, 294, 21, LMTStatus::Impossible},
                                   {18, 10, 294, 168, LMTStatus::Impossible},
                                   {18, 9, 297, 151, LMTStatus::Certain},
                                   {18, 12, 297, 194, LMTStatus::Impossible},
                                   {18, 7, 297, 112, LMTStatus::Impossible},
                                   {18, 13, 289, 211, LMTStatus::Certain},
                                   {19, 27, 300, 407, LMTStatus::Impossible},
                                   {19, 25, 300, 371, LMTStatus::Impossible},
                                   {19, 22, 301, 335, LMTStatus::Certain},
                                   {19, 20, 302, 314, LMTStatus::CH},
                                   {19, 2, 303, 35, LMTStatus::CH},
                                   {19, 1, 305, 23, LMTStatus::Certain},
                                   {19, 10, 305, 169, LMTStatus::Impossible},
                                   {19, 9, 306, 152, LMTStatus::Certain},
                                   {19, 18, 310, 291, LMTStatus::Certain},
                                   {19, 16, 310, 264, LMTStatus::Impossible},
                                   {19, 14, 310, 228, LMTStatus::Impossible},
                                   {19, 15, 310, 245, LMTStatus::Impossible},
                                   {19, 17, 300, 277, LMTStatus::Certain},
                                   {20, 24, 313, 356, LMTStatus::Impossible},
                                   {20, 23, 313, 344, LMTStatus::Impossible},
                                   {20, 21, 314, 324, LMTStatus::CH},
                                   {20, 19, 316, 301, LMTStatus::CH},
                                   {20, 17, 316, 276, LMTStatus::Impossible},
                                   {20, 22, 313, 334, LMTStatus::Certain},
                                   {20, 15, 313, 242, LMTStatus::Impossible},
                                   {20, 14, 313, 224, LMTStatus::Impossible},
                                   {20, 16, 313, 259, LMTStatus::Impossible},
                                   {20, 27, 313, 405, LMTStatus::Impossible},
                                   {20, 25, 313, 369, LMTStatus::Impossible},
                                   {21, 23, 323, 343, LMTStatus::Certain},
                                   {21, 31, 324, 478, LMTStatus::CH},
                                   {21, 20, 325, 313, LMTStatus::CH},
                                   {21, 22, 326, 333, LMTStatus::Certain},
                                   {21, 15, 330, 241, LMTStatus::Possible},
                                   {21, 14, 330, 223, LMTStatus::Impossible},
                                   {21, 16, 330, 258, LMTStatus::Impossible},
                                   {21, 27, 330, 404, LMTStatus::Impossible},
                                   {21, 25, 331, 368, LMTStatus::Possible},
                                   {21, 24, 322, 355, LMTStatus::Certain},
                                   {22, 24, 333, 357, LMTStatus::Possible},
                                   {22, 21, 334, 325, LMTStatus::Certain},
                                   {22, 20, 335, 316, LMTStatus::Certain},
                                   {22, 19, 336, 300, LMTStatus::Certain},
                                   {22, 17, 337, 275, LMTStatus::Certain},
                                   {22, 15, 341, 243, LMTStatus::Certain},
                                   {22, 14, 341, 225, LMTStatus::Impossible},
                                   {22, 16, 341, 260, LMTStatus::Impossible},
                                   {22, 27, 341, 406, LMTStatus::Impossible},
                                   {22, 25, 332, 370, LMTStatus::Possible},
                                   {23, 31, 343, 479, LMTStatus::Certain},
                                   {23, 21, 345, 322, LMTStatus::Certain},
                                   {23, 20, 345, 312, LMTStatus::Impossible},
                                   {23, 24, 346, 354, LMTStatus::Certain},
                                   {23, 25, 350, 366, LMTStatus::Possible},
                                   {23, 27, 350, 401, LMTStatus::Impossible},
                                   {23, 29, 350, 441, LMTStatus::Impossible},
                                   {23, 28, 350, 420, LMTStatus::Impossible},
                                   {23, 26, 351, 385, LMTStatus::Certain},
                                   {23, 35, 352, 528, LMTStatus::Certain},
                                   {23, 34, 342, 513, LMTStatus::Certain},
                                   {23, 33, 342, 501, LMTStatus::Impossible},
                                   {24, 23, 355, 345, LMTStatus::Certain},
                                   {24, 21, 357, 331, LMTStatus::Certain},
                                   {24, 20, 357, 311, LMTStatus::Impossible},
                                   {24, 22, 362, 332, LMTStatus::Possible},
                                   {24, 17, 362, 274, LMTStatus::Impossible},
                                   {24, 15, 362, 240, LMTStatus::Impossible},
                                   {24, 14, 362, 222, LMTStatus::Impossible},
                                   {24, 16, 362, 257, LMTStatus::Impossible},
                                   {24, 25, 363, 367, LMTStatus::Certain},
                                   {24, 27, 364, 402, LMTStatus::Possible},
                                   {24, 26, 354, 386, LMTStatus::Possible},
                                   {24, 35, 354, 529, LMTStatus::Impossible},
                                   {25, 23, 367, 346, LMTStatus::Possible},
                                   {25, 24, 368, 362, LMTStatus::Certain},
                                   {25, 21, 370, 330, LMTStatus::Possible},
                                   {25, 20, 370, 321, LMTStatus::Impossible},
                                   {25, 22, 373, 341, LMTStatus::Possible},
                                   {25, 19, 373, 299, LMTStatus::Impossible},
                                   {25, 17, 373, 273, LMTStatus::Impossible},
                                   {25, 15, 376, 239, LMTStatus::Certain},
                                   {25, 14, 376, 221, LMTStatus::Impossible},
                                   {25, 16, 376, 256, LMTStatus::Impossible},
                                   {25, 27, 379, 403, LMTStatus::Certain},
                                   {25, 29, 379, 442, LMTStatus::Impossible},
                                   {25, 28, 379, 421, LMTStatus::Impossible},
                                   {25, 26, 366, 387, LMTStatus::Possible},
                                   {25, 35, 366, 530, LMTStatus::Impossible},
                                   {26, 35, 385, 532, LMTStatus::Certain},
                                   {26, 34, 385, 514, LMTStatus::Impossible},
                                   {26, 33, 385, 502, LMTStatus::Impossible},
                                   {26, 31, 385, 480, LMTStatus::Impossible},
                                   {26, 23, 386, 350, LMTStatus::Certain},
                                   {26, 24, 387, 364, LMTStatus::Possible},
                                   {26, 25, 388, 379, LMTStatus::Possible},
                                   {26, 27, 393, 415, LMTStatus::Certain},
                                   {26, 15, 393, 237, LMTStatus::Impossible},
                                   {26, 14, 393, 219, LMTStatus::Impossible},
                                   {26, 16, 393, 254, LMTStatus::Impossible},
                                   {26, 7, 393, 108, LMTStatus::Impossible},
                                   {26, 29, 394, 440, LMTStatus::Certain},
                                   {26, 28, 381, 419, LMTStatus::Certain},
                                   {26, 38, 381, 587, LMTStatus::Impossible},
                                   {26, 46, 381, 706, LMTStatus::Impossible},
                                   {26, 36, 381, 546, LMTStatus::Impossible},
                                   {26, 40, 381, 614, LMTStatus::Impossible},
                                   {26, 42, 381, 647, LMTStatus::Impossible},
                                   {26, 43, 381, 663, LMTStatus::Impossible},
                                   {27, 23, 402, 347, LMTStatus::Impossible},
                                   {27, 24, 403, 363, LMTStatus::Possible},
                                   {27, 25, 409, 376, LMTStatus::Certain},
                                   {27, 21, 409, 329, LMTStatus::Impossible},
                                   {27, 20, 409, 320, LMTStatus::Impossible},
                                   {27, 22, 409, 340, LMTStatus::Impossible},
                                   {27, 19, 409, 298, LMTStatus::Impossible},
                                   {27, 17, 409, 272, LMTStatus::Impossible},
                                   {27, 15, 410, 238, LMTStatus::Certain},
                                   {27, 14, 413, 220, LMTStatus::Certain},
                                   {27, 16, 413, 255, LMTStatus::Impossible},
                                   {27, 30, 413, 463, LMTStatus::Impossible},
                                   {27, 29, 415, 443, LMTStatus::Certain},
                                   {27, 28, 415, 422, LMTStatus::Impossible},
                                   {27, 26, 402, 388, LMTStatus::Certain},
                                   {27, 35, 402, 531, LMTStatus::Impossible},
                                   {28, 34, 418, 516, LMTStatus::Impossible},
                                   {28, 35, 419, 533, LMTStatus::Certain},
                                   {28, 26, 426, 394, LMTStatus::Certain},
                                   {28, 23, 426, 349, LMTStatus::Impossible},
                                   {28, 25, 426, 378, LMTStatus::Impossible},
                                   {28, 27, 426, 414, LMTStatus::Impossible},
                                   {28, 17, 426, 286, LMTStatus::Impossible},
                                   {28, 15, 426, 251, LMTStatus::Impossible},
                                   {28, 14, 426, 236, LMTStatus::Impossible},
                                   {28, 29, 429, 457, LMTStatus::Certain},
                                   {28, 16, 429, 252, LMTStatus::Impossible},
                                   {28, 7, 429, 106, LMTStatus::Impossible},
                                   {28, 30, 436, 461, LMTStatus::Certain},
                                   {28, 8, 436, 129, LMTStatus::Impossible},
                                   {28, 4, 436, 62, LMTStatus::Impossible},
                                   {28, 44, 436, 679, LMTStatus::Impossible},
                                   {28, 37, 436, 568, LMTStatus::Impossible},
                                   {28, 49, 436, 750, LMTStatus::Impossible},
                                   {28, 36, 436, 548, LMTStatus::Impossible},
                                   {28, 46, 437, 707, LMTStatus::Certain},
                                   {28, 40, 418, 615, LMTStatus::Certain},
                                   {28, 43, 418, 664, LMTStatus::Impossible},
                                   {28, 32, 418, 491, LMTStatus::Impossible},
                                   {29, 26, 443, 393, LMTStatus::Certain},
                                   {29, 23, 443, 348, LMTStatus::Impossible},
                                   {29, 25, 443, 377, LMTStatus::Impossible},
                                   {29, 27, 446, 413, LMTStatus::Certain},
                                   {29, 17, 446, 285, LMTStatus::Impossible},
                                   {29, 15, 446, 250, LMTStatus::Impossible},
                                   {29, 14, 447, 235, LMTStatus::Certain},
                                   {29, 16, 450, 253, LMTStatus::Certain},
                                   {29, 13, 450, 206, LMTStatus::Impossible},
                                   {29, 7, 450, 107, LMTStatus::Impossible},
                                   {29, 30, 457, 462, LMTStatus::Certain},
                                   {29, 8, 457, 130, LMTStatus::Impossible},
                                   {29, 4, 457, 63, LMTStatus::Impossible},
                                   {29, 37, 457, 569, LMTStatus::Impossible},
                                   {29, 49, 457, 751, LMTStatus::Impossible},
                                   {29, 36, 457, 549, LMTStatus::Impossible},
                                   {29, 46, 457, 708, LMTStatus::Impossible},
                                   {29, 28, 440, 426, LMTStatus::Certain},
                                   {30, 32, 461, 492, LMTStatus::Impossible},
                                   {30, 34, 461, 517, LMTStatus::Impossible},
                                   {30, 35, 461, 534, LMTStatus::Impossible},
                                   {30, 28, 462, 429, LMTStatus::Certain},
                                   {30, 29, 465, 450, LMTStatus::Certain},
                                   {30, 27, 465, 412, LMTStatus::Impossible},
                                   {30, 14, 465, 234, LMTStatus::Impossible},
                                   {30, 16, 466, 271, LMTStatus::Certain},
                                   {30, 13, 467, 205, LMTStatus::Certain},
                                   {30, 7, 474, 105, LMTStatus::Certain},
                                   {30, 12, 474, 191, LMTStatus::Impossible},
                                   {30, 10, 474, 164, LMTStatus::Impossible},
                                   {30, 8, 474, 131, LMTStatus::Impossible},
                                   {30, 4, 474, 64, LMTStatus::Impossible},
                                   {30, 3, 474, 45, LMTStatus::Impossible},
                                   {30, 37, 474, 570, LMTStatus::Impossible},
                                   {30, 49, 476, 752, LMTStatus::Certain},
                                   {30, 36, 476, 550, LMTStatus::Impossible},
                                   {30, 46, 461, 709, LMTStatus::Certain},
                                   {30, 40, 461, 616, LMTStatus::Impossible},
                                   {31, 21, 479, 323, LMTStatus::CH},
                                   {31, 23, 482, 342, LMTStatus::Certain},
                                   {31, 26, 482, 384, LMTStatus::Impossible},
                                   {31, 35, 482, 527, LMTStatus::Impossible},
                                   {31, 34, 483, 512, LMTStatus::Certain},
                                   {31, 33, 485, 500, LMTStatus::Certain},
                                   {31, 46, 485, 702, LMTStatus::Impossible},
                                   {31, 32, 486, 487, LMTStatus::Certain},
                                   {31, 47, 478, 721, LMTStatus::CH},
                                   {32, 31, 488, 485, LMTStatus::Certain},
                                   {32, 33, 489, 510, LMTStatus::Certain},
                                   {32, 34, 490, 524, LMTStatus::Certain},
                                   {32, 35, 494, 542, LMTStatus::Certain},
                                   {32, 28, 494, 439, LMTStatus::Impossible},
                                   {32, 30, 494, 458, LMTStatus::Impossible},
                                   {32, 46, 494, 701, LMTStatus::Impossible},
                                   {32, 40, 496, 610, LMTStatus::Certain},
                                   {32, 41, 496, 627, LMTStatus::Impossible},
                                   {32, 43, 497, 659, LMTStatus::Certain},
                                   {32, 42, 499, 642, LMTStatus::Certain},
                                   {32, 39, 499, 600, LMTStatus::Impossible},
                                   {32, 47, 487, 722, LMTStatus::Certain},
                                   {33, 31, 504, 483, LMTStatus::Certain},
                                   {33, 23, 504, 353, LMTStatus::Impossible},
                                   {33, 26, 504, 383, LMTStatus::Impossible},
                                   {33, 35, 504, 526, LMTStatus::Impossible},
                                   {33, 34, 510, 511, LMTStatus::Certain},
                                   {33, 46, 510, 703, LMTStatus::Impossible},
                                   {33, 40, 510, 611, LMTStatus::Impossible},
                                   {33, 38, 510, 582, LMTStatus::Impossible},
                                   {33, 42, 510, 644, LMTStatus::Impossible},
                                   {33, 43, 510, 660, LMTStatus::Impossible},
                                   {33, 32, 500, 488, LMTStatus::Certain},
                                   {34, 33, 512, 504, LMTStatus::Certain},
                                   {34, 31, 513, 482, LMTStatus::Certain},
                                   {34, 23, 515, 352, LMTStatus::Certain},
                                   {34, 26, 515, 382, LMTStatus::Impossible},
                                   {34, 35, 524, 525, LMTStatus::Certain},
                                   {34, 28, 524, 417, LMTStatus::Impossible},
                                   {34, 30, 524, 459, LMTStatus::Impossible},
                                   {34, 46, 524, 704, LMTStatus::Impossible},
                                   {34, 36, 524, 544, LMTStatus::Impossible},
                                   {34, 40, 524, 612, LMTStatus::Impossible},
                                   {34, 38, 524, 583, LMTStatus::Impossible},
                                   {34, 42, 524, 645, LMTStatus::Impossible},
                                   {34, 43, 524, 661, LMTStatus::Impossible},
                                   {34, 32, 511, 489, LMTStatus::Certain},
                                   {35, 34, 528, 515, LMTStatus::Certain},
                                   {35, 33, 528, 503, LMTStatus::Impossible},
                                   {35, 31, 528, 481, LMTStatus::Impossible},
                                   {35, 23, 532, 351, LMTStatus::Certain},
                                   {35, 24, 532, 365, LMTStatus::Impossible},
                                   {35, 25, 532, 380, LMTStatus::Impossible},
                                   {35, 27, 532, 416, LMTStatus::Impossible},
                                   {35, 26, 533, 381, LMTStatus::Certain},
                                   {35, 28, 539, 418, LMTStatus::Certain},
                                   {35, 30, 539, 460, LMTStatus::Impossible},
                                   {35, 7, 539, 104, LMTStatus::Impossible},
                                   {35, 49, 539, 749, LMTStatus::Impossible},
                                   {35, 46, 539, 705, LMTStatus::Impossible},
                                   {35, 36, 539, 545, LMTStatus::Impossible},
                                   {35, 40, 542, 613, LMTStatus::Certain},
                                   {35, 42, 542, 646, LMTStatus::Impossible},
                                   {35, 43, 542, 662, LMTStatus::Impossible},
                                   {35, 32, 525, 490, LMTStatus::Certain},
                                   {36, 40, 547, 619, LMTStatus::Certain},
                                   {36, 34, 547, 519, LMTStatus::Impossible},
                                   {36, 35, 547, 538, LMTStatus::Impossible},
                                   {36, 26, 547, 397, LMTStatus::Impossible},
                                   {36, 46, 553, 720, LMTStatus::Certain},
                                   {36, 28, 553, 435, LMTStatus::Impossible},
                                   {36, 29, 553, 455, LMTStatus::Impossible},
                                   {36, 30, 553, 475, LMTStatus::Impossible},
                                   {36, 7, 553, 123, LMTStatus::Impossible},
                                   {36, 12, 553, 189, LMTStatus::Impossible},
                                   {36, 49, 559, 746, LMTStatus::Certain},
                                   {36, 4, 559, 59, LMTStatus::Impossible},
                                   {36, 3, 559, 42, LMTStatus::Impossible},
                                   {36, 0, 559, 3, LMTStatus::Impossible},
                                   {36, 44, 559, 677, LMTStatus::Impossible},
                                   {36, 37, 559, 566, LMTStatus::Impossible},
                                   {36, 38, 543, 585, LMTStatus::Certain},
                                   {36, 41, 543, 629, LMTStatus::Impossible},
                                   {36, 42, 543, 649, LMTStatus::Impossible},
                                   {36, 43, 543, 666, LMTStatus::Impossible},
                                   {37, 42, 564, 653, LMTStatus::Impossible},
                                   {37, 41, 567, 638, LMTStatus::Certain},
                                   {37, 38, 567, 595, LMTStatus::Impossible},
                                   {37, 36, 567, 558, LMTStatus::Impossible},
                                   {37, 49, 575, 764, LMTStatus::Possible},
                                   {37, 28, 575, 433, LMTStatus::Impossible},
                                   {37, 29, 575, 453, LMTStatus::Impossible},
                                   {37, 30, 575, 473, LMTStatus::Impossible},
                                   {37, 9, 575, 158, LMTStatus::Impossible},
                                   {37, 12, 575, 204, LMTStatus::Impossible},
                                   {37, 4, 575, 78, LMTStatus::Impossible},
                                   {37, 3, 575, 39, LMTStatus::Impossible},
                                   {37, 0, 576, 0, LMTStatus::Possible},
                                   {37, 44, 578, 674, LMTStatus::Certain},
                                   {37, 48, 578, 736, LMTStatus::Impossible},
                                   {37, 39, 564, 604, LMTStatus::Certain},
                                   {37, 45, 564, 693, LMTStatus::Impossible},
                                   {37, 47, 564, 727, LMTStatus::Impossible},
                                   {38, 43, 584, 667, LMTStatus::Impossible},
                                   {38, 33, 584, 507, LMTStatus::Impossible},
                                   {38, 34, 584, 521, LMTStatus::Impossible},
                                   {38, 40, 585, 620, LMTStatus::Certain},
                                   {38, 36, 588, 559, LMTStatus::Certain},
                                   {38, 46, 588, 719, LMTStatus::Impossible},
                                   {38, 26, 588, 395, LMTStatus::Impossible},
                                   {38, 49, 597, 766, LMTStatus::Certain},
                                   {38, 7, 597, 121, LMTStatus::Impossible},
                                   {38, 8, 597, 126, LMTStatus::Impossible},
                                   {38, 4, 597, 58, LMTStatus::Impossible},
                                   {38, 3, 597, 41, LMTStatus::Impossible},
                                   {38, 0, 597, 2, LMTStatus::Impossible},
                                   {38, 44, 597, 676, LMTStatus::Impossible},
                                   {38, 37, 597, 565, LMTStatus::Impossible},
                                   {38, 39, 597, 603, LMTStatus::Impossible},
                                   {38, 41, 584, 630, LMTStatus::Certain},
                                   {38, 42, 584, 651, LMTStatus::Impossible},
                                   {39, 47, 601, 730, LMTStatus::Impossible},
                                   {39, 32, 601, 498, LMTStatus::Impossible},
                                   {39, 42, 602, 656, LMTStatus::Certain},
                                   {39, 41, 604, 639, LMTStatus::Certain},
                                   {39, 38, 604, 596, LMTStatus::Impossible},
                                   {39, 37, 606, 578, LMTStatus::Certain},
                                   {39, 49, 606, 763, LMTStatus::Impossible},
                                   {39, 44, 607, 687, LMTStatus::Certain},
                                   {39, 48, 608, 745, LMTStatus::Certain},
                                   {39, 5, 609, 80, LMTStatus::Certain},
                                   {39, 45, 601, 695, LMTStatus::Certain},
                                   {40, 32, 613, 494, LMTStatus::Certain},
                                   {40, 33, 613, 506, LMTStatus::Impossible},
                                   {40, 34, 613, 520, LMTStatus::Impossible},
                                   {40, 35, 615, 539, LMTStatus::Certain},
                                   {40, 26, 615, 398, LMTStatus::Impossible},
                                   {40, 28, 617, 437, LMTStatus::Certain},
                                   {40, 30, 617, 477, LMTStatus::Impossible},
                                   {40, 46, 619, 700, LMTStatus::Certain},
                                   {40, 49, 619, 747, LMTStatus::Impossible},
                                   {40, 36, 620, 543, LMTStatus::Certain},
                                   {40, 38, 621, 584, LMTStatus::Certain},
                                   {40, 41, 624, 628, LMTStatus::Certain},
                                   {40, 47, 624, 725, LMTStatus::Impossible},
                                   {40, 42, 624, 648, LMTStatus::Impossible},
                                   {40, 43, 610, 665, LMTStatus::Certain},
                                   {41, 42, 626, 652, LMTStatus::Certain},
                                   {41, 43, 628, 668, LMTStatus::Certain},
                                   {41, 32, 628, 495, LMTStatus::Impossible},
                                   {41, 40, 630, 621, LMTStatus::Certain},
                                   {41, 36, 630, 560, LMTStatus::Impossible},
                                   {41, 38, 631, 597, LMTStatus::Certain},
                                   {41, 49, 637, 765, LMTStatus::Certain},
                                   {41, 7, 637, 120, LMTStatus::Impossible},
                                   {41, 10, 637, 163, LMTStatus::Impossible},
                                   {41, 4, 637, 57, LMTStatus::Impossible},
                                   {41, 3, 637, 40, LMTStatus::Impossible},
                                   {41, 0, 637, 1, LMTStatus::Impossible},
                                   {41, 44, 638, 675, LMTStatus::Possible},
                                   {41, 37, 639, 564, LMTStatus::Certain},
                                   {41, 39, 625, 602, LMTStatus::Certain},
                                   {41, 45, 625, 692, LMTStatus::Impossible},
                                   {41, 47, 625, 726, LMTStatus::Impossible},
                                   {42, 32, 643, 497, LMTStatus::Certain},
                                   {42, 43, 652, 669, LMTStatus::Certain},
                                   {42, 33, 652, 508, LMTStatus::Impossible},
                                   {42, 34, 652, 522, LMTStatus::Impossible},
                                   {42, 35, 652, 540, LMTStatus::Impossible},
                                   {42, 26, 652, 399, LMTStatus::Impossible},
                                   {42, 40, 652, 623, LMTStatus::Impossible},
                                   {42, 36, 652, 561, LMTStatus::Impossible},
                                   {42, 49, 652, 767, LMTStatus::Impossible},
                                   {42, 38, 652, 598, LMTStatus::Impossible},
                                   {42, 41, 656, 625, LMTStatus::Certain},
                                   {42, 37, 656, 563, LMTStatus::Impossible},
                                   {42, 44, 656, 673, LMTStatus::Impossible},
                                   {42, 48, 656, 735, LMTStatus::Impossible},
                                   {42, 39, 657, 601, LMTStatus::Certain},
                                   {42, 45, 658, 691, LMTStatus::Certain},
                                   {42, 47, 642, 724, LMTStatus::Certain},
                                   {43, 32, 665, 496, LMTStatus::Certain},
                                   {43, 33, 665, 509, LMTStatus::Impossible},
                                   {43, 34, 665, 523, LMTStatus::Impossible},
                                   {43, 35, 665, 541, LMTStatus::Impossible},
                                   {43, 26, 665, 400, LMTStatus::Impossible},
                                   {43, 28, 665, 438, LMTStatus::Impossible},
                                   {43, 40, 668, 624, LMTStatus::Certain},
                                   {43, 36, 668, 562, LMTStatus::Impossible},
                                   {43, 38, 668, 581, LMTStatus::Impossible},
                                   {43, 41, 669, 626, LMTStatus::Certain},
                                   {43, 42, 659, 643, LMTStatus::Certain},
                                   {43, 45, 659, 690, LMTStatus::Impossible},
                                   {43, 47, 659, 723, LMTStatus::Impossible},
                                   {44, 47, 674, 728, LMTStatus::Impossible},
                                   {44, 42, 674, 654, LMTStatus::Impossible},
                                   {44, 37, 675, 576, LMTStatus::Certain},
                                   {44, 41, 678, 637, LMTStatus::Possible},
                                   {44, 38, 678, 594, LMTStatus::Impossible},
                                   {44, 36, 678, 557, LMTStatus::Impossible},
                                   {44, 49, 680, 762, LMTStatus::Possible},
                                   {44, 28, 680, 432, LMTStatus::Impossible},
                                   {44, 4, 681, 77, LMTStatus::Possible},
                                   {44, 0, 686, 16, LMTStatus::Certain},
                                   {44, 3, 686, 55, LMTStatus::Impossible},
                                   {44, 8, 686, 143, LMTStatus::Impossible},
                                   {44, 11, 686, 179, LMTStatus::Impossible},
                                   {44, 5, 686, 82, LMTStatus::Impossible},
                                   {44, 48, 687, 737, LMTStatus::Certain},
                                   {44, 39, 674, 606, LMTStatus::Certain},
                                   {44, 45, 674, 694, LMTStatus::Impossible},
                                   {45, 47, 691, 733, LMTStatus::CH},
                                   {45, 43, 691, 670, LMTStatus::Impossible},
                                   {45, 42, 695, 657, LMTStatus::Certain},
                                   {45, 41, 695, 640, LMTStatus::Impossible},
                                   {45, 37, 695, 579, LMTStatus::Impossible},
                                   {45, 44, 695, 688, LMTStatus::Impossible},
                                   {45, 39, 698, 609, LMTStatus::Certain},
                                   {45, 48, 698, 744, LMTStatus::Impossible},
                                   {45, 8, 698, 141, LMTStatus::Impossible},
                                   {45, 5, 699, 93, LMTStatus::Certain},
                                   {45, 11, 689, 175, LMTStatus::CH},
                                   {46, 40, 707, 617, LMTStatus::Certain},
                                   {46, 32, 707, 493, LMTStatus::Impossible},
                                   {46, 31, 707, 484, LMTStatus::Impossible},
                                   {46, 33, 707, 505, LMTStatus::Impossible},
                                   {46, 34, 707, 518, LMTStatus::Impossible},
                                   {46, 35, 707, 537, LMTStatus::Impossible},
                                   {46, 26, 707, 396, LMTStatus::Impossible},
                                   {46, 28, 709, 436, LMTStatus::Certain},
                                   {46, 29, 709, 456, LMTStatus::Impossible},
                                   {46, 30, 718, 476, LMTStatus::Certain},
                                   {46, 13, 718, 218, LMTStatus::Impossible},
                                   {46, 7, 718, 103, LMTStatus::Impossible},
                                   {46, 9, 718, 146, LMTStatus::Impossible},
                                   {46, 12, 718, 190, LMTStatus::Impossible},
                                   {46, 8, 718, 128, LMTStatus::Impossible},
                                   {46, 4, 718, 61, LMTStatus::Impossible},
                                   {46, 3, 718, 44, LMTStatus::Impossible},
                                   {46, 0, 718, 5, LMTStatus::Impossible},
                                   {46, 49, 720, 748, LMTStatus::Certain},
                                   {46, 38, 720, 586, LMTStatus::Impossible},
                                   {46, 36, 700, 547, LMTStatus::Certain},
                                   {47, 31, 722, 486, LMTStatus::CH},
                                   {47, 32, 724, 499, LMTStatus::Certain},
                                   {47, 43, 724, 671, LMTStatus::Impossible},
                                   {47, 42, 733, 658, LMTStatus::Certain},
                                   {47, 40, 733, 622, LMTStatus::Impossible},
                                   {47, 41, 733, 641, LMTStatus::Impossible},
                                   {47, 37, 733, 580, LMTStatus::Impossible},
                                   {47, 44, 733, 672, LMTStatus::Impossible},
                                   {47, 48, 733, 734, LMTStatus::Impossible},
                                   {47, 39, 733, 599, LMTStatus::Impossible},
                                   {47, 5, 733, 79, LMTStatus::Impossible},
                                   {47, 11, 733, 176, LMTStatus::Impossible},
                                   {47, 45, 721, 689, LMTStatus::CH},
                                   {48, 47, 737, 729, LMTStatus::Impossible},
                                   {48, 42, 737, 655, LMTStatus::Impossible},
                                   {48, 37, 737, 577, LMTStatus::Impossible},
                                   {48, 44, 739, 686, LMTStatus::Certain},
                                   {48, 49, 739, 761, LMTStatus::Impossible},
                                   {48, 0, 740, 15, LMTStatus::Certain},
                                   {48, 3, 743, 54, LMTStatus::Certain},
                                   {48, 8, 743, 142, LMTStatus::Impossible},
                                   {48, 11, 743, 178, LMTStatus::Impossible},
                                   {48, 5, 745, 81, LMTStatus::Certain},
                                   {48, 45, 745, 696, LMTStatus::Impossible},
                                   {48, 39, 737, 607, LMTStatus::Certain},
                                   {49, 36, 748, 553, LMTStatus::Certain},
                                   {49, 40, 748, 618, LMTStatus::Impossible},
                                   {49, 46, 752, 718, LMTStatus::Certain},
                                   {49, 35, 752, 536, LMTStatus::Impossible},
                                   {49, 28, 752, 434, LMTStatus::Impossible},
                                   {49, 29, 752, 454, LMTStatus::Impossible},
                                   {49, 30, 754, 474, LMTStatus::Certain},
                                   {49, 16, 754, 270, LMTStatus::Impossible},
                                   {49, 7, 758, 122, LMTStatus::Certain},
                                   {49, 9, 758, 145, LMTStatus::Impossible},
                                   {49, 12, 758, 188, LMTStatus::Impossible},
                                   {49, 8, 758, 127, LMTStatus::Impossible},
                                   {49, 4, 760, 60, LMTStatus::Certain},
                                   {49, 3, 760, 43, LMTStatus::Impossible},
                                   {49, 0, 762, 4, LMTStatus::Possible},
                                   {49, 48, 762, 738, LMTStatus::Impossible},
                                   {49, 44, 764, 678, LMTStatus::Possible},
                                   {49, 39, 764, 605, LMTStatus::Impossible},
                                   {49, 37, 765, 567, LMTStatus::Possible},
                                   {49, 41, 766, 631, LMTStatus::Certain},
                                   {49, 38, 746, 588, LMTStatus::Certain},
                                   {49, 42, 746, 650, LMTStatus::Impossible}};

    return {ordered_points, halfedges_from_info<Halfedge>(ordered_points, info)};
}

template<typename OutputIterator, typename Exclude>
void generate_complete(OutputIterator out, const std::vector<std::size_t> &indices, const Exclude &exclude) {
    std::set<std::pair<std::size_t, std::size_t>> exclude_set;
    for(auto e : exclude) {
        exclude_set.emplace(e.first, e.second);
        exclude_set.emplace(e.second, e.first);
    }
    for(std::size_t i = 0; i < indices.size(); ++i) {
        for(std::size_t j = i + 1; j < indices.size(); ++j) {
            auto pair = std::make_pair(indices[i], indices[j]);
            if(exclude_set.count(pair))
                continue;
            *out = pair;
            ++out;
        }
    }
}

std::pair<Points, std::vector<Halfedge>> test_case2() {
    Points points = fake_lmt_test1_points();
    std::vector<std::pair<std::size_t, std::size_t>> skeleton{{4, 7}, {7, 9}, {9, 10},  {10, 11}, {11, 8}, {8, 6},
                                                              {6, 5}, {5, 4}, {14, 15}, {2, 17},  {12, 17}};
    std::vector<std::pair<std::size_t, std::size_t>> possible{
        {0, 4},   {0, 5},   {0, 6},  {0, 7},   {1, 4},   {1, 5},   {1, 6},  {1, 8},  {2, 4},
        {2, 5},   {2, 7},   {2, 9},  {3, 5},   {3, 6},   {3, 8},   {3, 11}, {12, 4}, {12, 7},
        {12, 9},  {12, 10}, {13, 6}, {13, 8},  {13, 10}, {13, 11}, {14, 7}, {14, 9}, {14, 10},
        {14, 11}, {15, 8},  {15, 9}, {15, 10}, {15, 11}, {17, 4},  {17, 7}, {17, 9}};
    std::vector<std::pair<std::size_t, std::size_t>> exclude{{4, 9}, {9, 11}, {4, 6}, {6, 11}};
    exclude.insert(exclude.end(), skeleton.begin(), skeleton.end());
    generate_complete(std::back_inserter(possible), {4, 5, 6, 7, 8, 9, 10, 11}, exclude);
    std::vector<Halfedge> halfedges = halfedges_from_graph(points, skeleton, possible);
    return {std::move(points), std::move(halfedges)};
}

TEST_CASE("Face collector - Known skeleton") {
    using LMTStatus = mwt::LMTStatus;
    CGAL::Protect_FPU_rounding rounder;
    auto [ordered_points, lmt_halfedges] = test_case1();
    std::vector<Halfedge *> lmt_skeleton;
    for(std::size_t i = 0; i < lmt_halfedges.size(); ++i) {
        if(lmt_halfedges[i].status == LMTStatus::Certain || lmt_halfedges[i].status == LMTStatus::CH) {
            lmt_skeleton.push_back(&lmt_halfedges[i]);
        }
    }
    mwt::Face_collector<Skeleton> face_collector(lmt_skeleton.begin(), lmt_skeleton.end());
    mwt::test::Face_collector_tester::verify_initial(face_collector);
    CHECK(face_collector.next());
    auto &f1 = face_collector.get_current_face_reference();
    CHECK(f1.is_simple());
    CHECK(f1.has_inner_edges());
    CHECK(!f1.inner_edges.empty());
    CHECK(f1.isolated_vertices.empty());
    CHECK(f1.hole_boundaries.empty());
    CHECK((f1.boundary.size() == 5 || f1.boundary.size() == 6));
    CHECK(face_collector.next());
    CHECK(f1.has_inner_edges());
    CHECK(!f1.inner_edges.empty());
    CHECK(f1.is_simple());
    CHECK(face_collector.next());
    CHECK(f1.has_inner_edges());
    CHECK(!f1.inner_edges.empty());
    CHECK(f1.is_simple());
    CHECK(!face_collector.next());
}

TEST_CASE("Face collector - Fixed test case 1") {
    CGAL::Protect_FPU_rounding rounder;
    Points points{
        {-263.531, -10.4604}, {-217.386, 69.4662},  {-125.095, -90.3853}, {-217.385, -90.3853}, {-125.095, 229.319},
        {-78.9513, 309.244},  {13.3396, 309.244},   {59.4843, 229.319},   {-263.531, 149.393},  {-217.386, 229.319},
        {151.777, 229.319},   {197.921, 149.393},   {151.777, 69.4662},   {197.921, -10.4604},  {151.777, -90.3853},
        {59.4855, -90.3853},  {13.3396, -170.31},   {-78.9502, -170.31},  {105.63, 109.669},    {109.828, 112.751},
        {101.432, 112.751},   {71.2287, 169.253},   {71.7987, 164.076},   {75.9967, 171.347},   {71.2287, 129.53},
        {75.9967, 127.436},   {71.7987, 134.707},   {140.03, 169.253},    {135.263, 171.347},   {139.46, 164.076},
        {140.03, 129.53},     {139.46, 134.707},    {135.263, 127.436},   {105.63, 189.115},    {101.432, 186.033},
        {109.828, 186.033},   {-32.8063, 189.597},  {-28.6083, 192.679},  {-37.0043, 192.679},  {-67.2069, 249.179},
        {-66.6368, 244.003},  {-62.4388, 251.276},  {-67.2069, 209.457},  {-62.4388, 207.362},  {-66.637, 214.634},
        {1.59427, 249.179},   {-3.17383, 251.276},  {1.02429, 244.003},   {1.59427, 209.457},   {1.02429, 214.634},
        {-3.17383, 207.362},  {-32.8063, 269.041},  {-37.0043, 265.959},  {-28.6083, 265.959},  {-171.243, 109.671},
        {-167.044, 112.753},  {-175.44, 112.753},   {-205.643, 169.255},  {-205.072, 164.078},  {-200.874, 171.35},
        {-205.643, 129.532},  {-200.874, 127.437},  {-205.072, 134.709},  {-136.841, 169.255},  {-141.61, 171.35},
        {-137.411, 164.078},  {-136.841, 129.532},  {-137.411, 134.709},  {-141.61, 127.437},   {-171.243, 189.117},
        {-175.44, 186.035},   {-167.044, 186.035},  {-171.243, -50.1824}, {-167.044, -47.1006}, {-175.44, -47.1006},
        {-205.643, 9.40096},  {-205.072, 4.22446},  {-200.874, 11.4957},  {-205.643, -30.3214}, {-200.874, -32.416},
        {-205.072, -25.1449}, {-136.841, 9.40105},  {-141.61, 11.4957},   {-137.411, 4.22446},  {-136.841, -30.3212},
        {-137.411, -25.1449}, {-141.61, -32.416},   {-171.243, 29.2622},  {-175.44, 26.1803},   {-167.044, 26.1803},
        {-32.8052, -130.109}, {-28.6072, -127.028}, {-37.0032, -127.028}, {-67.2058, -70.5253}, {-66.6357, -75.7023},
        {-62.4377, -68.4314}, {-67.2057, -110.247}, {-62.4377, -112.344}, {-66.6358, -105.071}, {1.59537, -70.5253},
        {-3.17273, -68.4314}, {1.02539, -75.7023},  {1.59537, -110.247},  {1.02539, -105.071},  {-3.17273, -112.344},
        {-32.8052, -50.6644}, {-37.0032, -53.7464}, {-28.6072, -53.7464}, {105.63, -50.1843},   {109.828, -47.1026},
        {101.432, -47.1026},  {71.2287, 9.39905},   {71.7987, 4.22275},   {75.9967, 11.4938},   {71.2287, -30.3231},
        {75.9967, -32.418},   {71.7987, -25.1467},  {140.03, 9.39935},    {135.263, 11.4938},   {139.46, 4.22275},
        {140.03, -30.3231},   {139.46, -25.1467},   {135.263, -32.418},   {105.63, 29.2605},    {101.432, 26.1785},
        {109.828, 26.1785},   {-125.095, 69.4662},  {-78.9513, 149.393},  {13.3396, -10.4604},  {-78.9502, -10.4604},
        {13.3396, 149.393},   {59.4855, 69.4662},   {-32.8052, 29.7438},  {-28.6072, 32.8257},  {-37.0032, 32.8257},
        {-67.2058, 89.3273},  {-66.6357, 84.1508},  {-62.4377, 91.4219},  {-67.2057, 49.6049},  {-62.4377, 47.5103},
        {-66.6358, 54.7817},  {1.59537, 89.3274},   {-3.17273, 91.4219},  {1.02539, 84.1508},   {1.59537, 49.605},
        {1.02539, 54.7817},   {-3.17273, 47.5103},  {-32.8052, 109.188},  {-37.0032, 106.107},  {-28.6072, 106.107},
        {-355.823, 149.393},  {-309.677, 29.7438},  {-305.479, 32.8257},  {-313.875, 32.8257},  {-344.078, 89.3273},
        {-343.508, 84.1508},  {-339.31, 91.4219},   {-344.078, 49.6049},  {-339.31, 47.5103},   {-343.508, 54.7817},
        {-275.276, 89.3274},  {-280.045, 91.4219},  {-275.846, 84.1508},  {-275.276, 49.605},   {-275.846, 54.7817},
        {-280.045, 47.5103},  {-309.677, 109.188},  {-313.875, 106.107},  {-305.479, 106.107},  {290.212, -10.4604},
        {290.212, 149.393},   {244.067, 29.7438},   {248.265, 32.8257},   {239.869, 32.8257},   {209.667, 89.3273},
        {210.237, 84.1508},   {214.435, 91.4219},   {209.667, 49.6049},   {214.435, 47.5103},   {210.237, 54.7817},
        {278.468, 89.3274},   {273.7, 91.4219},     {277.898, 84.1508},   {278.468, 49.605},    {277.898, 54.7817},
        {273.7, 47.5103},     {244.067, 109.188},   {239.869, 106.107},   {248.265, 106.107},   {197.92, 309.246},
        {105.63, 269.524},    {109.828, 272.606},   {101.432, 272.606},   {71.2287, 289.385},   {75.9967, 287.29},
        {71.7987, 294.562},   {140.03, 289.385},    {139.46, 294.562},    {135.262, 287.29},    {-171.24, 269.524},
        {-167.042, 272.606},  {-175.438, 272.606},  {-205.641, 289.385},  {-200.873, 287.29},   {-205.071, 294.562},
        {-136.84, 289.385},   {-137.41, 294.562},   {-141.608, 287.29},   {-355.822, 309.247},  {-263.531, 309.247},
        {-309.676, 189.598},  {-305.478, 192.68},   {-313.874, 192.68},   {-344.077, 249.181},  {-343.507, 244.005},
        {-339.309, 251.276},  {-344.077, 209.459},  {-339.309, 207.364},  {-343.507, 214.636},  {-275.276, 249.181},
        {-280.044, 251.276},  {-275.846, 244.005},  {-275.276, 209.459},  {-275.846, 214.636},  {-280.044, 207.364},
        {-309.676, 269.042},  {-313.874, 265.961},  {-305.478, 265.961},  {-355.822, -10.4604}, {-355.821, -170.313},
        {-309.676, -130.109}, {-305.478, -127.027}, {-313.874, -127.027}, {-344.077, -70.5253}, {-343.507, -75.7023},
        {-339.309, -68.4314}, {-344.077, -110.248}, {-339.309, -112.342}, {-343.507, -105.071}, {-275.276, -70.5253},
        {-280.044, -68.4314}, {-275.846, -75.7023}, {-275.276, -110.248}, {-275.846, -105.071}, {-280.044, -112.342},
        {-309.676, -50.6644}, {-313.874, -53.7464}, {-305.478, -53.7464}, {-263.531, -170.312}, {-205.641, -150.451},
        {-205.071, -155.628}, {-200.873, -148.357}, {-136.84, -150.451},  {-141.608, -148.357}, {-137.41, -155.628},
        {-171.24, -130.59},   {-175.438, -133.672}, {-167.042, -133.672}, {197.923, -170.312},  {71.2307, -150.451},
        {71.8007, -155.628},  {75.9987, -148.356},  {140.033, -150.451},  {135.264, -148.356},  {139.463, -155.628},
        {105.632, -130.59},   {101.434, -133.672},  {109.83, -133.672},   {290.212, -170.313},  {244.067, -130.109},
        {248.265, -127.027},  {239.869, -127.027},  {209.666, -70.5253},  {210.236, -75.7023},  {214.434, -68.4314},
        {209.666, -110.248},  {214.434, -112.342},  {210.236, -105.071},  {278.467, -70.5253},  {273.699, -68.4314},
        {277.897, -75.7023},  {278.467, -110.248},  {277.897, -105.071},  {273.699, -112.342},  {244.067, -50.6644},
        {239.869, -53.7464},  {248.265, -53.7464},  {290.212, 309.246},   {244.067, 189.597},   {248.265, 192.679},
        {239.869, 192.679},   {209.667, 249.181},   {210.237, 244.004},   {214.435, 251.275},   {209.667, 209.458},
        {214.435, 207.364},   {210.237, 214.635},   {278.468, 249.181},   {273.7, 251.275},     {277.898, 244.004},
        {278.468, 209.458},   {277.898, 214.635},   {273.7, 207.364},     {244.067, 269.042},   {239.869, 265.96},
        {248.265, 265.96}};
    Tree tree(points.begin(), points.end());
    Filter filter(&tree);
    filter.compute_remaining_edges();
    mwt::LMTSkeleton<Traits, Tree> skeleton_builder(&tree, std::move(filter.get_edges()));
    skeleton_builder.lmt_main_loop();
    skeleton_builder.advanced_lmt_loop();
    mwt::Face_collector<Skeleton> face_collector(skeleton_builder);
    std::size_t count_simple = 0, count_nonsimple = 0;
    while(face_collector.next()) {
        auto &count = (face_collector.get_current_face_reference().is_simple() ? count_simple : count_nonsimple);
        ++count;
    }
    CHECK(count_simple == 23);
    CHECK(count_nonsimple == 1);
}

TEST_CASE("Face collector - Fixed test case 2") {
    CGAL::Protect_FPU_rounding rounder;
    std::pair<Points, std::vector<Halfedge>> data = test_case2();
    auto &points = data.first;
    auto &halfedges = data.second;
    std::vector<Halfedge *> skeleton;
    for(Halfedge &he : halfedges) {
        if(he.status != mwt::LMTStatus::Possible)
            skeleton.push_back(&he);
    }
    std::mt19937_64 rng(151617);
    for(std::size_t i = 0; i < 100; ++i) {
        std::shuffle(skeleton.begin(), skeleton.end(), rng);
        mwt::Face_collector<Skeleton> face_collector(skeleton.begin(), skeleton.end());
        bool got_simple = false, got_non_simple = false;
        for(int i = 0; i < 2; ++i) {
            CHECK(face_collector.next());
            auto &f1 = face_collector.get_current_face_reference();
            CHECK(f1.has_inner_edges());
            CHECK(f1.isolated_vertices.empty());
            if(f1.is_simple()) {
                got_simple = true;
                CHECK(f1.isolated_vertices.empty());
                CHECK(f1.hole_boundaries.empty());
                CHECK(f1.boundary.size() == 8);
                CHECK(f1.inner_edges.size() == 32);
            } else {
                got_non_simple = true;
                CHECK(f1.boundary.size() == 9);
                CHECK(f1.hole_boundaries.size() == 8);
                CHECK(f1.inner_edges.size() == 70);
            }
        }
        CHECK(!face_collector.next());
        CHECK(got_simple);
        CHECK(got_non_simple);
        for(Halfedge &he : halfedges) {
            he.flagged = false;
        }
    }
}

TEST_CASE("Face collector - Fixed test case 3") {
    CGAL::Protect_FPU_rounding rounder;
    Points points{{128, 720}, {128, 352}, {160, 672}, {160, 576}, {192, 480}, {208, 432}, {224, 400}, {240, 480},
                  {272, 464}, {368, 448}, {384, 704}, {352, 640}, {528, 544}, {544, 720}, {544, 352}, {256, 576},
                  {256, 544}, {272, 528}, {288, 544}, {320, 528}, {336, 496}, {368, 560}};

    std::vector<std::pair<std::size_t, std::size_t>> skeleton_edges{
        {0, 2},   {0, 3},   {0, 10},  {1, 3},   {1, 4},   {1, 5},   {1, 6},   {1, 9},   {2, 3},   {2, 10},
        {3, 4},   {4, 5},   {5, 6},   {6, 9},   {7, 8},   {9, 12},  {9, 14},  {10, 11}, {10, 12}, {10, 13},
        {12, 13}, {12, 14}, {15, 16}, {15, 21}, {16, 17}, {17, 18}, {17, 19}, {19, 20}, {19, 21}};

    std::vector<std::pair<std::size_t, std::size_t>> possible_edges{
        {2, 6},   {2, 7},   {2, 15},  {2, 21},  {2, 12},  {2, 11},  {3, 7},   {3, 15},  {3, 10},  {4, 8},   {4, 7},
        {4, 16},  {4, 15},  {5, 9},   {5, 8},   {5, 7},   {5, 16},  {5, 15},  {6, 12},  {6, 20},  {6, 8},   {6, 7},
        {7, 20},  {7, 19},  {7, 17},  {7, 16},  {7, 15},  {8, 9},   {8, 12},  {8, 20},  {8, 16},  {9, 10},  {9, 21},
        {9, 20},  {9, 17},  {10, 21}, {10, 15}, {11, 12}, {11, 21}, {11, 15}, {12, 21}, {12, 19}, {12, 20}, {20, 21},
        {15, 17}, {15, 18}, {15, 19}, {16, 18}, {16, 21}, {17, 21}, {18, 19}, {18, 21}};

    std::vector<std::size_t> nonsimple_face_boundary{2, 3, 4, 5, 6, 9, 12, 10, 11, 10};
    std::vector<std::size_t> nonsimple_face_holes{7, 8, 15, 21, 19, 20, 19, 17, 16};
    std::vector<std::size_t> simple_face_boundary{15, 16, 17, 18, 17, 19, 21};

    CHECK(points.size() == 22);

    std::vector<Halfedge> halfedges = halfedges_from_graph(points, skeleton_edges, possible_edges);
    std::vector<Halfedge *> skeleton;
    for(Halfedge &he : halfedges) {
        if(he.status != mwt::LMTStatus::Possible)
            skeleton.push_back(&he);
    }
    mwt::Face_collector<Skeleton> face_collector(skeleton.begin(), skeleton.end());
    CHECK(face_collector.next());
    auto &f1 = face_collector.get_current_face_reference();

    auto point_index = [&](const Point &p) { return std::find(points.begin(), points.end(), p) - points.begin(); };

    std::vector<std::size_t> boundary_order;
    std::vector<std::size_t> hole_order;
    auto collect_and_check_boundaries = [&]() {
        boundary_order.clear();
        hole_order.clear();
        for(auto *he : f1.boundary) {
            boundary_order.push_back(point_index(he->source()));
        }
        for(auto *he : f1.hole_boundaries) {
            hole_order.push_back(point_index(he->source()));
        }
        if(f1.is_simple()) {
            CHECK(same_cyclic_order(boundary_order.begin(), boundary_order.end(), simple_face_boundary.begin(),
                                    simple_face_boundary.end()));
            CHECK(f1.hole_boundaries.empty());
        } else {
            CHECK(same_cyclic_order(boundary_order.begin(), boundary_order.end(), nonsimple_face_boundary.begin(),
                                    nonsimple_face_boundary.end()));
            CHECK(same_cyclic_order(hole_order.begin(), hole_order.end(), nonsimple_face_holes.begin(),
                                    nonsimple_face_holes.end()));
            mwt::Face_analyzer<Skeleton> analyzer;
            analyzer.analyze_face(f1);
            const auto &edge_info = analyzer.edge_info();
            const auto &vertex_info = analyzer.empty_triangles();
            for(auto *he : f1.boundary) {
                CHECK(edge_info.is_defined(he));
            }
            for(auto *he : f1.hole_boundaries) {
                CHECK(edge_info.is_defined(he));
            }
            for(auto *he : f1.inner_edges) {
                CHECK(edge_info.is_defined(he));
            }
        }
        CHECK(f1.isolated_vertices.empty());
    };

    collect_and_check_boundaries();
    CHECK(face_collector.next());
    auto &f2 = face_collector.get_current_face_reference();
    CHECK(&f1 == &f2);
    collect_and_check_boundaries();
}
