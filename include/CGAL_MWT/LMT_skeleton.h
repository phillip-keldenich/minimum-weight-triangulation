#ifndef CGAL_MWT_LMT_SKELETON_H_INCLUDED_
#define CGAL_MWT_LMT_SKELETON_H_INCLUDED_

#include "LMT_halfedge.h"
#include <CGAL/number_utils.h>
#include <boost/range/iterator_range_core.hpp>
#include <cmath>
#include <cstdint>
#include <tbb/combinable.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_invoke.h>

namespace mwt {

namespace test {
class LMTSkeletonTester;
}

namespace detail {

template<typename PointType, typename HalfedgeIterator>
void sort_halfedges_around_point(PointType origin, HalfedgeIterator begin, HalfedgeIterator end) {
    using Halfedge = typename std::iterator_traits<HalfedgeIterator>::value_type;
    using Traits = typename Halfedge::Traits;
    using Left_turn_2 = typename Traits::Left_turn_2;

    // partition into quadrants first and order based on that;
    // then use the much simpler left_turn_2 predicate instead of
    // one slow complex comparator
    auto is_above = [&origin](const Halfedge &he) -> bool { return origin.y() < he.target().y(); };
    auto is_right = [&origin](const Halfedge &he) -> bool { return origin.x() < he.target().x(); };
    auto is_left = [&origin](const Halfedge &he) -> bool { return origin.x() > he.target().x(); };
    auto is_left_turn = [&origin](const Halfedge &he1, const Halfedge &he2) -> bool {
        return Left_turn_2{}(origin, he1.target(), he2.target());
    };
    auto m2 = std::partition(begin, end, is_above);
    auto m1 = std::partition(begin, m2, is_right);
    auto m3 = std::partition(m2, end, is_left);
    std::sort(begin, m1, is_left_turn);
    std::sort(m1, m2, is_left_turn);
    std::sort(m2, m3, is_left_turn);
    std::sort(m3, end, is_left_turn);

    // fix twin pointers (use the stored and unchanged twins),
    // and initialize next pointers
    for(auto it = begin; it != end; ++it) {
        it->twin()->set_twin(&*it);
        it->next = &*it + 1;
    }
    // fix the past-the-end next pointer to correctly point to the begin
    if(begin != end)
        std::prev(end)->next = &*begin;
}

template<typename HalfedgeIterator> std::size_t init_start_pointers(HalfedgeIterator begin, HalfedgeIterator end) {
    auto e = &*begin;
    std::size_t convex_hull_size = 0;
    for(auto edge = begin; edge != end; ++edge) {
        while(!edge->has_on_right_side(e->next->target())) {
            e = e->next;
            if(e->next == &*edge) {
                ++convex_hull_size;
                edge->status = edge->twin()->status = LMTStatus::CH;
                edge->flagged = edge->twin()->flagged = false;
                break;
            }
        }
        if(e->next != &*edge) {
            edge->twin()->j_start = e;
        }
    }
    return convex_hull_size;
}

} // namespace detail

template<typename Traits_, typename Tree_> class LMTSkeleton {
  public:
    using Tree = Tree_;
    using Traits = Traits_;
    using Kernel = typename Traits::Kernel;
    using Halfedge = LMTHalfedge<Traits>;
    using Point_2 = typename Traits::Point_2;
    using Segment_2 = typename Traits::Segment_2;
    using Left_turn_2 = typename Traits::Left_turn_2;
    using PointIterator = typename Tree::Iterator;

    /**
     * Construct an LMT skeleton builder from a tree (acting as
     * a set of points) and a set of edges that remain after
     * the diamond filter test.
     * If the edges are moved into this constructor, their container
     * is cleared after construction of this LMTSkeleton object.
     */
    template<typename DiamondFilterEdges>
    LMTSkeleton(const Tree *tree, DiamondFilterEdges &&edges) : m_tree(tree), m_offsets(m_tree->size() + 1, 0) {
        p_halfedges_from_diamond_filter(std::forward<DiamondFilterEdges>(edges));
        clear_if_move(std::forward<DiamondFilterEdges>(edges));
        CGAL_expensive_postcondition(p_validate_initial_invariants());
    }

    /**
     * Get the number of edges/points on the convex hull.
     */
    std::size_t num_convex_hull_edges() const noexcept { return m_convex_hull_size; }

    /**
     * Get the number of edges in the triangulation, based on
     * the number of points and the convex hull size.
     */
    std::size_t expected_total_triangulation_edges() const noexcept {
        return 3 * num_points() - 3 - num_convex_hull_edges();
    }

    /**
     * Run the LMT main loop, i.e., the first stage of the LMT
     * skeleton construction.
     * If run_parallel is true, the algorithm is parallelized
     * using the TBB library.
     */
    void lmt_main_loop(bool run_parallel = true) {
        if(run_parallel) {
            p_lmt_main_loop_parallel();
            p_mark_certain_edges_parallel();
        } else {
            p_lmt_main_loop_sequential();
            p_mark_certain_edges_sequential();
        }
    }

    /**
     * Run the second stage of the LMT skeleton construction.
     */
    void advanced_lmt_loop() {
        std::vector<Halfedge *> stack;
        stack.reserve(m_candidates.size());
        for(Halfedge *e : m_candidates) {
            stack.push_back(e);
            e->flagged = e->twin()->flagged = true;
        }

        p_advanced_lmt_loop(stack);

        for(Halfedge *e : m_candidates) {
            if(e->status == LMTStatus::Certain) {
                m_skeleton.push_back(e);
            } else if(e->status == LMTStatus::Possible && !e->has_crossing_edges()) {
                e->status = e->twin()->status = LMTStatus::Certain;
                m_skeleton.push_back(e);
            }
        }

        m_candidates.erase(std::remove_if(m_candidates.begin(), m_candidates.end(),
                                          [](const Halfedge *e) { return e->status != LMTStatus::Possible; }),
                           m_candidates.end());
    }

    std::size_t num_points() const noexcept { return m_offsets.size() - 1; }

    std::vector<Halfedge> &get_all_halfedges() noexcept { return m_lmt_halfedges; }

    const std::vector<Halfedge> &get_all_halfedges() const noexcept { return m_lmt_halfedges; }

    void print_as_json(std::ostream &out) {
        auto to_i = [](const auto &x) { return CGAL::to_double(x); };
        auto print_point = [&](const auto &p) { out << "[" << to_i(p.x()) << ", " << to_i(p.y()) << "]"; };
        out << std::setprecision(17) << "{ \"points\": [";
        bool first = true;
        for(auto it = m_tree->points_begin(); it != m_tree->points_end(); ++it) {
            if(first) {
                first = false;
            } else {
                out << ", ";
            }
            print_point(*it);
        }
        out << "],\n\"edges\": [";
        first = true;
        for(const auto &e : m_lmt_halfedges) {
            if(!e.is_primary())
                continue;
            if(e.status == LMTStatus::Impossible)
                continue;
            if(first) {
                first = false;
            } else {
                out << ", ";
            }
            out << "{ \"source\": ";
            print_point(*e.source_handle());
            out << ", \"target\": ";
            print_point(*e.target_handle());
            out << ", \"status\": \"" << e.status << "\" }";
        }
        out << "],\n\"all_halfedges\": [";
        first = true;
        for(const auto &e : m_lmt_halfedges) {
            if(first) {
                first = false;
            } else {
                out << ", ";
            }
            out << "{ \"id\": " << (&e - m_lmt_halfedges.data())
                << ", \"source_id\": " << (e.source_handle() - &*m_tree->points_begin())
                << ", \"target_id\": " << (e.target_handle() - &*m_tree->points_begin()) << ", \"status\": \""
                << e.status << "\", \"twin_id\": " << (e.twin() - m_lmt_halfedges.data())
                << ", \"next_id\": " << (e.next_edge() - m_lmt_halfedges.data()) << " }";
        }
        out << "] }\n";
    }

    std::vector<Halfedge *> &candidates() noexcept { return m_candidates; }
    std::vector<Halfedge *> &skeleton() noexcept { return m_skeleton; }

    std::size_t total_halfedges() const noexcept { return m_lmt_halfedges.size(); }
    Halfedge *halfedges_begin() noexcept { return m_lmt_halfedges.data(); }
    Halfedge *halfedges_end() noexcept { return m_lmt_halfedges.data() + m_lmt_halfedges.size(); }
    const Halfedge *halfedges_begin() const noexcept { return m_lmt_halfedges.data(); }
    const Halfedge *halfedges_end() const noexcept { return m_lmt_halfedges.data() + m_lmt_halfedges.size(); }

    const Tree &get_tree() const noexcept { return *m_tree; }

    boost::iterator_range<Halfedge *> halfedges_around(const Point_2 *point_in_tree) noexcept {
        std::size_t point_index = p_point_index_of(point_in_tree);
        return {&m_lmt_halfedges[m_offsets[point_index]], &m_lmt_halfedges[m_offsets[point_index + 1]]};
    }

    boost::iterator_range<const Halfedge *> halfedges_around(const Point_2 *point_in_tree) const noexcept {
        std::size_t point_index = p_point_index_of(point_in_tree);
        return {&m_lmt_halfedges[m_offsets[point_index]], &m_lmt_halfedges[m_offsets[point_index + 1]]};
    }

  private:
    void p_advanced_lmt_loop(std::vector<Halfedge *> &stack) {
        while(!stack.empty()) {
            Halfedge *e = stack.back();
            stack.pop_back();
            e->flagged = e->twin()->flagged = false;
            if(!e->find_local_minimal_certificate()) {
                e->restack_neighborhood(stack);
                e->status = e->twin()->status = LMTStatus::Impossible;
                p_lazy_delete(e);
            }
        }
    }

    /**
     * The actual LMT main loop.
     */
    void p_lmt_loop(std::vector<Halfedge *> &stack) {
        while(!stack.empty()) {
            auto e = stack.back();
            stack.pop_back();
            e->flagged = e->twin()->flagged = false;
            if(!e->find_certificate()) {
                e->restack_edges(stack);
                e->status = e->twin()->status = LMTStatus::Impossible;
                p_lazy_delete(e);
            }
        }
    }

    /**
     * Lazy edge deletion.
     */
    void p_lazy_delete(Halfedge *e) {
        auto s = p_point_index_of(e->source_handle());
        auto t = p_point_index_of(e->target_handle());
        e->unlink(&m_lmt_halfedges[m_offsets[s]], &m_lmt_halfedges[m_offsets[s + 1]]);
        e->twin()->unlink(&m_lmt_halfedges[m_offsets[t]], &m_lmt_halfedges[m_offsets[t + 1]]);
    }

    /**
     * Run the LMT main loop sequentially.
     */
    void p_lmt_main_loop_sequential() {
        std::vector<Halfedge *> stack;
        stack.reserve(64);
        std::size_t count = 0;
        for(auto &he : m_lmt_halfedges) {
            if(!he.is_primary())
                continue;
            if(he.status == LMTStatus::CH)
                continue;
            stack.push_back(&he);
            p_lmt_loop(stack);
        }
    }

    /*
     * Run the LMT main loop in parallel.
     * The parallelization is based on data partitioning. We split the points recursively in two disjoint
     * sets and process them in parallel. Each thread can safely process edges that connect vertices in
     * their respective set. Working concurrently on edges between different point sets leads to a data race and
     * has to be avoided. To increase concurrency, the partitioning is done at the leave level and the unsafe
     * edges are passed back to the caller who will process them in serial. It is important that all edges are
     * marked as being processed (i.e. on_stack = true) initially.
     */
    void p_lmt_main_loop_parallel() {
        auto pbegin = m_tree->points_begin();
        auto pend = m_tree->points_end();
        auto mid = pbegin + std::distance(pbegin, pend) / 2;
        std::vector<Halfedge *> left;
        std::vector<Halfedge *> right;
        tbb::parallel_invoke([&] { p_partition_and_loop(pbegin, mid, left); },
                             [&] { p_partition_and_loop(mid, pend, right); });
        p_lmt_loop(left);
        p_lmt_loop(right);
    }

    /**
     * The recursive data partitioning and parallel processing.
     * Unsafe edges cannot be dealt with in a thread-safe way and
     * are passed back to the caller who can deal with them.
     */
    void p_partition_and_loop(PointIterator begin, PointIterator end, std::vector<Halfedge *> &unsafe_edges) {
        std::vector<Halfedge *> safe_edges;
        if(std::distance(begin, end) < 1000) {
            p_partition_and_loop_small(begin, end, safe_edges, unsafe_edges);
        } else {
            p_partition_and_loop_large(begin, end, safe_edges, unsafe_edges);
        }
    }

    /**
     * Partitioning and looping in case there is still enough points.
     */
    void p_partition_and_loop_large(PointIterator begin, PointIterator end, std::vector<Halfedge *> &safe_edges,
                                    std::vector<Halfedge *> &unsafe_edges) {
        auto mid = begin + (end - begin) / 2;
        std::vector<Halfedge *> left;
        std::vector<Halfedge *> right;
        tbb::parallel_invoke([&] { p_partition_and_loop(begin, mid, left); },
                             [&] { p_partition_and_loop(mid, end, right); });

        const Point_2 *beg_point = &*begin;
        const Point_2 *end_point = &*end;
        auto handle_edge = [&](Halfedge *e) {
            if(e->target_handle() < beg_point || e->target_handle() >= end_point) {
                unsafe_edges.push_back(e);
            } else {
                safe_edges.push_back(e);
                p_lmt_loop(safe_edges);
            }
        };

        for(Halfedge *e : left)
            handle_edge(e);
        for(Halfedge *e : right)
            handle_edge(e);
    }

    /**
     * Looping in case there is not enough points.
     */
    void p_partition_and_loop_small(PointIterator begin, PointIterator end, std::vector<Halfedge *> &safe_edges,
                                    std::vector<Halfedge *> &unsafe_edges) {
        const Point_2 *beg_point = &*begin;
        const Point_2 *end_point = &*end;
        std::size_t point_begin_index = p_point_index_of(begin);
        std::size_t point_end_index = p_point_index_of(end);
        auto he_begin = m_lmt_halfedges.begin() + m_offsets[point_begin_index];
        auto he_end = m_lmt_halfedges.begin() + m_offsets[point_end_index];
        for(auto e = he_begin; e != he_end; ++e) {
            if(!e->is_primary() || e->status == LMTStatus::CH)
                continue;
            if(e->target_handle() < beg_point || e->target_handle() >= end_point) {
                unsafe_edges.push_back(&*e);
            } else {
                safe_edges.push_back(&*e);
                p_lmt_loop(safe_edges);
            }
        }
    }

    std::size_t p_point_index_of(PointIterator iter) const noexcept {
        return std::size_t(iter - m_tree->points_begin());
    }

    std::size_t p_point_index_of(const Point_2 *point) const noexcept {
        return std::size_t(point - &*m_tree->points_begin());
    }

    friend class test::LMTSkeletonTester;

    void p_validate_twin_relationship() {
        std::size_t npoints = m_tree->size();
        auto points = m_tree->points_begin();
        std::size_t correct_count = 0;
        for(std::size_t i = 0; i < npoints; ++i) {
            auto source = points[i];
            auto nbegin = m_lmt_halfedges.begin() + m_offsets[i];
            auto nend = m_lmt_halfedges.begin() + m_offsets[i + 1];
            for(auto it = nbegin; it != nend; ++it) {
                Halfedge *self = &*it;
                Halfedge *twin = self->twin();
                if(twin->twin() != self) {
                    throw std::logic_error("twin relationship is not symmetric");
                }
                if(twin->target() != source) {
                    throw std::logic_error("twin has wrong target");
                }
                if(twin->flagged != self->flagged) {
                    throw std::logic_error("self and twin have different flags");
                }
                if(twin->status != self->status) {
                    throw std::logic_error("self and twin have different status");
                }
                bool primal = self->is_primary();
                bool dual = twin->is_primary();
                if(primal == dual) {
                    throw std::logic_error("self and twin are both primary or both non-primary");
                }
            }
        }
    }

    void p_validate_candidates_and_skeleton() {
        for(Halfedge *e : m_skeleton) {
            if(e->status != LMTStatus::CH && e->status != LMTStatus::Certain) {
                throw std::logic_error("skeleton edge has wrong status");
            }
            if(!e->is_primary()) {
                throw std::logic_error("skeleton edge is not primary");
            }
        }
        for(Halfedge *e : m_candidates) {
            if(e->status != LMTStatus::Possible) {
                throw std::logic_error("candidate edge has wrong status");
            }
            if(!e->is_primary()) {
                throw std::logic_error("candidate edge is not primary");
            }
        }
        std::vector<Halfedge *> check;
        for(std::vector<Halfedge *> *hev : {&m_skeleton, &m_candidates}) {
            check = *hev;
            std::sort(check.begin(), check.end());
            auto last = std::unique(check.begin(), check.end());
            if(last != check.end()) {
                throw std::logic_error("duplicate edge in skeleton or candidates");
            }
        }
    }

    void p_validate_initial_invariants() {
        p_validate_twin_relationship();
        for(const auto &he : m_lmt_halfedges) {
            if(he.status == LMTStatus::CH) {
                if(he.flagged) {
                    throw std::logic_error("CH edge is flagged");
                }
                if(he.twin()->flagged) {
                    throw std::logic_error("CH edge's twin is flagged");
                }
            } else {
                if(!he.flagged) {
                    throw std::logic_error("non-CH edge is not flagged");
                }
                if(he.status != LMTStatus::Possible) {
                    throw std::logic_error("non-CH edge has wrong status");
                }
            }
        }
    }

    template<typename DiamondFilterEdges> void p_halfedges_from_diamond_filter(const DiamondFilterEdges &init_edges) {
        auto points = m_tree->points_begin();
        m_offsets.back() = 2 * init_edges.size();
        for(auto [u, v] : init_edges) {
            ++m_offsets[u];
            ++m_offsets[v];
        }
        for(auto i = m_tree->size(); i > 0; --i) {
            m_offsets[i - 1] = m_offsets[i] - m_offsets[i - 1];
        }
        std::vector<std::size_t> counts(m_tree->size(), 0);
        m_lmt_halfedges.resize(2 * init_edges.size());
        for(auto [u, v] : init_edges) {
            std::size_t i = m_offsets[u] + counts[u]++;
            std::size_t j = m_offsets[v] + counts[v]++;
            m_lmt_halfedges[i] = Halfedge{&points[v], &m_lmt_halfedges[j]};
            m_lmt_halfedges[j] = Halfedge{&points[u], &m_lmt_halfedges[i]};
        }
        p_sort_around_points();
    }

    template<typename HEIter> void p_sort_around_point(std::size_t point_index, HEIter begin, HEIter end) {
        Point_2 origin = m_tree->points_begin()[point_index];
        detail::sort_halfedges_around_point(origin, begin, end);
        m_convex_hull_size += detail::init_start_pointers(begin, end);
    }

    void p_sort_around_points() {
        const std::size_t npoints = m_tree->size();
        for(std::size_t pi = 0; pi < npoints; ++pi) {
            auto begin = m_lmt_halfedges.begin() + m_offsets[pi];
            auto end = m_lmt_halfedges.begin() + m_offsets[pi + 1];
            p_sort_around_point(pi, begin, end);
        }
    }

    template<typename DiamondFilterEdges> static void clear_if_move(DiamondFilterEdges &&edges) noexcept {
        ClearIfMove<DiamondFilterEdges>::clear(std::forward<DiamondFilterEdges>(edges));
    }

    template<typename DEL_> struct ClearIfMove {
        using DEL = std::remove_cv_t<std::remove_reference_t<DEL_>>;

        static void clear(const DEL &) noexcept {}

        static void clear(DEL &&del) noexcept {
            del.clear();
            del.shrink_to_fit();
        }
    };

    void p_mark_certain_edges_sequential() {
        const auto npoints = m_tree->size();
        m_skeleton.reserve(std::size_t(npoints * 3));
        m_candidates.reserve(std::size_t(npoints * 1.5));

        // Mark all the edges that have a certificate but no crossing edges as certain.
        for(auto &e : m_lmt_halfedges) {
            if(!e.is_primary() || e.status == LMTStatus::Impossible)
                continue;
            if(e.status == LMTStatus::CH) {
                m_skeleton.push_back(&e);
                continue;
            }
            if(e.status == LMTStatus::Possible) {
                if(!e.has_crossing_edges()) {
                    e.status = e.twin()->status = LMTStatus::Certain;
                    m_skeleton.push_back(&e);
                } else {
                    m_candidates.push_back(&e);
                }
            }
        }
    }

    void p_mark_certain_edges_parallel() {
        tbb::combinable<std::vector<Halfedge *>> certainized;
        Halfedge *begin = m_lmt_halfedges.data();
        Halfedge *end = begin + m_lmt_halfedges.size();
        tbb::parallel_for(tbb::blocked_range<Halfedge *>(begin, end, 2048), [&](const auto &range) {
            std::vector<Halfedge *> &local = certainized.local();
            for(Halfedge &e : range) {
                if(!e.is_primary() || e.status != LMTStatus::Possible)
                    continue;
                if(!e.has_crossing_edges()) {
                    local.push_back(&e);
                }
            }
        });
        certainized.combine_each([&](const auto &vec) {
            for(Halfedge *e : vec) {
                e->status = e->twin()->status = LMTStatus::Certain;
                m_skeleton.push_back(e);
            }
        });
        for(auto &e : m_lmt_halfedges) {
            if(!e.is_primary())
                continue;
            if(e.status == LMTStatus::CH) {
                m_skeleton.push_back(&e);
            } else if(e.status == LMTStatus::Possible) {
                m_candidates.push_back(&e);
            }
        }
    }

    /**
     * The tree (used purely as point container).
     */
    const Tree *m_tree;

    /**
     * All halfedges are stored in a single vector.
     */
    std::vector<Halfedge> m_lmt_halfedges;

    /**
     * The halfedges of vertex i start at m_offsets[i] and end at m_offsets[i+1].
     * They are sorted around the vertex in CCW order.
     */
    std::vector<std::size_t> m_offsets;

    /**
     * The edges in the LMT skeleton.
     */
    std::vector<Halfedge *> m_skeleton;

    /**
     * The candidate edges after the LMT skeleton computation,
     * i.e., edges that are still possibly in the MWT then.
     */
    std::vector<Halfedge *> m_candidates;

    /**
     * The number of edges on the convex hull.
     */
    std::size_t m_convex_hull_size{0};
};

} // namespace mwt

#endif
