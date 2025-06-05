#ifndef CGAL_MWT_EXACT_SIMPLE_FACE_TRIANGULATOR_H_INCLUDED_
#define CGAL_MWT_EXACT_SIMPLE_FACE_TRIANGULATOR_H_INCLUDED_

#include "Compare_weights.h"
#include "Dynamic_program_utils.h"
#include <CGAL/FPU.h>
#include <boost/unordered/unordered_flat_set.hpp>
#include <deque>
#include <variant>

namespace mwt {

/**
 * Interval-based simple face triangulator.
 */
template<typename Skeleton_> class Interval_simple_face_triangulator {
    using Kernel = typename Skeleton_::Kernel;
    using FT = typename Kernel::FT;
    using Skeleton = Skeleton_;
    using Face = mwt::Face<Skeleton>;
    using Halfedge = typename Skeleton::Halfedge;
    using Handle = Halfedge *;
    using Interval = CGAL::Interval_nt_advanced;
    using Iterator = typename std::vector<Handle>::iterator;
    struct EdgeData;
    using EdgeMap = mwt::Unique_hash_map<const Handle, EdgeData>;

  public:
    explicit Interval_simple_face_triangulator() noexcept {};

    /**
     * Call the interval-based simple face triangulator.
     */
    void operator()(Face &&face) {
        if(!face.is_simple()) {
            throw std::runtime_error("Interval_simple_face_triangulator: face is not simple");
        }
        p_reset(std::move(face));
        p_dp_sort_inner();
        m_inner_edges.push_back(m_boundary.front());
        m_real_inner_begin = p_init_lengths_for_boundary_distance_2();
        m_real_inner_begin = p_init_lengths_for_boundary_distance_3();
        p_run_dp();
        p_check_is_unique();
        if(!m_is_unique)
            return;
        p_reconstruct_triangulation_unique();
    }

    /**
     * Check if an optimal MWT was uniquely identified
     * by the dynamic programming algorithm using interval arithmetic.
     * In that case, the triangulation is already reconstructed
     * and all edges are marked as certain.
     */
    bool is_unique() const noexcept { return m_is_unique; }

    /**
     * Set the unique flag (indicating that we have identified a truly optimal solution).
     * If the original solution using intervals was not unique, this is usually done
     * by an exact method following up on the interval DP.
     */
    void set_unique(bool is_unique) noexcept { m_is_unique = is_unique; }

    /**
     * Reconstruct the triangulation from the DP information;
     * is either implicitly done by operator() if the interval routine was successful,
     * or can be called explicitly if the triangulation was made unique by an exact method.
     */
    void reconstruct_triangulation() {
        if(!m_is_unique) {
            throw std::logic_error(
                "Interval_simple_face_triangulator: reconstruct_triangulation called on non-unique solution");
        }
        p_check_is_unique();
        if(!m_is_unique) {
            throw std::logic_error("Interval_simple_face_triangulator: reconstruct_triangulation called on non-unique "
                                   "solution incorrectly flagged as unique");
        }
        p_reconstruct_triangulation_unique();
    }

  private:
    void p_run_dp() {
        DPCore core(this);
        core();
    }

    /**
     * Check if the best solution found by dynamic programming
     * is 'unique' in the sense that interval arithmetic could
     * prove that no better potential solution is possible.
     */
    void p_check_is_unique() {
        if(!m_edge_data.is_defined(m_boundary.front())) {
            throw std::logic_error("Interval_simple_face_triangulator: boundary halfedge properties not defined");
        }
        const auto &outer_edge_data = m_edge_data[m_boundary.front()];
        if(!outer_edge_data.is_unique()) {
            m_is_unique = false;
            return;
        }
        auto [i, j] = outer_edge_data.unique_pair();
        if(i == nullptr || j == nullptr) {
            throw std::logic_error("Interval_simple_face_triangulator: boundary halfedge i or j are null");
        }
        m_is_unique = true;
        m_reconstruction_stack.push_back(i);
        m_reconstruction_stack.push_back(j);
        m_inner_edges.clear();
        while(!m_reconstruction_stack.empty()) {
            Handle e = m_reconstruction_stack.back();
            m_reconstruction_stack.pop_back();
            if(e->status != LMTStatus::Possible)
                continue;
            m_inner_edges.push_back(e);
            auto &e_data = m_edge_data[e];
            if(!e_data.is_unique()) {
                m_reconstruction_stack.clear();
                m_is_unique = false;
                return;
            }
            auto [in, jn] = e_data.unique_pair();
            if(in != nullptr)
                m_reconstruction_stack.push_back(in);
            if(jn != nullptr)
                m_reconstruction_stack.push_back(jn);
        }
    }

    /**
     * Use the information gathered during a successful
     * p_check_is_unique call to mark the edges of the triangulation.
     */
    void p_reconstruct_triangulation_unique() {
        // only called if m_is_unique is true;
        // the edges are already identified
        for(Handle e : m_inner_edges) {
            e->status = e->twin()->status = LMTStatus::Certain;
        }
        detail::mark_remaining_as_impossible(m_face_inner);
    }

    /**
     * Sort the inner edges by increasing boundary distance.
     * This is needed for the DP to work correctly (since
     * it constructs the triangulation for edges with longer
     * boundary distance based on the shorter ones).
     */
    void p_dp_sort_inner() { std::sort(m_inner_edges.begin(), m_inner_edges.end(), m_boundary_distance); }

    /**
     * Reset the triangulator with a new face.
     */
    void p_reset(Face &&face) {
        m_nonunique_list.clear();
        m_reconstruction_stack.clear();
        m_boundary = std::move(face.boundary);
        m_face_inner = std::move(face.inner_edges);
        std::uint32_t bsize = detail::init_edge_data_and_indices(m_edge_data, m_inner_edges, m_boundary);
        m_boundary_distance = detail::Boundary_distance<EdgeMap>(m_edge_data, bsize);
    }

    /**
     * Init lengths of all edges with boundary distance 2,
     * i.e., the shortest internal lengths.
     */
    Iterator p_init_lengths_for_boundary_distance_2() {
        Iterator it = m_inner_edges.begin();
        while(it != m_inner_edges.end()) {
            Handle e = *it;
            auto &edata = m_edge_data[e];
            if(m_boundary_distance.compute(edata.src, edata.tar) != 2)
                break;
            edata.min_weight = p_weight(e);
            ++it;
        }
        return it;
    }

    /**
     * Init lengths of all edges with boundary distance 3,
     * i.e., edges where there needs to be exactly one internal edge
     * to the left of the edge. This case is special because we
     * can do the exact comparison without using sqrt.
     */
    Iterator p_init_lengths_for_boundary_distance_3() {
        Iterator it = m_real_inner_begin;
        while(it != m_inner_edges.end()) {
            Handle e = *it;
            auto &edata = m_edge_data[e];
            if(m_boundary_distance.compute(edata.src, edata.tar) != 3)
                break;
            p_handle_boundary_distance_3(e, edata);
            ++it;
        }
        return it;
    }

    /**
     * Compare by squared distance if the underlying FT is exact.
     */
    bool p_is_lower_squared_distance_exact_ft(Handle old_i, Handle old_j, Handle new_i, Handle new_j) {
        auto compute_sqdist = [](Handle i, Handle j) -> FT {
            if(i->status == LMTStatus::Possible) {
                return CGAL::squared_distance(i->source(), i->target());
            } else {
                return CGAL::squared_distance(j->source(), j->target());
            }
        };
        return compute_sqdist(new_i, new_j) < compute_sqdist(old_i, old_j);
    }

    /**
     * Exactly compare by squared distance even if the underlying FT is inexact.
     */
    bool p_is_lower_squared_distance_inexact_ft(Handle old_i, Handle old_j, Handle new_i, Handle new_j) {
        using ExactFT = CGAL::Exact_predicates_exact_constructions_kernel::FT;
        using ExactP = CGAL::Exact_predicates_exact_constructions_kernel::Point_2;
        auto squared_distance = [](Handle e) -> ExactFT {
            ExactP s(e->source().x(), e->source().y());
            ExactP t(e->target().x(), e->target().y());
            return CGAL::squared_distance(s, t);
        };
        auto compute_sqdist = [&](Handle i, Handle j) -> ExactFT {
            if(i->status == LMTStatus::Possible) {
                return squared_distance(i);
            } else {
                return squared_distance(j);
            }
        };
        return compute_sqdist(new_i, new_j) < compute_sqdist(old_i, old_j);
    }

    /**
     * Compare the two triangles on the left of an edge
     * with boundary distance 3 by their (squared) distance.
     */
    bool p_is_lower_squared_distance(Handle old_i, Handle old_j, Handle new_i, Handle new_j) {
        CGAL_assertion(old_i->status != LMTStatus::Possible || old_j->status != LMTStatus::Possible);
        CGAL_assertion(old_i->status == LMTStatus::Possible || old_j->status == LMTStatus::Possible);
        CGAL_assertion(new_i->status != LMTStatus::Possible || new_j->status != LMTStatus::Possible);
        CGAL_assertion(new_i->status == LMTStatus::Possible || new_j->status == LMTStatus::Possible);
        if constexpr(CGAL::Algebraic_structure_traits<FT>::Is_exact::value) {
            return p_is_lower_squared_distance_exact_ft(old_i, old_j, new_i, new_j);
        } else {
            return p_is_lower_squared_distance_inexact_ft(old_i, old_j, new_i, new_j);
        }
    }

    /**
     * Compute the best way to triangulate the left-hand side of an
     * edge with boundary distance 3. Uses squared distance comparison
     * to break ties to avoid using any sqrt.
     * Helps, in particular, with triangulation of empty (near-)rectangles in the input.
     */
    void p_handle_boundary_distance_3(Handle e, EdgeData &edata) {
        e->reset();
        Interval current_min_weight(std::numeric_limits<double>::max(), std::numeric_limits<double>::infinity());
        Handle current_candidate_i = nullptr, current_candidate_j = nullptr;
        while(e->next_triangle()) {
            if(detail::is_bad_triangle(e, edata, m_edge_data))
                continue;
            Interval weight = 0;
            Handle i = e->i, j = e->j->twin();
            if(i->status == LMTStatus::Possible) {
                weight += m_edge_data[i].min_weight;
            }
            if(j->status == LMTStatus::Possible) {
                weight += m_edge_data[j].min_weight;
            }
            if(current_candidate_i == nullptr || weight.sup() < current_min_weight.inf()) {
                current_min_weight = weight;
                current_candidate_i = i;
                current_candidate_j = j;
            } else if(weight.inf() < current_min_weight.sup()) {
                if(p_is_lower_squared_distance(current_candidate_i, current_candidate_j, i, j)) {
                    current_min_weight = weight;
                    current_candidate_i = i;
                    current_candidate_j = j;
                }
            }
        }
        edata.candidate = CandidateEntry{current_candidate_i, current_candidate_j};
        edata.min_weight = current_min_weight + p_weight(e);
    }

    /**
     * Weight of an edge by handle.
     */
    Interval p_weight(Handle e) {
        auto s = e->source();
        auto t = e->target();
        Interval wx = s.x(), wy = s.y();
        wx -= t.x();
        wy -= t.y();
        wx *= wx;
        wy *= wy;
        wx += wy;
        return CGAL::sqrt(wx);
    }

  public:
    /**
     * Entry in the candidate list.
     * The candidate list is used to store multiple
     * candidates in case interval arithmetic cannot
     * decide which one is the best.
     */
    struct CandidateEntry {
        Handle i, j;

        CandidateEntry() noexcept : i(nullptr), j(nullptr) {}
        CandidateEntry(Handle i, Handle j) noexcept : i(i), j(j) {}
    };

    using CandidateIterator = typename std::vector<CandidateEntry>::const_iterator;

  private:
    /**
     * Datastructure capturing where multiple candidates
     * are stored in the m_nonunique_list.
     */
    struct CandidateSublist {
        std::size_t candidate_begin, num_candidates;
    };

    /**
     * EdgeData structure storing information on the edges
     * used in the dynamic programming algorithm.
     */
    struct EdgeData {
        EdgeData() noexcept
            : src(0), tar(0), min_weight(std::numeric_limits<double>::max(), std::numeric_limits<double>::infinity()),
              candidate(CandidateEntry{}), antenna(false) {}

        EdgeData(std::uint32_t src, std::uint32_t tar) noexcept
            : src(src), tar(tar),
              min_weight(std::numeric_limits<double>::max(), std::numeric_limits<double>::infinity()),
              candidate(CandidateEntry{}), antenna(false) {}

        /**
         * After DP with interval arithmetic, check whether
         * the solution is unique, i.e., there is exactly one
         * candidate triangle on the left of this edge which
         * has been uniquely identified by IA to lead to the
         * best triangulation of the left side of this edge.
         * This is a *SHALLOW* uniqueness check, i.e., some
         * decendant edges may have multiple candidates.
         */
        bool is_unique() const noexcept { return std::holds_alternative<CandidateEntry>(candidate); }

        std::pair<Handle, Handle> unique_pair() const noexcept {
            const CandidateEntry *ce = std::get_if<CandidateEntry>(&candidate);
            return std::make_pair(ce->i, ce->j);
        }

        void store_multiple(std::size_t begin, std::size_t num) { candidate = CandidateSublist{begin, num}; }

        std::uint32_t src, tar;
        Interval min_weight;
        std::variant<CandidateEntry, CandidateSublist> candidate;
        bool antenna; //< marks antenna edges
    };

    /**
     * Core implementation of the dynamic programming algorithm.
     */
    class DPCore {
      public:
        explicit DPCore(Interval_simple_face_triangulator *that) : that(that) {}

        void operator()() {
            Iterator begin = that->m_real_inner_begin;
            Iterator end = that->m_inner_edges.end();
            std::for_each(begin, end, [&](Halfedge *e) { p_handle_halfedge(e); });
        }

      private:
        void p_handle_halfedge(Halfedge *e) {
            e->reset();
            current_candidates.clear();
            current_min_weight = Interval(std::numeric_limits<double>::max(), std::numeric_limits<double>::infinity());
            auto &edata = that->m_edge_data[e];
            while(e->next_triangle()) {
                p_handle_triangle(e, edata);
            }
            current_min_weight += that->p_weight(e);
            if(current_candidates.size() == 1) {
                edata.candidate = current_candidates.front();
                edata.min_weight = current_min_weight;
            } else {
                auto &nul = that->m_nonunique_list;
                edata.store_multiple(nul.size(), current_candidates.size());
                edata.min_weight = current_min_weight;
                nul.insert(nul.end(), current_candidates.begin(), current_candidates.end());
            }
        }

        void p_handle_triangle(Halfedge *e, EdgeData &edata) {
            if(detail::is_bad_triangle(e, edata, that->m_edge_data))
                return;
            Interval weight = p_calculate_weight(e);
            CGAL::Uncertain<bool> comp = (weight < current_min_weight);
            if(CGAL::certainly_not(comp))
                return;
            Handle i = e->i, j = e->j->twin();
            if(CGAL::certainly(comp)) {
                current_candidates.assign(1, CandidateEntry{i, j});
                current_min_weight = weight;
                return;
            }
            current_candidates.push_back(CandidateEntry{i, j});
            auto new_lb = (std::min)(weight.inf(), current_min_weight.inf());
            auto new_ub = (std::max)(weight.sup(), current_min_weight.sup());
            current_min_weight = Interval{new_lb, new_ub};
        }

        Interval p_calculate_weight(Halfedge *e) {
            auto &edge_map = that->m_edge_data;
            Interval weight(0);
            Handle i = e->i, j = e->j->twin();
            if(i->status == LMTStatus::Possible) {
                weight += edge_map[i].min_weight;
            }
            if(j->status == LMTStatus::Possible) {
                weight += edge_map[j].min_weight;
            }
            return weight;
        }

        Interval_simple_face_triangulator *that;
        std::vector<CandidateEntry> current_candidates;
        Interval current_min_weight;
    };

    std::vector<CandidateEntry> m_nonunique_list;
    std::vector<Handle> m_inner_edges;
    std::vector<Handle> m_boundary;
    std::vector<Handle> m_reconstruction_stack;
    std::vector<Handle> m_face_inner;
    detail::Boundary_distance<EdgeMap> m_boundary_distance;
    Iterator m_real_inner_begin;
    EdgeMap m_edge_data;
    bool m_is_unique;
    template<typename S> friend class NonUniqueDPTreeBFS;

  public:
    class DPTree {
        Interval_simple_face_triangulator *that;
        Halfedge *node;
        EdgeData *node_data;

      public:
        explicit DPTree(Interval_simple_face_triangulator *that)
            : that(that), node(that->m_boundary.front()), node_data(&that->m_edge_data[node]) {}

        explicit DPTree(Interval_simple_face_triangulator *that, Halfedge *node, EdgeData *node_data) noexcept
            : that(that), node(node), node_data(node_data) {}

        EdgeData &data() noexcept { return *node_data; }
        const EdgeData &data() const noexcept { return *node_data; }

        bool shallow_unique() const noexcept { return node_data->is_unique(); }

        Halfedge *halfedge() const noexcept { return node; }

        std::pair<Handle, Handle> shallow_unique_pair() const noexcept { return node_data->unique_pair(); }

        std::pair<CandidateIterator, CandidateIterator> candidate_subrange() const noexcept {
            CGAL_precondition(!shallow_unique());
            CandidateSublist sublist = std::get<CandidateSublist>(node_data->candidate);
            auto slbegin = that->m_nonunique_list.begin() + sublist.candidate_begin;
            auto slend = slbegin + sublist.num_candidates;
            return {slbegin, slend};
        }

        std::vector<DPTree> children() const {
            std::vector<DPTree> c;
            children_into(c);
        }

        void child_halfedges_into(std::vector<Halfedge *> &c) const {
            auto push_if = [&](Halfedge *h) {
                if(h != nullptr && h->status == LMTStatus::Possible) {
                    c.push_back(h);
                }
            };
            c.clear();
            if(shallow_unique()) {
                auto [i, j] = shallow_unique_pair();
                push_if(i);
                push_if(j);
            } else {
                CandidateSublist sublist = std::get<CandidateSublist>(node_data->candidate);
                auto slbegin = that->m_nonunique_list.begin() + sublist.candidate_begin;
                auto slend = slbegin + sublist.num_candidates;
                std::for_each(slbegin, slend, [&](const CandidateEntry &ce) {
                    push_if(ce.i);
                    push_if(ce.j);
                });
            }
        }

        void children_into(std::vector<DPTree> &c) const {
            c.clear();
            if(shallow_unique()) {
                auto [i, j] = shallow_unique_pair();
                p_push_child(c, i);
                p_push_child(c, j);
            } else {
                CandidateSublist sublist = std::get<CandidateSublist>(node_data->candidate);
                auto slbegin = that->m_nonunique_list.begin() + sublist.candidate_begin;
                auto slend = slbegin + sublist.num_candidates;
                std::for_each(slbegin, slend, [&](const CandidateEntry &ce) {
                    p_push_child(c, ce.i);
                    p_push_child(c, ce.j);
                });
            }
        }

        std::vector<Halfedge *> ancestor_edges() {
            std::vector<Halfedge *> result;
            add_ancestor_edges_to(result);
            return result;
        }

        void ancestor_edges_into(std::vector<Halfedge *> &result) {
            CGAL_precondition(node->status == LMTStatus::Possible);
            CGAL_precondition(node_data->is_unique());
            result.clear();
            ancestor_edges_into(result, *std::get_if<CandidateEntry>(&node_data->candidate));
        }

        void ancestor_edges_into(std::vector<Halfedge *> &result, CandidateEntry assumed_entry) {
            result.clear();
            auto &edge_map = that->m_edge_data;
            auto i = assumed_entry.i;
            auto j = assumed_entry.j;
            if(i != nullptr && i->status == LMTStatus::Possible && edge_map.lookup(i)) {
                result.push_back(i);
            }
            if(j != nullptr && j->status == LMTStatus::Possible && edge_map.lookup(j)) {
                result.push_back(j);
            }
            p_bfs(result);
        }

      private:
        void p_bfs(std::vector<Halfedge *> &result) {
            auto &edge_map = that->m_edge_data;
            std::size_t index = 0;
            while(index < result.size()) {
                Halfedge *current = result[index++];
                CGAL_assertion(current->status == LMTStatus::Possible);
                EdgeData *current_node_data = &edge_map[current];
                CGAL_assertion(current_node_data->is_unique());
                auto [i, j] = current_node_data->unique_pair();
                if(i != nullptr && i->status == LMTStatus::Possible) {
                    result.push_back(i);
                }
                if(j != nullptr && j->status == LMTStatus::Possible) {
                    result.push_back(j);
                }
            }
        }

        void p_push_child(std::vector<DPTree> &c, Handle child) const {
            if(child == nullptr)
                return;
            if(child->status != LMTStatus::Possible)
                return;
            auto lookup = that->m_edge_data.lookup(child);
            if(!lookup)
                return;
            c.emplace_back(that, child, &that->m_edge_data[lookup]);
        }
    };
};

template<typename Skeleton_> class NonUniqueDPTreeBFS {
    using Skeleton = Skeleton_;
    using ISFT = mwt::Interval_simple_face_triangulator<Skeleton>;
    using Halfedge = typename Skeleton::Halfedge;
    using Kernel = typename Skeleton::Kernel;
    using DPTree = typename ISFT::DPTree;
    using CandidateEntry = typename ISFT::CandidateEntry;

    ISFT *m_dp;
    std::vector<Halfedge *> m_current_children;
    std::vector<DPTree> m_queue;
    std::vector<std::pair<DPTree, std::size_t>> m_nonunique_edges;
    std::vector<Halfedge *> m_ancestor_edges;
    std::vector<Halfedge *> m_best_ancestor_edges;
    boost::unordered_flat_set<Halfedge *> m_visited_edges;
    CompareWeights<Kernel> m_compare;

  public:
    explicit NonUniqueDPTreeBFS(ISFT *dp) : m_dp(dp) {}

    void dfs_collect_potentially_used_halfedges() {
        p_reset();
        auto &edge_data_map = m_dp->m_edge_data;
        const auto &boundary_distance = m_dp->m_boundary_distance;
        m_queue.emplace_back(m_dp);
        m_visited_edges.emplace(m_queue.back().halfedge());
        if(!m_queue.back().shallow_unique()) {
            m_nonunique_edges.emplace_back(m_queue.back(), std::numeric_limits<std::size_t>::max());
        }
        while(!m_queue.empty()) {
            DPTree current = m_queue.back();
            m_queue.pop_back();
            current.child_halfedges_into(m_current_children);
            for(Halfedge *e : m_current_children) {
                auto [it, inserted] = m_visited_edges.insert(e);
                if(!inserted)
                    continue;
                auto &edge_data = edge_data_map[e];
                m_queue.emplace_back(m_dp, e, &edge_data);
                auto &back = m_queue.back();
                if(!back.shallow_unique()) {
                    m_nonunique_edges.emplace_back(back, boundary_distance.compute(edge_data.src, edge_data.tar));
                }
            }
        }
    }

    void make_edges_unique() {
        std::sort(m_nonunique_edges.begin(), m_nonunique_edges.end(),
                  [&](const auto &e1, const auto &e2) { return e1.second < e2.second; });
        for(const auto &entry : m_nonunique_edges) {
            CGAL_assertion(!entry.first.shallow_unique());
            p_make_unique_with_children_unique(entry.first);
            CGAL_assertion(entry.first.shallow_unique());
        }
    }

  private:
    void p_reset() {
        m_nonunique_edges.clear();
        m_queue.clear();
        m_visited_edges.clear();
    }

    /**
     * Compute a unique candidate for the edge, based on the
     * precondition that all children (and their ancestors) have
     * previously been made unique already.
     */
    void p_make_unique_with_children_unique(DPTree t) {
        auto [beg, end] = t.candidate_subrange();
        CGAL_precondition(beg != end);
        CGAL_precondition(std::next(beg) != end);
        t.ancestor_edges_into(m_best_ancestor_edges, *beg);
        CandidateEntry best = *beg;
        for(++beg; beg != end; ++beg) {
            t.ancestor_edges_into(m_ancestor_edges, *beg);
            m_compare.reset();
            // set the current best as LHS, and *beg as RHS
            for(Halfedge *e : m_best_ancestor_edges) {
                m_compare.add_lhs(e->source(), e->target());
            }
            for(Halfedge *e : m_ancestor_edges) {
                m_compare.add_rhs(e->source(), e->target());
            }
            // if Sign(LHS - RHS) > 0, LHS > RHS, so we update the current best
            if(m_compare.is_positive()) {
                best = *beg;
                m_best_ancestor_edges.swap(m_ancestor_edges);
            }
        }
        t.data().candidate = best;
        CGAL_postcondition(t.shallow_unique());
    }
};

template<typename Skeleton> void fix_nonunique_interval_dp(Interval_simple_face_triangulator<Skeleton> &dp) {
    if(dp.is_unique())
        return;
    NonUniqueDPTreeBFS<Skeleton> bfs(&dp);
    bfs.dfs_collect_potentially_used_halfedges();
    bfs.make_edges_unique();
    dp.set_unique(true);
    dp.reconstruct_triangulation();
}

template<typename Skeleton_> struct Exact_simple_face_triangulator {
    using Skeleton = Skeleton_;
    using Kernel = typename Skeleton::Kernel;
    using FT = typename Kernel::FT;
    using Face = mwt::Face<Skeleton>;

    Exact_simple_face_triangulator(FaceTriangulatorOptions) {}

    void operator()(Face &&face) {
        CGAL::Protect_FPU_rounding protect_rounding;
        Interval_simple_face_triangulator<Skeleton> interval_triangulator;
        interval_triangulator(std::move(face));
        if(interval_triangulator.is_unique())
            return;
        fix_nonunique_interval_dp(interval_triangulator);
    }
};

} // namespace mwt

#endif
