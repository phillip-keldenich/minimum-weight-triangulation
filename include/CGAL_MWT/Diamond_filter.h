#ifndef CGAL_MWT_DIAMOND_FILTER_H_INCLUDED_
#define CGAL_MWT_DIAMOND_FILTER_H_INCLUDED_

#include "Exact_diamond_filtered_search_driver.h"
#include <CGAL/Interval_arithmetic.h>
#include <cstdint>
#include <tbb/blocked_range.h>
#include <tbb/concurrent_vector.h>
#include <tbb/parallel_for.h>
#include <type_traits>
#include <utility>
#include <vector>

namespace mwt {

template<typename Tree_, bool Parallel_, typename PointIndex_ = std::uint32_t> class Diamond_edge_filter {
  public:
    static constexpr bool parallel_search = Parallel_;
    static constexpr std::size_t expected_edges_per_node = 13;

    using Tree = Tree_;
    using PointIndex = PointIndex_;
    using Edge = std::pair<PointIndex, PointIndex>;
    using Edges = std::conditional_t<parallel_search, tbb::concurrent_vector<Edge>, std::vector<Edge>>;

    explicit Diamond_edge_filter(const Tree *tree) noexcept : m_tree(tree) {}

    /**
     * Compute a (as small as possible) superset of the edges
     * that satisfy the diamond property and might thus be part of the MWT.
     */
    void compute_remaining_edges() {
        m_edges.clear();
        m_edges.reserve(expected_edges_per_node * m_tree->size());

        if constexpr(parallel_search) {
            compute_remaining_edges_parallel(m_tree->points_begin(), m_tree->points_end());
        } else {
            compute_remaining_edges_sequential(m_tree->points_begin(), m_tree->points_end(), m_edges);
        }
    }

    const Edges &get_edges() const noexcept { return m_edges; }

    Edges &get_edges() noexcept { return m_edges; }

  private:
    template<typename PointConstIterator>
    void compute_remaining_edges_parallel(PointConstIterator begin, PointConstIterator end) {
        tbb::parallel_for(tbb::blocked_range<PointConstIterator>(begin, end, 1024), [&](const auto &range) {
            if(CGAL::FPU_get_cw() != CGAL_FE_UPWARD) {
                CGAL::FPU_set_cw(CGAL_FE_UPWARD);
            }
            std::vector<Edge> edges;
            edges.reserve(expected_edges_per_node * range.size());
            compute_remaining_edges_sequential(range.begin(), range.end(), edges);
            m_edges.grow_by(edges.begin(), edges.end());
        });
    }

    template<typename PointConstIterator>
    void compute_remaining_edges_sequential(PointConstIterator begin, PointConstIterator end,
                                            std::vector<Edge> &edges) {
        using Traits = mwt::Mwt_traits_2<typename Tree::Kernel>;
        Exact_diamond_filtered_search_driver<Traits, typename Tree::Iterator> driver(m_tree);
        const auto pbegin = m_tree->points_begin();
        for(auto it = begin; it != end; ++it) {
            const std::uint32_t qindex(it - pbegin);
            driver.enumerate_filtered_neighbors_of(*it, [&](auto niter) {
                std::uint32_t nindex(niter - pbegin);
                edges.emplace_back(qindex, nindex);
            });
        }
    }

    const Tree *m_tree;
    Edges m_edges;
};

} // namespace mwt

#endif
