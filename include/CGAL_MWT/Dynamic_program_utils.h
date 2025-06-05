#ifndef CGAL_MWT_DYNAMIC_PROGRAM_UTILS_H_INCLUDED_
#define CGAL_MWT_DYNAMIC_PROGRAM_UTILS_H_INCLUDED_

#include "Compare_weights.h"
#include "Face_analyzer.h"
#include "Face_collector.h"
#include "LMT_halfedge.h"
#include "Unique_hash_map.h"
#include <cstdint>

namespace mwt {

namespace detail {

/**
 * Mark edges still marked as possible as impossible;
 * usually done after the dynamic programming is done.
 */
template<typename InnerEdgeIterable> void mark_remaining_as_impossible(const InnerEdgeIterable &inner) {
    for(auto *h : inner) {
        if(h->status == LMTStatus::Possible) {
            h->status = LMTStatus::Impossible;
            h->twin()->status = LMTStatus::Impossible;
        }
    }
}

/**
 * Enumerate the inner edges around the source of the
 * given halfedge e, assigning the given index as source
 * index (and target index of the twins).
 * The edges are added to inner_edges.
 */
template<typename EdgeMapType, typename InnerEdgeContainer, typename EdgeHandleType>
void enumerate_inner_from(EdgeMapType &edge_data, InnerEdgeContainer &inner_edges, EdgeHandleType e,
                          std::uint32_t assign_index) {
    auto *next = e->next_edge();
    while(next->status == LMTStatus::Possible) {
        inner_edges.push_back(next);
        edge_data[next].src = assign_index;
        edge_data[next->twin()].tar = assign_index;
        next = next->next_edge();
    }
}

/**
 * Initialize the edge data for all boundary and inner edges,
 * enumerating all inner edges in the process.
 * Returns the boundary size.
 * Halfedges sharing an endpoint can have different
 * indices for that point in degenerate cases, e.g.,
 * if the polygon contains an antenna and thus has
 * the same vertex twice on the boundary.
 */
template<typename EdgeMapType, typename InnerEdgeContainer, typename BoundaryContainer>
std::uint32_t init_edge_data_and_indices(EdgeMapType &edge_data, InnerEdgeContainer &inner_edges,
                                         const BoundaryContainer &boundary) {
    using EdgeDataType = typename EdgeMapType::data_type;
    inner_edges.clear();
    edge_data.clear();

    std::uint32_t index = 0;
    for(auto *e : boundary) {
        auto *twin = e->twin();
        auto &new_edata = edge_data[e];
        new_edata.src = index;
        new_edata.tar = index + 1;
        auto twin_edata = edge_data.lookup(twin);
        if(twin_edata) {
            new_edata.antenna = true;
            edge_data[twin_edata].antenna = true;
        } else {
            edge_data[twin] = EdgeDataType(index + 1, index);
        }
        enumerate_inner_from(edge_data, inner_edges, e, index);
        ++index;
    }
    auto replace_overshoot = [index](std::uint32_t &val) {
        if(val == index) {
            val = 0;
        }
    };
    replace_overshoot(edge_data[boundary.back()].tar);
    replace_overshoot(edge_data[boundary.back()->twin()].src);
    CGAL_postcondition(index == boundary.size());
    return index;
}

/**
 * Detect triangles that are non-empty or otherwise unsuitable,
 * e.g., because of situations arising from degenerate cases
 * with faces that have antennas.
 */
template<typename EdgeMapType, typename EdgeData, typename Halfedge>
bool is_bad_triangle(Halfedge *e, EdgeData &edge_data, EdgeMapType &edge_map) {
    // check presence of both halfedges in our map.
    // if some edges are not inside the face, we must not use the triangle.
    Halfedge *i = e->i;
    auto iitem = edge_map.lookup(i);
    if(iitem == nullptr) {
        return true;
    }
    Halfedge *j = e->j;
    auto jitem = edge_map.lookup(j);
    if(jitem == nullptr) {
        return true;
    }
    const auto &idata = edge_map[iitem];
    const auto &jdata = edge_map[jitem];

    /**
     * Compute the (boundary index-wise) twin of the given edge;
     * takes the halfedge twin and returns its boundary source and target (swapped).
     */
    auto conceptual_twin = [](Halfedge *edge, EdgeMapType &edge_map) -> std::pair<std::uint32_t, std::uint32_t> {
        Halfedge *twin = edge->twin();
        const auto &data = edge_map[twin];
        return {data.tar, data.src};
    };

    /**
     * Check for any match between 1 or 2 values on each side.
     */
    struct MatchFinder {
        MatchFinder(std::uint32_t ind1, std::uint32_t ind2) noexcept {
            indices1[0] = ind1;
            indices2[0] = ind2;
            indices1[1] = std::numeric_limits<std::uint32_t>::max();
            indices2[1] = std::numeric_limits<std::uint32_t>::max() - 1;
        }

        void add1(std::uint32_t ind1) { indices1[1] = ind1; }

        void add2(std::uint32_t ind2) { indices2[1] = ind2; }

        bool have_match() const {
            return (indices1[0] == indices2[1]) | (indices1[1] == indices2[0]) | (indices1[1] == indices2[1]);
        }

        std::uint32_t indices1[2];
        std::uint32_t indices2[2];
    };

    // check for invalid configurations causing
    // triangles to be misidentified as empty
    if(jdata.src != edge_data.tar) {
        MatchFinder match(jdata.src, edge_data.tar);
        if(jdata.antenna) {
            auto [ctsrc, _] = conceptual_twin(j, edge_map);
            match.add1(ctsrc);
        }
        if(edge_data.antenna) {
            auto [_, cttgt] = conceptual_twin(e, edge_map);
            match.add2(cttgt);
        }
        if(!match.have_match()) {
            return true;
        }
    }
    if(idata.tar != jdata.tar) {
        MatchFinder match(idata.tar, jdata.tar);
        if(idata.antenna) {
            auto [_, cttgt] = conceptual_twin(i, edge_map);
            match.add1(cttgt);
        }
        if(jdata.antenna) {
            auto [_, cttgt] = conceptual_twin(j, edge_map);
            match.add2(cttgt);
        }
        if(!match.have_match()) {
            return true;
        }
    }
    if(idata.src != edge_data.src) {
        MatchFinder match(idata.src, edge_data.src);
        if(idata.antenna) {
            auto [ctsrc, _] = conceptual_twin(i, edge_map);
            match.add1(ctsrc);
        }
        if(edge_data.antenna) {
            auto [ctsrc, _] = conceptual_twin(e, edge_map);
            match.add2(ctsrc);
        }
        if(!match.have_match()) {
            return true;
        }
    }
    return false;
}

/**
 * Compare or compute the hop distance
 * between points on the boundary given their boundary index.
 */
template<typename EdgeDataMap> struct Boundary_distance {
    EdgeDataMap *edge_data;
    std::uint32_t boundary_size;

    Boundary_distance() = default;
    ~Boundary_distance() = default;
    Boundary_distance(const Boundary_distance &) = default;
    Boundary_distance &operator=(const Boundary_distance &) = default;

    explicit Boundary_distance(EdgeDataMap &edge_data, std::uint32_t boundary_size)
        : edge_data(&edge_data), boundary_size(boundary_size) {}

    std::uint32_t compute(std::uint32_t src, std::uint32_t tar) const {
        return (src - tar + boundary_size) % boundary_size;
    }

    template<typename Handle> bool operator()(const Handle e1, const Handle e2) const {
        const auto &edata1 = (*edge_data)[e1];
        const auto &edata2 = (*edge_data)[e2];
        std::uint32_t bdist1 = compute(edata1.src, edata1.tar);
        std::uint32_t bdist2 = compute(edata2.src, edata2.tar);
        return bdist1 < bdist2;
    }
};

} // namespace detail

} // namespace mwt

#endif
