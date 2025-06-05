#ifndef CGAL_MWT_FACE_COLLECTOR_H_INCLUDED_
#define CGAL_MWT_FACE_COLLECTOR_H_INCLUDED_

#include "LMT_skeleton.h"
#include <CGAL/Kernel/global_functions.h>

namespace mwt {

/**
 * Check the orientation of the given boundary (CW or CCW).
 * CW corresponds to looking at the outer face of the boundary,
 * CCW to the inner face, as is usual under standard assumptions.
 */
template<typename HalfedgeHandle> CGAL::Orientation boundary_orientation(const std::vector<HalfedgeHandle> &boundary) {
    using IterTraits = std::iterator_traits<HalfedgeHandle>;
    using IterValue = typename IterTraits::value_type;
    using Kernel = typename IterValue::Traits;
    using Less_xy_2 = typename Kernel::Less_xy_2;

    if(boundary.empty())
        return CGAL::CLOCKWISE; // default empty boundaries to clockwise
    auto lex_min = std::min_element(boundary.begin(), boundary.end(), [](HalfedgeHandle lhs, HalfedgeHandle rhs) {
        return Less_xy_2{}(lhs->target(), rhs->target());
    });
    auto next = std::next(lex_min);
    if(next == boundary.end())
        next = boundary.begin();
    HalfedgeHandle min_target = *lex_min;
    HalfedgeHandle next_target = *next;
    if(min_target->has_on_left_side(next_target->target())) {
        return CGAL::COUNTERCLOCKWISE;
    }
    // we cannot have collinearity here
    return CGAL::CLOCKWISE;
}

/**
 * Represents a face in the skeleton.
 */
template<typename Skeleton_> struct Face {
    using Skeleton = Skeleton_;
    using Tree = typename Skeleton::Tree;
    using Traits = typename Skeleton::Traits;
    using Kernel = typename Traits::Kernel;
    using Halfedge = typename Skeleton::Halfedge;
    using Point_2 = typename Skeleton::Point_2;
    using Left_turn_2 = typename Skeleton::Left_turn_2;
    using Less_xy_2 = typename Kernel::Less_xy_2;

    std::vector<Halfedge *> boundary;
    std::vector<Halfedge *> inner_edges;
    std::vector<Halfedge *> hole_boundaries;
    std::vector<const Point_2 *> isolated_vertices;

    void clear() noexcept {
        boundary.clear();
        inner_edges.clear();
        hole_boundaries.clear();
        isolated_vertices.clear();
    }

    void reserve() {
        boundary.reserve(64);
        inner_edges.reserve(128);
        hole_boundaries.reserve(64);
        isolated_vertices.reserve(16);
    }

    bool has_inner_edges() const noexcept { return !inner_edges.empty(); }

    bool is_simple() const noexcept { return hole_boundaries.empty() && isolated_vertices.empty(); }

    void print_as_json(std::ostream &output) const {
        output << "{\"boundary\": ";
        print_halfedge_list(output, boundary);
        output << ", \"hole_boundaries\": ";
        print_halfedge_list(output, hole_boundaries);
        output << ", \"inner_edges\": ";
        print_halfedge_list(output, inner_edges);
        output << "}";
    }

    void print_halfedge_list(std::ostream &output, const std::vector<Halfedge *> &list) const {
        output << "[";
        bool first = true;
        for(auto edge : list) {
            if(first) {
                first = false;
            } else {
                output << ", ";
            }
            auto s = edge->source();
            auto t = edge->target();
            output << std::setprecision(19);
            output << "{\"source\": [" << s.x() << ", " << s.y() << "], \"target\": [" << t.x() << ", " << t.y()
                   << "]}";
        }
        output << "]";
    }
};

namespace test {
class Face_collector_tester;
}

/**
 * Given halfedges that are certainly in the MWT,
 * collect faces that are not yet triangulated.
 */
template<typename Skeleton_> class Face_collector {
  public:
    using Skeleton = Skeleton_;

  private:
    using Tree = typename Skeleton::Tree;
    using Traits = typename Skeleton::Traits;
    using Kernel = typename Traits::Kernel;
    using Halfedge = typename Skeleton::Halfedge;
    using Point_2 = typename Skeleton::Point_2;
    using Left_turn_2 = typename Skeleton::Left_turn_2;
    using Less_xy_2 = typename Kernel::Less_xy_2;

    friend class test::Face_collector_tester;

  public:
    using Face = mwt::Face<Skeleton>;

    explicit Face_collector(Skeleton &skeleton)
        : Face_collector(skeleton.skeleton().begin(), skeleton.skeleton().end()) {}

    template<typename HalfedgeIterator>
    Face_collector(HalfedgeIterator skeleton_begin, HalfedgeIterator skeleton_end)
        : m_skeleton_begin(&*skeleton_begin), m_skeleton_end(&*skeleton_end), m_skeleton_current(m_skeleton_begin) {
        m_current_face.reserve();
        CGAL_expensive_precondition(p_vertify_initial());
    }

    /**
     * Switch to the next face that is not already fully triangulated
     * in the LMT skeleton.
     * Returns true if we found such a face, false otherwise.
     */
    bool next() {
        m_current_face.clear();
        while(m_skeleton_current < m_skeleton_end) {
            Halfedge *current = *m_skeleton_current++;
            if(current->flagged)
                continue;
            if(p_collect_face(current)) {
                return true;
            }
            m_current_face.clear();
        }
        return false;
    }

    /**
     * Get a copy of the current face.
     */
    Face get_current_face_copy() const { return m_current_face; }

    /**
     * Get a reference to the current face.
     */
    Face &get_current_face_reference() { return m_current_face; }

  private:
    /**
     * Debugging/testing method: Check that the LMT skeleton is valid.
     */
    void p_verify_initial() {
        for(auto it = m_skeleton_begin; it != m_skeleton_end; ++it) {
            Halfedge *edge = *it;
            if(edge->flagged)
                throw std::logic_error("Face_collector: initial skeleton has flagged edges");
            if(edge->status != LMTStatus::Certain && edge->status != LMTStatus::CH) {
                throw std::logic_error("Face_collector: initial skeleton has uncertain edges");
            }
            p_verify_next(edge, true);
        }
    }

    /**
     * Debugging/testing method: Check the next pointers of the given edge.
     */
    void p_verify_next(Halfedge *edge, bool must_be_unvisited) {
        Halfedge *current = edge->next_edge();
        while(current != edge) {
            if(must_be_unvisited && current->flagged) {
                throw std::logic_error("Face_collector: initial possible edge is flagged");
            }
            if(current->status == LMTStatus::Impossible) {
                throw std::logic_error("Face_collector: initial graph structure has impossible edges");
            }
            current = current->next_edge();
        }
    }

    /**
     * Collect a face bounded by the given halfedge.
     */
    bool p_collect_face(Halfedge *current) {
        p_traverse_boundary_and_collect_edges(current, m_current_face.boundary, m_current_face.inner_edges);
        if(!m_current_face.has_inner_edges())
            return false;
        p_check_for_nonsimple_face();
        return true;
    }

    /**
     * After finding a non-triangulated face, check whether it is simple.
     */
    void p_check_for_nonsimple_face() {
        // boundary orientation value
        CGAL::Orientation orientation = CGAL::ZERO;
        auto &current_inner = m_current_face.inner_edges;
        for(std::size_t i = 0; i < current_inner.size(); ++i) {
            // loop is index-based because the vector can
            // be changed during the loop
            Halfedge *inner = current_inner[i];
            if(!inner->twin()->flagged) {
                // we have found an inner (possible) edge which we only saw from one direction
                // this means that we have encountered a non-simple face.
                // we could be in different situations:
                //  - we traversed the outer boundary of a face with a hole or
                //    an isolated vertex, with the target of the inner edge isolated or on the hole
                //  - we traversed a hole inside a face, and the target of the inner edge is on the face
                //    boundary, or isolated, or on another hole in the same face
                if(orientation == CGAL::ZERO) {
                    // compute the orientation
                    orientation = boundary_orientation(m_current_face.boundary);
                }

                if(orientation == CGAL::COUNTERCLOCKWISE) {
                    // we already traversed the (outer) boundary of the face
                    // traverse the hole (or identify the point as isolated vertex)
                    std::size_t old_size = m_current_face.hole_boundaries.size();
                    p_traverse_boundary_and_collect_edges(inner->twin(), m_current_face.hole_boundaries, current_inner);
                    std::size_t new_size = m_current_face.hole_boundaries.size();
                    if(new_size == old_size) {
                        // we have found an isolated vertex
                        m_current_face.isolated_vertices.push_back(inner->target_handle());
                    }
                } else {
                    // we traversed a hole; we must traverse the outer boundary eventually,
                    // because one of the inner edges must lead us there. what we collected
                    // is the boundary of a hole inside the face (which we have not yet found so far)
                    m_current_face.hole_boundaries.insert(m_current_face.hole_boundaries.end(),
                                                          m_current_face.boundary.begin(),
                                                          m_current_face.boundary.end());
                    m_current_face.boundary.clear();
                    orientation = CGAL::ZERO;
                    p_traverse_boundary_and_collect_edges(inner->twin(), m_current_face.boundary, current_inner);
                    if(m_current_face.boundary.empty()) {
                        // we have found an isolated vertex
                        m_current_face.isolated_vertices.push_back(inner->target_handle());
                    }
                }
            }
        }
    }

    /**
     * Traverse the boundary of a face bounded by the halfedge start.
     * Collect the boundary edges in collect_boundary and the inner edges in collect_inner.
     * The face is traversed by following the Halfedge::next pointers, which take us
     * (counterclockwise) around each vertex, enumerating all the possible and certain edges;
     * 'inner edges' are possible edges with status == LMTStatus::Possible, whereas boundary
     * edges are always certain edges of the MWT (status == LMTStatus::Certain or LMTStatus::CH).
     * If the boundary order returned by this method is CCW, then the inner edges are on the
     * inside of the face; if the boundary order is CW, then the inner edges are on the outside.
     */
    void p_traverse_boundary_and_collect_edges(Halfedge *start, std::vector<Halfedge *> &collect_boundary,
                                               std::vector<Halfedge *> &collect_inner) {
        std::size_t old_size = collect_boundary.size();
        start->flagged = true;
        if(start->status == LMTStatus::Possible) {
            collect_inner.push_back(start);
        }
        Halfedge *next = start;
        do {
            next = next->next_edge();

            // skip over possible edges (collect them in collect_inner)
            while(next->status == LMTStatus::Possible && next != start) {
                collect_inner.push_back(next);
                next->flagged = true;
                next = next->next_edge();
            }

            // either all edges are possible (isolated vertex), or we
            // have found the (twin of the) next boundary edge
            if(next->status != LMTStatus::Possible) {
                next = next->twin();
                next->flagged = true;
                collect_boundary.push_back(next);
            }
        } while(next != start);

        // change to standard orientation: in our mode,
        // a CCW boundary traversal corresponds to looking
        // at 'inner edges' on the outside of the face,
        // and a CW traversal corresponds to looking at
        // 'inner edges' on the inside of the face.
        // this goes against standard assumptions, so we
        // reverse the order here
        std::reverse(collect_boundary.begin() + old_size, collect_boundary.end());
    }

    Halfedge **m_skeleton_begin;
    Halfedge **m_skeleton_end;
    Halfedge **m_skeleton_current;
    Face m_current_face;
};

} // namespace mwt

#endif
