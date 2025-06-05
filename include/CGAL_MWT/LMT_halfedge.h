#ifndef CGAL_MWT_LMT_HALFEDGE_H
#define CGAL_MWT_LMT_HALFEDGE_H

#include <CGAL/Algebraic_structure_traits.h>
#include <CGAL/assertions.h>
#include <CGAL/enum.h>
#include <CGAL/number_utils.h>
#include <cmath>
#include <iostream>

namespace mwt {

/**
 * Enum that describes the status of an edge in the LMT.
 */
enum class LMTStatus : char {
    Possible = 0,   //< the edge is currently still possible
    Certain = 1,    //< the edge is certain (i.e., definitely in the MWT)
    Impossible = 2, //< the edge is impossible (i.e., definitely not in the MWT)
    CH = 3          //< the edge is on the convex hull (i.e., definitely in the MWT)
};

inline std::ostream &operator<<(std::ostream &output, LMTStatus status) {
    switch(status) {
    case LMTStatus::Possible:
        return output << "possible";
    case LMTStatus::Certain:
        return output << "certain";
    case LMTStatus::Impossible:
        return output << "impossible";
    case LMTStatus::CH:
        return output << "ch";
    default:
        return output << "unknown";
    }
}

/**
 * Halfedge struct.
 */
template<typename Traits_> struct LMTHalfedge {
    using Traits = Traits_;
    using FT = typename Traits::FT;
    using Point_2 = typename Traits::Point_2;
    using Left_turn_2 = typename Traits::Left_turn_2;
    using Compute_squared_distance_2 = typename Traits::Compute_squared_distance_2;
    using Less_distance_2 = typename Traits::Less_distance_2;
    using Quadrilateral_convex = typename Traits::Quadrilateral_convex;

    /**
     * Default halfedge constructor, constructs an invalid halfedge.
     */
    LMTHalfedge() = default;

    /**
     * Create a new halfedge targeting a point.
     */
    LMTHalfedge(const Point_2 *target, LMTHalfedge *twn) noexcept
        : tar{target}, twn{twn}, next{nullptr}, i{nullptr}, j{nullptr}, j_start{nullptr}, status{LMTStatus::Possible},
          flagged{true} {}

    /**
     * Find the next (most likely empty) triangle on the left side of this halfedge.
     * In a full edge set the triangle would be guaranteed to be empty. However, the candidate edges are prefiltered
     * by the diamond test and there is a chance (due to a missing edge to a nearby point) that we find a non-empty
     * triangle. The overall algorithm is still correct in the sense that it produces a subset of the MWT.
     * Additionally, it seems that all of the non-empty triangles get eliminated anyway by the end of the algorithm.
     */
    bool next_triangle() {
        CGAL_precondition(status != LMTStatus::Impossible);

        if(j == nullptr) {
            // We already iterated through all triangles in an earlier call.
            return false;
        }

        // Break the current triangle.
        // (On the first call after a reset it might actually create the very first valid triangle)
        j = j->next_edge();

        while(j->target_handle() != source_handle()) {
            // Rotate the left edge until we hit the same target (triangle found) or pass beyond it.
            while(i->target_handle() != j->target_handle() && !j->has_on_left_side(i->target())) {
                i = i->next_edge();
            }

            // Now do the same with the right edge.
            // If we hit the source vertex, we went so far that there are no more triangles.
            while(j->target_handle() != i->target_handle() && j->target_handle() != source_handle() &&
                  j->has_on_left_side(i->target())) {
                j = j->next_edge();
            }

            // Repeat the above process until a triangle is found or one of edges went too far
            if(j->target_handle() == i->target_handle()) {
                CGAL_assertion(has_on_left_side(i->target()));
                return true;
            } else if(!has_on_left_side(i->target())) {
                j = nullptr;
                return false;
            }
        }

        // Right edge went too far
        j = nullptr;
        return false;
    }

    void reset() {
        CGAL_precondition(j_start != nullptr);
        CGAL_precondition(status != LMTStatus::Impossible);
        i = next_edge();
        CGAL_postcondition(i != this);
        j = j_start;
    }

    // Tries to find the next valid certificate for this edge.
    // A certificate is a pair of triangles that form a quadrilateral q such that q is either non-convex or this
    // edge is the shorter diagonal in q.
    bool find_certificate() {
        CGAL_precondition(is_primary());
        CGAL_precondition(status != LMTStatus::Impossible);

        if(!has_certificate()) { // Called for the first time for this edge
            reset();
            twin()->reset();
            if(!twin()->next_triangle())
                return false;
        }

        // After restacking this edge, at least one of the triangles is invalid.
        // Special treatment is needed if the triangle on the twin side is invalid.
        if(!twin()->has_valid_triangle()) {
            if(twin()->next_triangle())
                reset();
            else
                return false;
        }

        // Try all quadrilaterals, i.e. all pairs of triangles on opposite sides of this edge.
        while(true) {
            if(!next_triangle()) {
                if(twin()->next_triangle()) {
                    reset();
                } else {
                    return false;
                }
            } else {
                CGAL_assertion(has_valid_triangle() && twin()->has_valid_triangle());
                // A quadrilateral is a certificate iff this edge is the shorter diagonal or it is not convex.
                if(!Less_distance_2{}(twin()->i->target(), i->target(), source(), target()) ||
                   !Quadrilateral_convex{}(source(), twin()->i->target(), target(), i->target())) {
                    return true;
                }
            }
        }
    }

    // TODO: Can probably be improved. Checks the same triangles multiple times for local minimality.
    bool find_local_minimal_certificate() {
        CGAL_precondition(is_primary());
        CGAL_precondition(status != LMTStatus::Impossible);
        while(!has_valid_triangle() || !twin()->has_valid_triangle() || !has_local_minimal_triangle() ||
              !twin()->has_local_minimal_triangle()) {
            if(!find_certificate()) {
                return false;
            }
        }

        return true;
    }

    // Restack those edges whose local minimal certificate might get invalidated when this edge becomes impossible.
    // It is somewhat difficult to figure out the exact set of edges that need to be restacked.
    // It is a subset of the neighborhood, i.e., the outgoing edges and their outgoing edges.
    // We restack them all.
    void restack_neighborhood(std::vector<LMTHalfedge *> &stack) const {
        CGAL_precondition(status != LMTStatus::Impossible);

        // Neighborhood of source vertex
        auto e = next_edge();
        while(e != this) {
            auto ne = e->twin();
            do {
                if(!ne->flagged && ne->status == LMTStatus::Possible) {
                    stack.push_back(ne->primary_edge());
                    ne->flagged = ne->twn->flagged = true;
                }
                ne = ne->next_edge();
            } while(ne != e->twin());
            e = e->next_edge();
        }

        // Neighborhood of target vertex
        e = twn->next_edge();
        while(e != twn) {
            auto ne = e->twin();
            do {
                if(!ne->flagged && ne->status == LMTStatus::Possible) {
                    stack.push_back(ne->primary_edge());
                    ne->flagged = ne->twn->flagged = true;
                }
                ne = ne->next_edge();
            } while(ne != e->twin());
            e = e->next_edge();
        }
    }

    // A triangle t is local minimal iff all three edges e of t meet either condition:
    // (i) e is on the convex hull
    // (ii) e is local minimal w.r.t a certificate that has t in common
    bool has_local_minimal_triangle() {
        CGAL_precondition(status != LMTStatus::Impossible); // This edge is local minimal
        CGAL_precondition(has_valid_triangle());

        // Make copies of i and j and try to find a certificate for them
        bool left_valid{false};
        bool right_valid{false};

        if(i->status == LMTStatus::CH) {
            left_valid = true;
        } else {
            LMTHalfedge left = *i;
            left.reset();

            while(!left_valid && left.next_triangle()) {
                left_valid |= !Less_distance_2{}(target(), left.i->target(), left.source(), left.target()) ||
                              !Quadrilateral_convex{}(left.source(), target(), left.target(), left.i->target());
            }
        }

        if(j->status == LMTStatus::CH) {
            right_valid = true;
        } else {
            LMTHalfedge right = *j->twin();
            right.reset();

            while(left_valid && !right_valid && right.next_triangle()) {
                right_valid |= !Less_distance_2{}(source(), right.i->target(), right.source(), right.target()) ||
                               !Quadrilateral_convex{}(right.source(), source(), right.target(), right.i->target());
            }
        }

        return left_valid && right_valid;
    }

    // Check all neighbors of the source vertex on the left side of this edge for an outgoing intersecting edge
    bool has_crossing_edges_on_left_side() const {
        auto n = next_edge();
        while(has_on_left_side(n->target())) {
            auto edge = n->twin()->next_edge();

            CGAL_assertion(has_on_left_side(edge->source()));

            while(edge->target_handle() != target_handle() && edge->has_on_left_side(target())) {
                // This holds because the edges are sorted by angle and the edge set
                // always contains a triangulation.
                // There should always be an edge that either connects to the target_handle or has the target on
                // the right side before there is an edge that wraps around so much that this would fail.
                CGAL_assertion(edge->has_on_right_side(source()));

                if(has_on_right_side(edge->target())) {
                    return true;
                }
                edge = edge->next_edge();

                // Same reason as above
                CGAL_assertion(edge->target_handle() != source_handle());
            }
            n = n->next_edge();
        }
        return false;
    }

    bool has_crossing_edges() const {
        // The candidate edges contain a valid triangulation at any point in time. Therefore, it suffices to
        // check the neighbors for outgoing intersecting edges. Indeed it suffices to check only the neighbors of
        // one side of one vertex.
        return has_crossing_edges_on_left_side();
    }

    // Restacks exactly those edges whose certificates become invalid when this edge becomes impossible.
    void restack_edges(std::vector<LMTHalfedge *> &stack) const {
        CGAL_precondition(status != LMTStatus::Impossible);

        auto e = next_edge();
        while(e != this) {
            if(!e->flagged && e->status == LMTStatus::Possible) {
                if(e->twn->j == this || e->i == this) {
                    stack.push_back(e->primary_edge());
                    e->flagged = e->twn->flagged = true;
                }
            }
            e = e->next_edge();
        }

        e = twn->next_edge();
        while(e != twn) {
            if(!e->flagged && e->status == LMTStatus::Possible) {
                if(e->twn->j == twn || e->i == twn) {
                    stack.push_back(e->primary_edge());
                    e->flagged = e->twn->flagged = true;
                }
            }
            e = e->next_edge();
        }
    }

    LMTHalfedge *next_edge() const { return next; }

    void unlink(LMTHalfedge *b, LMTHalfedge *e) {
        CGAL_precondition(status == LMTStatus::Impossible);
        CGAL_precondition(b != e);
        CGAL_precondition(next != this);

        --e;
        CGAL_assertion(b->source_handle() == source_handle());
        CGAL_assertion(e->source_handle() == source_handle());

        auto prev = this - 1;
        if(prev < b)
            prev = e;

        do {
            prev->next = next;
            --prev;
            if(prev < b)
                prev = e;
        } while(prev->next == this);
    }

    LMTHalfedge *primary_edge() { return is_primary() ? this : twin(); }
    const LMTHalfedge *primary_edge() const { return is_primary() ? this : twin(); }

    // Define one of the halfedges of each pair as the primary one. The choice is more or less arbitrary.
    bool is_primary() const { return source_handle() < target_handle(); }

    const Point_2 &target() const { return *tar; }

    const Point_2 &source() const { return twn->target(); }

    const Point_2 *target_handle() const { return tar; }

    const Point_2 *source_handle() const { return twn->tar; }

    const LMTHalfedge *twin() const {
        CGAL_assertion(twn != nullptr);
        return twn;
    }

    LMTHalfedge *twin() {
        CGAL_assertion(twn != nullptr);
        return twn;
    }

    void set_twin(LMTHalfedge *twin) { twn = twin; }

    bool has_certificate() const {
        bool b = ((i != nullptr) & (j != nullptr));
        CGAL_assertion(!b || (twin()->i != nullptr && twin()->j != nullptr));
        return b;
    }

    bool has_valid_triangle() const {
        return i != nullptr && j != nullptr && i->status != LMTStatus::Impossible &&
               j->status != LMTStatus::Impossible && i->target_handle() == j->target_handle();
    }

    bool has_on_left_side(const Point_2 &p) const { return Left_turn_2{}(source(), target(), p); }

    bool has_on_right_side(const Point_2 &p) const { return Left_turn_2{}(target(), source(), p); }

    bool is_certain() const { return this->status == LMTStatus::Certain; }

    FT length() const {
        return typename CGAL::Algebraic_structure_traits<FT>::Sqrt{}(Compute_squared_distance_2{}(source(), target()));
    }

    template<typename NumType> NumType weight() const {
        const auto &s = source();
        const auto &t = target();
        NumType dx = s.x();
        dx -= t.x();
        NumType dy = s.y();
        dy -= t.y();
        dx *= dx;
        dy *= dy;
        dx += dy;
        if constexpr(std::is_floating_point<NumType>::value) {
            return std::sqrt(dx);
        } else {
            return CGAL::sqrt(dx);
        }
    }

    const Point_2 *tar;   //< Target of this halfedge.
    LMTHalfedge *twn;     //< Twin of this halfedge.
    LMTHalfedge *next;    //< The next outgoing halfedge of our source vertex that is not impossible in ccw order.
    LMTHalfedge *i;       //< The left halfedge of the certificate.
    LMTHalfedge *j;       //< The right halfedge of the certificate.
    LMTHalfedge *j_start; //< Entry point for the triangle search.
    LMTStatus status;     //< The status of this halfedge.
    bool flagged;         //< Whether this halfedge is flagged (e.g., on the stack or visited).
};

} // namespace mwt

#endif // MWT_LMT_HALFEDGE_H
