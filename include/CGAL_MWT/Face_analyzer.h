#ifndef CGAL_MWT_FACE_ANALYZER_H_INCLUDED_
#define CGAL_MWT_FACE_ANALYZER_H_INCLUDED_

#include "Face_collector.h"
#include "Unique_hash_map.h"

namespace mwt {

/**
 * Option struct passed to all face triangulators.
 */
struct FaceTriangulatorOptions {
    bool lp_verbose;
    bool use_cuts;
};

/**
 * Simple struct representing an empty triangle in a face.
 */
template<typename Halfedge_> struct Face_triangle {
    using Halfedge = Halfedge_;
    using EdgeHandle = const Halfedge *;
    using MutableEdgeHandle = Halfedge *;

    /**
     * The edges making up this triangle,
     * in CCW order around the triangle (interior of the triangle to the left).
     * The first edge is the one with the smallest pointer.
     */
    MutableEdgeHandle edges[3];

    /**
     * Compute the total weight of the possible
     * (not certain/convex hull) edges of this triangle.
     */
    template<typename FT> FT possible_weight() const noexcept {
        FT result = 0;
        for(EdgeHandle edge : edges) {
            if(edge->status == LMTStatus::Possible) {
                result += edge->template weight<FT>();
            }
        }
        return result;
    }
};

/**
 * Class that analyzes a face before it is triangulated
 * using either DP or LP/MIP.
 */
template<typename Skeleton_> class Face_analyzer {
  public:
    using Skeleton = Skeleton_;
    using Halfedge = typename Skeleton::Halfedge;
    using EdgeHandle = const Halfedge *;
    using MutableEdgeHandle = Halfedge *;
    using Face = mwt::Face<Skeleton>;

    struct EdgeInfo {
        /**
         * true iff a triangle to the left of this edge
         * has to be selected, i.e., if the edge is a boundary
         * edge, hole edge or an antenna reaching into the face.
         */
        // bool must_have_left_triangle;

        /**
         * Indices of empty (or almost-always-empty) triangles
         * to the left of this edge.
         */
        std::vector<std::size_t> left_triangles;
    };

    using Triangle = Face_triangle<Halfedge>;

    Face_analyzer() {}

    void clear() {
        m_edge_info.clear();
        m_empty_triangles.clear();
    }

    Triangle &get_triangle(std::size_t triangle_id) noexcept { return m_empty_triangles[triangle_id]; }

    const Triangle &get_triangle(std::size_t triangle_id) const noexcept { return m_empty_triangles[triangle_id]; }

    void analyze_face(const Face &face) {
        clear();
        m_face = &face;
        p_enter_edges(face);
        p_find_triangles(face);
    }

    const std::vector<Triangle> &empty_triangles() const noexcept { return m_empty_triangles; }

    const Unique_hash_map<EdgeHandle, EdgeInfo> &edge_info() const noexcept { return m_edge_info; }

    const Face &face() const noexcept { return *m_face; }

  private:
    /**
     * Enter all halfedges in this face into the edge_info map.
     */
    void p_enter_edges(const Face &face) {
        m_edge_info.reserve(face.boundary.size() + face.hole_boundaries.size() + face.inner_edges.size());
        for(EdgeHandle edge : face.boundary) {
            m_edge_info[edge] = EdgeInfo{};
        }
        for(EdgeHandle edge : face.hole_boundaries) {
            m_edge_info[edge] = EdgeInfo{};
        }
        for(EdgeHandle edge : face.inner_edges) {
            m_edge_info[edge] = EdgeInfo{};
        }
    }

    void p_find_triangles(const Face &face) {
        auto treat_edge = [&](auto *edge) {
            edge->reset();
            while(edge->next_triangle()) {
                p_add_triangle_if_primary(edge);
            }
        };
        for(auto *edge : face.boundary) {
            treat_edge(edge);
        }
        for(auto *edge : face.hole_boundaries) {
            treat_edge(edge);
        }
        for(auto *edge : face.inner_edges) {
            treat_edge(edge);
        }
    }

    /**
     * Add the triangle to the left of the given edge that is
     * currently stored in the halfedge datastructure to our
     * triangle datastructure if this edge is the primary edge
     * of that triangle.
     */
    void p_add_triangle_if_primary(MutableEdgeHandle edge) {
        MutableEdgeHandle e2 = edge->j;
        MutableEdgeHandle e3 = edge->i->twin();
        if(edge < e2 && edge < e3) { // check that edge is the primary edge of the triangle
            auto i1 = m_edge_info.lookup(edge);
            auto i2 = m_edge_info.lookup(e2);
            auto i3 = m_edge_info.lookup(e3);
            if(!i1 || !i2 || !i3)
                return; // not all edges are part of the face/possible in the face
            std::size_t index = m_empty_triangles.size();
            m_empty_triangles.emplace_back(Triangle{edge, e2, e3});
            m_edge_info[i1].left_triangles.push_back(index);
            m_edge_info[i2].left_triangles.push_back(index);
            m_edge_info[i3].left_triangles.push_back(index);
        }
    }

    const Face *m_face{nullptr};
    Unique_hash_map<EdgeHandle, EdgeInfo> m_edge_info;
    std::vector<Triangle> m_empty_triangles;
};

} // namespace mwt

#endif
