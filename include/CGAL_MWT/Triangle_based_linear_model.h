#ifndef CGAL_MWT_TRIANGLE_BASED_LINEAR_MODEL_H_INCLUDED_
#define CGAL_MWT_TRIANGLE_BASED_LINEAR_MODEL_H_INCLUDED_

#include "Face_analyzer.h"

namespace mwt {

template<typename Scalar, typename Skeleton, typename AddCallback /*(Scalar coeff, std::size_t var_index)*/,
         typename FinalizeCallback /*(Scalar rhs)*/>
inline void create_possible_halfedge_constraint(const Face_analyzer<Skeleton> &analyzer,
                                                const typename Skeleton::Halfedge *halfedge, AddCallback &&add_callback,
                                                FinalizeCallback &&finalize_callback) {
    if(!halfedge->is_primary()) {
        // for possible (as opposed to certain) halfedges,
        // we only need a constraint for one direction; so
        // we only add a constraint for the 'primary' direction
        return;
    }

    auto *t = halfedge->twin();
    const auto &analyzer_info = analyzer.edge_info();
    const auto &p_info = analyzer_info[halfedge];
    const auto &t_info = analyzer_info[t];

    for(std::size_t tind : p_info.left_triangles) {
        add_callback(1, tind);
    }
    for(std::size_t tind : t_info.left_triangles) {
        add_callback(-1, tind);
    }
    finalize_callback(0);
}

template<typename Scalar, typename Skeleton, typename AddCallback /*(Scalar coeff, std::size_t var_index)*/,
         typename FinalizeCallback /*(Scalar rhs)*/>
inline void create_certain_halfedge_constraint(const Face_analyzer<Skeleton> &analyzer,
                                               const typename Skeleton::Halfedge *halfedge, AddCallback &&add_callback,
                                               FinalizeCallback &&finalize_callback) {
    const auto &analyzer_info = analyzer.edge_info();
    const auto &h_info = analyzer_info[halfedge];
    for(std::size_t tind : h_info.left_triangles) {
        add_callback(1, tind);
    }
    finalize_callback(1);
}

template<typename Scalar, typename Skeleton, typename AddCallback /*(Scalar coeff, std::size_t var_index)*/,
         typename FinalizeCallback /*(Scalar rhs)*/>
inline void create_halfedge_constraint(const Face_analyzer<Skeleton> &analyzer,
                                       const typename Skeleton::Halfedge *halfedge, AddCallback &&add_callback,
                                       FinalizeCallback &&finalize_callback) {
    if(halfedge->status == LMTStatus::Possible) {
        create_possible_halfedge_constraint<Scalar>(analyzer, halfedge, add_callback, finalize_callback);
        return;
    }
    create_certain_halfedge_constraint<Scalar>(analyzer, halfedge, add_callback, finalize_callback);
}

} // namespace mwt

#endif
