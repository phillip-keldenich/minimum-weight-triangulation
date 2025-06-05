#ifndef CGAL_MWT_TOTAL_WEIGHT_H_INCLUDED_
#define CGAL_MWT_TOTAL_WEIGHT_H_INCLUDED_

#include "LMT_skeleton.h"
#include <CGAL/Delaunay_triangulation_2.h>

namespace mwt {

/**
 * Compute the total weight of the given Delaunay triangulation,
 * using the given NumberType to calculate each individual edge weight
 * and the total edge weight.
 */
template<typename NumberType, typename DT, typename TDS>
NumberType total_weight(const CGAL::Delaunay_triangulation_2<DT, TDS> &triangulation) {
    NumberType total(0);
    for(auto it = triangulation.finite_edges_begin(); it != triangulation.finite_edges_end(); ++it) {
        auto seg = triangulation.segment(it);
        NumberType x1(seg.source().x());
        NumberType x2(seg.target().x());
        NumberType y1(seg.source().y());
        NumberType y2(seg.target().y());
        x1 -= x2;
        y1 -= y2;
        x1 *= x1;
        y1 *= y1;
        x1 += y1;
        total += CGAL::sqrt(x1);
    }
    return total;
}

template<typename NumberType, typename Traits, typename Tree>
NumberType total_weight(const LMTSkeleton<Traits, Tree> &skeleton) {
    NumberType total(0);
    for(auto it = skeleton.halfedges_begin(); it != skeleton.halfedges_end(); ++it) {
        if(!it->is_primary() || it->status == LMTStatus::Impossible) {
            continue;
        }
        if(it->status == LMTStatus::Possible) {
            throw std::logic_error("There should not be any possible edges left!");
        }
        auto s = it->source();
        auto t = it->target();
        NumberType x1(s.x()), x2(t.x());
        NumberType y1(s.y()), y2(t.y());
        x1 -= x2;
        y1 -= y2;
        x1 *= x1;
        y1 *= y1;
        x1 += y1;
        total += CGAL::sqrt(x1);
    }
    return total;
}

} // namespace mwt

#endif
