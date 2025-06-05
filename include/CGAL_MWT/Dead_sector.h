#ifndef CGAL_MWT_DEAD_SECTOR_H_INCLUDED_
#define CGAL_MWT_DEAD_SECTOR_H_INCLUDED_

#include "Is_interval.h"
#include "Mwt_traits.h"
#include <CGAL/Uncertain.h>
#include <tuple>

namespace mwt {

/**
 * Given the two values, if FT is an interval type,
 * return the interval [lb, ub], otherwise return lb as FT.
 */
template<typename FT> FT interval_or_lb(double lb, double ub) {
    if constexpr(is_interval_number_type_v<FT>) {
        return FT(lb, ub);
    } else {
        return FT(lb);
    }
}

/**
 * Given the two values, if FT is an interval type,
 * return the interval [lb, ub], otherwise return ub as FT.
 */
template<typename FT> FT interval_or_ub(double lb, double ub) {
    if constexpr(is_interval_number_type_v<FT>) {
        return FT(lb, ub);
    } else {
        return FT(ub);
    }
}

template<typename IntervalType, std::enable_if_t<is_interval_number_type_v<IntervalType>, int> = 0>
static inline auto inf_or_value(IntervalType nt) noexcept {
    return nt.inf();
}

template<typename IntervalType, std::enable_if_t<is_interval_number_type_v<IntervalType>, int> = 0>
static inline auto sup_or_value(IntervalType nt) noexcept {
    return nt.sup();
}

template<typename NonIntervalType, std::enable_if_t<!is_interval_number_type_v<NonIntervalType>, int> = 0>
static inline const NonIntervalType &inf_or_value(const NonIntervalType &nt) noexcept {
    return nt;
}

template<typename NonIntervalType, std::enable_if_t<!is_interval_number_type_v<NonIntervalType>, int> = 0>
static inline const NonIntervalType &sup_or_value(const NonIntervalType &nt) noexcept {
    return nt;
}

/**
 * Base function to compute dead sectors.
 * The result is an interval of pseudoangles
 * in which points can safely be filtered out.
 * In case there is no sector, the LB will be above the UB.
 * The resulting tuple has a third entry (bool),
 * which describes whether the sector is bounded by the two vertices.
 */
template<typename FT>
inline std::tuple<FT, FT, bool> compute_sector(const FT &px, const FT &py, const FT &lx, const FT &ly, const FT &rx,
                                               const FT &ry) {
    std::tuple<FT, FT, bool> result{FT(1), FT(0), false};
    const FT vlx = lx - px, vly = ly - py;
    const FT vrx = rx - px, vry = ry - py;
    const FT dot = vlx * vrx + vly * vry;
    if(CGAL::possibly(dot <= 0))
        return result;

    // cos^2( <(l,r) ) * (||l|| * ||r||)^2 = (l * r)^2 = dot^2
    const FT dot_squared = dot * dot;
    const FT nl = vlx * vlx + vly * vly; // ||l||
    const FT nr = vrx * vrx + vry * vry; // ||r||
    const FT denom = (nl * nr);          // ||l|| * ||r||
    if(CGAL::possibly(dot_squared <= detail::SQUARED_COS_2DPA * denom)) {
        // in this case, we cannot be certain a sector exists.
        // it is fine to return an empty sector here.
        return result;
    }

    // check whether the angle is <= alpha
    Compute_pseudo_angle<FT> compute_angle{};
    auto comp = (dot_squared >= detail::SQUARED_COS_DPA * denom);
    if(CGAL::certainly(comp)) {
        // angle <= alpha: The dead sector is bounded by the two vertices
        std::get<2>(result) = true;
        std::get<0>(result) = compute_angle(vrx, vry);
        std::get<1>(result) = compute_angle(vlx, vly);
    } else {
        // alpha < angle <= 2*alpha: We have to rotate
        // the vector vl by alpha to the right to get the upper bound
        // and the vector vr by alpha to the left to get the lower bound.
        // Rotate right vector to the left (ccw)
        const auto cos_dpa = interval_or_ub<FT>(detail::COS_DPA_LB, detail::COS_DPA);
        const auto sin_dpa = interval_or_lb<FT>(detail::SIN_DPA, detail::SIN_DPA_UB);
        const FT slx = cos_dpa * vrx - sin_dpa * vry;
        const FT sly = sin_dpa * vrx + cos_dpa * vry;
        const FT srx = cos_dpa * vlx + sin_dpa * vly;
        const FT sry = cos_dpa * vly - sin_dpa * vlx;
        std::get<0>(result) = compute_angle(srx, sry);
        std::get<1>(result) = compute_angle(slx, sly);
        if constexpr(is_interval_number_type_v<FT>) {
            if(!CGAL::certainly(!comp)) {
                // we are unsure whether the angle is <= alpha or > alpha;
                // make sure to bound the above result by the two vertices' angles
                const FT lb = compute_angle(vrx, vry);
                const FT ub = compute_angle(vlx, vly);
                std::get<0>(result) = FT((std::max)(std::get<0>(result).inf(), lb.inf()),
                                         (std::max)(std::get<0>(result).sup(), lb.sup()));
                std::get<1>(result) = FT((std::min)(std::get<1>(result).inf(), ub.inf()),
                                         (std::min)(std::get<1>(result).sup(), ub.sup()));
            }
        }
    }

    if constexpr(!CGAL::Algebraic_structure_traits<FT>::Is_exact::value) {
        // check that what we are returning is an interval with LB <= UB;
        // otherwise, clear the result and return an empty sector.
        if(!CGAL::certainly(std::get<0>(result) <= std::get<1>(result))) {
            std::get<0>(result) = FT(1);
            std::get<1>(result) = FT(0);
            std::get<2>(result) = false;
        }
    }
    return result;
}

} // namespace mwt

#endif
