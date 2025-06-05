#ifndef CGAL_MWT_ARE_ALL_COLLINEAR_H_INCLUDED_
#define CGAL_MWT_ARE_ALL_COLLINEAR_H_INCLUDED_

#include <CGAL/Kernel/global_functions.h>
#include <algorithm>
#include <utility>

namespace mwt {

/**
 * Check if all points pointed to by the given iterators are collinear.
 */
template<typename PointIterator> bool are_all_points_collinear(PointIterator begin, PointIterator end) {
    if(std::distance(begin, end) < 3) {
        return true;
    }
    auto p0 = *begin++;
    auto p1 = *begin++;
    for(; begin != end; ++begin) {
        if(!CGAL::collinear(p0, p1, *begin)) {
            return false;
        }
    }
    return true;
}

} // namespace mwt

#endif
