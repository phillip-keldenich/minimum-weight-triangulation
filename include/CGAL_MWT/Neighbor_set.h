#ifndef CGAL_MWT_NEIGHBOR_SET_H_INCLUDED_
#define CGAL_MWT_NEIGHBOR_SET_H_INCLUDED_

#include <boost/container/flat_set.hpp>

namespace mwt {

template<typename Kernel_, typename PointIterator_> struct Point_angle_pair {
    using Kernel = Kernel_;
    using Point_const_iterator = PointIterator_;
    using FT = typename Kernel::FT;

    Point_const_iterator point;
    FT polar_angle;

    Point_angle_pair() = default;
    Point_angle_pair(Point_const_iterator p, FT angle) : point{p}, polar_angle{angle} {}

    bool operator<(const Point_angle_pair &rhs) const { return polar_angle < rhs.polar_angle; }
    bool operator<=(const Point_angle_pair &rhs) const { return polar_angle <= rhs.polar_angle; }
};

template<typename Kernel, typename PointIterator>
using Point_angle_set = boost::container::flat_set<Point_angle_pair<Kernel, PointIterator>>;

} // namespace mwt

#endif
