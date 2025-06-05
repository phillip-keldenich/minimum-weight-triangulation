#ifndef CGAL_MWT_QUADTREE_AUX_H_INCLUDED_
#define CGAL_MWT_QUADTREE_AUX_H_INCLUDED_

#include <CGAL/Iso_rectangle_2.h>
#include <CGAL/Point_2.h>

namespace mwt {

//--------------------------------------------------------------------------------------------------------------
// Some helper enums and structs that are used to Hilbert sort the point set while building the quadtree
//--------------------------------------------------------------------------------------------------------------
enum class Order { XY, YX };
enum class Dimension { X, Y };

template<Order> struct Split;
template<> struct Split<Order::XY> {
    static constexpr Order same = Order::XY;
    static constexpr Order other = Order::YX;
};
template<> struct Split<Order::YX> {
    static constexpr Order same = Order::YX;
    static constexpr Order other = Order::XY;
};

/**
 * Comparison predicate (compare dimension (x or y) of
 * points to some fixed value), either < (reverse_order false)
 * or >= (reverse_order true).
 */
template<class K, Dimension d, bool reverse_order> struct Fixed_value_cmp_2;

/**
 * Reversed-order implementation.
 */
template<class K, Dimension d> struct Fixed_value_cmp_2<K, d, true> {
    using FT = typename K::FT;
    using Point_2 = typename K::Point_2;

    FT value;

    Fixed_value_cmp_2(const FT &v) : value(v) {}

    bool operator()(const Point_2 &p) const { return !Fixed_value_cmp_2<K, d, false>{value}(p); }
};

/**
 * Non-reversed-order comparison on x-coordinate.
 */
template<class K> struct Fixed_value_cmp_2<K, Dimension::X, false> {
    using FT = typename K::FT;
    using Point_2 = typename K::Point_2;

    FT value;

    Fixed_value_cmp_2(const FT &v) : value(v) {}
    bool operator()(const Point_2 &p) const { return p.x() < value; }
};

/**
 * Non-reversed-order comparison on y-coordinate.
 */
template<class K> struct Fixed_value_cmp_2<K, Dimension::Y, false> {
    using FT = typename K::FT;
    using Point_2 = typename K::Point_2;

    FT value;

    Fixed_value_cmp_2(const FT &v) : value(v) {}

    bool operator()(const Point_2 &p) const { return p.y() < value; }
};

/**
 * Compare points by coordinate.
 */
template<class K, Dimension d, bool reverse_order> struct Coordinate_cmp_2;
template<class K, Dimension d> struct Coordinate_cmp_2<K, d, true> {
    using FT = typename K::FT;
    using Point_2 = typename K::Point_2;

    bool operator()(const Point_2 &p, const Point_2 &q) const { return Coordinate_cmp_2<K, d, false>{}(q, p); }
};

/**
 * Compare points by x-coordinate.
 */
template<class K> struct Coordinate_cmp_2<K, Dimension::X, false> {
    using FT = typename K::FT;
    using Point_2 = typename K::Point_2;

    bool operator()(const Point_2 &p, const Point_2 &q) const { return p.x() < q.x(); }
};

/**
 * Compare points by y-coordinate.
 */
template<class K> struct Coordinate_cmp_2<K, Dimension::Y, false> {
    using FT = typename K::FT;
    using Point_2 = typename K::Point_2;

    bool operator()(const Point_2 &p, const Point_2 &q) const { return p.y() < q.y(); }
};

/**
 * Compute the squared distance from a point to a rectangle.
 */
template<typename Kernel>
auto squared_distance(const CGAL::Point_2<Kernel> &p, const CGAL::Iso_rectangle_2<Kernel> &q) {
    using FT = typename Kernel::FT;
    FT result{0};
    if(p.x() < q.xmin()) {
        FT dx{q.xmin() - p.x()};
        result += dx * dx;
    } else if(p.x() > q.xmax()) {
        FT dx{p.x() - q.xmax()};
        result += dx * dx;
    }
    if(p.y() < q.ymin()) {
        FT dy{q.ymin() - p.y()};
        result += dy * dy;
    } else if(p.y() > q.ymax()) {
        FT dy{p.y() - q.ymax()};
        result += dy * dy;
    }
    return result;
}

} // namespace mwt

#endif
