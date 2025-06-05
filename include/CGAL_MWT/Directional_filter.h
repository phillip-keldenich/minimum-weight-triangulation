#ifndef CGAL_MWT_DIRECTIONAL_FILTER_H_INCLUDED_
#define CGAL_MWT_DIRECTIONAL_FILTER_H_INCLUDED_

#include "Dead_sector.h"
#include "Interval_set.h"
#include "Is_interval.h"
#include "Static_quadtree.h"
#include <CGAL/assertions.h>
#include <boost/container/flat_set.hpp>
#include <cassert>

namespace mwt {

namespace detail {

/**
 * The value for which we insert [-2, -DIRECTIONAL_FILTER_PANGLE_CUTOFF]
 * and [DIRECTIONAL_FILTER_PANGLE_CUTOFF, 2] into the dead sector set
 * from the start to remove:
 *  - the discontinuity at pi
 *  - duplicate reporting of points that we have already seen before
 */
static constexpr double DIRECTIONAL_FILTER_PANGLE_CUTOFF = 1.25;

/**
 * The pseudoangle value from which we stop considering points even
 * for the purpose of updating our filter.
 */
static constexpr double DIRECTIONAL_FILTER_UPDATE_CUTOFF = 1.45;

} // namespace detail

namespace test {
template<typename Kernel_> class Test_exact_diamond_filter;
}

template<typename Kernel_> class Exact_diamond_filter {
  public:
    using Kernel = Kernel_;
    using FT = typename Kernel::FT;
    using Point_2 = typename Kernel::Point_2;
    using Iso_rectangle_2 = typename Kernel::Iso_rectangle_2;
    using IntervalNT = CGAL::Interval_nt_advanced;

    static constexpr bool has_exact_nt = CGAL::Algebraic_structure_traits<FT>::Is_exact::value;

    Exact_diamond_filter() : m_query{}, m_dead_sectors{} { m_dead_sectors.reserve(24); }

    /**
     * Reset the filter to be reused for
     * the next query point.
     */
    void reset(const Point_2 &query) {
        m_query = query;
        m_dead_sectors.clear();
        m_dead_sectors.emplace_back(-1.7976931348623157e+308, -detail::DIRECTIONAL_FILTER_PANGLE_CUTOFF);
        m_dead_sectors.emplace_back(detail::DIRECTIONAL_FILTER_PANGLE_CUTOFF, 1.7976931348623157e+308);
    }

    /**
     * Print the current dead sectors to the given output stream.
     */
    void print(std::ostream &output) const {
        for(std::size_t i = 0; i < m_dead_sectors.size(); ++i) {
            output << "Sector #" << i << ": [" << m_dead_sectors[i].pseudoangle_lb << ", "
                   << m_dead_sectors[i].pseudoangle_ub << "]\n    low = ";
            if(m_dead_sectors[i].low_point) {
                output << "(" << *m_dead_sectors[i].low_point << "), high = ";
            } else {
                output << "null, high = ";
            }
            if(m_dead_sectors[i].high_point) {
                output << "(" << *m_dead_sectors[i].high_point << ")\n";
            } else {
                output << "null\n";
            }
        }
    }

    /**
     * Print the current dead sectors to the output stream in JSON format.
     */
    void print_as_json(std::ostream &output) const {
        output << "{ \"sectors\": [";
        output << std::setprecision(19);
        bool first = true;
        for(const auto &sector : m_dead_sectors) {
            if(first) {
                first = false;
            } else {
                output << ", ";
            }
            output << "{\"lb\": " << sector.pseudoangle_lb << ", \"ub\": " << sector.pseudoangle_ub;
            if(sector.low_point) {
                output << ", \"low\": [" << sector.low_point->x() << ", " << sector.low_point->y() << "]";
            } else {
                output << ", \"low\": null";
            }
            if(sector.high_point) {
                output << ", \"high\": [" << sector.high_point->x() << ", " << sector.high_point->y() << "]";
            } else {
                output << ", \"high\": null";
            }
            output << "}";
        }
        output << "]}";
    }

    /**
     * Return the number of dead sectors.
     */
    std::size_t size() const noexcept { return m_dead_sectors.size(); }

    /**
     * Insert a new dead sector spanned by two points.
     */
    template<typename PointIterator> void insert_sector(PointIterator l, PointIterator r) {
        if constexpr(has_exact_nt) {
            auto [lb, ub, point_angle] = compute_sector(m_query.x(), m_query.y(), l->x(), l->y(), r->x(), r->y());
            if(ub <= lb)
                return;
            if(point_angle) {
                CGAL_precondition(!CGAL::left_turn(m_query, *l, *r));
                Sector_inserter inserter{this, lb, ub, &*r, &*l};
                inserter.insert_with_points();
            } else {
                Sector_inserter inserter{this, lb, ub, nullptr, nullptr};
                inserter.insert_without_points();
            }
        } else {
            IntervalNT qx(m_query.x()), qy(m_query.y());
            IntervalNT lx(l->x()), ly(l->y());
            IntervalNT rx(r->x()), ry(r->y());
            auto [lbi, ubi, point_angle] = compute_sector(qx, qy, lx, ly, rx, ry);
            if(CGAL::possibly(ubi <= lbi))
                return;
            if(point_angle) {
                Sector_inserter inserter{this, lbi.sup(), ubi.inf(), &*r, &*l};
                inserter.insert_with_points();
            } else {
                Sector_inserter inserter{this, lbi.sup(), ubi.inf(), nullptr, nullptr};
                inserter.insert_without_points();
            }
        }
    }

    /**
     * Insert a precomputed sector.
     */
    void insert_sector(FT pseudoangle_lb, FT pseudoangle_ub, const Point_2 *low, const Point_2 *high) {
        if(low) {
            Sector_inserter inserter{this, pseudoangle_lb, pseudoangle_ub, low, high};
            inserter.insert_with_points();
        } else {
            Sector_inserter inserter{this, pseudoangle_lb, pseudoangle_ub, nullptr, nullptr};
            inserter.insert_without_points();
        }
    }

    /**
     * Return true if all remaining points can safely be filtered out.
     */
    bool operator()() const noexcept { return m_dead_sectors.size() == 1; }

    /**
     * Check whether the given bounding box can safely be filtered out.
     */
    bool operator()(const Iso_rectangle_2 &rect) const { return p_check_rect(rect); }

    /**
     * Return true if the point can be safely filtered out.
     */
    bool operator()(const Point_2 &point) const noexcept {
        if constexpr(has_exact_nt) {
            FT angle = Compute_pseudo_angle<FT>{}(point.x() - m_query.x(), point.y() - m_query.y());
            auto it = m_dead_sectors.begin();
            while(it->pseudoangle_ub < angle)
                ++it;
            return it->pseudoangle_lb <= angle;
        } else {
            IntervalNT vx{point.x()};
            vx -= m_query.x();
            IntervalNT vy{point.y()};
            vy -= m_query.y();
            IntervalNT angle = Compute_pseudo_angle<IntervalNT>{}(vx, vy);
            auto it = m_dead_sectors.begin();
            double ainf = angle.inf(), asup = angle.sup();
            while(it->pseudoangle_ub < ainf)
                ++it;
            return (it->pseudoangle_lb <= ainf) & (it->pseudoangle_ub >= asup);
        }
    }

    /**
     * Return true if all points in the given pseudoangle interval
     * can safely be filtered out.
     */
    bool contains_pseudo_angle_interval(const FT &lb, const FT &ub) const {
        CGAL_precondition(lb <= ub);
        auto it = m_dead_sectors.begin();
        while(it->pseudoangle_ub < lb)
            ++it;
        return (it->pseudoangle_lb <= lb) & (it->pseudoangle_ub >= ub);
    }

    template<typename NT, std::enable_if_t<is_interval_number_type_v<NT>, int> = 0>
    bool contains_pseudo_angle_interval(NT iv) const {
        return contains_pseudo_angle_interval(iv.inf(), iv.sup());
    }

    template<typename NT, std::enable_if_t<!is_interval_number_type_v<NT>, int> = 0>
    bool contains_pseudo_angle_interval(NT x) const {
        return contains_pseudo_angle_interval(x, x);
    }

  private:
    template<typename K_> friend class test::Test_exact_diamond_filter;

    template<typename = void> FT p_angle_inexact_nt(FT x1, FT y1, bool ub) const {
        Compute_pseudo_angle<IntervalNT> compute_angle{};
        IntervalNT vx{x1}, vy{y1};
        vx -= m_query.x();
        vy -= m_query.y();
        IntervalNT a = compute_angle(vx, vy);
        return ub ? a.sup() : a.inf();
    }

    template<typename = void> FT p_angle_exact_nt(const FT &x1, const FT &y1, bool /*ub*/) const {
        Compute_pseudo_angle<FT> compute_angle{};
        return compute_angle(x1 - m_query.x(), y1 - m_query.y());
    }

    bool p_check_rect(const Iso_rectangle_2 &rect) const {
        auto angle = [&](const FT &x1, const FT &y1, bool ub) -> FT {
            if constexpr(has_exact_nt) {
                return p_angle_exact_nt(x1, y1, ub);
            } else {
                return p_angle_inexact_nt(x1, y1, ub);
            }
        };

        if(rect.xmax() <= m_query.x()) { // completely to the left
            if(rect.ymax() <= m_query.y()) {
                if(rect.ymax() == m_query.y()) {
                    // we would get 2 instead of -2 for a1 here
                    FT a2 = angle(rect.xmax(), rect.ymin(), true);
                    return a2 <= m_dead_sectors.front().pseudoangle_ub;
                }
                // completely below
                FT a1 = angle(rect.xmin(), rect.ymax(), false);
                FT a2 = angle(rect.xmax(), rect.ymin(), true);
                return contains_pseudo_angle_interval(a1, a2);
            } else if(rect.ymin() < m_query.y()) {
                // intersects the horizontal line, i.e., our
                // discontinuity, requiring special handling
                FT a1 = angle(rect.xmax(), rect.ymax(), false);
                FT a2 = angle(rect.xmax(), rect.ymin(), true);
                return (a1 >= m_dead_sectors.back().pseudoangle_lb) & (a2 <= m_dead_sectors.front().pseudoangle_ub);
            } else {
                // completely above
                FT a1 = angle(rect.xmax(), rect.ymax(), false);
                FT a2 = angle(rect.xmin(), rect.ymin(), true);
                return contains_pseudo_angle_interval(a1, a2);
            }
        }

        if(rect.xmin() < m_query.x()) { // intersects the vertical line
            // intersects the vertical line
            if(rect.ymax() < m_query.y()) {
                // completely below
                FT a1 = angle(rect.xmin(), rect.ymax(), false);
                FT a2 = angle(rect.xmax(), rect.ymax(), true);
                return contains_pseudo_angle_interval(a1, a2);
            } else if(rect.ymin() < m_query.y()) {
                // intersects the rectangle
                return false;
            } else {
                // completely above
                FT a1 = angle(rect.xmax(), rect.ymin(), false);
                FT a2 = angle(rect.xmin(), rect.ymin(), true);
                return contains_pseudo_angle_interval(a1, a2);
            }
        }

        // completely to the right
        if(rect.ymax() <= m_query.y()) {
            // completely below
            FT a1 = angle(rect.xmin(), rect.ymin(), false);
            FT a2 = angle(rect.xmax(), rect.ymax(), true);
            return contains_pseudo_angle_interval(a1, a2);
        } else if(rect.ymin() < m_query.y()) {
            // intersects the horizontal line
            FT a1 = angle(rect.xmin(), rect.ymin(), false);
            FT a2 = angle(rect.xmin(), rect.ymax(), true);
            return contains_pseudo_angle_interval(a1, a2);
        } else {
            // completely above
            FT a1 = angle(rect.xmax(), rect.ymin(), false);
            FT a2 = angle(rect.xmin(), rect.ymax(), true);
            return contains_pseudo_angle_interval(a1, a2);
        }
    }

    struct Dead_sector {
        /**
         * Construct a dead sector given two double values.
         * This is usually done only for the initial dead sectors.
         */
        Dead_sector(double lb, double ub)
            : pseudoangle_lb(lb), pseudoangle_ub(ub), high_point(nullptr), low_point(nullptr) {}

        /**
         * Construct a dead sector given two FT values indicating
         * a covered interval, and optionally two pointers to points.
         */
        Dead_sector(const FT &lb, const FT &ub, const Point_2 *low, const Point_2 *high)
            : pseudoangle_lb(lb), pseudoangle_ub(ub), high_point(high), low_point(low) {}

        // an interval of pseudo-angles;
        // all points within this interval
        // relative to the query point
        // are covered by this dead sector
        FT pseudoangle_lb, pseudoangle_ub;

        // the point with the highest angle
        // that we know is covered by this sector
        const Point_2 *high_point;

        // the point with the lowest angle
        // that we know is covered by this sector
        const Point_2 *low_point;
    };

    using DSContainer = std::vector<Dead_sector>;
    using DSIter = typename DSContainer::iterator;

    struct Sector_inserter {
        explicit Sector_inserter(Exact_diamond_filter *that, FT pseudoangle_lb, FT pseudoangle_ub, const Point_2 *low,
                                 const Point_2 *high)
            : that(that), new_lb(pseudoangle_lb), new_ub(pseudoangle_ub), new_low(low), new_high(high) {
            assert(!new_low || new_high);
            assert(new_low || !new_high);
        }

        void check_and_merge_upwards(DSIter from) {
            DSIter removed_begin = std::next(from), removed_end = removed_begin;
            const Point_2 *new_low = from->low_point;
            const Point_2 *new_high = from->high_point;
            FT new_ub = from->pseudoangle_ub;

            auto should_merge = [&](DSIter it) {
                if(it->pseudoangle_lb <= new_ub) {
                    if(!new_low)
                        new_low = it->low_point;
                    return true;
                }
                if(it->low_point) {
                    if(new_high && new_high == it->low_point) {
                        if(!new_low)
                            new_low = it->low_point;
                        return true;
                    }
                }
                return false;
            };

            while(removed_end != that->m_dead_sectors.end() && should_merge(removed_end)) {
                if(removed_end->pseudoangle_ub > new_ub) {
                    new_ub = removed_end->pseudoangle_ub;
                    if(removed_end->high_point)
                        new_high = removed_end->high_point;
                }
                ++removed_end;
            }

            if(removed_begin == removed_end)
                return;
            from->low_point = new_low;
            from->high_point = new_high;
            from->pseudoangle_ub = new_ub;
            that->m_dead_sectors.erase(removed_begin, removed_end);
        }

        void merge_into_first_no_points() {
            auto &merge_with = that->m_dead_sectors.front();
            if(merge_with.pseudoangle_ub >= new_ub)
                return;
            merge_with.pseudoangle_ub = new_ub;
            check_and_merge_upwards(that->m_dead_sectors.begin());
        }

        void merge_into_last_no_points() {
            auto it = std::prev(that->m_dead_sectors.end());
            if(it->pseudoangle_lb <= new_lb)
                return;
            auto beg = that->m_dead_sectors.begin();
            const Point_2 *new_low = nullptr;
            while(it != beg) {
                auto next = std::prev(it);
                if(next->pseudoangle_ub < new_lb)
                    break;
                it = next;
                if(next->low_point)
                    new_low = next->low_point;
            }
            it->pseudoangle_lb = (std::min)(new_lb, it->pseudoangle_lb);
            it->pseudoangle_ub = 3.0;
            it->high_point = nullptr;
            it->low_point = new_low;
            that->m_dead_sectors.erase(it + 1, that->m_dead_sectors.end());
        }

        void insert_without_points() {
            /*if(new_lb <= -1.0) {
                merge_into_first_no_points();
                return;
            }
            if(new_ub > 1.0) {
                merge_into_last_no_points();
                return;
            }*/

            first_ub_above = find_first_ub_above();
            // we always have first_ub_above->pseudoangle_ub >= new_lb;
            // the sector before first_ub_above does not exist or it has
            // pseudoangle_ub < new_lb, meaning that we cannot possibly
            // merge with it without a point
            if(first_ub_above->pseudoangle_ub > new_ub) {
                // cases 2, 4, 5: the upper bound of first_ub_above is above our upper bound
                if(first_ub_above->pseudoangle_lb > new_ub) {
                    // case 2: we are completely disjoint from first_ub_above
                    that->m_dead_sectors.emplace(first_ub_above, new_lb, new_ub, nullptr, nullptr);
                    // we can stop early here
                    return;
                }
                if(first_ub_above->pseudoangle_lb <= new_lb) {
                    // case 4: we are completely contained in first_ub_above;
                    // since we have no points, we need not do anything here
                    return;
                }
                // case 5: we overlap with first_ub_above on our right side
                first_ub_above->pseudoangle_lb = new_lb;
                return;
            }
            // cases 1, 3: the upper bound of first_ub_above is at most our upper bound
            if(first_ub_above->pseudoangle_lb < new_lb) {
                // case 1: we overlap with first_ub_above on our right side
                first_ub_above->pseudoangle_ub = new_ub;
                check_and_merge_upwards(first_ub_above);
                return;
            }
            // case 3: we completely contain first_ub_above
            first_ub_above->pseudoangle_lb = new_lb;
            first_ub_above->pseudoangle_ub = new_ub;
            check_and_merge_upwards(first_ub_above);
        }

        bool check_merge_with_previous(DSIter &first_ub_above) {
            if(first_ub_above == that->m_dead_sectors.begin())
                return false;
            DSIter prev = std::prev(first_ub_above);
            if(prev->high_point == new_low) {
                // we must merge with the previous sector
                prev->pseudoangle_ub = new_ub;
                prev->high_point = new_high;
                first_ub_above = prev;
                return true;
            }
            return false;
        }

        void insert_with_points() {
            first_ub_above = find_first_ub_above();
            if(check_merge_with_previous(first_ub_above)) {
                check_and_merge_upwards(first_ub_above);
                return;
            }

            if(first_ub_above->pseudoangle_ub > new_ub) {
                // cases 2, 4, 5: the upper bound of first_ub_above is above our upper bound
                if(first_ub_above->pseudoangle_lb > new_ub && first_ub_above->low_point != new_high) {
                    // case 2: we are completely disjoint from first_ub_above
                    that->m_dead_sectors.emplace(first_ub_above, new_lb, new_ub, new_low, new_high);
                    // we can stop early here
                    return;
                } else if(first_ub_above->low_point == new_high) {
                    // we overlap with first_ub_above
                    if(new_lb < first_ub_above->pseudoangle_lb) {
                        first_ub_above->pseudoangle_lb = new_lb;
                        first_ub_above->low_point = new_low;
                    }
                    // no upward merging necessary
                    return;
                }
                if(first_ub_above->pseudoangle_lb <= new_lb) {
                    // case 4: we are completely contained in first_ub_above
                    if(!first_ub_above->low_point)
                        first_ub_above->low_point = new_low;
                    if(!first_ub_above->high_point)
                        first_ub_above->high_point = new_high;
                    return;
                }
                // case 5: we overlap with first_ub_above on our right side
                first_ub_above->pseudoangle_lb = new_lb;
                first_ub_above->low_point = new_low;
                if(!first_ub_above->high_point)
                    first_ub_above->high_point = new_high;
                return;
            } else {
                // cases 1, 3: the upper bound of first_ub_above is at most our upper bound
                if(first_ub_above->pseudoangle_lb < new_lb) {
                    // case 1: we overlap with first_ub_above on our left side
                    first_ub_above->pseudoangle_ub = new_ub;
                    first_ub_above->high_point = new_high;
                    if(!first_ub_above->low_point)
                        first_ub_above->low_point = new_low;
                    check_and_merge_upwards(first_ub_above);
                    return;
                }
                // case 3: we completely contain first_ub_above
                first_ub_above->pseudoangle_lb = new_lb;
                first_ub_above->pseudoangle_ub = new_ub;
                first_ub_above->low_point = new_low;
                first_ub_above->high_point = new_high;
                check_and_merge_upwards(first_ub_above);
                return;
            }
        }

        DSIter find_first_ub_above() {
            // extend the interval if it reaches into the initial pseudo-sectors
            if(new_lb > detail::DIRECTIONAL_FILTER_PANGLE_CUTOFF)
                new_lb = detail::DIRECTIONAL_FILTER_PANGLE_CUTOFF;
            if(new_ub < -detail::DIRECTIONAL_FILTER_PANGLE_CUTOFF)
                new_ub = -detail::DIRECTIONAL_FILTER_PANGLE_CUTOFF;
            // this ensures that we will always return a valid iterator here
            auto it = that->m_dead_sectors.begin();
            while(it->pseudoangle_ub < new_lb)
                ++it;
            // now it->pseudoangle_ub >= new_lb
            return it;
        }

        Exact_diamond_filter *that;
        FT new_lb, new_ub;
        const Point_2 *new_low;
        const Point_2 *new_high;
        DSIter first_ub_above;
    };

    Point_2 m_query;
    DSContainer m_dead_sectors;
};

template<typename Traits>
inline std::ostream &operator<<(std::ostream &output, const Exact_diamond_filter<Traits> &filter) {
    filter.print(output);
    return output;
}

} // namespace mwt

#endif
