#ifndef CGAL_MWT_EXACT_DIAMOND_FILTERED_SEARCH_DRIVER_H_INCLUDED_
#define CGAL_MWT_EXACT_DIAMOND_FILTERED_SEARCH_DRIVER_H_INCLUDED_

#include "Directional_filter.h"

namespace mwt {

namespace test {
class Test_exact_diamond_filtered_search_driver;
}

template<typename Traits_, typename PointIterator_, typename ExpansionCallback = void>
class Exact_diamond_filtered_search_driver {
  public:
    friend class test::Test_exact_diamond_filtered_search_driver;

    using Traits = Traits_;
    using Kernel = typename Traits::Kernel;
    using PointIterator = PointIterator_;
    using FT = typename Kernel::FT;
    using Tree = Point_quadtree<Kernel, PointIterator>;
    using Filter = Exact_diamond_filter<Kernel>;
    using Point = typename Traits::Point_2;
    using Interval = CGAL::Interval_nt_advanced;
    static constexpr bool has_exact_nt = CGAL::Algebraic_structure_traits<FT>::Is_exact::value;

    /**
     * Driver for exact incremental filtered search.
     * Constructor version without expansion callback.
     */
    explicit Exact_diamond_filtered_search_driver(const Tree *tree)
        : m_tree(tree), m_query{0, 0}, m_points{PointEntryLess{m_query}}, m_search{tree, m_query} {}

    /**
     * Driver for exact incremental filtered search.
     * Constructor version with expansion callback.
     */
    template<typename CallbackType>
    explicit Exact_diamond_filtered_search_driver(const Tree *tree, CallbackType &&exp_cb)
        : m_tree(tree), m_query{0, 0}, m_points{PointEntryLess{m_query}},
          m_search{tree, m_query, std::forward<CallbackType>(exp_cb)} {}

    /**
     * Call the given callable with all neighbors
     * that we cannot filter out by the diamond property.
     * We only report points that are lexicographically
     * greater than query, assuming that all points are
     * used as query point at some point.
     */
    template<typename Callable /*(PointIerator)*/>
    void enumerate_filtered_neighbors_of(const Point &query, Callable &&callback) {
        //        auto before = std::chrono::steady_clock::now();
        p_reset(query);
        CGAL_assertion(*m_search.current().handle() == query);

        auto &filter = m_search.filter();
        auto check_filter = [&](const Point &p, const auto &angle) {
            if(!typename Kernel::Less_xy_2{}(query, p)) {
                return true;
            }
            return filter.contains_pseudo_angle_interval(angle);
        };

        while(m_search.next()) {
            p_apply_waiting_sectors();

            if(filter()) {
                // early abort: 360 degree
                // angle around query is covered
                break;
            }

            const auto cur_handle = m_search.current().handle();
            PointEntry entry{cur_handle, query};
            if(CGAL::certainly(CGAL::abs(entry.angle) >= detail::DIRECTIONAL_FILTER_UPDATE_CUTOFF)) {
                // we can be certain that this point fails the
                // Less_xy_2 test, and furthermore it is very unlikely
                // that it does contribute to the filter
                continue;
            }

            auto [it, did_insert] = m_points.insert(entry);
            if(!did_insert) {
                // the point is collinear to an earlier point and the query,
                // so we can safely ignore it.
                continue;
            }

            const auto &point = *cur_handle;
            if(!check_filter(point, entry.angle)) {
                // the filter cannot immediately filter out this point
                if(!p_filtered_by_extended_diamond_check(it, point)) {
                    std::invoke(callback, cur_handle);
                }
            }

            p_construct_dead_sectors(it, point);
        }
    }

  private:
    template<typename PointSetIter>
    bool p_filtered_by_extended_diamond_check(PointSetIter inserted, const Point &point) {
        Traits traits{};
        auto diamond_test = traits.diamond_test_2_object();

        {
            auto i = std::next(inserted);
            auto e = m_points.end();
            for(;;) {
                if(i == e || inf_or_value(i->angle) - sup_or_value(inserted->angle) >= 0.525) {
                    return false;
                }
                if(diamond_test(m_query, point, *i->point))
                    break;
                ++i;
            }
        }

        {
            auto i = inserted;
            auto e = m_points.begin();
            for(;;) {
                if(i == e)
                    return false;
                --i;
                if(inf_or_value(inserted->angle) - sup_or_value(i->angle) >= 0.525) {
                    return false;
                }
                if(diamond_test(point, m_query, *i->point))
                    return true;
            }
        }
    }

    template<typename PointSetIter> void p_construct_dead_sectors(PointSetIter inserted, const Point &point) {
        auto handle_dead_sector = [&](const std::tuple<FT, FT, bool> &sector, const Point *low, const Point *high) {
            const auto &lb = std::get<0>(sector);
            const auto &ub = std::get<1>(sector);
            if(ub <= lb)
                return;
            if(m_search.filter().contains_pseudo_angle_interval(lb, ub))
                return;
            if(!std::get<2>(sector)) {
                low = high = nullptr;
            }
            FT activation_distance =
                sup_or_value(m_search.current().squared_distance_interval_or_exact() * detail::SQUARED_2COS_DPA);
            m_waiting_sectors.emplace_back(sup_or_value(lb), inf_or_value(ub), low, high, activation_distance);
        };

        auto next = std::next(inserted);
        if(next != m_points.end()) {
            auto sec = p_compute_sector(*next->point, point);
            handle_dead_sector(sec, &*inserted->point, &*next->point);
        }
        if(inserted != m_points.begin()) {
            auto prev = std::prev(inserted);
            auto sec = p_compute_sector(point, *prev->point);
            handle_dead_sector(sec, &*prev->point, &*inserted->point);
        }
    }

    std::tuple<FT, FT, bool> p_compute_sector(const Point &left, const Point &right) {
        if constexpr(has_exact_nt) {
            return compute_sector(m_query.x(), m_query.y(), left.x(), left.y(), right.x(), right.y());
        } else {
            Interval qx{m_query.x()}, qy{m_query.y()};
            Interval lx{left.x()}, ly{left.y()};
            Interval rx{right.x()}, ry{right.y()};
            auto sec = compute_sector(qx, qy, lx, ly, rx, ry);
            return {std::get<0>(sec).sup(), std::get<1>(sec).inf(), std::get<2>(sec)};
        }
    }

    void p_apply_waiting_sectors() {
        const auto &cpoint = m_search.current();
        while(!m_waiting_sectors.empty() &&
              cpoint.squared_distance_certainly_at_least(m_waiting_sectors.front().squared_activation_distance)) {
            const auto &sector = m_waiting_sectors.front();
            m_search.filter().insert_sector(sector.pseudoangle_lb, sector.pseudoangle_ub, sector.low_point,
                                            sector.high_point);
            m_waiting_sectors.pop_front();
        }
    }

    /**
     * Reset the search to a new query point.
     */
    void p_reset(const Point &query) {
        m_query = query;
        auto found = m_search.start_search(query);
        CGAL_assertion(found);
        static_cast<void>(found);
        m_points = PointSet{PointEntryLess{query}};
        m_points.reserve(128);
        m_waiting_sectors.clear();
    }

    /**
     * AngleType is either FT (for exact number types) or CGAL::Interval_nt_advanced.
     */
    using AngleType = std::conditional_t<has_exact_nt, FT, CGAL::Interval_nt_advanced>;

    /**
     * A point entry in the PointSet of points
     * angularly sorted around query point.
     */
    struct PointEntry {
        PointIterator point;
        AngleType angle;

        static AngleType compute_angle(const Point &query, const Point &point) {
            if constexpr(has_exact_nt) {
                return Compute_pseudo_angle<FT>{}(point.x() - query.x(), point.y() - query.y());
            } else {
                Interval qx{point.x()}, qy{point.y()};
                qx -= query.x();
                qy -= query.y();
                return Compute_pseudo_angle<Interval>{}(qx, qy);
            }
        }

        PointEntry(PointIterator point, const Point &query) : point(point), angle(compute_angle(query, *point)) {}
    };

    /**
     * Comparator for point entries.
     */
    struct PointEntryLess {
        /**
         * Comparison by stored pseudoangle interval first
         * and by quadrant/orientation predicate second.
         */
        bool operator()(const PointEntry &e1, const PointEntry &e2) const {
            if(CGAL::certainly(e1.angle < e2.angle))
                return true;
            if(CGAL::certainly(e1.angle >= e2.angle))
                return false;
            if(e1.point == e2.point)
                return false;
            bool e1_below = (e1.point->y() <= query.y());
            bool e2_below = (e2.point->y() <= query.y());
            if(e1_below != e2_below)
                return e1_below;
            CGAL::Orientation o = CGAL::orientation(query, *e1.point, *e2.point);
            if(o == CGAL::LEFT_TURN)
                return true;
            return false;
        }

        Point query;
    };

    using PointSet = boost::container::flat_set<PointEntry, PointEntryLess>;
    using Search = typename Tree::template Incremental_search<Filter, Traits, ExpansionCallback>;

    struct SectorWithDistance {
        FT pseudoangle_lb, pseudoangle_ub;
        const Point *low_point, *high_point;
        FT squared_activation_distance;

        SectorWithDistance(FT pseudoangle_lb, FT pseudoangle_ub, const Point *low_point, const Point *high_point,
                           FT squared_activation_distance)
            : pseudoangle_lb(std::move(pseudoangle_lb)), pseudoangle_ub(std::move(pseudoangle_ub)),
              low_point(low_point), high_point(high_point),
              squared_activation_distance(std::move(squared_activation_distance)) {}
    };

    const Tree *m_tree;
    Point m_query;
    PointSet m_points;
    Search m_search;
    std::deque<SectorWithDistance> m_waiting_sectors;
};

} // namespace mwt

#endif
