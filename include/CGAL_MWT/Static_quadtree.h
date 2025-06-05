#ifndef CGAL_MWT_STATIC_QUADTREE_H_INCLUDED_
#define CGAL_MWT_STATIC_QUADTREE_H_INCLUDED_

#include <algorithm>
#include <array>
#include <boost/heap/d_ary_heap.hpp>
#include <iterator>
#include <memory>
#include <queue>
#include <vector>

#include "Is_interval.h"
#include "Quadtree_aux.h"
#include "Search_aux.h"
#include <CGAL/Bbox_2.h>
#include <CGAL/Exact_rational.h>
#include <CGAL/Interval_nt.h>
#include <CGAL/Iso_rectangle_2.h>
#include <CGAL/enum.h>
#include <CGAL/squared_distance_2.h>

namespace mwt {

/**
 * Node order in the quad tree:
 * _____________
 * |     |     |
 * |  1  |  2  |
 * |_____|_____|
 * |     |     |
 * |  0  |  3  |
 * |_____|_____|
 */
template<typename Traits_, typename Iterator_> class Quadtree_node {
  public:
    using Traits = Traits_;
    using Iterator = Iterator_;

    using FT = typename Traits::FT;
    using Point_2 = typename Traits::Point_2;
    using Iso_rectangle_2 = typename Traits::Iso_rectangle_2;
    using Bbox_2 = CGAL::Bbox_2;

    const Iso_rectangle_2 bbox;
    Quadtree_node *children;
    const Iterator begin;
    const Iterator end;

    /**
     * Create a node from a bounding box rectangle and a range of points.
     */
    Quadtree_node(const Iso_rectangle_2 &box, Iterator begin, Iterator end)
        : bbox{box}, children{nullptr}, begin{begin}, end{end} {
        CGAL_expensive_assertion(
            std::none_of(begin, end, [&](const Point_2 &p) { return bbox.has_on_unbounded_side(p); }));
    }

    /**
     * Create a node from a Bbox_2 and a range of points.
     */
    Quadtree_node(const Bbox_2 &box, Iterator begin, Iterator end) : Quadtree_node(Iso_rectangle_2(box), begin, end) {}

    /**
     * Create a node from min and max coordinates and a range of points.
     */
    Quadtree_node(const FT &xmin, const FT &ymin, const FT &xmax, const FT &ymax, Iterator begin, Iterator end)
        : Quadtree_node(Iso_rectangle_2{xmin, ymin, xmax, ymax}, begin, end) {}

    /**
     * Check whether the node has any children.
     */
    bool is_leaf() const { return children == nullptr; }

    /**
     * Check whether the node is empty, i.e., has any points.
     */
    bool is_empty() const { return begin == end; }

    /**
     * Get the bounding box of this node.
     */
    const Iso_rectangle_2 &bounding_box() const { return bbox; }

    /**
     * Get the begin of the point range.
     */
    Iterator points_begin() const { return begin; }

    /**
     * Get the end of the point range.
     */
    Iterator points_end() const { return end; }

    /**
     * Get the begin of the node's children.
     */
    const Quadtree_node *children_begin() const { return children; }

    /**
     * Get the end of the node's children.
     */
    const Quadtree_node *children_end() const { return children + 4; }

    /**
     * Compute the squared distance from a
     * point to the bounding box of this node.
     */
    FT squared_distance(const Point_2 &p) const {
        FT distance{0};
        if(p.x() < bbox.xmin()) {
            FT dx = bbox.xmin() - p.x();
            distance += dx * dx;
        } else if(p.x() > bbox.xmax()) {
            FT dx = p.x() - bbox.xmax();
            distance += dx * dx;
        }

        if(p.y() < bbox.ymin()) {
            FT dy = bbox.ymin() - p.y();
            distance += dy * dy;
        } else if(p.y() > bbox.ymax()) {
            FT dy = p.y() - bbox.ymax();
            distance += dy * dy;
        }
        return distance;
    }

    /**
     * Compute the squared distance from a point to our children.
     */
    std::array<FT, 4> squared_distance_to_children(const Point_2 &p) const {
        return squared_distance_to_children(p, false);
    };

    /**
     * Compute the square distance from a point to our children,
     * with a flag specifying whether the point is inside our bounding box.
     */
    std::array<FT, 4> squared_distance_to_children(const Point_2 &p, bool p_inside_bbox) const {
        CGAL_precondition(children != nullptr);

        if(p_inside_bbox) {
            return squared_distance_to_children_inside_bbox(p);
        } else {
            std::array<FT, 4> distance{0, 0, 0, 0};

            distance[0] = children[0].squared_distance(p);
            distance[1] = children[1].squared_distance(p);
            distance[2] = children[2].squared_distance(p);
            distance[3] = children[3].squared_distance(p);

            return distance;
        }
    }

  private:
    std::array<FT, 4> squared_distance_to_children_inside_bbox(const Point_2 &p) const {
        CGAL_precondition(!bbox.has_on_unbounded_side(p));

        std::array<FT, 4> distance{0, 0, 0, 0};
        const FT mx = FT{0.5} * (bbox.xmin() + bbox.xmax());
        const FT my = FT{0.5} * (bbox.ymin() + bbox.ymax());
        const FT dx = (p.x() - mx) * (p.x() - mx);
        const FT dy = (p.y() - my) * (p.y() - my);
        const FT dxdy = dx + dy;

        if(p.x() < mx) {
            if(p.y() < my) {
                distance[1] = dy;
                distance[2] = dxdy;
                distance[3] = dx;
            } else {
                distance[0] = dy;
                distance[2] = dx;
                distance[3] = dxdy;
            }
        } else {
            if(p.y() < my) {
                distance[0] = dx;
                distance[1] = dxdy;
                distance[2] = dy;
            } else {
                distance[0] = dxdy;
                distance[1] = dx;
                distance[3] = dy;
            }
        }

        CGAL_expensive_postcondition(distance[0] == children[0].squared_distance(p));
        CGAL_expensive_postcondition(distance[1] == children[1].squared_distance(p));
        CGAL_expensive_postcondition(distance[2] == children[2].squared_distance(p));
        CGAL_expensive_postcondition(distance[3] == children[3].squared_distance(p));
        return distance;
    }
};

/**
 * Only split a quadtree node if it has more than SPLIT_THRESHOLD points.
 */
static constexpr int SPLIT_THRESHOLD = 16;

/**
 * Make an effort to use a more balanced split if
 * more than BALANCE_THRESHOLD * num_points points
 * are on one side of the split.
 */
static constexpr double BALANCE_THRESHOLD = 0.95;

struct EmptyBase {
    template<typename Node> void report_internal_expansion(const Node *, int) const noexcept {}
    template<typename Node> void report_leaf_expansion(const Node *) const noexcept {}
};

template<typename Callback> struct ExpansionCallbackContainer {
    Callback callback;

    ExpansionCallbackContainer(Callback callback) : callback{std::move(callback)} {}

    template<typename Node> void report_internal_expansion(const Node *node, int num_children) const noexcept {
        callback.report_internal_expansion(node, num_children);
    }

    template<typename Node> void report_leaf_expansion(const Node *node) const noexcept {
        callback.report_leaf_expansion(node);
    }

    template<typename Callback_> void set_expansion_callback(Callback_ &&callback) {
        this->callback = std::forward<Callback_>(callback);
    }

    Callback &get_expansion_callback() noexcept { return callback; }

    const Callback &get_expansion_callback() const noexcept { return callback; }
};

template<typename Kernel_, typename PointIterator_, typename NodeAllocator = void> class Point_quadtree {
  public:
    using Kernel = Kernel_;
    using Iterator = PointIterator_;

    static_assert(std::is_convertible_v<typename std::iterator_traits<Iterator>::iterator_category,
                                        std::random_access_iterator_tag>,
                  "Point_quadtree: PointIterator must be a random access iterator!");

    using FT = typename Kernel::FT;
    using Point_2 = typename Kernel::Point_2;
    using Iso_rectangle_2 = typename Kernel::Iso_rectangle_2;
    using Segment_2 = typename Kernel::Segment_2;

    class Node;
    class Node {
        const Iterator m_points_begin;
        const Iterator m_points_end;
        Node *m_children;
        FT xmin, xmax, ymin, ymax;

        using FTCRef = typename std::conditional_t<!std::is_floating_point_v<FT>, const FT &, FT>;

        friend Point_quadtree;

        template<typename NodeAlloc> void destroy(NodeAlloc &alloc) {
            using Alloc_traits = std::allocator_traits<NodeAlloc>;
            if(m_children != nullptr) {
                for(int i = 0; i < 4; ++i) {
                    m_children[i].destroy(alloc);
                    Alloc_traits::destroy(alloc, m_children + i);
                }
                Alloc_traits::deallocate(alloc, m_children, 4);
                m_children = nullptr;
            }
        }

      public:
        bool empty() const noexcept { return m_points_begin == m_points_end; }
        bool leaf() const noexcept { return m_children == nullptr; }
        Iterator points_begin() const noexcept { return m_points_begin; }
        Iterator points_end() const noexcept { return m_points_end; }
        Node *children_begin() noexcept { return m_children; }
        Node *children_end() noexcept { return m_children + 4; }
        const Node *children_begin() const noexcept { return m_children; }
        const Node *children_end() const noexcept { return m_children + 4; }
        Iso_rectangle_2 bbox() const noexcept { return Iso_rectangle_2{xmin, ymin, xmax, ymax}; }
        FTCRef get_xmin() const noexcept { return xmin; }
        FTCRef get_xmax() const noexcept { return xmax; }
        FTCRef get_ymin() const noexcept { return ymin; }
        FTCRef get_ymax() const noexcept { return ymax; }
        std::size_t size() const noexcept {
            return static_cast<std::size_t>(std::distance(points_begin(), points_end()));
        }

        /**
         * Check if the bounding box contains the given point.
         */
        bool contains(const Point_2 &p) const noexcept {
            return xmin <= p.x() && xmax >= p.x() && ymin <= p.y() && ymax >= p.y();
        }

        Node(Iterator points_begin, Iterator points_end, FT xmin, FT ymin, FT xmax, FT ymax)
            : m_points_begin{points_begin}, m_points_end{points_end}, m_children{nullptr}, xmin{xmin}, xmax{xmax},
              ymin{ymin}, ymax{ymax} {}

        /**
         * Compute the squared distance from a
         * point to the bounding box of this node.
         */
        FT squared_distance(const Point_2 &p) const {
            FT distance{0};
            const FT px = p.x(), py = p.y();
            if(px < xmin) {
                FT dx = xmin - px;
                distance += dx * dx;
            } else if(px > xmax) {
                FT dx = px - xmax;
                distance += dx * dx;
            }
            if(py < ymin) {
                FT dy = ymin - py;
                distance += dy * dy;
            } else if(py > ymax) {
                FT dy = py - ymax;
                distance += dy * dy;
            }
            return distance;
        }

        /**
         * Expand this non-leaf node for a query point inside its bounding box.
         * Children with zero distance (i.e., containing the query) are inserted
         * into zero_node_buffer by calling push_back, other nodes are inserted
         * into node_queue by calling emplace.
         */
        template<typename ZeroNodeBuffer, typename NodeQueue>
        void expand_children_query_contained(const Point_2 &query, ZeroNodeBuffer &zero_node_buffer,
                                             NodeQueue &node_queue) const {
            CGAL_precondition(contains(query));
            CGAL_precondition(!leaf());
            const Node *c0 = m_children;
            const FT &mx = c0->get_xmax();
            const FT &my = c0->get_ymax();

            // assert that we have the right corner
            // picked for splitting coordinates
            CGAL_assertion(mx != get_xmin() || c0->get_xmin() == c0->get_xmax());
            CGAL_assertion(mx != get_xmax() || m_children[3].get_xmin() == m_children[3].get_xmax());
            CGAL_assertion(my != get_ymin() || c0->get_ymin() == c0->get_ymax());
            CGAL_assertion(my != get_ymax() || m_children[1].get_ymin() == m_children[1].get_ymax());

            auto queue_emplace = [&](const Node *child, bool on_border) {
                if(child->empty())
                    return;
                if(on_border)
                    zero_node_buffer.push_back(child);
                else
                    node_queue.emplace(query, child);
            };

            if(query.x() < mx) {
                if(query.y() < my) {
                    if(!c0->empty())
                        zero_node_buffer.push_back(c0);
                    queue_emplace(c0 + 1, false);
                    queue_emplace(c0 + 2, false);
                    queue_emplace(c0 + 3, false);
                } else {
                    if(!c0[1].empty())
                        zero_node_buffer.push_back(c0 + 1);
                    queue_emplace(c0, query.y() == my);
                    queue_emplace(c0 + 2, false);
                    queue_emplace(c0 + 3, false);
                }
            } else {
                if(query.y() < my) {
                    if(!c0[3].empty())
                        zero_node_buffer.push_back(c0 + 3);
                    queue_emplace(c0, query.x() == mx);
                    queue_emplace(c0 + 1, false);
                    queue_emplace(c0 + 2, false);
                } else {
                    if(!c0[2].empty())
                        zero_node_buffer.push_back(c0 + 2);
                    queue_emplace(c0, (query.x() == mx) & (query.y() == my));
                    queue_emplace(c0 + 1, query.x() == mx);
                    queue_emplace(c0 + 3, query.y() == my);
                }
            }
        }

        /**
         * Verify the given node and its invariants.
         */
        void verify(int level = 0) const {
            auto oob_point = [&](const Point_2 &p) { return bbox().has_on_unbounded_side(p); };
            if(std::any_of(m_points_begin, m_points_end, oob_point)) {
                throw std::logic_error("Incorrect tree node: Point outside bounding box!");
            }
            if(leaf()) {
                if(std::distance(m_points_begin, m_points_end) >= SPLIT_THRESHOLD) {
                    throw std::logic_error("Incorrect tree node: Leaf node with too many points!");
                }
            } else {
                for(int i{0}; i < 4; ++i) {
                    m_children[i].verify(level + 1);
                }
            }
        }
    };

    void print_node_to_json(std::ostream &out, bool &first, const Node *node, int level) const {
        bool is_leaf = node->leaf();
        if(!first) {
            out << ",\n\t";
        } else {
            first = false;
        }
        out << "{ \"xmin\": " << int(node->get_xmin()) << ", \"ymin\": " << int(node->get_ymin())
            << ", \"xmax\": " << int(node->get_xmax()) << ", \"ymax\": " << int(node->get_ymax())
            << ", \"leaf\": " << (is_leaf ? "true" : "false") << ", \"level\": " << level << " }";
        if(!is_leaf) {
            for(int i = 0; i < 4; ++i) {
                print_node_to_json(out, first, node->children_begin() + i, level + 1);
            }
        }
    }

    void print_nodelist_to_json(std::ostream &out) const {
        if(!m_root) {
            out << "{\"nodes\": []}";
        } else {
            out << "{ \"nodes\": [\n\t";
            bool first = true;
            print_node_to_json(out, first, m_root, 0);
            out << "\t\n]\n}";
        }
    }

    Iso_rectangle_2 bbox() const noexcept {
        if(m_root) {
            return m_root->bbox();
        } else {
            return Iso_rectangle_2{0, 0, 0, 0};
        }
    }

    /**
     * Implementation of an incremental search with a filtering predicate.
     * The incremental search starts from a given query point and lists
     * all points of the point set stored in the tree in order
     * of ascending distance to the query point; points for which the filtering
     * predicate cannot return false may be skipped in the enumeration.
     * For more guarantees on filtering, call the filtering predicate on each
     * reported point and skip those that are filtered out.
     * Care is taken that the search finds points in exactly the right order
     * even if the number type of the given Traits/Kernel is not exact, at
     * least when the input points have double coordinates.
     *
     * The filtering predicate has to support the following interface:
     *  void reset(const Point_2& query)
     *      - reset the filter's query point to another point
     *      - needed to reuse the incremental search structure
     *        for more than one query
     *  bool operator()(const Point_2& p) const;
     *      - can be called with any point in the point set
     *      - shall return true if the point should not be
     *        returned by the incremental search
     *      - if exactness in either direction (i.e., w.r.t. false negative
     *        and false positive points) is required,
     *        the filter predicate should be exact in that direction.
     *  bool operator()(const Iso_rectangle_2& box) const;
     *      - can be called with any bounding box of a node
     *      - shall return true if no point in the bounding box
     *        should be returned by the incremental search
     *      - if the search must not suffer from falsely filtered
     *        out points, the predicate must not return true for
     *        a box that might contain unfiltered points.
     *  bool operator()() const;
     *      - shall return true if all points of interest have been
     *        found and the search can be aborted without missing
     *        any potentially unfiltered points.
     *      - a conservative implementation can just return false.
     */
    template<class Filter_, typename Traits_ = Kernel, typename ExpansionCallback = void>
    class Incremental_search : public std::conditional_t<!std::is_void_v<ExpansionCallback>,
                                                         ExpansionCallbackContainer<ExpansionCallback>, EmptyBase> {
      public:
        using Super = std::conditional_t<!std::is_void_v<ExpansionCallback>,
                                         ExpansionCallbackContainer<ExpansionCallback>, EmptyBase>;
        using Filter = Filter_;
        using Traits = Traits_;
        using Node_handle = const Node *;
        using Point_handle = Iterator;
        using FT = typename Traits::FT;

        static constexpr bool has_exact_nt = CGAL::Algebraic_structure_traits<FT>::Is_exact::value;
        static constexpr bool has_expansion_callback = !std::is_void_v<ExpansionCallback>;

      private:
        template<typename AvoidNeedlessInstantiation_> class Node_with_distance_FT_nonexact;

        /**
         * We have two datastructures for storing a point
         * with a distance. One depends on an exact number
         * type as FT, the other does not and uses interval
         * arithmetic and only computes the exact distance
         * if necessary.
         */
        template<typename AvoidNeedlessInstantiation_> class Point_with_distance_FT_nonexact {
          public:
            using ExactNT = CGAL::Exact_rational;

            Point_with_distance_FT_nonexact(const Point_2 &query, Point_handle distance_to)
                : m_distance_interval(CGAL::squared_distanceC2(
                      CGAL::Interval_nt_advanced(query.x()), CGAL::Interval_nt_advanced(query.y()),
                      CGAL::Interval_nt_advanced(distance_to->x()), CGAL::Interval_nt_advanced(distance_to->y()))),
                  m_distance_to(distance_to) {}

            Point_handle handle() const noexcept { return m_distance_to; }

            /**
             * Cheaply check if the squared distance from the query
             * to *this is certainly at least the given squared distance.
             */
            bool squared_distance_certainly_at_least(const FT &sqdist) const {
                return sqdist <= m_distance_interval.inf();
            }

            CGAL::Interval_nt_advanced squared_distance_interval_or_exact() const noexcept {
                return m_distance_interval;
            }

            class CompareGreater {
              public:
                explicit CompareGreater(const Point_2 &query) : query(query) {}

                CompareGreater(const CompareGreater &other) noexcept : query(other.query) {}

                CompareGreater &operator=(const CompareGreater &other) noexcept {
                    query = other.query;
                    qx.reset();
                    qy.reset();
                    return *this;
                }

                CompareGreater(CompareGreater &&) = default;
                CompareGreater &operator=(CompareGreater &&) = default;

                /**
                 * Compare two distances. Compute the exact distance
                 * only if necessary.
                 */
                bool operator()(const Point_with_distance_FT_nonexact &d1,
                                const Point_with_distance_FT_nonexact &d2) const {
                    const auto &d1i = d1.m_distance_interval;
                    const auto &d2i = d2.m_distance_interval;
                    if(d1i.inf() > d2i.sup())
                        return true;
                    if(d1i.sup() <= d2i.inf())
                        return false;
                    if(d1.m_distance_to == d2.m_distance_to)
                        return false;
                    return d1.get_exact_distance(get_qx(), get_qy()) > d2.get_exact_distance(get_qx(), get_qy());
                }

              private:
                const ExactNT &get_qx() const {
                    if(!qx) {
                        qx = std::make_unique<ExactNT>(query.x());
                    }
                    return *qx;
                }

                const ExactNT &get_qy() const {
                    if(!qy) {
                        qy = std::make_unique<ExactNT>(query.y());
                    }
                    return *qy;
                }

                Point_2 query;
                mutable std::unique_ptr<ExactNT> qx;
                mutable std::unique_ptr<ExactNT> qy;
            };

          private:
            const ExactNT &get_exact_distance(const ExactNT &qx, const ExactNT &qy) const {
                if(!m_exact_distance) {
                    m_exact_distance = std::make_unique<ExactNT>(
                        CGAL::squared_distanceC2(qx, qy, ExactNT(m_distance_to->x()), ExactNT(m_distance_to->y())));
                }
                return *m_exact_distance;
            }

            friend class Node_with_distance_FT_nonexact<AvoidNeedlessInstantiation_>;

            CGAL::Interval_nt_advanced m_distance_interval;
            Point_handle m_distance_to;
            mutable std::unique_ptr<ExactNT> m_exact_distance;
        };

        /**
         * Store a tree node and its distance.
         */
        template<typename AvoidNeedlessInstantiation_> class Node_with_distance_FT_nonexact {
          public:
            Node_with_distance_FT_nonexact(const Point_2 &query, const Node *node)
                : m_node(node), m_distance_to(compute_closest(query)),
                  m_distance_interval(CGAL::squared_distanceC2(
                      CGAL::Interval_nt_advanced(query.x()), CGAL::Interval_nt_advanced(query.y()),
                      CGAL::Interval_nt_advanced(m_distance_to.x()), CGAL::Interval_nt_advanced(m_distance_to.y()))) {}

            const Node *node() const noexcept { return m_node; }

            class CompareGreater {
              public:
                explicit CompareGreater(const Point_2 &query) : query(query) {}

                /**
                 * Compare two distances. Compute the exact distance
                 * only if necessary.
                 */
                bool operator()(const Node_with_distance_FT_nonexact &d1,
                                const Node_with_distance_FT_nonexact &d2) const {
                    const auto &d1i = d1.m_distance_interval;
                    const auto &d2i = d2.m_distance_interval;
                    if(d1i.inf() > d2i.sup())
                        return true;
                    if(d1i.sup() <= d2i.inf())
                        return false;
                    if(d1.m_distance_to == d2.m_distance_to)
                        return false;
                    return d1.get_exact_distance(get_qx(), get_qy()) > d2.get_exact_distance(get_qx(), get_qy());
                }

                template<typename X_>
                bool operator()(const Node_with_distance_FT_nonexact &d1,
                                const Point_with_distance_FT_nonexact<X_> &d2) const {
                    const auto &d1i = d1.m_distance_interval;
                    const auto &d2i = d2.m_distance_interval;
                    if(d1i.inf() > d2i.sup())
                        return true;
                    if(d1i.sup() <= d2i.inf())
                        return false;
                    return d1.get_exact_distance(get_qx(), get_qy()) > d2.get_exact_distance(get_qx(), get_qy());
                }

                CompareGreater(const CompareGreater &other) noexcept : query(other.query) {}

                CompareGreater &operator=(const CompareGreater &other) noexcept {
                    query = other.query;
                    qx.reset();
                    qy.reset();
                    return *this;
                }

                CompareGreater(CompareGreater &&) = default;
                CompareGreater &operator=(CompareGreater &&) = default;

              private:
                using ExactNT = CGAL::Exact_rational;

                const ExactNT &get_qx() const {
                    if(!qx) {
                        qx = std::make_unique<ExactNT>(query.x());
                    }
                    return *qx;
                }

                const ExactNT &get_qy() const {
                    if(!qy) {
                        qy = std::make_unique<ExactNT>(query.y());
                    }
                    return *qy;
                }

                Point_2 query;
                mutable std::unique_ptr<ExactNT> qx;
                mutable std::unique_ptr<ExactNT> qy;
            };

            using ExactNT = CGAL::Exact_rational;

            Point_2 compute_closest(const Point_2 &query) {
                const FT &bx = m_node->get_xmin();
                const FT &by = m_node->get_ymin();
                const FT &cx = m_node->get_xmax();
                const FT &cy = m_node->get_ymax();
                const FT &qx = query.x();
                const FT &qy = query.y();
                FT close_x, close_y;
                if(qx < bx) {
                    close_x = bx;
                } else if(qx > cx) {
                    close_x = cx;
                } else {
                    close_x = qx;
                }
                if(qy < by) {
                    close_y = by;
                } else if(qy > cy) {
                    close_y = cy;
                } else {
                    close_y = qy;
                }
                return Point_2{close_x, close_y};
            }

            const ExactNT &get_exact_distance(const ExactNT &qx, const ExactNT &qy) const {
                if(!m_exact_distance) {
                    m_exact_distance = std::make_unique<ExactNT>(
                        CGAL::squared_distanceC2(qx, qy, ExactNT(m_distance_to.x()), ExactNT(m_distance_to.y())));
                }
                return *m_exact_distance;
            }

            const Node *m_node;
            Point_2 m_distance_to;
            CGAL::Interval_nt_advanced m_distance_interval;
            mutable std::unique_ptr<ExactNT> m_exact_distance;
        };

        /**
         * Datastructure to store a point (or rather, point handle)
         * with its distance to the query point (relies on an exact
         * number type as FT).
         */
        template<typename AvoidNeedlessInstantiation_> class Point_with_distance_FT_exact {
          public:
            Point_with_distance_FT_exact(const Point_2 &query, Point_handle distance_to)
                : m_distance(CGAL::squared_distance(query, *distance_to)), m_distance_to(distance_to) {}

            Point_handle handle() const noexcept { return m_distance_to; }

            /**
             * Cheaply check if the squared distance from the query
             * to *this is certainly at least the given squared distance.
             */
            bool squared_distance_certainly_at_least(const FT &sqdist) const { return sqdist <= m_distance; }

            FT squared_distance_interval_or_exact() const { return m_distance; }

            class CompareGreater {
              public:
                CompareGreater(const Point_2 & /*query*/) noexcept {}

                /**
                 * Compare two distances.
                 */
                bool operator()(const Point_with_distance_FT_exact &d1, const Point_with_distance_FT_exact &d2) const {
                    return d1.m_distance > d2.m_distance;
                }
            };

            FT m_distance;
            Point_handle m_distance_to;
        };

        /**
         * Datastructure to store a tree node with
         * the distance from its bounding box to a query point.
         */
        template<typename AvoidNeedlessInstantiation_> class Node_with_distance_FT_exact {
          public:
            Node_with_distance_FT_exact(const Point_2 &query, const Node *node)
                : m_distance(mwt::squared_distance(query, node->bbox())), m_node(node) {}

            const Node *node() const noexcept { return m_node; }

            class CompareGreater {
              public:
                CompareGreater(const Point_2 & /*query*/) noexcept {}

                bool operator()(const Node_with_distance_FT_exact &d1, const Node_with_distance_FT_exact &d2) const {
                    return d1.m_distance > d2.m_distance;
                }

                template<typename X_>
                bool operator()(const Node_with_distance_FT_exact &d1,
                                const Point_with_distance_FT_exact<X_> &d2) const {
                    return d1.m_distance > d2.m_distance;
                }
            };

          private:
            FT m_distance;
            const Node *m_node;
        };

        /**
         * Select the type used to store points and distances
         * based on our detection of whether the number type
         * is exact or not.
         */
        template<typename X>
        using Point_with_distance_template =
            std::conditional_t<has_exact_nt, Point_with_distance_FT_exact<X>, Point_with_distance_FT_nonexact<X>>;
        using Point_with_distance = Point_with_distance_template<void>;

        /**
         * Select the type used to store nodes and distances
         * based on our detection of whether the number type
         * is exact or not.
         */
        template<typename X>
        using Node_with_distance_template =
            std::conditional_t<has_exact_nt, Node_with_distance_FT_exact<X>, Node_with_distance_FT_nonexact<X>>;
        using Node_with_distance = Node_with_distance_template<void>;

        /**
         * The query point of our search (i.e., the point
         * from which we want to list the nearest neighbors
         * in increasing order).
         */
        Point_2 m_query;

        /**
         * The search tree.
         */
        const Point_quadtree *m_tree;

        /**
         * The filter.
         */
        Filter m_filter;

        /**
         * The queue of points that we have encountered in
         * expanded nodes.
         */
        Resettable_priority_queue<Point_with_distance, typename Point_with_distance::CompareGreater> m_point_queue;

        /**
         * The queue of nodes that we still need to expand.
         */
        Resettable_priority_queue<Node_with_distance, typename Node_with_distance::CompareGreater> m_node_queue;

        /**
         * Buffer for expanding nodes that contain the query point.
         */
        std::vector<Node_handle> m_zero_node_buffer;

      public:
        /**
         * Create a new incremental search instance.
         */
        template<typename... FilterArgs,
                 std::enable_if_t<sizeof...(FilterArgs) >= 0 && !has_expansion_callback, int> = 0>
        explicit Incremental_search(const Point_quadtree *tree, const Point_2 &query, FilterArgs &&...filter_args)
            : m_query(query), m_tree(tree), m_filter(std::forward<FilterArgs>(filter_args)...),
              m_point_queue(typename Point_with_distance::CompareGreater(query)),
              m_node_queue(typename Node_with_distance::CompareGreater(query)) {
            m_point_queue.reserve(200);
            m_node_queue.reserve(200);
        }

        /**
         * Create a new incremental search instance (constructor for the case with expansion callback).
         */
        template<typename CallbackType, typename... FilterArgs,
                 std::enable_if_t<sizeof...(FilterArgs) >= 0 && has_expansion_callback, int> = 0>
        explicit Incremental_search(const Point_quadtree *tree, const Point_2 &query, CallbackType &&exp_cb,
                                    FilterArgs &&...filter_args)
            : Super(std::forward<CallbackType>(exp_cb)), m_query(query), m_tree(tree),
              m_filter(std::forward<FilterArgs>(filter_args)...),
              m_point_queue(typename Point_with_distance::CompareGreater(query)),
              m_node_queue(typename Node_with_distance::CompareGreater(query)) {
            m_point_queue.reserve(200);
            m_node_queue.reserve(200);
        }

        /**
         * Reset the search to a new query point.
         */
        bool start_search(const Point_2 &point) {
            m_query = point;
            m_point_queue.reset(point);
            m_node_queue.reset(point);
            m_filter.reset(point);
            m_zero_node_buffer.clear();
            p_init_search();
            return !m_point_queue.empty();
        }

        /**
         * Return the current point in the search.
         */
        const Point_with_distance &current() const {
            CGAL_precondition(!m_point_queue.empty());
            return m_point_queue.top();
        }

        /**
         * Move to the next point in the search;
         * if there is not next point, return false.
         */
        bool next() {
            auto &points = m_point_queue;
            auto &nodes = m_node_queue;
            if(m_filter()) {
                points.clear();
                nodes.clear();
                return false;
            }

            const auto &node_greater = nodes.get_comparator();
            if(!points.empty()) {
                // go to next point
                points.pop();
            }

            if(points.empty()) {
                // no points, expand nodes until we have some
                if(!p_expand_until_point()) {
                    // no more unfiltered points
                    return false;
                }
            }

            while(!nodes.empty() && !node_greater(nodes.top(), points.top())) {
                p_expand_top();
            }
            return true;
        }

        /**
         * Invoke callable for each point in the search.
         * Replaces calls to start_search, current, and next.
         */
        template<typename Callable> void for_each(const Point_2 &point, Callable &&callable) {
            for_each(point, false, std::forward<Callable>(callable));
        }

        /**
         * Invoke callable for each point in the search.
         * Replaces calls to start_search, current, and next.
         * If skip_first is true, the first point is not returned;
         * this is used mostly to exclude the query point itself,
         * which has distance 0 to itself (if it is in the point),
         * and will thus always be the first point to be reported.
         * The callable can return void or a type convertible to bool;
         * if it returns bool, its return value is checked and
         * if it is false, the search stops.
         */
        template<typename Callable> void for_each(const Point_2 &point, bool skip_first, Callable &&callable) {
            if(!start_search(point))
                return;
            if(skip_first) {
                if(!next())
                    return;
            }
            for(;;) {
                Point_handle cur = current().handle();
                if(!m_filter(*cur)) {
                    using CallableResult = std::remove_reference_t<decltype(std::invoke(callable, current().handle()))>;
                    if constexpr(std::is_convertible_v<CallableResult, bool>) {
                        auto b = static_cast<bool>(std::invoke(callable, cur));
                        if(!b)
                            return;
                    } else {
                        static_assert(std::is_void_v<std::remove_cv_t<CallableResult>>,
                                      "Callable must return void or convertible to bool!");
                        std::invoke(callable, cur);
                    }
                }
                if(!next())
                    return;
            }
        }

        Filter &filter() { return m_filter; }

        const Filter &filter() const { return m_filter; }

        /**
         * Print the search state as JSON.
         */
        void print_as_json(std::ostream &output) const {
            output << "{ \"open_nodes\": [";
            output << std::setprecision(19);
            bool first = true;
            for(const auto &node : m_node_queue.get_container()) {
                if(!first) {
                    output << ", ";
                } else {
                    first = false;
                }
                output << "{ \"xmin\": " << node.node()->get_xmin() << ", \"ymin\": " << node.node()->get_ymin()
                       << ", \"xmax\": " << node.node()->get_xmax() << ", \"ymax\": " << node.node()->get_ymax()
                       << " }";
            }
            output << "]}";
        }

      private:
        void p_init_search() {
            if(m_tree->empty())
                return;
            if(!m_tree->m_root->contains(m_query)) {
                throw std::invalid_argument("Query point lies outside of the bounding box of the point set!");
            }
            m_zero_node_buffer.push_back(m_tree->m_root);
            p_expand_zero_buffer();
            CGAL_postcondition(!m_point_queue.empty());
            CGAL_postcondition((m_node_queue.empty() || !m_node_queue.top().node()->contains(m_query)));
        }

        void p_expand_zero_buffer() {
            while(!m_zero_node_buffer.empty()) {
                Node_handle node = m_zero_node_buffer.back();
                m_zero_node_buffer.pop_back();
                if(node->empty() || m_filter(node->bbox()))
                    continue;
                if(node->leaf()) {
                    p_expand_leaf(node);
                } else {
                    Super::report_internal_expansion(node, 4);
                    node->expand_children_query_contained(m_query, m_zero_node_buffer, m_node_queue);
                }
            }
        }

        void p_expand_leaf(const Node *node) {
            CGAL_precondition(node->leaf());
            Super::report_leaf_expansion(node);
            for(Point_handle it = node->points_begin(); it != node->points_end(); ++it) {
                m_point_queue.push(Point_with_distance(m_query, it));
            }
        }

        bool p_expand_until_point() {
            while(!m_node_queue.empty()) {
                const Node *node = m_node_queue.top().node();
                m_node_queue.pop();
                if(node->empty())
                    continue;
                if(node->leaf()) {
                    p_expand_leaf(node);
                    return true;
                } else {
                    p_expand_internal(node);
                }
            }
            return false;
        }

        void p_expand_top() {
            const Node *node = m_node_queue.top().node();
            m_node_queue.pop();
            if(node->leaf()) {
                p_expand_leaf(node);
            } else {
                p_expand_internal(node);
            }
        }

        void p_expand_internal(const Node *node) {
            int count = 0;
            const Node *chld = node->children_begin();
            for(int i = 0; i < 4; ++i) {
                if(chld[i].empty() || m_filter(chld[i].bbox()))
                    continue;
                m_node_queue.emplace(m_query, chld + i);
                ++count;
            }
            Super::report_internal_expansion(node, count);
        }
    };

    using Node_allocator = std::conditional_t<std::is_void_v<NodeAllocator>, std::allocator<Node>, NodeAllocator>;
    using Alloc_traits = std::allocator_traits<Node_allocator>;

    Point_quadtree(Iterator points_begin, Iterator points_end, Node_allocator alloc = {})
        : m_points_begin{points_begin}, m_points_end{points_end}, m_root{nullptr}, m_allocator{alloc} {
        p_build_tree();
    }

    ~Point_quadtree() {
        if(m_root) {
            m_root->destroy(m_allocator);
            Alloc_traits::destroy(m_allocator, m_root);
            Alloc_traits::deallocate(m_allocator, m_root, 1);
        }
    }

    bool empty() const noexcept { return !m_root; }

    void verify() const {
        if(m_root)
            m_root->verify();
    }

    Iterator points_begin() const noexcept { return m_points_begin; }
    Iterator points_end() const noexcept { return m_points_end; }
    std::size_t size() const noexcept { return std::distance(m_points_begin, m_points_end); }

  private:
    static Iso_rectangle_2 p_rect(Iterator pbegin, Iterator pend) {
        if(pbegin == pend)
            return Iso_rectangle_2{};
        FT xmin = pbegin->x(), xmax = xmin, ymin = pbegin->y(), ymax = ymin;
        for(auto it = pbegin + 1; it != pend; ++it) {
            const FT x = it->x();
            const FT y = it->y();
            if(x < xmin) {
                xmin = x;
            } else if(x > xmax) {
                xmax = x;
            }
            if(y < ymin) {
                ymin = y;
            } else if(y > ymax) {
                ymax = y;
            }
        }
        return Iso_rectangle_2{xmin, ymin, xmax, ymax};
    }

    void p_build_tree() {
        if(m_points_begin == m_points_end)
            return;
        auto bbox = p_rect(m_points_begin, m_points_end);
        m_root = Alloc_traits::allocate(m_allocator, 1);
        Alloc_traits::construct(m_allocator, m_root, m_points_begin, m_points_end, bbox.xmin(), bbox.ymin(),
                                bbox.xmax(), bbox.ymax());
        p_hilbert_split_node<Order::XY, false>(m_root);
    }

    template<Order order, bool reverse> void p_hilbert_split_node(Node *node) {
        if(std::distance(node->points_begin(), node->points_end()) < SPLIT_THRESHOLD) {
            return;
        }
        if(p_hilbert_split_normal<order, reverse>(node)) {
            return;
        }
        p_hilbert_split_balanced<order, reverse>(node);
    }

    template<Order order, bool reverse> bool p_hilbert_split_normal(Node *node) {
        const FT xmid = FT{0.5} * (node->get_xmin() + node->get_xmax());
        const FT ymid = FT{0.5} * (node->get_ymin() + node->get_ymax());

        Iterator m[5];
        m[0] = node->points_begin();
        m[4] = node->points_end();
        if constexpr(order == Order::XY) {
            m[2] = std::partition(m[0], m[4], Fixed_value_cmp_2<Kernel, Dimension::X, reverse>{xmid});
            m[1] = std::partition(m[0], m[2], Fixed_value_cmp_2<Kernel, Dimension::Y, reverse>{ymid});
            m[3] = std::partition(m[2], m[4], Fixed_value_cmp_2<Kernel, Dimension::Y, !reverse>{ymid});
        } else {
            m[2] = std::partition(m[0], m[4], Fixed_value_cmp_2<Kernel, Dimension::Y, reverse>{ymid});
            m[1] = std::partition(m[0], m[2], Fixed_value_cmp_2<Kernel, Dimension::X, reverse>{xmid});
            m[3] = std::partition(m[2], m[4], Fixed_value_cmp_2<Kernel, Dimension::X, !reverse>{xmid});
        }

        std::size_t max_size = 0;
        std::size_t total_size(std::distance(node->points_begin(), node->points_end()));
        for(int i = 0; i < 4; ++i) {
            std::size_t size = std::distance(m[i], m[i + 1]);
            if(size > max_size)
                max_size = size;
        }
        if(total_size * BALANCE_THRESHOLD <= max_size)
            return false;

        CreateChildren<Split<order>::same, reverse>{this}.create(node, +m, xmid, ymid);
        return true;
    }

    template<Order order, bool reverse> void p_hilbert_split_balanced(Node *node) {
        Iterator m[5];
        m[0] = node->points_begin();
        m[4] = node->points_end();
        Iterator mid = m[0] + std::distance(m[0], m[4]) / 2;
        FT xmid;
        FT ymid;
        if constexpr(order == Order::XY) {
            std::nth_element(m[0], mid, m[4], Coordinate_cmp_2<Kernel, Dimension::X, reverse>{});
            xmid = mid->x();
            m[2] = std::partition(m[0], m[4], Fixed_value_cmp_2<Kernel, Dimension::X, reverse>{xmid});
            if(std::distance(m[0], m[2]) >= std::distance(m[2], m[4])) {
                mid = m[0] + std::distance(m[0], m[2]) / 2;
                std::nth_element(m[0], mid, m[2], Coordinate_cmp_2<Kernel, Dimension::Y, reverse>{});
                ymid = mid->y();
            } else {
                mid = m[2] + std::distance(m[2], m[4]) / 2;
                std::nth_element(m[2], mid, m[4], Coordinate_cmp_2<Kernel, Dimension::Y, reverse>{});
                ymid = mid->y();
            }
            m[1] = std::partition(m[0], m[2], Fixed_value_cmp_2<Kernel, Dimension::Y, reverse>{ymid});
            m[3] = std::partition(m[2], m[4], Fixed_value_cmp_2<Kernel, Dimension::Y, !reverse>{ymid});
        } else {
            std::nth_element(m[0], mid, m[4], Coordinate_cmp_2<Kernel, Dimension::X, reverse>{});
            ymid = mid->y();
            m[2] = std::partition(m[0], m[4], Fixed_value_cmp_2<Kernel, Dimension::Y, reverse>{ymid});
            if(std::distance(m[0], m[2]) >= std::distance(m[2], m[4])) {
                mid = m[0] + std::distance(m[0], m[2]) / 2;
                std::nth_element(m[0], mid, m[2], Coordinate_cmp_2<Kernel, Dimension::X, reverse>{});
                xmid = mid->x();
            } else {
                mid = m[2] + std::distance(m[2], m[4]) / 2;
                std::nth_element(m[2], mid, m[4], Coordinate_cmp_2<Kernel, Dimension::X, reverse>{});
                xmid = mid->x();
            }
            m[1] = std::partition(m[0], m[2], Fixed_value_cmp_2<Kernel, Dimension::X, reverse>{xmid});
            m[3] = std::partition(m[2], m[4], Fixed_value_cmp_2<Kernel, Dimension::X, !reverse>{xmid});
        }

        CreateChildren<Split<order>::same, reverse>{this}.create(node, +m, xmid, ymid);
    }

    template<Order order, bool reverse> struct CreateChildren {
        Point_quadtree *that;

        void create(Node *node, Iterator *m, FT xmid, FT ymid) {
            Node_allocator &alloc = that->m_allocator;
            node->m_children = alloc.allocate(4);
            const FT xmin = node->get_xmin();
            const FT xmax = node->get_xmax();
            const FT ymin = node->get_ymin();
            const FT ymax = node->get_ymax();
            auto *children = node->m_children;
            if constexpr(order == Order::XY && !reverse) {
                Alloc_traits::construct(alloc, children + 0, m[0], m[1], xmin, ymin, xmid, ymid);
                Alloc_traits::construct(alloc, children + 1, m[1], m[2], xmin, ymid, xmid, ymax);
                Alloc_traits::construct(alloc, children + 2, m[2], m[3], xmid, ymid, xmax, ymax);
                Alloc_traits::construct(alloc, children + 3, m[3], m[4], xmid, ymin, xmax, ymid);

                that->p_hilbert_split_node<Split<order>::other, reverse>(children + 0);
                that->p_hilbert_split_node<Split<order>::same, reverse>(children + 1);
                that->p_hilbert_split_node<Split<order>::same, reverse>(children + 2);
                that->p_hilbert_split_node<Split<order>::other, !reverse>(children + 3);
            } else if constexpr(order == Order::XY && reverse) {
                Alloc_traits::construct(alloc, children + 0, m[2], m[3], xmin, ymin, xmid, ymid);
                Alloc_traits::construct(alloc, children + 1, m[3], m[4], xmin, ymid, xmid, ymax);
                Alloc_traits::construct(alloc, children + 2, m[0], m[1], xmid, ymid, xmax, ymax);
                Alloc_traits::construct(alloc, children + 3, m[1], m[2], xmid, ymin, xmax, ymid);

                that->p_hilbert_split_node<Split<order>::same, reverse>(children + 0);
                that->p_hilbert_split_node<Split<order>::other, !reverse>(children + 1);
                that->p_hilbert_split_node<Split<order>::other, reverse>(children + 2);
                that->p_hilbert_split_node<Split<order>::same, reverse>(children + 3);
            } else if constexpr(order == Order::YX && !reverse) {
                Alloc_traits::construct(alloc, children + 0, m[0], m[1], xmin, ymin, xmid, ymid);
                Alloc_traits::construct(alloc, children + 1, m[3], m[4], xmin, ymid, xmid, ymax);
                Alloc_traits::construct(alloc, children + 2, m[2], m[3], xmid, ymid, xmax, ymax);
                Alloc_traits::construct(alloc, children + 3, m[1], m[2], xmid, ymin, xmax, ymid);

                that->p_hilbert_split_node<Split<order>::other, reverse>(children + 0);
                that->p_hilbert_split_node<Split<order>::other, !reverse>(children + 1);
                that->p_hilbert_split_node<Split<order>::same, reverse>(children + 2);
                that->p_hilbert_split_node<Split<order>::same, reverse>(children + 3);
            } else {
                Alloc_traits::construct(alloc, children + 0, m[2], m[3], xmin, ymin, xmid, ymid);
                Alloc_traits::construct(alloc, children + 1, m[1], m[2], xmin, ymid, xmid, ymax);
                Alloc_traits::construct(alloc, children + 2, m[0], m[1], xmid, ymid, xmax, ymax);
                Alloc_traits::construct(alloc, children + 3, m[3], m[4], xmid, ymin, xmax, ymid);

                that->p_hilbert_split_node<Split<order>::same, reverse>(children + 0);
                that->p_hilbert_split_node<Split<order>::same, reverse>(children + 1);
                that->p_hilbert_split_node<Split<order>::other, reverse>(children + 2);
                that->p_hilbert_split_node<Split<order>::other, !reverse>(children + 3);
            }
        }
    };

    Iterator m_points_begin, m_points_end;
    Node *m_root;
    Node_allocator m_allocator;
};

template<typename Traits_, typename RandomAccessIterator> class Static_quadtree {
  public:
    using Traits = Traits_;
    using Iterator = RandomAccessIterator;

    using FT = typename Traits::FT;
    using Point_2 = typename Traits::Point_2;

    using Node = Quadtree_node<Traits, Iterator>;
    using Node_allocator = std::allocator<Node>;
    using Alloc_traits = std::allocator_traits<Node_allocator>;

    Static_quadtree(Iterator begin, Iterator end)
        : points_begin{begin}, points_end{end}, node_allocator{}, root_node{nullptr} {
        build_tree(begin, end);
    }

    ~Static_quadtree() {
        if(root_node != nullptr) {
            deallocate_children(root_node);
            Alloc_traits::destroy(node_allocator, root_node);
            node_allocator.deallocate(root_node, 1);
        }
    }

    void verify_tree() { verify_node(root_node); }

    void verify_node(const Node *node) {
        CGAL_assertion(std::none_of(node->begin, node->end,
                                    [&](const Point_2 &p) { return node->bbox.has_on_unbounded_side(p); }));

        if(node->is_leaf()) {
            CGAL_assertion(std::distance(node->begin, node->end) < SPLIT_THRESHOLD);
        } else {
            for(int i{0}; i < 4; ++i) {
                verify_node(node->children + i);
            }
        }
    }

    bool is_empty() const { return root_node == nullptr; }

    const Node *root() const { return root_node; }

    Iterator begin() const { return points_begin; }

    Iterator end() const { return points_end; }

  private:
    void deallocate_children(Node *node) {
        if(node->children == nullptr) {
            return;
        }

        deallocate_children(node->children + 0);
        deallocate_children(node->children + 1);
        deallocate_children(node->children + 2);
        deallocate_children(node->children + 3);

        Alloc_traits::destroy(node_allocator, node->children + 0);
        Alloc_traits::destroy(node_allocator, node->children + 1);
        Alloc_traits::destroy(node_allocator, node->children + 2);
        Alloc_traits::destroy(node_allocator, node->children + 3);

        node_allocator.deallocate(node->children, 4);
    }

    void build_tree(Iterator begin, Iterator end) {
        if(begin == end)
            return;

        auto bbox = bbox_2(begin, end);

        root_node = node_allocator.allocate(1);
        Alloc_traits::construct(node_allocator, root_node, bbox, begin, end);

        hilbert_split_node<Order::XY, false>(root_node);
    }

    template<Order order, bool reverse> void hilbert_split_node(Node *node) {
        if(std::distance(node->begin, node->end) < SPLIT_THRESHOLD) {
            return;
        }

        const FT xmin = node->bounding_box().xmin();
        const FT ymin = node->bounding_box().ymin();
        const FT xmax = node->bounding_box().xmax();
        const FT ymax = node->bounding_box().ymax();

        const FT xmid = FT{0.5} * (xmin + xmax);
        const FT ymid = FT{0.5} * (ymin + ymax);

        Iterator m0 = node->begin;
        Iterator m4 = node->end;
        Iterator m1, m2, m3;

        if(order == Order::XY) {
            m2 = std::partition(m0, m4, Fixed_value_cmp_2<Traits, Dimension::X, reverse>{xmid});
            m1 = std::partition(m0, m2, Fixed_value_cmp_2<Traits, Dimension::Y, reverse>{ymid});
            m3 = std::partition(m2, m4, Fixed_value_cmp_2<Traits, Dimension::Y, !reverse>{ymid});
        }
        if(order == Order::YX) {
            m2 = std::partition(m0, m4, Fixed_value_cmp_2<Traits, Dimension::Y, reverse>{ymid});
            m1 = std::partition(m0, m2, Fixed_value_cmp_2<Traits, Dimension::X, reverse>{xmid});
            m3 = std::partition(m2, m4, Fixed_value_cmp_2<Traits, Dimension::X, !reverse>{xmid});
        }

        node->children = node_allocator.allocate(4);

        if(order == Order::XY && !reverse) {
            Alloc_traits::construct(node_allocator, node->children + 0, xmin, ymin, xmid, ymid, m0, m1);
            Alloc_traits::construct(node_allocator, node->children + 1, xmin, ymid, xmid, ymax, m1, m2);
            Alloc_traits::construct(node_allocator, node->children + 2, xmid, ymid, xmax, ymax, m2, m3);
            Alloc_traits::construct(node_allocator, node->children + 3, xmid, ymin, xmax, ymid, m3, m4);

            hilbert_split_node<Split<order>::other, reverse>(node->children + 0);
            hilbert_split_node<Split<order>::same, reverse>(node->children + 1);
            hilbert_split_node<Split<order>::same, reverse>(node->children + 2);
            hilbert_split_node<Split<order>::other, !reverse>(node->children + 3);
        }
        if(order == Order::XY && reverse) {
            Alloc_traits::construct(node_allocator, node->children + 0, xmin, ymin, xmid, ymid, m2, m3);
            Alloc_traits::construct(node_allocator, node->children + 1, xmin, ymid, xmid, ymax, m3, m4);
            Alloc_traits::construct(node_allocator, node->children + 2, xmid, ymid, xmax, ymax, m0, m1);
            Alloc_traits::construct(node_allocator, node->children + 3, xmid, ymin, xmax, ymid, m1, m2);

            hilbert_split_node<Split<order>::same, reverse>(node->children + 0);
            hilbert_split_node<Split<order>::other, !reverse>(node->children + 1);
            hilbert_split_node<Split<order>::other, reverse>(node->children + 2);
            hilbert_split_node<Split<order>::same, reverse>(node->children + 3);
        }
        if(order == Order::YX && !reverse) {
            Alloc_traits::construct(node_allocator, node->children + 0, xmin, ymin, xmid, ymid, m0, m1);
            Alloc_traits::construct(node_allocator, node->children + 1, xmin, ymid, xmid, ymax, m3, m4);
            Alloc_traits::construct(node_allocator, node->children + 2, xmid, ymid, xmax, ymax, m2, m3);
            Alloc_traits::construct(node_allocator, node->children + 3, xmid, ymin, xmax, ymid, m1, m2);

            hilbert_split_node<Split<order>::other, reverse>(node->children + 0);
            hilbert_split_node<Split<order>::other, !reverse>(node->children + 1);
            hilbert_split_node<Split<order>::same, reverse>(node->children + 2);
            hilbert_split_node<Split<order>::same, reverse>(node->children + 3);
        }
        if(order == Order::YX && reverse) {
            Alloc_traits::construct(node_allocator, node->children + 0, xmin, ymin, xmid, ymid, m2, m3);
            Alloc_traits::construct(node_allocator, node->children + 1, xmin, ymid, xmid, ymax, m1, m2);
            Alloc_traits::construct(node_allocator, node->children + 2, xmid, ymid, xmax, ymax, m0, m1);
            Alloc_traits::construct(node_allocator, node->children + 3, xmid, ymin, xmax, ymid, m3, m4);

            hilbert_split_node<Split<order>::same, reverse>(node->children + 0);
            hilbert_split_node<Split<order>::same, reverse>(node->children + 1);
            hilbert_split_node<Split<order>::other, reverse>(node->children + 2);
            hilbert_split_node<Split<order>::other, !reverse>(node->children + 3);
        }
    }

    const Iterator points_begin;
    const Iterator points_end;

    Node_allocator node_allocator;

    Node *root_node;
};

} // namespace mwt

#endif
