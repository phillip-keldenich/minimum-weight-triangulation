#ifndef MWT_FILTERED_INCREMENTAL_SEARCH_H
#define MWT_FILTERED_INCREMENTAL_SEARCH_H

#include <array>
#include <boost/container/flat_set.hpp>
#include <boost/heap/d_ary_heap.hpp>
#include <type_traits>
#include <utility>
#include <vector>

#include "Static_quadtree.h"

namespace mwt {

template<typename Tree, template<typename T> class Filter_> class Filtered_incremental_search {
  public:
    using Traits = typename Tree::Traits;
    using Point_iterator = typename Tree::Iterator;
    using Node = typename Tree::Node;
    using Filter = Filter_<Traits>;

    using Node_handle = const Node *;
    using Point_handle = Point_iterator;

    using FT = typename Traits::FT;
    using Point_2 = typename Traits::Point_2;

    template<typename Handle> struct Handle_with_distance {
        Handle handle;
        FT distance;

        Handle_with_distance(Handle h, FT d) : handle{h}, distance{d} {};

        bool operator<(const Handle_with_distance &rhs) const { return distance < rhs.distance; }
        bool operator<=(const Handle_with_distance &rhs) const { return distance <= rhs.distance; }
        bool operator>(const Handle_with_distance &rhs) const { return distance > rhs.distance; }
        bool operator>=(const Handle_with_distance &rhs) const { return distance >= rhs.distance; }
    };

    using Node_with_distance = Handle_with_distance<Node_handle>;
    using Point_with_distance = Handle_with_distance<Point_handle>;

    using Node_queue = boost::heap::d_ary_heap<Node_with_distance, boost::heap::arity<4>,
                                               boost::heap::compare<std::greater<Node_with_distance>>>;

    using Point_queue = boost::heap::d_ary_heap<Point_with_distance, boost::heap::arity<4>,
                                                boost::heap::compare<std::greater<Point_with_distance>>>;

    template<bool b> using Contains_query = std::integral_constant<bool, b>;

    Filtered_incremental_search(const Tree &tree) : tree{tree}, filter{} {
        nodes.reserve(200);
        points.reserve(250);
    }

    bool search(const Point_2 &q) {
        query = q;
        filter.reset(q);
        nodes.clear();
        points.clear();
        init_search();
        return !points.empty();
    }

    bool next() {
        if(!points.empty()) {
            points.pop();
        }
        while(!nodes.empty() && points.empty()) {
            expand_top_node(Contains_query<false>{});
        }
        if(!nodes.empty() && !points.empty()) {
            while(!nodes.empty() && nodes.top().distance <= points.top().distance) {
                expand_top_node(Contains_query<false>{});
            }
        }
        return !points.empty();
    }

    const Point_with_distance &current() const {
        CGAL_precondition(!points.empty());
        return points.top();
    }

    void init_search() {
        if(tree.is_empty())
            return;

        nodes.emplace(tree.root(), FT{0});

        // Expand all nodes that contain the query
        while(!nodes.empty() && nodes.top().distance == FT{0}) {
            expand_top_node(Contains_query<true>{});
        }
    }

    void expand_leaf_node() {
        CGAL_precondition(!nodes.empty());
        Node_handle node = nodes.top().handle;
        CGAL_precondition(node->is_leaf());

        nodes.pop();

        for(auto p = node->begin; p != node->end; ++p) {
            const FT distance = typename Traits::Compute_squared_distance_2{}(*p, query);
            points.emplace(p, distance);
        }
    }

    void expand_top_node(Contains_query<true>) {
        CGAL_precondition(!nodes.empty());

        while(!nodes.top().handle->is_leaf()) {
            Node_handle node = nodes.top().handle;
            bool top_reused = false;
            const auto distance = node->squared_distance_to_children(query, true);
            for(int i{0}; i < 4; ++i) {
                Node_handle child = node->children + i;
                if(!child->is_empty()) {
                    if(!top_reused && distance[i] == 0) {
                        // There is always at least one child with distance 0.
                        // Thus, we can replace the top node and save one pop and push in each iteration
                        Node_with_distance &top = const_cast<Node_with_distance &>(nodes.top());
                        top.handle = child;
                        top_reused = true;
                    } else {
                        nodes.emplace(child, distance[i]);
                    }
                }
            }
        }
        expand_leaf_node();
    }

    void expand_top_node(Contains_query<false>) {
        CGAL_precondition(!nodes.empty());

        while(!nodes.top().handle->is_leaf()) {
            Node_handle node = nodes.top().handle;
            nodes.pop();

            const auto distance = node->squared_distance_to_children(query, false);
            for(int i{0}; i < 4; ++i) {
                Node_handle child = node->children + i;
                if(!child->is_empty() && !filter(child->bbox)) {
                    nodes.emplace(child, distance[i]);
                }
            }

            if(nodes.empty())
                return;
        }
        expand_leaf_node();
    }

    const Tree &tree;
    Point_2 query;
    Filter filter;
    Node_queue nodes;
    Point_queue points;
};

} // namespace mwt

#endif // MWT_FILTERED_INCREMENTAL_SEARCH_H
