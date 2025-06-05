#ifndef MWT_INTERVAL_SET_H
#define MWT_INTERVAL_SET_H

#include <algorithm>
#include <iostream>
#include <iterator>
#include <utility>
#include <vector>

#include <CGAL/assertions.h>

namespace mwt {

/**
 * A simple Interval template with ordering operators (<,<=)
 * defined by the LB (left side) value.
 */
template<class T> struct Interval {
    T left;
    T right;

    Interval(const T &l, const T &r) : left{l}, right{r} {}

    bool operator<(const Interval &rhs) const { return left < rhs.left; }
    bool operator<=(const Interval &rhs) const { return left <= rhs.left; }

    friend inline std::ostream &operator<<(std::ostream &stream, const Interval &interval) {
        return stream << "[" << interval.left << ", " << interval.right << "]";
    }
};

template<class T> inline bool operator<(const Interval<T> &lhs, const T &val) { return lhs.left < val; }
template<class T> inline bool operator<(const T &val, const Interval<T> &rhs) { return val < rhs.left; }

template<class T> inline bool operator<=(const Interval<T> &lhs, const T &val) { return lhs.left <= val; }
template<class T> inline bool operator<=(const T &val, const Interval<T> &rhs) { return val <= rhs.left; }

/**
 * A set of intervals.
 */
template<class T> class Interval_set {
  public:
    using Interval_type = Interval<T>;
    using Interval_container = std::vector<Interval_type>;
    using Iterator = typename Interval_container::iterator;
    using Const_iterator = typename Interval_container::const_iterator;

    /**
     * Create an interval set covering a certain total interval.
     */
    Interval_set(T min, T max) : min{min}, max{max} {
        CGAL_precondition(min < max);
        intervals.reserve(16);
    }

    Const_iterator begin() const { return intervals.cbegin(); }

    Const_iterator end() const { return intervals.cend(); }

    void insert(const T &left, const T &right) {
        if(intervals.empty()) {
            intervals.emplace_back(left, right);
            return;
        }

        auto lb = priv_lower_bound(intervals.begin(), intervals.end(), left);
        if(lb != intervals.begin() && left <= std::prev(lb)->right) {
            //|__|  |_lb_|
            //  |___....
            auto prev = std::prev(lb);

            if(prev->right < right) {
                prev->right = right;

                auto i = lb;
                for(; i != intervals.end() && i->left <= right; ++i) {
                    prev->right = std::max(i->right, right);
                }
                intervals.erase(lb, i);
            }
        } else {
            //|__|     |_lb_|
            //     |__....
            if(lb == intervals.end() || right < lb->left) {
                //|__|     |_lb_|
                //     |__|
                intervals.emplace(lb, left, right);
            } else {
                //|__|     |_lb_|
                //     |_____......
                lb->left = left;

                auto i = lb;
                for(; i != intervals.end() && i->left <= right; ++i) {
                    lb->right = std::max(i->right, right);
                }

                if(i > ++lb) {
                    intervals.erase(lb, i);
                }
            }
        }
    }

    bool empty() const { return intervals.empty(); }

    void clear() { intervals.clear(); }

    bool completely_covered() const {
        return !intervals.empty() && ((intervals.front().left <= min) & (intervals.front().right >= max));
    }

    bool contains(const T &left, const T &right) const {
        if(intervals.empty())
            return false;

        auto ub = priv_upper_bound(intervals.begin(), intervals.end(), left);
        if(ub != intervals.begin()) {
            --ub;
            return right <= ub->right;
        } else {
            return false;
        }
    }

    // Special function if the interval wraps around
    bool contains(const T &left, const T &right, bool) const {
        if(intervals.empty()) {
            return false;
        }
        const auto &b = intervals.back();
        const auto &f = intervals.front();
        return (b.left <= left) & (b.right >= max) & (f.left <= right) & (f.right >= right);
    }

    bool contains(const T &val) const {
        if(intervals.empty()) {
            return false;
        }
        auto ub = priv_upper_bound(intervals.begin(), intervals.end(), val);
        if(ub != intervals.begin()) {
            --ub;
            return val <= ub->right;
        } else {
            return false;
        }
    }

    void print() const {
        for(auto &i : intervals) {
            std::cout << i << " ";
        }
        std::cout << std::endl;
    }

    size_t size() const { return intervals.size(); }

    T get_range_min() const { return min; }
    T get_range_max() const { return max; }

  private:
    template<typename It, typename V> static It priv_lower_bound(It begin, It end, const V &val) {
        auto len = std::distance(begin, end);
        while(len > 1) {
            const auto half = len >> 1;
            auto mid = begin + half;

            begin = (*mid < val) ? mid : begin;
            len -= half;
        }
        return begin + (*begin < val);
    }

    template<typename It, typename V> static It priv_upper_bound(It begin, It end, const V &val) {
        auto len = std::distance(begin, end);
        while(len > 1) {
            const auto half = len >> 1;
            auto mid = begin + half;

            begin = (*mid <= val) ? mid : begin;
            len -= half;
        }
        return begin + (*begin <= val);
    }

    const T min;
    const T max;
    Interval_container intervals;
};

template<typename T> class Circular_interval_set {
  public:
    Circular_interval_set(const T min, const T max) : set{min, max} {}

    void insert(const T &left, const T &right) {
        if(left <= right) {
            set.insert(left, right);
        } else {
            // Interval wraps around.
            set.insert(left, get_range_max());
            set.insert(get_range_min(), right);
        }
    }

    bool contains(const T &left, const T &right) const {
        if(left <= right) {
            return set.contains(left, right);
        } else {
            return set.contains(left, right, true);
        }
    }

    bool contains(const std::pair<T, T> &interval) const { return contains(interval.first, interval.second); }

    bool contains(const T &value) const { return set.contains(value); }

    void clear() { set.clear(); }

    bool completely_covered() const { return set.completely_covered(); }

    size_t size() const { return set.size(); }

    void print() const { set.print(); }

    T get_range_min() const { return set.get_range_min(); }

    T get_range_max() const { return set.get_range_max(); }

  private:
    using Interval_set = mwt::Interval_set<T>;
    Interval_set set;
};

} // namespace mwt

#endif // MWT_INTERVAL_SET_H
