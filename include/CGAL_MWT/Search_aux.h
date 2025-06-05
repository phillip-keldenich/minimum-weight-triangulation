#ifndef CGAL_MWT_SEARCH_AUX_H_INCLUDED_
#define CGAL_MWT_SEARCH_AUX_H_INCLUDED_

#include <queue>
#include <vector>

namespace mwt {

/**
 * Extension of std::priority_queue that allows
 * resetting the comparison predicate (emptying the queue in the process).
 * Essentially allows us to avoid re-allocating the container all the time.
 */
template<typename T, typename Comparison>
class Resettable_priority_queue : public std::priority_queue<T, std::vector<T>, Comparison> {
    using Base = std::priority_queue<T, std::vector<T>, Comparison>;

  public:
    using Base::Base;

    /**
     * Clear the queue and construct a new comparison predicate
     * from the given arguments.
     */
    template<typename... ComparisonArgs> void reset(ComparisonArgs &&...compargs) {
        this->clear();
        this->comp = Comparison{std::forward<ComparisonArgs>(compargs)...};
    }

    void clear() { this->c.clear(); }

    void reserve(std::size_t n) { this->c.reserve(n); }

    const Comparison &get_comparator() const noexcept { return this->comp; }

    const std::vector<T> &get_container() const noexcept { return this->c; }
};

} // namespace mwt

#endif
