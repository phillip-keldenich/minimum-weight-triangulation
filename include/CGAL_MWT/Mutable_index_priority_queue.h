#ifndef CGAL_MWT_MUTABLE_INDEX_PRIORITY_QUEUE_H_INCLUDED_
#define CGAL_MWT_MUTABLE_INDEX_PRIORITY_QUEUE_H_INCLUDED_

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <exception>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <vector>

#include <iostream>

namespace mipq {

struct DefaultParameters {
    template<typename KeyType> using KeyCompare = std::less<KeyType>;
    static constexpr bool cache_key_in_heap = true;
};

/**
 * Implements a mutable priority queue of indices.
 * These indices are usually indices of elements in some container
 * that need to be kept track of in a priority queue;
 * the mapping from index -> item should not usually change while using the queue.
 *
 * The number of indices can not cheaply be increased (essentially,
 * requires resetting the container).
 *
 * The priority queue needs to have a 'get key' functor that returns
 * the comparison key of a given index (similar to, e.g., python's sort(..., key=...)).
 * Optionally, the key can be cached in the queue (this is useful if the key is relatively expensive to compute).
 *
 * Internally, the priority queue keeps track of the position in the heap of each item index.
 * This allows the priority queue to efficiently support the following operations:
 *  - initialize the queue with a given set of indices from 0 to n-1 (O(n)).
 *  - determine whether a given index is in the queue (O(1), essentially one vector access)
 *  - determine the index with the smallest key (O(1), accessing the first element of the vector)
 *  - pop the index with the smallest key (O(log n)).
 *  - erase any item index from the queue O(log n).
 *  - notify the queue that the key of an item has changed (O(log n)).
 */
template<typename GetKeyFunctor_ /*(std::size_t index)*/, typename Parameters_ = DefaultParameters>
class MutableIndexPriorityQueue {
  public:
    using GetKeyFunctor = GetKeyFunctor_;
    using Parameters = Parameters_;
    using KeyType = std::decay_t<decltype(std::declval<const GetKeyFunctor &>()(std::declval<std::size_t>()))>;
    using KeyCompare = typename Parameters::template KeyCompare<KeyType>;

    explicit MutableIndexPriorityQueue(const GetKeyFunctor &get_key, std::size_t n = 0) : m_get_key(get_key) {
        reset_queue(n);
    }

    explicit MutableIndexPriorityQueue(GetKeyFunctor &&get_key, std::size_t n = 0) : m_get_key(std::move(get_key)) {
        reset_queue(n);
    }

    /**
     * Reset the queue to contain the indices from 0 to n-1.
     */
    void reset_queue(std::size_t n) {
        m_n = n;
        m_element_index_to_heap_index.resize(n);
        std::iota(m_element_index_to_heap_index.begin(), m_element_index_to_heap_index.end(), std::size_t(0));
        m_heap_elements.clear();
        m_heap_elements.reserve(n);
        for(std::size_t i = 0; i < n; ++i) {
            if constexpr(Parameters::cache_key_in_heap) {
                m_heap_elements.push_back({i, m_get_key(i)});
            } else {
                m_heap_elements.push_back({i});
            }
        }
        p_make_heap();
    }

    /**
     * Reset the queue to contain the indices from 0 to n-1;
     * also reset the get key functor.
     */
    void reset_queue(std::size_t n, const GetKeyFunctor &get_key) {
        m_get_key = get_key;
        reset_queue(n);
    }

    /**
     * Reset the queue to contain the indices from 0 to n-1;
     * also reset the get key functor.
     */
    void reset_queue(std::size_t n, GetKeyFunctor &&get_key) {
        m_get_key = std::move(get_key);
        reset_queue(n);
    }

    /**
     * Check if the queue is empty.
     */
    bool empty() const noexcept { return m_n == 0; }

    /**
     * Get the number of elements in the queue.
     */
    std::size_t size() const noexcept { return m_n; }

    /**
     * Get the maximum index that can be in the queue.
     */
    std::size_t index_bound() const noexcept { return m_element_index_to_heap_index.size(); }

    /**
     * Check if the given index is in the queue.
     */
    bool is_in_queue(std::size_t index) const noexcept { return m_element_index_to_heap_index[index] < m_n; }

    /**
     * Re-insert a deleted index from 0, ..., n-1
     * into the queue.
     */
    void reinsert(std::size_t index) noexcept {
        m_element_index_to_heap_index[index] = m_n;
        if constexpr(Parameters::cache_key_in_heap) {
            m_heap_elements.push_back({index, m_get_key(index)});
        } else {
            m_heap_elements.push_back({index});
        }
        ++m_n;
        p_bubble_up(m_n - 1);
    }

    /**
     * Erase an index from the queue.
     * Returns true if the index was in the queue.
     */
    bool erase(std::size_t s) {
        std::size_t heap_index = m_element_index_to_heap_index[s];
        if(heap_index >= m_n) {
            return false;
        }
        p_swap_heap_indices(heap_index, --m_n);
        m_heap_elements.pop_back();
        m_element_index_to_heap_index[s] = std::numeric_limits<std::size_t>::max();
        p_sift_down(heap_index);
        return true;
    }

    std::size_t top() const noexcept { return m_heap_elements[0].index; }

    std::size_t pop() noexcept {
        assert(!empty());
        std::size_t result = m_heap_elements[0].index;
        m_element_index_to_heap_index[result] = std::numeric_limits<std::size_t>::max();
        p_swap_heap_indices(0, --m_n);
        m_heap_elements.pop_back();
        p_sift_down(0);
        return result;
    }

    void notify_key_increased(std::size_t index) noexcept {
        std::size_t heap_index = m_element_index_to_heap_index[index];
        p_recache_key(heap_index, index);
        assert(heap_index < m_n);
        p_sift_down(heap_index);
    }

    void notify_key_decreased(std::size_t index) noexcept {
        std::size_t heap_index = m_element_index_to_heap_index[index];
        p_recache_key(heap_index, index);
        assert(heap_index < m_n);
        p_bubble_up(heap_index);
    }

    void write_heap(std::ostream &out) {
        std::size_t layer_size = 1;
        std::size_t current_in_layer = 0;
        for(std::size_t i = 0; i < m_n; ++i) {
            out << "(" << m_heap_elements[i].index << ", " << p_get_key(m_heap_elements[i]) << ") ";
            if(++current_in_layer == layer_size) {
                std::cout << std::endl;
                layer_size *= 4;
                current_in_layer = 0;
            }
        }
        if(current_in_layer) {
            out << std::endl;
        }
    }

    /**
     * Use this method to notify the queue that the key of the given index has changed.
     * Like the other notify_key_... methods, this method will
     * not work to update more than one key at a time.
     */
    void notify_key_changed(std::size_t index) noexcept {
        std::size_t heap_index = m_element_index_to_heap_index[index];
        assert(heap_index < m_n);
        if constexpr(Parameters::cache_key_in_heap) {
            KeyType old_key = p_get_key(m_heap_elements[heap_index]);
            KeyType new_key = m_get_key(index);
            m_heap_elements[heap_index].cached_key = new_key;
            if(!KeyCompare{}(old_key, new_key)) {
                // new_key <= old_key
                p_bubble_up(heap_index);
                return;
            } else {
                // new_key > old_key
                p_sift_down(heap_index);
                return;
            }
        } else {
            p_bubble_up(heap_index);
            p_sift_down(heap_index);
        }
    }

    /**
     * Testing/debugging method that checks that the
     * classes' invariants are satisfied.
     */
    void verify_invariants() const {
        for(std::size_t i = 0; i < m_n; ++i) {
            std::size_t c0 = p_child3(i) - 3;
            std::size_t cmax = (std::min)(c0 + 4, m_n);
            KeyType pval = p_get_key(m_heap_elements[i]);
            if constexpr(Parameters::cache_key_in_heap) {
                KeyType noncached = m_get_key(m_heap_elements[i].index);
                if(KeyCompare{}(pval, noncached) || KeyCompare{}(noncached, pval)) {
                    throw std::runtime_error("Heap invariant violated: cached key does not match non-cached key");
                }
            }
            for(std::size_t c = c0; c < cmax; ++c) {
                KeyType cval = p_get_key(m_heap_elements[c]);
                if(KeyCompare{}(p_get_key(m_heap_elements[c]), p_get_key(m_heap_elements[i]))) {
                    throw std::runtime_error("Heap invariant violated: child has smaller key than parent");
                }
            }
        }
        for(std::size_t original_index = 0; original_index < index_bound(); ++original_index) {
            std::size_t heap_index = m_element_index_to_heap_index[original_index];
            if(heap_index >= m_n) {
                continue;
            }
            if(m_heap_elements[heap_index].index != original_index) {
                throw std::runtime_error("Heap invariant violated: element at heap index does not match element index");
            }
        }
    }

  private:
    void p_recache_key(std::size_t heap_index, std::size_t element_index) {
        if constexpr(Parameters::cache_key_in_heap) {
            m_heap_elements[heap_index].cached_key = m_get_key(element_index);
        }
    }

    void p_make_heap() {
        if(m_n <= 1)
            return;

        std::size_t max_with_children = (m_n - 1) / 4;
        for(std::size_t i = max_with_children; i < m_n; --i) {
            p_sift_down(i);
        }
    }

    void p_sift_down(std::size_t i) noexcept {
        std::size_t c3 = p_child3(i);
        while(c3 < m_n) {
            KeyType best_val = p_get_key(m_heap_elements[i]);
            std::size_t best_child = i;
            for(std::size_t j = 0; j < 4; ++j) {
                std::size_t c = c3 - j;
                if(KeyCompare{}(p_get_key(m_heap_elements[c]), best_val)) {
                    best_val = p_get_key(m_heap_elements[c]);
                    best_child = c;
                }
            }
            if(best_child == i) {
                return;
            }
            p_swap_heap_indices(i, best_child);
            i = best_child;
            c3 = p_child3(i);
        }
        if(c3 - 3 < m_n) {
            KeyType best_val = p_get_key(m_heap_elements[i]);
            std::size_t best_child = i;
            for(std::size_t c = c3 - 3; c < m_n; ++c) {
                if(KeyCompare{}(p_get_key(m_heap_elements[c]), best_val)) {
                    best_val = p_get_key(m_heap_elements[c]);
                    best_child = c;
                }
            }
            if(best_child != i) {
                p_swap_heap_indices(i, best_child);
            }
        }
    }

    void p_bubble_up(std::size_t i) noexcept {
        KeyType current_val = p_get_key(m_heap_elements[i]);
        while(i > 0) {
            std::size_t p = p_parent(i);
            if(KeyCompare{}(current_val, p_get_key(m_heap_elements[p]))) {
                p_swap_heap_indices(i, p);
                i = p;
            } else {
                return;
            }
        }
    }

    void p_swap_heap_indices(std::size_t i1, std::size_t i2) noexcept {
        std::size_t item1 = m_heap_elements[i1].index;
        std::size_t item2 = m_heap_elements[i2].index;
        m_element_index_to_heap_index[item1] = i2;
        m_element_index_to_heap_index[item2] = i1;
        std::swap(m_heap_elements[i1], m_heap_elements[i2]);
    }

    static constexpr std::size_t p_parent(std::size_t i) noexcept { return (i - 1) / 4; }

    static constexpr std::size_t p_child3(std::size_t i) noexcept { return 4 * i + 4; }

    struct HeapElementWithKey {
        std::size_t index;
        KeyType cached_key;
    };

    struct HeapElementWithoutKey {
        std::size_t index;
    };

    using HeapElementType =
        std::conditional_t<Parameters::cache_key_in_heap, HeapElementWithKey, HeapElementWithoutKey>;

    KeyType p_get_key(const HeapElementType &element) const {
        if constexpr(Parameters::cache_key_in_heap) {
            return element.cached_key;
        } else {
            return m_get_key(element.index);
        }
    }

    GetKeyFunctor m_get_key;
    std::vector<std::size_t> m_element_index_to_heap_index;
    std::vector<HeapElementType> m_heap_elements;
    std::size_t m_n;
};

template<typename GetKeyFunctor>
explicit MutableIndexPriorityQueue(GetKeyFunctor &&get_key, std::size_t n = 0)
    -> MutableIndexPriorityQueue<std::decay_t<GetKeyFunctor>, DefaultParameters>;

} // namespace mipq

#endif
