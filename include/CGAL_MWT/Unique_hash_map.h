#ifndef CGAL_MWT_UNIQUE_HASH_MAP_H_INCLUDED_
#define CGAL_MWT_UNIQUE_HASH_MAP_H_INCLUDED_

/**
 * The code in this file borrows code from
 * CGAL::Unique_hash_map and is based on the same
 * underlying implementation.
 * The only change is that we add a method 'lookup' to
 * do a lookup by key which can be used to avoid
 * double lookups in situations like this:
 *  - check if key is in map with is_defined
 *  - if it is not, return
 *  - otherwise, read/modify the value
 * In this map, we can use item = map.lookup(key), check
 * if item is nullptr (key not in map) and
 * otherwise read/modify the value using map[item].
 */

#include <CGAL/Unique_hash_map.h>
#include <CGAL/disable_warnings.h>

namespace mwt {

template<class Key_, class Data_, class UniqueHashFunction = CGAL::Handle_hash_function,
         class Allocator_ = CGAL_ALLOCATOR(Data_)>
class Unique_hash_map {
  public:
    typedef Key_ Key;
    typedef Data_ Data;
    typedef UniqueHashFunction Hash_function;
    typedef Allocator_ Allocator;

    // STL compliance
    typedef Key_ key_type;
    typedef Data_ data_type;
    typedef UniqueHashFunction hasher;

    typedef Unique_hash_map<Key, Data, Hash_function, Allocator> Self;

    typedef CGAL::internal::chained_map<Data, Allocator> Map;
    typedef typename Map::Item Item;

  private:
    Hash_function m_hash_function;
    Map m_map;

    template<class It, class Iterator_category> void reserve_impl(It, It, Iterator_category) {}

    template<class It> void reserve_impl(It b, It e, std::forward_iterator_tag) { m_map.reserve(std::distance(b, e)); }

  public:
    Unique_hash_map() = default;

    Unique_hash_map(const Data &deflt, std::size_t table_size = Map::default_size) : m_map(table_size, deflt) {}

    Unique_hash_map(const Data &deflt, std::size_t table_size, const Hash_function &fct)
        : m_hash_function(fct), m_map(table_size, deflt) {}

    Unique_hash_map(Key first1, Key beyond1, Data first2) { insert(first1, beyond1, first2); }
    Unique_hash_map(Key first1, Key beyond1, Data first2, const Data &deflt, std::size_t table_size = 1,
                    const Hash_function &fct = Hash_function())
        : m_hash_function(fct), m_map(table_size, deflt) {
        insert(first1, beyond1, first2);
    }

    void reserve(std::size_t n) { m_map.reserve(n); }

    Data default_value() const { return m_map.cxdef(); }

    Hash_function hash_function() const { return m_hash_function; }

    void clear() { m_map.clear(); }

    void clear(const Data &deflt) {
        m_map.clear();
        m_map.xdef() = deflt;
    }

    bool is_defined(const Key &key) const { return m_map.lookup(m_hash_function(key)) != 0; }

    Item lookup(const Key &key) const { return m_map.lookup(m_hash_function(key)); }

    const Data &operator[](const Key &key) const {
        Item p = m_map.lookup(m_hash_function(key));
        if(p != 0)
            return m_map.inf(p);
        return m_map.cxdef();
    }

    Data &operator[](const Key &key) { return m_map.access(m_hash_function(key)); }

    Data &operator[](Item item) noexcept {
        CGAL_assertion(item != nullptr);
        return m_map.inf(item);
    }

    Data insert(Key first1, Key beyond1, Data first2) {
        reserve_impl(first1, beyond1, typename std::iterator_traits<Key>::iterator_category());
        for(; first1 != beyond1; (++first1, ++first2)) {
            operator[](first1) = first2;
        }
        return first2;
    }
};

} // namespace mwt

#include <CGAL/enable_warnings.h>

#endif
