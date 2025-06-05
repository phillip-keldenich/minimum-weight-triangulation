#ifndef CGAL_MWT_CROSSING_FREE_H_INCLUDED_
#define CGAL_MWT_CROSSING_FREE_H_INCLUDED_

#include <iostream>
#include <random>
#include <vector>

namespace mwt {

/**
 * Check a set of edges for intersections in the interior;
 * endpoints shared by two edges do not count,
 * but an endpoint in the interior of another edge does.
 *
 * CGAL::do_curves_intersect requires exact constructions
 * to produce meaningful results for this task and requires way too much
 * memory and time for large point sets (often more than computing the MWT).
 * That is unfortunate (since the guarantees required by
 * do_curves_intersect could, for segments at least, be handled by code that
 * does _not_ require exact constructions).
 */
template<typename Kernel_> class CrossingFreenessVerifier {
  public:
    using Kernel = Kernel_;
    using FT = typename Kernel::FT;
    using Point_2 = typename Kernel::Point_2;
    using Segment_2 = typename Kernel::Segment_2;

    explicit CrossingFreenessVerifier(std::vector<Segment_2> segments)
        : m_rng(std::random_device{}()), m_segments(std::move(segments)),
          m_segment_locations(m_segments.size(), nullptr) {
        p_fix_segments();
        p_setup_events();
        p_skiplist_init();
    }

    ~CrossingFreenessVerifier() noexcept { p_skiplist_destroy(); }

    bool operator()() {
        if(m_events.empty()) {
            return m_intersection.has_value();
        }
        for(const Event &event : m_events) {
            if(event.entering) {
                if(p_handle_entering(event.segment_index))
                    break;
            } else {
                if(p_handle_leaving(event.segment_index))
                    break;
            }
        }
        m_events.clear();
        return m_intersection.has_value();
    }

    const std::vector<Segment_2> &get_segments() const noexcept { return m_segments; }

    std::pair<std::size_t, std::size_t> get_intersecting_segment_indices() const noexcept { return *m_intersection; }

    std::pair<Segment_2, Segment_2> get_intersecting_segments() const noexcept {
        return {m_segments[m_intersection->first], m_segments[m_intersection->second]};
    }

  private:
    /**
     * Event data.
     */
    struct Event {
        std::size_t segment_index;
        Point_2 point;
        bool entering;

        bool operator<(const Event &other) const noexcept {
            auto cres = CGAL::compare_xy(point, other.point);
            if(cres == CGAL::SMALLER)
                return true;
            if(cres == CGAL::LARGER)
                return false;
            return !entering && other.entering;
        }
    };

    /**
     * Node in the skiplist datastructure (used for sweep line).
     */
    struct SkipNode {
        std::size_t segment_index; //< index of the segment
        SkipNode *prev, *next;     //< prev and next pointer at this level
        std::int32_t level;        //< the level at which this node is
        std::int32_t height;       //< the total height of this node, i.e., its maximum level
    };

    /**
     * Fix segments so that the source point is always
     * the xy-lexicographically smaller of the two points.
     */
    void p_fix_segments() {
        for(auto &seg : m_segments) {
            if(seg.source() > seg.target()) {
                seg = Segment_2(seg.target(), seg.source());
            }
        }
    }

    /**
     * Setup the event queue (sorted vector of events).
     */
    void p_setup_events() {
        m_events.clear();
        m_events.reserve(2 * m_segments.size());
        for(std::size_t i = 0, nsegs = m_segments.size(); i < nsegs; ++i) {
            const Segment_2 &seg = m_segments[i];
            m_events.push_back({i, seg.source(), true});
            m_events.push_back({i, seg.target(), false});
        }
        std::sort(m_events.begin(), m_events.end());
    }

    /**
     * Handle an entering segment.
     */
    bool p_handle_entering(std::size_t segment_index) {
        SkipNode *new_node = p_skiplist_insert(segment_index);
        if(!new_node)
            return true;
        SkipNode *prev = new_node->prev, *next = new_node->next;
        if(prev->prev) {
            if(p_check_segment_against(segment_index, prev))
                return true;
        }
        if(next->next) {
            if(p_check_segment_against(segment_index, next))
                return true;
        }
        return false;
    }

    void p_print_skiplist() {
        for(SkipNode *node = m_current_head; node; node = node->next) {
            if(node->segment_index >= m_segments.size()) {
                std::cout << "NIL ";
            } else {
                const auto &seg = m_segments[node->segment_index];
                std::cout << "[" << node->segment_index << " (" << seg.source() << " -> " << seg.target() << ")] ";
            }
        }
        std::cout << std::endl;
    }

    /**
     * Handle a leaving segment.
     */
    bool p_handle_leaving(std::size_t segment_index) {
        SkipNode *rem_node = m_segment_locations[segment_index];
        SkipNode *prev = rem_node->prev, *next = rem_node->next;
        if(prev->prev && next->next) {
            if(p_check_segment_against(prev->segment_index, next))
                return true;
        }
        p_skiplist_erase(rem_node);
        return false;
    }

    /**
     * Check two segments for intersection.
     */
    bool p_check_segment_against(std::size_t segment_index, SkipNode *adjacent_node) {
        std::size_t other_index = adjacent_node->segment_index;
        Segment_2 s1 = m_segments[segment_index];
        Segment_2 s2 = m_segments[other_index];
        if(!CGAL::do_intersect(s1, s2))
            return false;
        if(s1.source() == s2.source() || s1.target() == s2.target() || s1.source() == s2.target() ||
           s1.target() == s2.source())
            return false;
        m_intersection.emplace(segment_index, other_index);
        return true;
    }

    /**
     * Search through the current level of the skip list until
     * callable returns false. The callable can also return
     * optional<bool> instead of bool to indicate situations in
     * which the callable may have to order segments that cannot be ordered,
     * e.g., a vertical segment intersecting another segment,
     * which always means an intersection; the offending segment will
     * populate m_intersection->first.
     * The function returns true iff it finds an intersection.
     */
    template<typename Callable>
    bool p_search_at_level(Callable &&callable, SkipNode *&current, SkipNode *&next, SkipNode *&next_next) {
        for(;;) {
            next_next = next->next;
            if(!next_next)
                return false;
            auto res = std::invoke(callable, next->segment_index);
            if constexpr(std::is_same_v<std::decay_t<decltype(res)>, bool>) {
                if(!res)
                    return false;
            } else {
                static_assert(std::is_same_v<std::decay_t<decltype(res)>, std::optional<bool>>);
                if(!res.has_value()) {
                    m_intersection.emplace(next->segment_index, std::numeric_limits<std::size_t>::max());
                    return true;
                }
                if(!*res)
                    return false;
            }
            current = next;
            next = next_next;
        }
    }

    /**
     * Destructor for skip nodes for use with std::unique_ptr.
     */
    struct SkipnodeDestructor {
        void operator()(SkipNode *node) const noexcept { that->p_skipnode_free(node); }
        CrossingFreenessVerifier *that;
    };

    /**
     * Free a skiplist node that has already been
     * unlinked from the skiplist. Small skiplist
     * nodes are placed into a freelist, while
     * large nodes are directly deallocated.
     */
    void p_skipnode_free(SkipNode *node_to_free) noexcept {
        auto h = node_to_free->height;
        if(h > 4) {
            delete[] node_to_free;
            return;
        }
        --h;
        node_to_free->next = m_small_free_lists[h];
        m_small_free_lists[h] = node_to_free;
    }

    /**
     * Allocate a skiplist node. Small skiplist nodes
     * are taken from a freelist. If the freelist is empty,
     * a growing block of nodes is allocated and put onto it.
     * Large skiplist nodes are allocated directly.
     * Initializes height, levels, and segment index.
     */
    SkipNode *p_skipnode_alloc(std::size_t segment_index, std::int32_t height) {
        SkipNode *result;
        if(height <= 4) {
            result = p_skipnode_alloc_small(height - 1);
        } else {
            result = new SkipNode[height];
        }
        p_skipnode_init_levels(result, segment_index, height);
        return result;
    }

    /**
     * Allocate a growing number of skip list nodes in a single
     * allocation and put them into the free list for the given height,
     * then return the first node.
     */
    SkipNode *p_skipnode_alloc_small_slow_path(std::int32_t flist_index) {
        std::size_t alloc_size = m_small_free_last_alloc[flist_index] * 2;
        const std::size_t height = flist_index + 1;
        std::unique_ptr<SkipNode[]> alloc_holder = std::make_unique<SkipNode[]>(alloc_size * height);
        m_small_allocs.push_back(alloc_holder.get());
        SkipNode *alloc = alloc_holder.release();
        for(std::size_t i = 1; i < alloc_size - 1; ++i) {
            alloc[i * height].next = alloc + (i + 1) * height;
        }
        alloc[(alloc_size - 1) * height].next = nullptr;
        m_small_free_lists[flist_index] = alloc + height;
        m_small_free_last_alloc[flist_index] = alloc_size;
        return alloc;
    }

    /**
     * Skipnode allocation for small heights.
     */
    SkipNode *p_skipnode_alloc_small(std::int32_t flist_index) {
        if(!m_small_free_lists[flist_index]) {
            return p_skipnode_alloc_small_slow_path(flist_index);
        }
        SkipNode *result = m_small_free_lists[flist_index];
        m_small_free_lists[flist_index] = result->next;
        return result;
    }

    /**
     * Initialize the levels of the given skiplist node.
     */
    void p_skipnode_init_levels(SkipNode *node, std::size_t segment_index, std::int32_t height) noexcept {
        for(std::int32_t i = 0; i < height; ++i) {
            node[i].segment_index = segment_index;
            node[i].height = height;
            node[i].level = i;
            node[i].prev = nullptr;
            node[i].next = nullptr;
        }
    }

    /**
     * Destroy and reset the skiplist and free
     * all memory allocated for the skiplist.
     */
    void p_skiplist_destroy() noexcept {
        SkipNode *current = m_current_head;
        while(current) {
            SkipNode *next = current->next;
            if(current->height > 4) {
                delete[] current;
            }
            current = next;
        }
        m_current_head = m_current_tail = nullptr;
        for(SkipNode *bulk : m_small_allocs) {
            delete[] bulk;
        }
        m_small_allocs.clear();
        std::fill_n(m_small_free_lists, 4, nullptr);
        m_small_free_last_alloc[0] = 2048;
        m_small_free_last_alloc[1] = 1024;
        m_small_free_last_alloc[2] = 512;
        m_small_free_last_alloc[3] = 256;
    }

    /**
     * Draw a single random bit (fair coin toss).
     */
    bool p_draw_bit() noexcept {
        if(m_rng_bits == 0) {
            m_rng_buf = m_rng();
            m_rng_bits = 64;
        }
        bool bit = m_rng_buf & 1;
        m_rng_buf >>= 1;
        --m_rng_bits;
        return bit;
    }

    /**
     * Draw bits until a 1 is drawn and return the number of bits drawn.
     */
    std::int32_t p_draw_bits_until_one() noexcept {
        std::int32_t bits_drawn = 0;
#if defined(__GNUC__) && __GNUC__ >= 4
        while(BOOST_UNLIKELY(!m_rng_buf)) {
            bits_drawn += m_rng_bits;
            m_rng_buf = m_rng();
            m_rng_bits = 64;
        }
        std::int32_t cdraw = __builtin_ctzll(m_rng_buf) + 1;
        bits_drawn += cdraw;
        m_rng_buf >>= cdraw;
        m_rng_bits -= cdraw;
#else
        while(!p_draw_bit()) {
            ++bits_drawn;
        }
        ++bits_drawn;
#endif
        constexpr std::size_t MAX_HEIGHT = sizeof(std::size_t) * CHAR_BIT;
        return BOOST_UNLIKELY(bits_drawn > MAX_HEIGHT) ? MAX_HEIGHT : bits_drawn;
    }

    /**
     * Allocate a pair of nodes.
     * Handles exceptions raised after an allocation.
     */
    std::pair<SkipNode *, SkipNode *> p_alloc_nodepair(std::size_t segment_index, std::int32_t height) {
        std::unique_ptr<SkipNode[], SkipnodeDestructor> n1{p_skipnode_alloc(segment_index, height), {this}};
        std::unique_ptr<SkipNode[], SkipnodeDestructor> n2{p_skipnode_alloc(segment_index, height), {this}};
        return {n1.release(), n2.release()};
    }

    /**
     * Initialize the skiplist datastructure.
     */
    void p_skiplist_init() {
        constexpr std::size_t nil_segment = std::numeric_limits<std::size_t>::max();
        auto [head, tail] = p_alloc_nodepair(nil_segment, 1);
        m_current_head = head;
        m_current_tail = tail;
        head->next = tail;
        tail->prev = head;
        m_current_max_level = 1;
    }

    /**
     * Create a skiplist node for the given segment index.
     * The node is created with a height chosen by coin tosses
     * until a 1 is drawn.
     */
    SkipNode *p_skipnode_create(std::size_t segment_index) {
        std::int32_t height = p_draw_bits_until_one();
        if(height > m_current_max_level) {
            height = std::min(height, m_current_max_level + 1);
            p_skiplist_realloc_head_and_tail(height);
            m_current_max_level = height;
        }
        SkipNode *node = p_skipnode_alloc(segment_index, height);
        return node;
    }

    /**
     * Reallocate the head and tail nodes whenever we create
     * a new node with a height greater than the current maximum level.
     */
    void p_skiplist_realloc_head_and_tail(std::int32_t new_height) {
        constexpr std::size_t nil_segment = std::numeric_limits<std::size_t>::max();
        auto [new_head, new_tail] = p_alloc_nodepair(nil_segment, new_height);
        new_head[new_height - 1].next = &new_tail[new_height - 1];
        new_tail[new_height - 1].prev = &new_head[new_height - 1];
        for(unsigned l = 0; l < new_height - 1; ++l) {
            SkipNode *head = &m_current_head[l];
            SkipNode *tail = &m_current_tail[l];
            SkipNode *nhead = &new_head[l];
            SkipNode *ntail = &new_tail[l];
            SkipNode *next = head->next;
            next->prev = nhead;
            nhead->next = next;
            SkipNode *prev = tail->prev;
            prev->next = ntail;
            ntail->prev = prev;
        }
        p_skipnode_free(m_current_head);
        p_skipnode_free(m_current_tail);
        m_current_head = new_head;
        m_current_tail = new_tail;
    }

    void p_set_inserting_segment(std::size_t segment_index) {
        m_inserting_segment_index = segment_index;
        m_inserting_segment = m_segments[segment_index];
        m_inserting_vertical = (m_inserting_segment.source().x() == m_inserting_segment.target().x());
    }

    static std::optional<bool> p_invert(std::optional<bool> x) noexcept {
        if(x)
            return !*x;
        return std::nullopt;
    }

    /**
     * Check if the first, vertical segment goes after the other segment.
     */
    std::optional<bool> p_first_vertical_goes_after(const Segment_2 &vertical, const Segment_2 &other) const {
        Point_2 vertical_lower = vertical.source();
        Point_2 other_source = other.source();
        if(vertical_lower.x() == other_source.x()) {
            if(vertical_lower.y() >= other_source.y()) {
                return true;
            }
            return std::nullopt;
        }
        auto o1 = CGAL::compare_y_at_x(vertical_lower, other);
        if(o1 == CGAL::EQUAL) {
            // lower point of vertical segment is inside the other segment
            return std::nullopt;
        }
        if(o1 == CGAL::LARGER) {
            // lower point of vertical segment lies above other segment
            return true;
        }
        auto o2 = CGAL::compare_y_at_x(vertical.target(), other);
        if(o2 == CGAL::SMALLER || vertical.target() == other.target()) {
            // upper point of vertical segment lies below other segment
            return false;
        }
        return std::nullopt;
    }

    /**
     * Check if the first of two non-vertical segments,
     * which is the segment we are just inserting,
     * goes after the other segment.
     */
    std::optional<bool> p_nonvertical_first_goes_after(const Segment_2 &seg1, const Segment_2 &seg2) const {
        const Point_2 &src1 = seg1.source();
        const Point_2 &src2 = seg2.source();
        if(src1.x() == src2.x()) {
            if(src1.y() > src2.y()) {
                return true;
            }
            // same source point - order by incline
            auto order = CGAL::orientation(src1, seg1.target(), seg2.target());
            if(order == CGAL::COLLINEAR) {
                return std::nullopt;
            }
            return order == CGAL::RIGHT_TURN;
        }
        auto comp = CGAL::compare_y_at_x(src1, seg2);
        if(comp == CGAL::EQUAL) {
            // src1 is in seg2's interior
            return std::nullopt;
        }
        return comp == CGAL::LARGER;
    }

    /**
     * Check if the inserting segment goes after the given other segment.
     */
    std::optional<bool> p_inserting_goes_after(std::size_t other_segment_index) {
        const Segment_2 &other_segment = m_segments[other_segment_index];
        bool other_vertical = (other_segment.source().x() == other_segment.target().x());
        if(m_inserting_vertical && other_vertical) {
            return std::nullopt;
        }
        if(m_inserting_vertical) {
            return p_first_vertical_goes_after(m_inserting_segment, other_segment);
        }
        if(other_vertical) {
            return p_invert(p_first_vertical_goes_after(other_segment, m_inserting_segment));
        }
        return p_nonvertical_first_goes_after(m_inserting_segment, other_segment);
    }

    /**
     * Insert a skiplist node. If an intersection is found,
     * the segments are recorded and nullptr is returned.
     * Otherwise, the inserted node is returned.
     */
    SkipNode *p_skiplist_insert(std::size_t segment_index) {
        std::unique_ptr<SkipNode[], SkipnodeDestructor> node{p_skipnode_create(segment_index), {this}};
        std::int32_t l = m_current_head->height - 1;
        const std::int32_t nl = node[0].height - 1;
        SkipNode *current = &m_current_head[l];
        SkipNode *next = current->next;
        SkipNode *next_next;
        p_set_inserting_segment(segment_index);

        /**
         * Operator to determine if the new node's segment goes after the given other segment.
         */
        auto new_goes_after = [this](std::size_t other_segment_index) -> std::optional<bool> {
            return this->p_inserting_goes_after(other_segment_index);
        };

        /**
         * Search above new node level.
         */
        while(l > nl) {
            if(p_search_at_level(new_goes_after, current, next, next_next)) {
                m_intersection->second = segment_index;
                return nullptr;
            }
            current -= 1;
            --l;
            next = current->next;
        }

        /**
         * Search below new node level;
         * record pointers to allow insertion.
         */
        for(;;) {
            if(p_search_at_level(new_goes_after, current, next, next_next)) {
                m_intersection->second = segment_index;
                return nullptr;
            }
            node[l].prev = current;
            node[l].next = next;
            if(l-- == 0)
                break;
            current -= 1;
            next = current->next;
        }

        /**
         * Change pointers of existing nodes to insert new node.
         */
        for(std::int32_t i = 0; i <= nl; ++i) {
            SkipNode *levelptr = &node[i];
            levelptr->prev->next = levelptr;
            levelptr->next->prev = levelptr;
        }

        /**
         * Record inserted position for later erasure.
         */
        m_segment_locations[segment_index] = node.get();
        return node.release();
    }

    /**
     * Erase the given node from the skiplist.
     */
    void p_skiplist_erase(SkipNode *node) {
        for(std::int32_t l = 0, h = node->height; l < h; ++l) {
            SkipNode *layer = &node[l];
            SkipNode *prev = layer->prev, *next = layer->next;
            prev->next = next;
            next->prev = prev;
        }
        p_skipnode_free(node);
    }

    // to avoid some memory allocations for the most common skip node heights
    SkipNode *m_small_free_lists[4] = //< pointers to free lists for heights 1, 2, 3, 4
        {nullptr, nullptr, nullptr, nullptr};
    std::size_t m_small_free_last_alloc[4] = //< last allocation multiple for heights 1, 2, 3, 4
        {2048, 1024, 512, 256};
    std::vector<SkipNode *> m_small_allocs; //< mass allocations for small heights to allow deletion
    SkipNode *m_current_head = nullptr;     //< the head of the skip list
    SkipNode *m_current_tail = nullptr;     //< the tail of the skip list
    std::mt19937_64 m_rng;
    std::uint64_t m_rng_buf = 0;
    std::size_t m_rng_bits = 0;
    std::int32_t m_current_max_level = 0;
    std::optional<std::pair<std::size_t, std::size_t>> m_intersection = std::nullopt;
    std::vector<Segment_2> m_segments;
    std::vector<SkipNode *> m_segment_locations;
    std::size_t m_inserting_segment_index = std::numeric_limits<std::size_t>::max();
    Segment_2 m_inserting_segment;
    bool m_inserting_vertical = false;
    std::vector<Event> m_events;
};

} // namespace mwt

#endif
