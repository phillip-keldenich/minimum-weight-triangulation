#ifndef CGAL_MWT_RATIONAL_SPARSE_GAUSSIAN_H_INCLUDED_
#define CGAL_MWT_RATIONAL_SPARSE_GAUSSIAN_H_INCLUDED_

#include "Mutable_index_priority_queue.h"
#include <boost/container/small_vector.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/iterator/transform_iterator.hpp>
#include <boost/range/iterator_range_core.hpp>
#include <vector>

namespace mwt {

template<typename NumberType_> class RationalSparseMatrix {
  public:
    using NumberType = NumberType_;

    struct Entry {
        std::size_t index;
        NumberType value;

        bool operator<(const Entry &other) const noexcept { return index < other.index; }

        bool operator==(const Entry &other) const noexcept { return index == other.index; }
    };

    template<typename VT> using ContainerType = boost::container::small_vector<VT, 16>;

    explicit RationalSparseMatrix(std::size_t num_columns)
        : m_num_cols(num_columns), m_rows_with_column(num_columns), m_column_scores(num_columns, 0),
          m_column_score_queue(GetColumnScoreKey{m_column_scores.data()}), m_col_is_eliminated(num_columns, false) {}

    void start_new_row() { m_tmp_row.clear(); }

    template<typename EntryNumberType> void new_row_add_entry(std::size_t index, EntryNumberType &&num_type) {
        m_tmp_row.push_back(Entry{index, std::forward<EntryNumberType>(num_type)});
    }

    template<typename EntryNumberType> void finish_row(EntryNumberType &&rhs) {
        p_new_row_from_temp();
        m_rhs.emplace_back(std::forward<EntryNumberType>(rhs));
    }

    template<typename RangeType, typename EntryNumberType>
    void new_row(const RangeType &entries, EntryNumberType &&rhs) {
        for(const auto &e : entries) {
            m_tmp_row.push_back(e);
        }
        p_new_row_from_temp();
        m_rhs.emplace_back(std::forward<EntryNumberType>(rhs));
    }

    /**
     * The number of rows.
     */
    std::size_t num_rows() const noexcept { return m_rows.size(); }

    /**
     * The number of columns.
     */
    std::size_t num_cols() const noexcept { return m_num_cols; }

    /**
     * Run the gaussian elimination; returns the rank of the matrix.
     */
    std::size_t eliminate() {
        m_current_rank = m_rows.size();
        p_recompute_rows_with_column_and_column_scores();
        p_resize_stamp_set();
        m_row_is_eliminated.assign(m_rows.size(), false);
        m_col_is_eliminated.assign(m_num_cols, false);
        m_row_elimination_order.clear();
        m_col_elimination_order.clear();
        for(;;) {
            std::size_t next_column = p_pick_column();
            if(next_column >= m_num_cols) {
                break;
            }
            auto [next_row, index_in_row] = p_pick_row(next_column);
            if(next_row >= m_rows.size()) {
                break;
            }
            p_perform_iteration(next_column, next_row, index_in_row);
        }
        return m_current_rank;
    }

    const ContainerType<Entry> &get_row(std::size_t index) const { return m_rows[index]; }

    ContainerType<Entry> &get_row(std::size_t index) { return m_rows[index]; }

    const NumberType &get_rhs(std::size_t index) const { return m_rhs[index]; }

    NumberType &get_rhs(std::size_t index) { return m_rhs[index]; }

    const std::vector<std::size_t> &row_elimination_order() const { return m_row_elimination_order; }

    const std::vector<std::size_t> &col_elimination_order() const { return m_col_elimination_order; }

    void remove_zero_rows() {
        std::vector<std::size_t> old_index_to_new;
        std::size_t count = 0;
        bool have_dropped = false;
        old_index_to_new.reserve(m_rows.size());
        for(std::size_t ri = 0, nr = m_rows.size(); ri < nr; ++ri) {
            auto &row = m_rows[ri];
            if(!row.empty()) {
                old_index_to_new.push_back(count);
                if(have_dropped) {
                    m_rows[count] = std::move(row);
                    m_rhs[count] = std::move(m_rhs[ri]);
                    m_row_is_integer[count] = m_row_is_integer[ri];
                }
                ++count;
            } else {
                old_index_to_new.push_back(std::numeric_limits<std::size_t>::max());
                have_dropped = true;
            }
        }
        m_rows.resize(count);
        m_rhs.resize(count);
        auto translate_to_new = [&](std::size_t old_row_index) { return old_index_to_new[old_row_index]; };
        std::transform(m_row_elimination_order.begin(), m_row_elimination_order.end(), m_row_elimination_order.begin(),
                       translate_to_new);
        for(auto &col : m_rows_with_column) {
            col.clear();
        }
        m_row_is_eliminated.assign(m_rows.size(), true);
    }

  private:
    void p_perform_iteration(std::size_t col, std::size_t row, std::size_t index_in_row) {
        m_row_is_eliminated[row] = true;
        m_col_is_eliminated[col] = true;
        m_row_elimination_order.push_back(row);
        m_col_elimination_order.push_back(col);
        m_column_score_updates.assign(m_rows[row].size(), 0);
        for(NonZeroRow r : m_nonzero_row_buffer) {
            if(r.row_index == row)
                continue;
            p_row_elimination(col, row, index_in_row, r.row_index, r.index_in_row);
        }

        // remove the eliminated row from the column score;
        // update the column score and the queue
        std::size_t update_index = 0;
        for(const auto &e : m_rows[row]) {
            std::size_t col_index = e.index;
            auto diff = m_column_score_updates[update_index++];
            diff -= 1;
            m_column_scores[col_index] += diff;
            /*if(diff > 0) {
                m_column_score_queue.notify_key_increased(col_index);
            } else if(diff < 0) {
                m_column_score_queue.notify_key_decreased(col_index);
            }*/
        }
    }

    /**
     * Add an appropriate multiple of the source row to the target row
     * to make the entry in column col of the target row 0.
     */
    void p_row_elimination(std::size_t col, std::size_t source_row, std::size_t col_in_source, std::size_t target_row,
                           std::size_t col_in_target) {
        const auto &source = m_rows[source_row];
        auto &target = m_rows[target_row];
        NumberType coeff = -target[col_in_target].value;
        coeff /= source[col_in_source].value;
        p_merge_add_to_temp_row(target_row, source, target, coeff);
        m_rhs[target_row] += m_rhs[source_row] * std::move(coeff);
        if(m_tmp_row.empty()) {
            // the resulting row has only zeros.
            m_row_is_eliminated[target_row] = true;
            if(m_rhs[target_row] == 0) {
                --m_current_rank;
            } else {
                throw std::logic_error("System of linear equations is infeasible!");
            }
        }
        p_row_from_temp(target_row);
    }

    void p_row_from_temp(std::size_t target_row) {
        auto &target = m_rows[target_row];
        target.assign(std::move_iterator(m_tmp_row.begin()), std::move_iterator(m_tmp_row.end()));
        m_tmp_row.clear();
    }

    void p_new_row_from_temp() {
        std::size_t last_index = std::numeric_limits<std::size_t>::max();
        std::sort(m_tmp_row.begin(), m_tmp_row.end());
        m_rows.emplace_back();
        auto &back_row = m_rows.back();
        for(auto &entry : m_tmp_row) {
            if(entry.index == last_index) {
                back_row.back().value += entry.value;
            } else {
                back_row.push_back(std::move(entry));
            }
            last_index = entry.index;
        }
        bool is_int =
            std::all_of(back_row.begin(), back_row.end(), [](const Entry &e) { return e.value.is_platform_int(); });
        m_tmp_row.clear();
        m_row_is_integer.push_back(is_int);
    }

    void p_merge_add_to_temp_row(std::size_t target_row, const ContainerType<Entry> &source,
                                 ContainerType<Entry> &target, const NumberType &coeff) {
        m_tmp_row.clear();
        auto source_it = source.begin();
        auto target_it = target.begin();
        auto source_end = source.end();
        auto target_end = target.end();
        bool new_is_integer = true;
        std::size_t col_score_update_index = 0;

        /**
         * Advance the pivot row iterator.
         */
        auto advance_source = [&](std::size_t sindex) -> bool {
            NumberType new_coeff = source_it->value * coeff;
            if(!m_col_is_eliminated[sindex]) {
                m_rows_with_column[sindex].push_back(target_row);
                m_column_score_updates[col_score_update_index] += 1;
            }
            new_is_integer &= new_coeff.is_platform_int();
            m_tmp_row.push_back(Entry{sindex, std::move(new_coeff)});
            ++col_score_update_index;
            return ++source_it == source_end;
        };

        /**
         * Advance the modified row iterator.
         */
        auto advance_target = [&](std::size_t tindex) -> bool {
            new_is_integer &= target_it->value.is_platform_int();
            m_tmp_row.push_back(std::move(*target_it));
            return ++target_it == target_end;
        };

        /*
         * Advance both iterators (same column index),
         * i.e., non-zero in pivot row hits non-zero in modified row.
         */
        auto advance_both = [&](std::size_t index) -> bool {
            NumberType &target_coeff = target_it->value;
            target_coeff += source_it->value * coeff;
            if(target_coeff == 0) {
                m_column_score_updates[col_score_update_index] -= 1;
            } else {
                new_is_integer &= target_coeff.is_platform_int();
                m_tmp_row.push_back(std::move(*target_it));
            }
            ++source_it, ++target_it, ++col_score_update_index;
            return source_it == source_end || target_it == target_end;
        };

        for(;;) {
            // similar idea to mergesort merge: the index sequences are ordered
            std::size_t sindex = source_it->index;
            std::size_t tindex = target_it->index;
            if(sindex < tindex) {
                if(advance_source(sindex))
                    break;
            } else if(sindex > tindex) {
                if(advance_target(tindex))
                    break;
            } else {
                if(advance_both(sindex))
                    break;
            }
        }
        // drain one potentially non-empty sequence;
        // at this point, the other sequence must be empty
        if(target_it != target_end) {
            while(!advance_target(target_it->index)) {
            }
        } else if(source_it != source_end) {
            while(!advance_source(source_it->index)) {
            }
        }
        m_row_is_integer[target_row] = new_is_integer;
    }

    /**
     * Select a column to eliminate.
     */
    std::size_t p_pick_column() {
        while(!m_column_score_queue.empty()) {
            std::size_t col = m_column_score_queue.pop();
            if(!m_col_is_eliminated[col] && m_column_scores[col] != 0) {
                return col;
            }
        }
        return std::numeric_limits<std::size_t>::max();

        // trivial queue-free version.
        // becomes a massive performance bottleneck for large matrices,
        // essentially taking up 99% of the runtime.
        /*std::size_t min_score_col = std::numeric_limits<std::size_t>::max();
        std::size_t min_score = min_score_col;
        for(std::size_t i = 0; i < m_num_cols; ++i) {
            // do not attempt to eliminate a column more than once
            if(m_col_is_eliminated[i]) continue;
            std::size_t score = m_column_scores[i];
            // do not attempt to eliminate using a column without
            // non-eliminated rows to choose from
            if(score == 0) continue;
            if(score < min_score) {
                min_score = score;
                min_score_col = i;
            }
        }
        return min_score_col;*/
    }

    /**
     * Pick a row as pivot row when eliminating column.
     * Fills m_nonzero_row_buffer with the rows containing
     * non-zeros in column and also includes the position
     * of the non-zero in the respective row.
     */
    std::pair<std::size_t, std::size_t> p_pick_row(std::size_t column) {
        m_nonzero_row_buffer.clear();
        std::size_t best_row_score = std::numeric_limits<std::size_t>::max();
        std::size_t index_in_best = 0;
        std::size_t best_row_index = std::numeric_limits<std::size_t>::max();
        bool best_is_integer = false;

        p_next_stamp();
        const auto &with_col = m_rows_with_column[column];
        for(std::size_t candidate_index : with_col) {
            if(!p_stamp_if_unstamped(candidate_index))
                continue;
            if(m_row_is_eliminated[candidate_index])
                continue;
            auto i_in_r = p_find_in_row(candidate_index, column);
            if(!i_in_r)
                continue;
            m_nonzero_row_buffer.push_back(NonZeroRow{candidate_index, *i_in_r});
            std::size_t score = m_rows[candidate_index].size();
            if(score < best_row_score ||
               (score == best_row_score && !best_is_integer && m_row_is_integer[candidate_index])) {
                best_is_integer = true;
                best_row_score = score;
                index_in_best = *i_in_r;
                best_row_index = candidate_index;
            }
        }
        return {best_row_index, index_in_best};
    }

    void p_recompute_rows_with_column() {
        for(auto &rwc : m_rows_with_column) {
            rwc.clear();
        }
        for(std::size_t i = 0, m = num_rows(); i < m; ++i) {
            for(const auto &e : m_rows[i]) {
                m_rows_with_column[e.index].push_back(i);
            }
        }
    }

    struct GetColumnScoreKey {
        explicit GetColumnScoreKey(const std::size_t *column_scores) : column_scores(column_scores) {}

        std::size_t operator()(std::size_t index) const noexcept {
            std::size_t raw_score = column_scores[index];
            return raw_score == 0 ? std::numeric_limits<std::size_t>::max() : raw_score;
        }

        const std::size_t *column_scores;
    };

    void p_recompute_rows_with_column_and_column_scores() {
        p_recompute_rows_with_column();
        m_column_scores.clear();
        for(const auto &rwc : m_rows_with_column) {
            m_column_scores.push_back(rwc.size());
        }
        m_column_score_queue.reset_queue(m_num_cols, GetColumnScoreKey{m_column_scores.data()});
    }

    std::optional<std::size_t> p_find_in_row(std::size_t row, std::size_t column) const {
        const auto &row_vec = m_rows[row];
        auto begin = row_vec.begin();
        auto end = row_vec.end();
        auto pos = std::partition_point(begin, end, [&](const Entry &e) { return e.index < column; });
        if(pos == end || pos->index != column)
            return std::nullopt;
        return std::size_t(pos - begin);
    }

    void p_next_stamp() {
        if(++m_stamp_value == 0) {
            ++m_stamp_value;
            std::fill(m_index_stamp.begin(), m_index_stamp.end(), std::uint16_t(0));
        }
    }

    void p_resize_stamp_set() { m_index_stamp.resize((std::max)(m_num_cols, m_rows.size()), std::uint16_t(0)); }

    bool p_is_stamped(std::size_t index) const noexcept { return m_index_stamp[index] == m_stamp_value; }

    bool p_stamp_if_unstamped(std::size_t index) noexcept {
        if(m_index_stamp[index] == m_stamp_value) {
            return false;
        }
        m_index_stamp[index] = m_stamp_value;
        return true;
    }

    // the total number of columns
    std::size_t m_num_cols;

    // the current rows of the matrix
    std::vector<ContainerType<Entry>> m_rows;

    // a temporary buffer for a row; used both during
    // row insertion and when doing row eliminations
    ContainerType<Entry> m_tmp_row;

    // the right-hand side values of the equalities
    std::vector<NumberType> m_rhs;

    // flags tracking whether a row is definitely integer;
    // it may happen that an all-integer row is not marked here
    std::vector<bool> m_row_is_integer;

    // a vector of vectors containing an over-approximation
    // of the rows that have a non-zero value for the given column
    // may contain more rows/duplicate rows than actually have a non-zero
    std::vector<ContainerType<std::size_t>> m_rows_with_column;

    // the row indices of the rows that were used to
    // eliminate, in the order they were used
    std::vector<std::size_t> m_row_elimination_order;
    // whether a row has been eliminated; rows that became
    // zero are marked as eliminated but do not appear in m_row_elimination_order
    std::vector<bool> m_row_is_eliminated;

    // the column indices of the columns that were used to
    // eliminate, in the order they were used
    std::vector<std::size_t> m_col_elimination_order;
    // whether a column has been eliminated
    std::vector<bool> m_col_is_eliminated;

    // column score: the number of non-eliminated rows
    // which have a non-zero value for the column
    std::vector<std::size_t> m_column_scores;
    mipq::MutableIndexPriorityQueue<GetColumnScoreKey> m_column_score_queue;
    std::vector<std::make_signed_t<std::size_t>> m_column_score_updates;

    // a 'stamp set' for temporarily marking columns/rows
    std::uint16_t m_stamp_value{0};
    std::vector<std::uint16_t> m_index_stamp;

    // buffer of nonzero rows for current pivot column
    // with the index of the column in the row
    struct NonZeroRow {
        std::size_t row_index;
        std::size_t index_in_row;
    };
    std::vector<NonZeroRow> m_nonzero_row_buffer;

    // current matrix rank
    std::size_t m_current_rank;
};

} // namespace mwt

#endif
