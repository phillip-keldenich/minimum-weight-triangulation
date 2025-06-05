#ifndef CGAL_MWT_RATIONAL_SPARSE_LU_H_INCLUDED_
#define CGAL_MWT_RATIONAL_SPARSE_LU_H_INCLUDED_

#include "CGAL_Rational_aux.h"
#include "Rational_or_int.h"
#include "Rational_sparse_gaussian.h"

namespace mwt {

template<typename NumberType_> class SymbolicCalcTable {
  public:
    using NumberType = NumberType_;

    struct SymbolicEntry {
        std::size_t symbol;
        NumberType value;
    };

    explicit SymbolicCalcTable(std::size_t num_symbols, std::size_t num_rows_cols)
        : m_num_symbols(num_symbols),
          m_output_index_to_index_map(num_rows_cols, std::numeric_limits<std::size_t>::max()),
          m_coefficient_table(num_symbols, NumberType(0)), m_is_nz(num_symbols, false), m_vector_entry_begins(1, 0) {}

    void clear() {
        m_symbolic_table.clear();
        m_vector_entry_begins.clear();
        m_vector_entry_begins.push_back(0);
    }

    template<typename SymbolRange> void read_symbolic_value(const SymbolRange &range) {
        for(const auto &entry : range) {
            m_coefficient_table[entry.symbol] = entry.value;
            if(!m_is_nz[entry.symbol]) {
                m_is_nz[entry.symbol] = true;
                m_current_nzs.push_back(entry.symbol);
            }
        }
    }

    void add_entry(std::optional<std::size_t> symbol) {
        if(!symbol)
            return;
        std::size_t s = *symbol;
        if(s >= m_coefficient_table.size())
            throw std::range_error("Symbol out of range");
        m_coefficient_table[s] = 1;
        if(!m_is_nz[s]) {
            m_is_nz[s] = true;
            m_current_nzs.push_back(s);
        }
    }

    void entry_from_row(std::size_t row_index) {
        std::size_t internal_index = m_output_index_to_index_map[row_index];
        std::size_t first_index = m_vector_entry_begins[internal_index];
        std::size_t last_index = m_vector_entry_begins[internal_index + 1];
        for(std::size_t i = first_index; i < last_index; ++i) {
            const auto &entry = m_symbolic_table[i];
            m_coefficient_table[entry.symbol] = entry.value;
            m_is_nz[entry.symbol] = true;
            m_current_nzs.push_back(entry.symbol);
        }
    }

    void scale_row(std::size_t output_index, const NumberType &factor) {
        std::size_t internal_index = m_output_index_to_index_map[output_index];
        std::size_t first_index = m_vector_entry_begins[internal_index];
        std::size_t last_index = m_vector_entry_begins[internal_index + 1];
        for(std::size_t i = first_index; i < last_index; ++i) {
            m_symbolic_table[i].value *= factor;
        }
    }

    void divide_coefficients(const NumberType &divisor) {
        for(std::size_t nz : m_current_nzs) {
            m_coefficient_table[nz] /= divisor;
        }
    }

    void add_scaled_to(std::vector<SymbolicEntry> &output, const NumberType &factor) {
        bool produced_zero = false;
        for(SymbolicEntry &e : output) {
            if(m_is_nz[e.symbol]) {
                m_is_nz[e.symbol] = false;
                e.value += factor * m_coefficient_table[e.symbol];
            }
        }
        for(std::size_t nz : m_current_nzs) {
            if(!m_is_nz[nz]) {
                m_is_nz[nz] = true;
            } else {
                if(m_coefficient_table[nz]) {
                    output.push_back(SymbolicEntry{nz, factor * m_coefficient_table[nz]});
                }
            }
        }
        auto new_end = std::remove_if(output.begin(), output.end(), [](const auto &val) { return !val.value; });
        output.erase(new_end, output.end());
    }

    void subtract_scaled_from(std::vector<SymbolicEntry> &output, const NumberType &factor) {
        for(SymbolicEntry &e : output) {
            if(m_is_nz[e.symbol]) {
                m_is_nz[e.symbol] = false;
                e.value -= factor * m_coefficient_table[e.symbol];
            }
        }
        for(std::size_t nz : m_current_nzs) {
            if(!m_is_nz[nz]) {
                m_is_nz[nz] = true;
            } else {
                if(m_coefficient_table[nz]) {
                    output.push_back(SymbolicEntry{nz, -factor * m_coefficient_table[nz]});
                }
            }
        }
        auto new_end = std::remove_if(output.begin(), output.end(), [](const auto &val) { return !val.value; });
        output.erase(new_end, output.end());
    }

    void permute_rows(const std::vector<std::size_t> &in_order, const std::vector<std::size_t> &out_order) {
        m_permute_tmp.resize(m_output_index_to_index_map.size());
        for(std::size_t i = 0, n = m_output_index_to_index_map.size(); i < n; ++i) {
            m_permute_tmp[out_order[i]] = m_output_index_to_index_map[in_order[i]];
        }
        m_permute_tmp.swap(m_output_index_to_index_map);
    }

    void subtract_scaled(std::size_t previous_output_index, const NumberType &factor) {
        std::size_t internal_index = m_output_index_to_index_map[previous_output_index];
        std::size_t first_index = m_vector_entry_begins[internal_index];
        std::size_t last_index = m_vector_entry_begins[internal_index + 1];
        for(std::size_t i = first_index; i < last_index; ++i) {
            const auto &entry = m_symbolic_table[i];
            if(!!(m_coefficient_table[entry.symbol] -= factor * entry.value)) {
                if(!m_is_nz[entry.symbol]) {
                    m_is_nz[entry.symbol] = true;
                    m_current_nzs.push_back(entry.symbol);
                }
            }
        }
    }

    void subtract_scaled(const std::vector<SymbolicEntry> &entries, const NumberType &factor) {
        for(const auto &entry : entries) {
            if(!!(m_coefficient_table[entry.symbol] -= factor * entry.value)) {
                if(!m_is_nz[entry.symbol]) {
                    m_is_nz[entry.symbol] = true;
                    m_current_nzs.push_back(entry.symbol);
                }
            }
        }
    }

    void divide_and_finalize(std::size_t output_index, const NumberType &divisor) {
        m_output_index_to_index_map[output_index] = m_vector_entry_begins.size() - 1;
        p_clear_table_to_vector_entry(m_symbolic_table, divisor);
        m_vector_entry_begins.push_back(m_symbolic_table.size());
    }

    void divide_and_finalize(std::vector<SymbolicEntry> &output, const NumberType &divisor) {
        output.clear();
        p_clear_table_to_vector_entry(output, divisor);
    }

    template<typename VectorType, typename Subtract = std::false_type>
    void scalar_product_into(const VectorType &nonsymbolic, const std::vector<std::vector<SymbolicEntry>> &symbolic,
                             std::vector<SymbolicEntry> &result, Subtract = {}) {
        for(const auto &entry : nonsymbolic) {
            for(const auto &sym : symbolic[entry.index]) {
                if constexpr(Subtract::value) {
                    m_coefficient_table[sym.symbol] -= entry.value * sym.value;
                } else {
                    m_coefficient_table[sym.symbol] += entry.value * sym.value;
                }
                if(!m_is_nz[sym.symbol]) {
                    m_is_nz[sym.symbol] = true;
                    m_current_nzs.push_back(sym.symbol);
                }
            }
        }
        result.clear();
        p_clear_table_to_vector_entry(result);
    }

    void print_nonzero_coefficients() const {
        for(std::size_t s : m_current_nzs) {
            std::cout << "s_" << s << ": " << m_coefficient_table[s] << "\n";
        }
    }

    void print_table() {
        for(std::size_t i = 0; i < m_output_index_to_index_map.size(); ++i) {
            std::size_t begin_table_index = m_output_index_to_index_map[i];
            if(begin_table_index >= m_vector_entry_begins.size()) {
                std::cout << "Row " << i << ": not present\n";
            } else {
                std::size_t first_index = m_vector_entry_begins[begin_table_index];
                std::size_t last_index = m_vector_entry_begins[begin_table_index + 1];
                std::cout << "Row " << i << ":";
                for(std::size_t j = first_index; j < last_index; ++j) {
                    const auto &entry = m_symbolic_table[j];
                    std::cout << " + " << entry.value << " * s_{" << entry.symbol << "}";
                }
                std::cout << "\n";
            }
        }
    }

    void output_coefficients_to_vector(std::vector<SymbolicEntry> &out, bool negate) {
        out.clear();
        if(negate) {
            p_clear_table_to_vector_entry_negate(out);
        } else {
            p_clear_table_to_vector_entry(out);
        }
    }

  private:
    void p_clear_table_to_vector_entry(std::vector<SymbolicEntry> &out) {
        for(std::size_t nzi : m_current_nzs) {
            auto &c = m_coefficient_table[nzi];
            if(!!c) {
                out.push_back(SymbolicEntry{nzi, std::move(c)});
            }
            c = NumberType(0);
            m_is_nz[nzi] = false;
        }
        m_current_nzs.clear();
    }

    void p_clear_table_to_vector_entry_negate(std::vector<SymbolicEntry> &out) {
        for(std::size_t nzi : m_current_nzs) {
            auto &c = m_coefficient_table[nzi];
            if(!!c) {
                c.negate();
                out.push_back(SymbolicEntry{nzi, std::move(c)});
            }
            c = NumberType(0);
            m_is_nz[nzi] = false;
        }
        m_current_nzs.clear();
    }

    void p_clear_table_to_vector_entry(std::vector<SymbolicEntry> &out, const NumberType &divisor) {
        for(std::size_t nzi : m_current_nzs) {
            auto &c = m_coefficient_table[nzi];
            if(!!c) {
                c /= divisor;
                out.push_back(SymbolicEntry{nzi, std::move(c)});
            }
            c = NumberType(0);
            m_is_nz[nzi] = false;
        }
        m_current_nzs.clear();
    }

    std::size_t m_num_symbols;
    std::vector<std::size_t> m_output_index_to_index_map;
    std::vector<std::size_t> m_permute_tmp;
    std::vector<NumberType> m_coefficient_table;
    std::vector<SymbolicEntry> m_symbolic_table;
    std::vector<std::size_t> m_current_nzs;
    std::vector<bool> m_is_nz;
    std::vector<std::size_t> m_vector_entry_begins;
};

template<typename NumberType, typename SymbolicEntryType>
NumberType substitute_symbols(const std::vector<SymbolicEntryType> &symbolic, const std::vector<NumberType> &values) {
    NumberType result(0);
    for(const auto &entry : symbolic) {
        result += entry.value * values[entry.symbol];
    }
    return result;
}

template<typename NumberType, typename SymbolicEntryType, typename CoefficientAndValueTransform>
auto substitute_symbols(const std::vector<SymbolicEntryType> &symbolic, const std::vector<NumberType> &values,
                        CoefficientAndValueTransform &&transform)
    -> std::decay_t<decltype(transform(symbolic[0].value, values[0]))> {
    using ResultType = std::decay_t<decltype(transform(symbolic[0].value, values[0]))>;
    ResultType result(0);
    for(const auto &entry : symbolic) {
        result += transform(entry.value, values[entry.symbol]);
    }
    return result;
}

/**
 * Class that manages storage of the
 * rows or columns of the L or U matrix
 * from our LU decomposition.
 */
template<typename NumberType_> class ContiguousFBSubstitutionStorage {
  public:
    using NumberType = NumberType_;
    using Entry = typename RationalSparseMatrix<NumberType>::Entry;
    using SymbolicTable = SymbolicCalcTable<NumberType>;
    using SymbolicEntry = typename SymbolicTable::SymbolicEntry;
    using Interval = CGAL::Interval_nt_advanced;

    /**
     * Create a ContiguousFBSubstitutionStorage instance.
     * @param rhs_order The elimination order of the rows or columns (whatever is associated with the RHS).
     * @param out_order The elimination order of the rows or columns (whatever is associated with the output).
     * @param reverse_order If true, the order is reversed, i.e., backward substitution is used.
     */
    ContiguousFBSubstitutionStorage(const std::vector<std::size_t> *rhs_order,
                                    const std::vector<std::size_t> *out_order, bool reverse_order)
        : m_rhs_order(rhs_order), m_out_order(out_order), m_reverse_order(reverse_order) {}

    /**
     * Set the matrix from the rows or columns given without transposing.
     */
    template<typename RowOrColCollectionType> void set_triangular_from_units(RowOrColCollectionType &&rows) {
        m_entries.clear();
        m_interval_entries.clear();
        m_unit_begins.clear();
        std::size_t current_begin = 0;
        auto handle_row = [&](std::size_t row_index, std::size_t col_index) {
            m_unit_begins.push_back(current_begin);
            const NumberType *pivot = nullptr;
            for(const auto &entry : rows[row_index]) {
                if(entry.index == col_index) {
                    pivot = &entry.value;
                } else {
                    if constexpr(!std::is_reference_v<RowOrColCollectionType>) {
                        m_entries.push_back(Entry{entry.index, std::move(entry.value)});
                    } else {
                        m_entries.push_back(Entry{entry.index, entry.value});
                    }
                }
            }
            assert(pivot);
            m_entries.push_back(Entry{row_index, *pivot});
            current_begin = m_entries.size();
        };
        const auto &rhso = *m_rhs_order;
        const auto &outo = *m_out_order;
        if(m_reverse_order) {
            for(std::size_t m = rows.size(), i = m - 1; i < m; --i) {
                handle_row(rhso[i], outo[i]);
            }
        } else {
            for(std::size_t i = 0, m = rows.size(); i < m; ++i) {
                handle_row(rhso[i], outo[i]);
            }
        }
        m_unit_begins.push_back(m_entries.size());
    }

    /**
     * Set the matrix from the columns given, transposing it in the process.
     */
    template<typename RowOrColCollectionType> void set_triangular_transpose(RowOrColCollectionType &&cols) {
        std::decay_t<RowOrColCollectionType> transpose(cols.size());
        for(std::size_t i = 0, m = cols.size(); i < m; ++i) {
            for(const auto &entry : cols[i]) {
                if constexpr(!std::is_reference_v<RowOrColCollectionType>) {
                    transpose[entry.index].push_back(Entry{i, std::move(entry.value)});
                } else {
                    transpose[entry.index].push_back(Entry{i, entry.value});
                }
            }
        }
        set_triangular_from_units(std::move(transpose));
    }

    void solve(const std::vector<NumberType> &rhs, std::vector<NumberType> &result) {
        const std::size_t m = m_rhs_order->size();
        result.resize(m);
        const auto &rhso = *m_rhs_order;
        const auto &outo = *m_out_order;
        if(m_reverse_order) {
            for(std::size_t i = m - 1; i < m; --i) {
                p_solve_iteration(m - 1 - i, rhs, result, rhso[i], outo[i]);
            }
        } else {
            for(std::size_t i = 0; i < m; ++i) {
                p_solve_iteration(i, rhs, result, rhso[i], outo[i]);
            }
        }
    }

    void solve_with_intervals(const std::vector<Interval> &rhs, std::vector<Interval> &result) {
        assert(CGAL::FPU_get_cw() == CGAL_FE_UPWARD);
        const std::size_t m = m_rhs_order->size();
        result.resize(m);
        if(m_interval_entries.empty()) {
            p_compute_interval_entries();
        }
        const auto &rhso = *m_rhs_order;
        const auto &outo = *m_out_order;
        if(m_reverse_order) {
            for(std::size_t i = m - 1; i < m; --i) {
                p_solve_intervals_iteration(m - 1 - i, rhs, result, rhso[i], outo[i]);
            }
        } else {
            for(std::size_t i = 0; i < m; ++i) {
                p_solve_intervals_iteration(i, rhs, result, rhso[i], outo[i]);
            }
        }
    }

    void solve_initial_with_symbols(SymbolicCalcTable<NumberType> &table,
                                    const std::vector<std::optional<std::size_t>> &rhs) {
        table.clear();
        const std::size_t m = m_rhs_order->size();
        const auto &rhso = *m_rhs_order;
        const auto &outo = *m_out_order;
        if(m_reverse_order) {
            for(std::size_t i = m - 1; i < m; --i) {
                p_solve_initial_iteration(m - 1 - i, rhs, table, rhso[i], outo[i]);
            }
        } else {
            for(std::size_t i = 0; i < m; ++i) {
                p_solve_initial_iteration(i, rhs, table, rhso[i], outo[i]);
            }
        }
    }

    void solve_final_with_symbols(SymbolicCalcTable<NumberType> &table,
                                  std::vector<std::vector<SymbolicEntry>> &output) {
        const std::size_t m = m_rhs_order->size();
        output.resize(m);
        const auto &rhso = *m_rhs_order;
        const auto &outo = *m_out_order;
        if(m_reverse_order) {
            for(std::size_t i = m - 1; i < m; --i) {
                p_solve_final_iteration(m - 1 - i, output, table, rhso[i], outo[i]);
            }
        } else {
            for(std::size_t i = 0; i < m; ++i) {
                p_solve_final_iteration(i, output, table, rhso[i], outo[i]);
            }
        }
    }

  private:
    void p_compute_interval_entries() {
        m_interval_entries.reserve(m_entries.size());
        std::transform(m_entries.begin(), m_entries.end(), std::back_inserter(m_interval_entries),
                       [](const auto &entry) { return mwt::exact_to_interval(entry.value); });
    }

    void p_solve_iteration(std::size_t unit_index, const std::vector<NumberType> &rhs, std::vector<NumberType> &result,
                           std::size_t rhs_index, std::size_t out_index) {
        auto &res = result[out_index];
        res = rhs[rhs_index];
        std::size_t unit_begin_i = m_unit_begins[unit_index];
        std::size_t unit_end_i = m_unit_begins[unit_index + 1] - 1;
        const Entry *unit_begin = m_entries.data() + unit_begin_i;
        const Entry *unit_end = m_entries.data() + unit_end_i;
        for(; unit_begin != unit_end; ++unit_begin) {
            res -= unit_begin->value * result[unit_begin->index];
        }
        res /= unit_end->value;
    }

    void p_solve_intervals_iteration(std::size_t unit_index, const std::vector<Interval> &rhs,
                                     std::vector<Interval> &result, std::size_t rhs_index, std::size_t out_index) {
        auto &res = result[out_index];
        res = rhs[rhs_index];
        std::size_t unit_begin_i = m_unit_begins[unit_index];
        std::size_t unit_end_i = m_unit_begins[unit_index + 1] - 1;
        const Entry *unit_begin = m_entries.data() + unit_begin_i;
        const Entry *unit_end = m_entries.data() + unit_end_i;
        const Interval *intervals_begin = m_interval_entries.data() + unit_begin_i;
        for(; unit_begin != unit_end; ++unit_begin, ++intervals_begin) {
            res -= *intervals_begin * result[unit_begin->index];
        }
        res /= *intervals_begin;
    }

    void p_solve_initial_iteration(std::size_t unit_index, const std::vector<std::optional<std::size_t>> &rhs,
                                   SymbolicCalcTable<NumberType> &table, std::size_t rhs_index, std::size_t out_index) {
        std::size_t unit_begin_i = m_unit_begins[unit_index];
        std::size_t unit_end_i = m_unit_begins[unit_index + 1] - 1;
        const Entry *unit_begin = m_entries.data() + unit_begin_i;
        const Entry *unit_end = m_entries.data() + unit_end_i;
        table.add_entry(rhs[rhs_index]);
        for(; unit_begin != unit_end; ++unit_begin) {
            table.subtract_scaled(unit_begin->index, unit_begin->value);
        }
        table.divide_and_finalize(out_index, unit_end->value);
    }

    void p_solve_final_iteration(std::size_t unit_index, std::vector<std::vector<SymbolicEntry>> &output,
                                 SymbolicCalcTable<NumberType> &table, std::size_t rhs_index, std::size_t out_index) {
        std::size_t unit_begin_i = m_unit_begins[unit_index];
        std::size_t unit_end_i = m_unit_begins[unit_index + 1] - 1;
        const Entry *unit_begin = m_entries.data() + unit_begin_i;
        const Entry *unit_end = m_entries.data() + unit_end_i;
        table.entry_from_row(rhs_index);
        for(; unit_begin != unit_end; ++unit_begin) {
            table.subtract_scaled(output[unit_begin->index], unit_begin->value);
        }
        table.divide_and_finalize(output[out_index], unit_end->value);
    }

    std::vector<std::size_t> m_unit_begins;
    std::vector<Entry> m_entries;
    std::vector<Interval> m_interval_entries;
    const std::vector<std::size_t> *m_rhs_order;
    const std::vector<std::size_t> *m_out_order;
    bool m_reverse_order;
};

template<typename NumberType_> class RationalSparseLU {
  public:
    using NumberType = NumberType_;
    using Interval = CGAL::Interval_nt_advanced;
    using Entry = typename RationalSparseMatrix<NumberType_>::Entry;
    using ContainerType = typename RationalSparseMatrix<NumberType_>::template ContainerType<Entry>;
    using IntervalContainerType = typename RationalSparseMatrix<NumberType_>::template ContainerType<Interval>;
    using IndexContainerType = typename RationalSparseMatrix<NumberType_>::template ContainerType<std::size_t>;

    struct ColumnMajor {};
    struct RowMajor {};

  private:
    void p_reset(bool is_new) {
        const std::size_t m = m_m;
        m_rows.resize(m);
        m_rows_with_column.resize(m);

        if(is_new) {
            m_row_elimination_order.reserve(m);
            m_col_elimination_order.reserve(m);
            m_diag_entries.reserve(m);
            m_col_score.reserve(m);
        }

        m_row_is_eliminated.assign(m, false);
        m_col_is_eliminated.assign(m, false);
        m_stamp_set.assign(m, std::uint8_t(0));
        m_stamp_value = 0;
        m_have_intervals = false;
        m_have_transposed_intervals = false;

        if(!is_new) {
            m_row_elimination_order.clear();
            m_col_elimination_order.clear();
            m_diag_entries.clear();
            m_l_entries.clear();
            m_col_score.clear();
            for(auto &x : m_rows_with_column) {
                x.clear();
            }
        }
    }

    explicit RationalSparseLU(std::size_t m)
        : m_m(m), m_row_queue(RowGetScore(this)),
          m_l_nontransposed(&m_row_elimination_order, &m_col_elimination_order, false),
          m_u_nontransposed(&m_row_elimination_order, &m_col_elimination_order, true),
          m_u_transposed(&m_col_elimination_order, &m_row_elimination_order, false),
          m_l_transposed(&m_col_elimination_order, &m_row_elimination_order, true) {
        p_reset(true);
    }

  public:
    template<typename ColumnRange>
    explicit RationalSparseLU(std::size_t m, ColumnMajor, const ColumnRange &col_range) : RationalSparseLU(m) {
        reset_columns(m, col_range, true);
    }

    template<typename ColumnRange>
    void reset_columns(std::size_t m, const ColumnRange &col_range, bool is_new = false) {
        if(!is_new) {
            p_reset(false);
        }
        std::size_t col_index = 0;
        for(const auto &col : col_range) {
            std::size_t col_size = 0;
            for(const auto &entry : col) {
                m_rows[entry.index].push_back(Entry{col_index, entry.value});
                m_rows_with_column[col_index].push_back(entry.index);
                ++col_size;
            }
            ++col_index;
            m_col_score.push_back(col_size);
        }
    }

    template<typename RowRange>
    explicit RationalSparseLU(std::size_t m, RowMajor, const RowRange &row_range) : RationalSparseLU(m) {
        m_col_score.assign(m, 0);
        std::size_t i = 0;
        for(const auto &row : row_range) {
            if(i == m) {
                ++i;
                break;
            }
            std::transform(row.begin(), row.end(), std::back_inserter(m_rows[i]), [&](const auto &e) {
                m_col_score[e.index] += 1;
                m_rows_with_column[e.index].push_back(i);
                return Entry{e.index, e.value};
            });
            assert(std::is_sorted(m_rows[i].begin(), m_rows[i].end()));
            ++i;
        }
        if(i != m) {
            throw std::logic_error("Row-major matrix has wrong number of rows");
        }
    }

    void compute_decomposition() {
        m_row_queue.reset_queue(m_m);
        for(std::size_t i = 0; i < m_m; ++i) {
            std::size_t row = p_pick_row();
            std::size_t col = p_pick_col(row);
            p_perform_iteration(row, col);
        }
        p_finalize_decomposition();
    }

    template<typename RHSNumberType>
    void solve_linear_system(const std::vector<RHSNumberType> &rhs, std::vector<RHSNumberType> &result,
                             std::vector<RHSNumberType> *tmp_storage = nullptr) {
        std::vector<RHSNumberType> tmp;
        if(!tmp_storage) {
            tmp_storage = &tmp;
        }
        m_l_nontransposed.solve(rhs, result);
        p_transform_with_d(result, *tmp_storage);
        m_u_nontransposed.solve(result, *tmp_storage);
        tmp_storage->swap(result);
    }

    using SymbolicTable = SymbolicCalcTable<NumberType>;
    using SymbolicEntry = typename SymbolicTable::SymbolicEntry;

    void solve_linear_system_with_symbols(SymbolicCalcTable<NumberType> &table,
                                          std::vector<std::vector<SymbolicEntry>> &output_values,
                                          const std::vector<std::optional<std::size_t>> &rhs) {
        m_l_nontransposed.solve_initial_with_symbols(table, rhs);
        p_transform_with_d(table);
        m_u_nontransposed.solve_final_with_symbols(table, output_values);
        // p_solve_with_l(table, rhs);
        // p_transform_with_d(table);
        // p_solve_with_u(table, output_values);
    }

    template<typename RHSNumberType>
    void solve_transposed_linear_system(const std::vector<RHSNumberType> &rhs, std::vector<RHSNumberType> &result,
                                        std::vector<RHSNumberType> *tmp_storage = nullptr) {
        std::vector<RHSNumberType> tmp;
        if(!tmp_storage) {
            tmp_storage = &tmp;
        }
        m_u_transposed.solve(rhs, result);
        p_transform_with_dtrans(result, *tmp_storage);
        m_l_transposed.solve(result, *tmp_storage);
        tmp_storage->swap(result);
        // p_ensure_have_transposed_decomposition();
        // p_solve_with_utrans(rhs, result);
        // p_transform_with_dtrans(result, *tmp_storage);
        // p_solve_with_ltrans(result, *tmp_storage);
    }

    void solve_transposed_with_intervals(const std::vector<Interval> &rhs, std::vector<Interval> &result,
                                         std::vector<Interval> *tmp_storage = nullptr) {
        std::vector<Interval> tmp;
        if(!tmp_storage) {
            tmp_storage = &tmp;
        }
        m_u_transposed.solve_with_intervals(rhs, result);
        p_transform_with_dtrans_intervals(result, *tmp_storage);
        m_l_transposed.solve_with_intervals(result, *tmp_storage);
        tmp_storage->swap(result);
        // p_ensure_have_transposed_intervals();
        // p_solve_with_utrans_intervals(rhs, result);
        // p_transform_with_dtrans_intervals(result, *tmp_storage);
        // p_solve_with_ltrans_intervals(result, *tmp_storage);
    }

    void solve_transposed_linear_system_with_symbols(SymbolicCalcTable<NumberType> &table,
                                                     std::vector<std::vector<SymbolicEntry>> &output_values,
                                                     const std::vector<std::optional<std::size_t>> &rhs) {
        m_u_transposed.solve_initial_with_symbols(table, rhs);
        p_transform_with_dtrans(table);
        m_l_transposed.solve_final_with_symbols(table, output_values);
        // p_ensure_have_transposed_decomposition();
        // p_solve_with_utrans(table, rhs);
        // p_transform_with_dtrans(table);
        // p_solve_with_ltrans(table, output_values);
    }

    void print_u() const {
        for(const auto &r : m_rows) {
            p_print_row(r);
        }
    }

    void print_l() const {
        for(const auto &r : m_l_rows) {
            p_print_row(r);
        }
    }

  private:
    // void p_ensure_have_transposed_intervals() {
    //     if(m_have_transposed_intervals) return;
    //     p_ensure_have_transposed_decomposition();
    //     p_transposed_intervals_from_decomposition();
    // }

    // void p_ensure_have_intervals() {
    //     if(m_have_intervals) return;
    //     p_intervals_from_decomposition();
    // }

    // void p_intervals_from_decomposition() {
    //     if(!m_have_transposed_intervals) {
    //         p_compute_diag_intervals();
    //     }
    //     m_l_rows_intervals.resize(m_m);
    //     m_rows_intervals.resize(m_m);
    //     for(std::size_t i = 0; i < m_m; ++i) {
    //         p_exact_to_interval_container(m_l_rows[i], m_l_rows_intervals[i]);
    //         p_exact_to_interval_container(m_rows[i], m_rows_intervals[i]);
    //     }
    //     m_have_intervals = true;
    // }

    // void p_transposed_intervals_from_decomposition() {
    //     if(!m_have_intervals) {
    //         p_compute_diag_intervals();
    //     }
    //     m_l_cols_intervals.resize(m_m);
    //     m_u_cols_intervals.resize(m_m);
    //     for(std::size_t i = 0; i < m_m; ++i) {
    //         p_exact_to_interval_container(m_l_cols[i], m_l_cols_intervals[i]);
    //         p_exact_to_interval_container(m_u_cols[i], m_u_cols_intervals[i]);
    //     }
    //     m_have_transposed_intervals = true;
    // }

    // template<typename ExactContainer, typename IntervalContainer>
    // void p_exact_to_interval_container(const ExactContainer& e, IntervalContainer& i) {
    //     i.clear();
    //     i.reserve(e.size());
    //     for(const auto& x : e) {
    //         i.emplace_back(mwt::to_interval(x.value));
    //     }
    // }

    void p_compute_diag_intervals() {
        m_diag_intervals.clear();
        m_diag_intervals.reserve(m_m);
        for(const auto &x : m_diag_entries) {
            m_diag_intervals.emplace_back(mwt::to_interval(x));
        }
    }

    template<typename RowType> void p_print_row(const RowType &row) const {
        auto iter = row.begin(), end = row.end();
        for(std::size_t i = 0; i < m_m; ++i) {
            if(iter != end && iter->index == i) {
                std::cout << iter->value << " \t";
                ++iter;
            } else {
                std::cout << "0 \t";
            }
        }
        std::cout << std::endl;
    }

    // template<typename RHSNumberType>
    // void p_solve_with_l(const std::vector<RHSNumberType>& rhs, std::vector<RHSNumberType>& result) {
    //     assert(m_m == rhs.size());
    //     result.resize(m_m);
    //     for(std::size_t index = 0; index < m_m; ++index) {
    //         std::size_t row_index = m_row_elimination_order[index];
    //         std::size_t col_index = m_col_elimination_order[index];
    //         auto& res = result[col_index];
    //         res = rhs[row_index];
    //         const NumberType* piv = nullptr;
    //         for(const auto& entry : m_l_rows[row_index]) {
    //             if(entry.index == col_index) {
    //                 piv = &entry.value;
    //             } else {
    //                 res -= entry.value * result[entry.index];
    //             }
    //         }
    //         res /= *piv;
    //     }
    // }

    // void p_solve_with_l(SymbolicCalcTable<NumberType>& table, const std::vector<std::optional<std::size_t>>& rhs) {
    //     table.clear();
    //     for(std::size_t index = 0; index < m_m; ++index) {
    //         std::size_t row_index = m_row_elimination_order[index];
    //         std::size_t col_index = m_col_elimination_order[index];
    //         std::size_t out_index = col_index;
    //         const NumberType* piv = nullptr;
    //         table.add_entry(rhs[row_index]);
    //         for(const auto& entry : m_l_rows[row_index]) {
    //             if(entry.index == col_index) {
    //                 piv = &entry.value;
    //             } else {
    //                 table.subtract_scaled(entry.index, entry.value);
    //             }
    //         }
    //         table.divide_and_finalize(out_index, *piv);
    //     }
    // }

    // void p_solve_with_utrans(SymbolicCalcTable<NumberType>& table, const std::vector<std::optional<std::size_t>>&
    // rhs) {
    //     table.clear();
    //     for(std::size_t index = 0; index < m_m; ++index) {
    //         std::size_t row_index = m_row_elimination_order[index];
    //         std::size_t col_index = m_col_elimination_order[index];
    //         std::size_t out_index = row_index;
    //         const NumberType* piv = nullptr;
    //         table.add_entry(rhs[col_index]);
    //         for(const auto& entry : m_u_cols[col_index]) {
    //             if(entry.index == row_index) {
    //                 piv = &entry.value;
    //             } else {
    //                 table.subtract_scaled(entry.index, entry.value);
    //             }
    //         }
    //         table.divide_and_finalize(out_index, *piv);
    //     }
    // }

    // template<typename RHSNumberType>
    // void p_solve_with_utrans(const std::vector<RHSNumberType>& rhs, std::vector<RHSNumberType>& result) {
    //     assert(m_m == rhs.size());
    //     result.resize(m_m);
    //     for(std::size_t index = 0; index < m_m; ++index) {
    //         std::size_t row_index = m_row_elimination_order[index];
    //         std::size_t col_index = m_col_elimination_order[index];
    //         auto& res = result[row_index];
    //         res = rhs[col_index];
    //         const NumberType* piv = nullptr;
    //         for(const auto& entry : m_u_cols[col_index]) {
    //             if(entry.index == row_index) {
    //                 piv = &entry.value;
    //             } else {
    //                 res -= entry.value * result[entry.index];
    //             }
    //         }
    //         res /= *piv;
    //     }
    // }

    // void p_solve_with_utrans_intervals(const std::vector<Interval>& rhs, std::vector<Interval>& result) {
    //     assert(m_m == rhs.size());
    //     result.resize(m_m);
    //     for(std::size_t index = 0; index < m_m; ++index) {
    //         std::size_t row_index = m_row_elimination_order[index];
    //         std::size_t col_index = m_col_elimination_order[index];
    //         auto& res = result[row_index];
    //         res = rhs[col_index];
    //         Interval pivot;
    //         const auto& exact = m_u_cols[col_index];
    //         const auto& intervals = m_u_cols_intervals[col_index];
    //         for(std::size_t i = 0, k = exact.size(); i < k; ++i) {
    //             Interval value = intervals[i];
    //             std::size_t eindex = exact[i].index;
    //             if(eindex == row_index) {
    //                 pivot = value;
    //             } else {
    //                 res -= value * result[eindex];
    //             }
    //         }
    //         res /= pivot;
    //     }
    // }

    template<typename RHSNumberType>
    void p_transform_with_dtrans(std::vector<RHSNumberType> &result, std::vector<RHSNumberType> &tmp) {
        tmp.resize(m_m);
        for(std::size_t i = 0; i < m_m; ++i) {
            std::size_t r = m_row_elimination_order[i];
            std::size_t c = m_col_elimination_order[i];
            const auto &x = m_diag_entries[i];
            tmp[c] = std::move(result[r]);
            tmp[c] *= x;
        }
        result.swap(tmp);
    }

    void p_transform_with_dtrans_intervals(std::vector<Interval> &result, std::vector<Interval> &tmp) {
        tmp.resize(m_m);
        for(std::size_t i = 0; i < m_m; ++i) {
            std::size_t r = m_row_elimination_order[i];
            std::size_t c = m_col_elimination_order[i];
            Interval x = m_diag_intervals[i];
            tmp[c] = result[r] * x;
        }
        result.swap(tmp);
    }

    void p_transform_with_dtrans(SymbolicCalcTable<NumberType> &table) {
        for(std::size_t i = 0; i < m_m; ++i) {
            std::size_t r = m_row_elimination_order[i];
            const auto &x = m_diag_entries[i];
            table.scale_row(r, x);
        }
        table.permute_rows(m_row_elimination_order, m_col_elimination_order);
    }

    void p_transform_with_d(SymbolicCalcTable<NumberType> &table) {
        for(std::size_t i = 0; i < m_m; ++i) {
            std::size_t c = m_col_elimination_order[i];
            const auto &x = m_diag_entries[i];
            table.scale_row(c, x);
        }
        table.permute_rows(m_row_elimination_order, m_col_elimination_order);
    }

    template<typename RHSNumberType>
    void p_transform_with_d(std::vector<RHSNumberType> &result, std::vector<RHSNumberType> &tmp) {
        tmp.resize(m_m);
        for(std::size_t i = 0; i < m_m; ++i) {
            std::size_t r = m_row_elimination_order[i];
            std::size_t c = m_col_elimination_order[i];
            const auto &x = m_diag_entries[i];
            tmp[r] = std::move(result[c]);
            tmp[r] *= x;
        }
        result.swap(tmp);
    }

    // template<typename RHSNumberType>
    // void p_solve_with_u(std::vector<RHSNumberType>& result, std::vector<RHSNumberType>& tmp) {
    //     for(std::size_t index = m_m - 1; index < m_m; --index) {
    //         std::size_t row_index = m_row_elimination_order[index];
    //         std::size_t col_index = m_col_elimination_order[index];

    //         auto& res = tmp[col_index];
    //         res = result[row_index];
    //         const NumberType* piv = nullptr;
    //         for(const auto& entry : m_rows[row_index]) {
    //             if(entry.index == col_index) {
    //                 piv = &entry.value;
    //             } else {
    //                 res -= entry.value * tmp[entry.index];
    //             }
    //         }
    //         res /= *piv;
    //     }
    //     result = std::move(tmp);
    // }

    // void p_solve_with_u(SymbolicCalcTable<NumberType>& table, std::vector<std::vector<SymbolicEntry>>& output)
    // {
    //     output.resize(m_m);
    //     for(std::size_t index = m_m - 1; index < m_m; --index) {
    //         std::size_t row_index = m_row_elimination_order[index];
    //         std::size_t col_index = m_col_elimination_order[index];
    //         std::size_t out_index = col_index;
    //         const NumberType* piv = nullptr;
    //         table.entry_from_row(row_index);
    //         for(const auto& entry : m_rows[row_index]) {
    //             if(entry.index == col_index) {
    //                 piv = &entry.value;
    //             } else {
    //                 table.subtract_scaled(output[entry.index], entry.value);
    //             }
    //         }
    //         table.divide_and_finalize(output[out_index], *piv);
    //     }
    // }

    // void p_solve_with_ltrans(SymbolicCalcTable<NumberType>& table, std::vector<std::vector<SymbolicEntry>>& output) {
    //     output.resize(m_m);
    //     for(std::size_t index = m_m - 1; index < m_m; --index) {
    //         std::size_t row_index = m_row_elimination_order[index];
    //         std::size_t col_index = m_col_elimination_order[index];
    //         std::size_t out_index = row_index;
    //         const NumberType* piv = nullptr;
    //         table.entry_from_row(col_index);
    //         for(const auto& entry : m_l_cols[col_index]) {
    //             if(entry.index == row_index) {
    //                 piv = &entry.value;
    //             } else {
    //                 table.subtract_scaled(output[entry.index], entry.value);
    //             }
    //         }
    //         table.divide_and_finalize(output[out_index], *piv);
    //     }
    // }

    // template<typename RHSNumberType>
    // void p_solve_with_ltrans(std::vector<RHSNumberType>& result, std::vector<RHSNumberType>& tmp) {
    //     for(std::size_t index = m_m - 1; index < m_m; --index) {
    //         std::size_t row_index = m_row_elimination_order[index];
    //         std::size_t col_index = m_col_elimination_order[index];

    //         auto& res = tmp[row_index];
    //         res = result[col_index];
    //         const NumberType* piv = nullptr;
    //         for(const auto& entry : m_l_cols[col_index]) {
    //             if(entry.index == row_index) {
    //                 piv = &entry.value;
    //             } else {
    //                 res -= entry.value * tmp[entry.index];
    //             }
    //         }
    //         res /= *piv;
    //     }
    //     result.swap(tmp);
    // }

    // void p_solve_with_ltrans_intervals(std::vector<Interval>& result, std::vector<Interval>& tmp) {
    //     for(std::size_t index = m_m - 1; index < m_m; --index) {
    //         std::size_t row_index = m_row_elimination_order[index];
    //         std::size_t col_index = m_col_elimination_order[index];

    //         auto& res = tmp[row_index];
    //         res = result[col_index];
    //         Interval pivot;
    //         const auto& exact = m_l_cols[col_index];
    //         const auto& intervals = m_l_cols_intervals[col_index];
    //         for(std::size_t i = 0, k = exact.size(); i < k; ++i) {
    //             std::size_t eindex = exact[i].index;
    //             Interval value = intervals[i];
    //             if(eindex == row_index) {
    //                 pivot = value;
    //             } else {
    //                 res -= value * tmp[eindex];
    //             }
    //         }
    //         res /= pivot;
    //     }
    //     result.swap(tmp);
    // }

    // void p_ensure_have_transposed_decomposition() {
    //     if(m_l_cols.empty() || m_l_cols[0].empty()) {
    //         p_compute_transposed_decomposition();
    //     }
    // }

    // void p_compute_transposed_decomposition() {
    //     m_l_cols.resize(m_m);
    //     m_u_cols.resize(m_m);
    //     for(std::size_t r = 0; r < m_m; ++r) {
    //         for(const auto& entry : m_l_rows[r]) {
    //             m_l_cols[entry.index].push_back(Entry{r, entry.value});
    //         }
    //         for(const auto& entry : m_rows[r]) {
    //             m_u_cols[entry.index].push_back(Entry{r, entry.value});
    //         }
    //     }
    // }

    void p_finalize_decomposition() {
        // turn L into a matrix with rows
        m_l_rows.resize(m_m);
        std::sort(m_l_entries.begin(), m_l_entries.end(),
                  [](const auto &a, const auto &b) { return a.row < b.row || (a.row == b.row && a.col < b.col); });
        for(auto &l_entry : m_l_entries) {
            m_l_rows[l_entry.row].push_back(Entry{l_entry.col, std::move(l_entry.value)});
        }
        m_l_entries.clear();
        m_l_transposed.set_triangular_transpose(m_l_rows);
        m_l_nontransposed.set_triangular_from_units(std::move(m_l_rows));
        m_u_transposed.set_triangular_transpose(m_rows);
        m_u_nontransposed.set_triangular_from_units(std::move(m_rows));
        p_compute_diag_intervals();
        for(auto &x : m_rows) {
            x.clear();
        }
        for(auto &x : m_l_rows) {
            x.clear();
        }
    }

    struct Triple {
        Triple(std::size_t row, Entry &&entry) noexcept : row(row), col(entry.index), value(std::move(entry.value)) {}

        std::size_t row;
        std::size_t col;
        NumberType value;
    };

    void p_perform_iteration(std::size_t row_index, std::size_t col_index) {
        const NumberType &pivot_x = m_diag_entries.back();
        const auto &pivot_row = m_rows[row_index];
        m_l_entries.emplace_back(row_index, Entry{col_index, pivot_x});
        p_mark_row_eliminated(row_index);
        m_col_is_eliminated[col_index] = true;
        p_bump_stamp();
        for(std::size_t orow : m_rows_with_column[col_index]) {
            if(m_row_is_eliminated[orow])
                continue;
            if(!p_stamp_row(orow))
                continue;
            std::size_t orow_index = p_find_in_row(orow, col_index);
            if(orow_index > m_m)
                continue;
            NumberType pivot_y = m_rows[orow][orow_index].value;
            pivot_y.negate();
            pivot_y /= pivot_x;
            p_row_elimination(row_index, col_index, orow, pivot_y);
        }
    }

    void p_row_elimination(std::size_t pivot_row, std::size_t pivot_col, std::size_t target_row,
                           const NumberType &factor) {
        const auto &source = m_rows[pivot_row];
        auto &target = m_rows[target_row];
        auto source_it = source.begin(), source_end = source.end();
        auto target_it = target.begin(), target_end = target.end();
        m_tmp_row.clear();

        auto advance_source = [&](std::size_t sindex) -> bool {
            m_tmp_row.push_back(Entry{sindex, source_it->value * factor});
            m_rows_with_column[sindex].push_back(target_row);
            ++m_col_score[sindex];
            return ++source_it == source_end;
        };

        auto advance_target = [&](std::size_t tindex) -> bool {
            m_tmp_row.emplace_back(std::move(*target_it));
            return ++target_it == target_end;
        };

        auto advance_both = [&](std::size_t index) -> bool {
            if(index == pivot_col) {
                m_l_entries.emplace_back(target_row, std::move(*target_it));
            } else {
                NumberType &target_coeff = target_it->value;
                target_coeff += source_it->value * factor;
                if(!target_coeff) {
                    --m_col_score[index];
                } else {
                    m_tmp_row.emplace_back(std::move(*target_it));
                }
            }
            ++source_it, ++target_it;
            return source_it == source_end || target_it == target_end;
        };

        for(;;) {
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
        if(target_it != target_end) {
            while(!advance_target(target_it->index)) {
            }
        } else if(source_it != source_end) {
            while(!advance_source(source_it->index)) {
            }
        }
        p_row_from_tmp(target_row);
    }

    void p_row_from_tmp(std::size_t target_row) {
        auto &target = m_rows[target_row];
        target.assign(std::move_iterator(m_tmp_row.begin()), std::move_iterator(m_tmp_row.end()));
        m_tmp_row.clear();
        if(!m_row_queue.is_in_queue(target_row)) {
            throw std::logic_error("Row not in queue");
        }
        m_row_queue.notify_key_changed(target_row);
    }

    bool p_stamp_row(std::size_t row) {
        auto &val = m_stamp_set[row];
        if(val == m_stamp_value) {
            return false;
        }
        val = m_stamp_value;
        return true;
    }

    void p_mark_row_eliminated(std::size_t row_index) {
        m_row_is_eliminated[row_index] = true;
        for(const auto &entry : m_rows[row_index]) {
            m_col_score[entry.index] -= 1;
        }
    }

    void p_bump_stamp() {
        if(++m_stamp_value == 0) {
            ++m_stamp_value;
            std::fill(m_stamp_set.begin(), m_stamp_set.end(), std::uint8_t(0));
        }
    }

    std::size_t p_find_in_row(std::size_t row_index, std::size_t col_index) {
        const auto &row = m_rows[row_index];
        auto pos = std::partition_point(row.begin(), row.end(), [&](const Entry &e) { return e.index < col_index; });
        if(pos == row.end() || pos->index != col_index) {
            return std::numeric_limits<std::size_t>::max();
        }
        return std::size_t(pos - row.begin());
    }

    struct RowGetScore {
        explicit RowGetScore(const RationalSparseLU *that) noexcept : that(that) {}

        std::size_t operator()(std::size_t row) const noexcept { return that->m_rows[row].size(); }

        const RationalSparseLU *that;
    };

    std::size_t p_pick_row() {
        while(!m_row_queue.empty()) {
            std::size_t row = m_row_queue.pop();
            if(!m_row_is_eliminated[row]) {
                return row;
            }
        }
        throw std::logic_error("Basis matrix is singular!");
    }

    std::size_t p_pick_col(std::size_t row_index) {
        std::size_t best_score = std::numeric_limits<std::size_t>::max();
        std::size_t best_col = std::numeric_limits<std::size_t>::max();
        const Entry *best_entry = nullptr;
        bool best_is_integer = false;
        bool best_is_unity = false;
        const auto &row = m_rows[row_index];
        for(const auto &entry : row) {
            std::size_t col_index = entry.index;
            std::size_t col_score = m_col_score[col_index];
            if(col_score == 0)
                throw std::logic_error("Basis matrix is singular!");
            if(col_score > best_score)
                continue;
            bool is_integer = entry.value.is_platform_int();
            bool is_unity = (is_integer && (entry.value == 1 || entry.value == -1));
            if(col_score < best_score || is_unity || (is_integer && !best_is_integer)) {
                best_score = col_score;
                best_col = col_index;
                best_is_integer = is_integer;
                best_is_unity = is_unity;
                best_entry = &entry;
            }
        }
        m_diag_entries.push_back(best_entry->value);
        m_row_elimination_order.push_back(row_index);
        m_col_elimination_order.push_back(best_col);
        return best_col;
    }

    std::size_t m_m;
    std::vector<std::size_t> m_row_elimination_order;
    std::vector<std::size_t> m_col_elimination_order;
    std::vector<IndexContainerType> m_rows_with_column;
    std::vector<Triple> m_l_entries;
    std::vector<NumberType> m_diag_entries;
    std::vector<Interval> m_diag_intervals;
    mipq::MutableIndexPriorityQueue<RowGetScore> m_row_queue;
    std::vector<ContainerType> m_rows;
    std::vector<ContainerType> m_l_rows;

    ContiguousFBSubstitutionStorage<NumberType> m_l_nontransposed;
    ContiguousFBSubstitutionStorage<NumberType> m_u_nontransposed;
    ContiguousFBSubstitutionStorage<NumberType> m_u_transposed;
    ContiguousFBSubstitutionStorage<NumberType> m_l_transposed;

    // std::vector<IntervalContainerType> m_rows_intervals;

    // std::vector<IntervalContainerType> m_l_rows_intervals;
    // std::vector<ContainerType> m_l_cols;
    // std::vector<IntervalContainerType> m_l_cols_intervals;
    // std::vector<ContainerType> m_u_cols;
    // std::vector<IntervalContainerType> m_u_cols_intervals;
    ContainerType m_tmp_row;
    std::vector<bool> m_row_is_eliminated;
    std::vector<bool> m_col_is_eliminated;
    std::vector<std::size_t> m_col_score;
    std::vector<std::uint8_t> m_stamp_set;
    std::uint8_t m_stamp_value{0};
    bool m_have_intervals{false};
    bool m_have_transposed_intervals{false};
};

} // namespace mwt

#endif
