#ifndef CGAL_MWT_RATIONAL_SPARSE_MATRIX_H_INCLUDED_
#define CGAL_MWT_RATIONAL_SPARSE_MATRIX_H_INCLUDED_

#include "CGAL_Rational_aux.h"
#include "LP_backends.h"
#include "Rational_or_int.h"
#include "Rational_sparse_LU.h"
#include "Rational_sparse_gaussian.h"
#include <CGAL/FPU.h>

namespace mwt {

template<typename NumberType_> class RationalSparseLPWithSymbolicObjective {
  public:
    using NumberType = NumberType_;

    template<typename VT> using ContainerType = boost::container::small_vector<VT, 16>;

    using Entry = typename RationalSparseMatrix<NumberType>::Entry;

    using ColumnType = ContainerType<Entry>;

    using ConstrIndex = std::size_t;
    using VarIndex = std::size_t;
    using BoundType = std::optional<NumberType>;
    using SymTable = SymbolicCalcTable<NumberType>;
    using SymEntry = typename SymTable::SymbolicEntry;
    using Interval = CGAL::Interval_nt_advanced;

    explicit RationalSparseLPWithSymbolicObjective(bool is_minimization) : m_is_minimization(is_minimization) {}

    // min c^Tx
    // Ax - Iw = 0
    // a <= w <= b, l <= x <= u
    // we know the non-basic values;
    // if we remove the related columns from (A|-I), we get a new RHS v, and a 'basis matrix'
    // where Bx_B = v must hold -> allows to get the values!
    // objective: (c_N - N^TB^{-T}c_B)x_N

    template<typename EntryRange>
    ConstrIndex add_constraint(BoundType lhs, BoundType rhs, const EntryRange &row_entries) {
        ConstrIndex index = m_num_constraints++;
        m_is_equality.push_back(lhs && rhs && *lhs == *rhs);
        m_lhs.emplace_back(std::move(lhs));
        m_rhs.emplace_back(std::move(rhs));
        for(const auto &entry : row_entries) {
            m_original_matrix[entry.index].emplace_back(Entry{index, entry.value});
        }
        return index;
    }

    VarIndex add_variable(BoundType lb, BoundType ub) {
        ++m_original_vars;
        m_var_lb.emplace_back(std::move(lb));
        m_var_ub.emplace_back(std::move(ub));
        m_original_matrix.emplace_back();
        return m_original_matrix.size() - 1;
    }

    VarIndex add_variables(std::size_t num_variables, const BoundType &lb, const BoundType &ub) {
        VarIndex result = m_original_vars;
        m_original_vars += num_variables;
        m_var_lb.insert(m_var_lb.end(), num_variables, lb);
        m_var_ub.insert(m_var_ub.end(), num_variables, ub);
        m_original_matrix.resize(m_original_vars);
        return result;
    }

    template<typename EntryRange> VarIndex add_variable(BoundType lb, BoundType ub, const EntryRange &col_entries) {
        VarIndex index = add_variable(std::move(lb), std::move(ub));
        m_original_matrix.back().assign(std::begin(col_entries), std::end(col_entries));
        std::sort(m_original_matrix.back().begin(), m_original_matrix.back().end());
        return index;
    }

    ConstrIndex add_equalities(RationalSparseMatrix<NumberType> &&matrix) {
        ConstrIndex result = m_num_constraints;
        for(std::size_t i = 0, m = matrix.num_rows(); i < m; ++i) {
            auto &row = matrix.get_row(i);
            auto &rhs = matrix.get_rhs(i);
            m_lhs.emplace_back(rhs);
            m_rhs.emplace_back(std::move(rhs));
            m_is_equality.push_back(true);
            for(auto &entry : row) {
                m_original_matrix[entry.index].emplace_back(Entry{m_num_constraints, std::move(entry.value)});
            }
            ++m_num_constraints;
        }
        return result;
    }

    std::size_t num_variables() const noexcept { return m_original_vars; }

    std::size_t num_constraints() const noexcept { return m_num_constraints; }

    template<typename LPBackend, typename ObjectiveCoefficientType> struct ImpreciseModel {
        using IModel = typename LPBackend::Model;
        using IVar = typename LPBackend::Var;
        using IConstr = typename LPBackend::Constr;
        using Interval = CGAL::Interval_nt_advanced;

        explicit ImpreciseModel(std::vector<ObjectiveCoefficientType> objective, bool verbose = false)
            : model{std::make_unique<IModel>(LPBackend::create_model(!verbose))},
              objective_exact(std::move(objective)) {
            p_exact_objective_to_approx();
        }

        explicit ImpreciseModel(std::vector<ObjectiveCoefficientType> objective, std::vector<Interval> objective_approx,
                                bool verbose = false)
            : model{std::make_unique<IModel>(LPBackend::create_model(!verbose))}, objective_exact(std::move(objective)),
              objective_approx(std::move(objective_approx)) {
            if(this->objective_approx.size() != this->objective_exact.size()) {
                if(this->objective_approx.empty()) {
                    p_exact_objective_to_approx();
                } else {
                    throw std::invalid_argument("Approximate objective must have the same number of coefficients as "
                                                "the objective, or be empty");
                }
            }
        }

        ImpreciseModel(const ImpreciseModel &) = delete;
        ImpreciseModel &operator=(const ImpreciseModel &) = delete;
        ImpreciseModel(ImpreciseModel &&) = default;
        ImpreciseModel &operator=(ImpreciseModel &&) = default;

        const std::vector<double> &get_variable_values() {
            if(variable_values.empty()) {
                LPBackend::get_solution(*model, vars, variable_values);
            }
            return variable_values;
        }

        double inexact_objective_value() const { return LPBackend::objective_value(*model); }

        void print_stats() const { LPBackend::print_model_stats(*model); }

        std::unique_ptr<IModel> model;
        std::vector<IVar> vars;
        std::vector<IConstr> constrs;
        std::vector<Interval> objective_approx;
        std::vector<ObjectiveCoefficientType> objective_exact;
        std::vector<double> variable_values;

      private:
        void p_exact_objective_to_approx() {
            CGAL::Protect_FPU_rounding protect;
            objective_approx.reserve(objective_exact.size());
            for(const auto &coeff : objective_exact) {
                objective_approx.push_back(mwt::to_interval(coeff));
            }
        }
    };

    template<typename LPBackend, typename ObjectiveCoefficientType>
    ImpreciseModel<LPBackend, ObjectiveCoefficientType>
    to_imprecise_model(std::vector<ObjectiveCoefficientType> objective, std::vector<Interval> objective_approx = {},
                       bool verbose = false) {
        ImpreciseModel<LPBackend, ObjectiveCoefficientType> result(std::move(objective), std::move(objective_approx),
                                                                   verbose);
        p_imprecise_create_vars(result);
        p_imprecise_create_constraints(result);
        p_imprecise_set_objective(result);
        LPBackend::model_construction_done(*result.model);
        return result;
    }

    struct ImpreciseModelError : public std::runtime_error {
        using std::runtime_error::runtime_error;
    };

    /**
     * Let the imprecise LP solver backend solve the imprecise model.
     * Throws an ImpreciseModelError if the LP is found to be infeasible or unbounded.
     */
    template<typename LPBackend, typename ObjectiveCoefficientType>
    void solve_imprecise_model(ImpreciseModel<LPBackend, ObjectiveCoefficientType> &imodel) {
        imodel.variable_values.clear();
        auto status = LPBackend::optimize(*imodel.model);
        switch(status) {
        case LPStatus::OPTIMAL:
            return;
        case LPStatus::INFEASIBLE:
            throw ImpreciseModelError("Infeasible LP");
        case LPStatus::UNBOUNDED:
            throw ImpreciseModelError("Unbounded LP");
        default:
            throw std::logic_error("Unknown LP status");
        }
    }

    /**
     * Extract the basis information from the imprecise model.
     */
    template<typename LPBackend, typename ObjectiveCoefficientType>
    void imprecise_load_basis(ImpreciseModel<LPBackend, ObjectiveCoefficientType> &imodel) {
        std::vector<BasicStatus> var_basis;
        std::vector<BasicStatus> con_basis;
        LPBackend::get_basis_status(*imodel.model, imodel.vars, var_basis);
        LPBackend::get_basis_status(*imodel.model, imodel.constrs, con_basis);
        load_basis(var_basis, con_basis);
    }

    /**
     * Load a given basis.
     */
    void load_basis(const std::vector<BasicStatus> &var_basis, const std::vector<BasicStatus> &con_basis) {
        m_basis.clear();
        m_nonbasis.clear();
        m_basis_index.resize(m_original_vars + m_num_constraints);
        m_nonbasic_index.resize(m_original_vars + m_num_constraints);

        std::size_t bindex = 0, nindex = 0;
        for(std::size_t i = 0; i < m_original_vars; ++i) {
            switch(var_basis[i]) {
            case BasicStatus::LOWER_BOUND:
                m_basis_index[i] = NONBASIC_LB;
                m_nonbasis.push_back(i);
                m_nonbasic_index[i] = nindex++;
                break;
            case BasicStatus::UPPER_BOUND:
                m_basis_index[i] = NONBASIC_UB;
                m_nonbasis.push_back(i);
                m_nonbasic_index[i] = nindex++;
                break;
            case BasicStatus::BASIC:
                m_basis_index[i] = bindex++;
                m_basis.push_back(i);
                m_nonbasic_index[i] = std::numeric_limits<std::size_t>::max();
                break;
            }
        }
        for(std::size_t i = 0; i < m_num_constraints; ++i) {
            bool is_leq = !m_lhs[i].has_value();
            switch(con_basis[i]) {
            case BasicStatus::LOWER_BOUND:
                m_basis_index[m_original_vars + i] = is_leq ? NONBASIC_UB : NONBASIC_LB;
                m_nonbasis.push_back(m_original_vars + i);
                m_nonbasic_index[m_original_vars + i] = nindex++;
                break;
            case BasicStatus::UPPER_BOUND:
                m_basis_index[m_original_vars + i] = is_leq ? NONBASIC_LB : NONBASIC_UB;
                m_nonbasis.push_back(m_original_vars + i);
                m_nonbasic_index[m_original_vars + i] = nindex++;
                break;
            case BasicStatus::BASIC:
                m_basis_index[m_original_vars + i] = bindex++;
                m_basis.push_back(m_original_vars + i);
                m_nonbasic_index[m_original_vars + i] = std::numeric_limits<std::size_t>::max();
                break;
            }
        }
    }

    /**
     * Re-factorize the basis matrix and compute
     * the precise values of the basic variables after
     * loading a new basis.
     */
    void compute_basic_values() {
        p_recompute_lu();
        m_tmp_rhs.assign(m_num_constraints, 0);
        for(std::size_t i = 0; i < m_original_vars; ++i) {
            const NumberType *value = nullptr;
            if(m_basis_index[i] == NONBASIC_LB) {
                value = &m_var_lb[i].value();
            } else if(m_basis_index[i] == NONBASIC_UB) {
                value = &m_var_ub[i].value();
            }
            if(value && !!*value) {
                for(const auto &entry : m_original_matrix[i]) {
                    m_tmp_rhs[entry.index] -= entry.value * *value;
                }
            }
        }
        for(std::size_t i = 0; i < m_num_constraints; ++i) {
            std::size_t bii = i + m_original_vars;
            const NumberType *value = nullptr;
            if(m_basis_index[bii] == NONBASIC_LB) {
                value = &m_lhs[i].value();
            } else if(m_basis_index[bii] == NONBASIC_UB) {
                value = &m_rhs[i].value();
            }
            if(value && !!*value) {
                m_tmp_rhs[i] += *value;
            }
        }
        m_lu->solve_linear_system(m_tmp_rhs, m_basic_values, &m_tmp_buffer);
    }

    bool is_primal_feasible() const {
        for(std::size_t index_in_basis = 0; index_in_basis < m_num_constraints; ++index_in_basis) {
            std::size_t col_index = m_basis[index_in_basis];
            const auto &bval = m_basic_values[index_in_basis];
            const BoundType *lb, *ub;
            if(col_index < m_original_vars) {
                lb = &m_var_lb[col_index];
                ub = &m_var_ub[col_index];
            } else {
                lb = &m_lhs[col_index - m_original_vars];
                ub = &m_rhs[col_index - m_original_vars];
            }
            if(*lb && bval < **lb) {
                return false;
            }
            if(*ub && **ub < bval) {
                return false;
            }
        }
        return true;
    }

    /**
     * Use interval arithmetic to compute the objective coefficients
     * of each non-basic variable.
     */
    template<typename LPBackend, typename ObjectiveCoefficientType>
    void compute_interval_objective_coefficients(const ImpreciseModel<LPBackend, ObjectiveCoefficientType> &imodel) {
        CGAL::Protect_FPU_rounding protect;
        p_compute_cb_interval(imodel.objective_approx);
        m_lu->solve_transposed_with_intervals(m_cb_interval, m_interval_buffer, &m_nonbasic_objective_intervals);
        m_nonbasic_objective_intervals.clear();
        for(std::size_t nbi : m_nonbasis) {
            if(nbi < m_original_vars) {
                Interval value = imodel.objective_approx[nbi];
                for(const auto &entry : m_original_matrix[nbi]) {
                    value -= mwt::to_interval(entry.value) * m_interval_buffer[entry.index];
                }
                m_nonbasic_objective_intervals.push_back(value);
            } else {
                m_nonbasic_objective_intervals.push_back(m_interval_buffer[nbi - m_original_vars]);
            }
        }
    }

    const std::vector<std::size_t> &index_in_nonbasis_of_potential_dual_infeasibilities() const {
        return m_potential_dual_infeasibilities;
    }

    /**
     * Use the intervals computed by compute_interval_objective_coefficients
     * to find nonbasis elements that may still lead to an improvement of the objective.
     * Return a vector of indices (into the list of nonbasic columns) which might be dual infeasible.
     */
    const std::vector<std::size_t> &compute_index_in_nonbasis_of_potential_dual_infeasibilities() {
        m_potential_dual_infeasibilities.clear();
        for(std::size_t i = 0, k = m_nonbasis.size(); i < k; ++i) {
            std::size_t nbi = m_nonbasis[i];
            // non-basic equality constraints are simultaneously at LB and UB.
            // we could simply declare them to be at the bound we desire.
            if(nbi >= m_original_vars && m_is_equality[nbi - m_original_vars])
                continue;
            bool is_at_ub = (m_basis_index[nbi] == NONBASIC_UB);
            auto interval = m_nonbasic_objective_intervals[i];
            bool possibly_negative = (interval.inf() < 0.0);
            bool possibly_positive = (interval.sup() > 0.0);
            if(possibly_negative == possibly_positive) {
                if(possibly_negative) {
                    m_potential_dual_infeasibilities.push_back(i);
                }
                continue;
            }
            if(m_is_minimization ^ (possibly_negative == is_at_ub)) {
                m_potential_dual_infeasibilities.push_back(i);
            }
        }
        return m_potential_dual_infeasibilities;
    }

    /**
     * After loading a new basis and computing the basic values,
     * compute the coefficients of the objective function symbolically,
     * i.e., assuming the ith original variable has objective coefficient
     * s_i (where s_i is a unique symbol). Constraint variables
     * have no (=0) objective coefficient.
     */
    void compute_symbolic_objective_coefficients() {
        p_ensure_have_table();
        p_compute_cb_symbolic();
        m_nonbasic_objective_coefficients.resize(m_nonbasis.size());
        m_lu->solve_transposed_linear_system_with_symbols(*m_symbolic_table, m_symbolic_buffer, m_cb_symbolic);
        std::size_t index = 0;
        for(std::size_t nbi : m_nonbasis) {
            if(nbi < m_original_vars) {
                m_symbolic_table->add_entry(nbi);
            }
            m_symbolic_table->scalar_product_into(p_column_by_index(nbi), m_symbolic_buffer,
                                                  m_nonbasic_objective_coefficients[index++], std::true_type{});
        }
    }

    const std::vector<NumberType> &basic_values() const noexcept { return m_basic_values; }

    const std::vector<std::size_t> &basis_column_indices() const noexcept { return m_basis; }

    const std::vector<std::size_t> &nonbasis_column_indices() const noexcept { return m_nonbasis; }

    const NumberType &solution_value(std::size_t var_index) {
        std::size_t bindex = m_basis_index[var_index];
        if(bindex < NONBASIC_UB) {
            // basic
            return m_basic_values[bindex];
        } else if(bindex == NONBASIC_UB) {
            // nonbasic at UB
            if(var_index < m_original_vars) {
                return m_var_ub[var_index].value();
            } else {
                return m_rhs[var_index - m_original_vars].value();
            }
        } else {
            // nonbasic at LB
            if(var_index < m_original_vars) {
                return m_var_lb[var_index].value();
            } else {
                return m_lhs[var_index - m_original_vars].value();
            }
        }
    }

    const std::vector<std::vector<SymEntry>> &nonbasic_objective_coefficients() const noexcept {
        return m_nonbasic_objective_coefficients;
    }

    const std::vector<Interval> &nonbasic_objective_intervals() const noexcept {
        return m_nonbasic_objective_intervals;
    }

    template<typename LPBackend, typename ObjectiveCoefficientType,
             typename CombineCallable /*(NumberType, ObjectiveCoefficientType) -> ObjectiveCoefficientType*/>
    std::vector<std::pair<std::size_t, ObjectiveCoefficientType>>
    compute_index_in_nonbasis_of_dual_infeasibilities_symbolically(
        const ImpreciseModel<LPBackend, ObjectiveCoefficientType> &inexact, CombineCallable &&callable) {
        std::vector<std::pair<std::size_t, ObjectiveCoefficientType>> result;
        for(std::size_t i_in_nb : m_potential_dual_infeasibilities) {
            std::size_t nbi = m_nonbasis[i_in_nb];
            auto &coefficients = m_nonbasic_objective_coefficients[i_in_nb];
            ObjectiveCoefficientType v = substitute_symbols(coefficients, inexact.objective_exact, callable);
            Interval interval = mwt::exact_to_interval(v);
            m_nonbasic_objective_intervals[i_in_nb] = interval;
            bool negative = false;
            if constexpr(is_real_embeddable<ObjectiveCoefficientType>()) {
                auto s = CGAL::sign(v);
                if(s == CGAL::ZERO)
                    continue;
                negative = (s == CGAL::NEGATIVE);
            } else {
                if(!v)
                    continue;
                negative = (v < 0);
            }
            bool is_at_ub = (m_basis_index[nbi] == NONBASIC_UB);
            if(m_is_minimization ^ (negative == is_at_ub)) {
                result.emplace_back(i_in_nb, std::move(v));
            }
        }
        return result;
    }

    template<typename LPBackend, typename ObjectiveCoefficientType,
             typename ComputeSignCallable /*(const std::vector<SymEntry>&) -> CGAL::Sign*/>
    std::vector<std::size_t> compute_index_in_nonbasis_of_dual_infeasibilities_with_custom_sign(
        const ImpreciseModel<LPBackend, ObjectiveCoefficientType> &inexact, ComputeSignCallable &&callable) {
        std::vector<std::size_t> result;
        for(std::size_t i_in_nb : m_potential_dual_infeasibilities) {
            std::size_t nbi = m_nonbasis[i_in_nb];
            const auto &coefficients = m_nonbasic_objective_coefficients[i_in_nb];
            auto sign = callable(coefficients);
            if(sign == CGAL::ZERO)
                continue;
            bool negative = (sign == CGAL::NEGATIVE);
            bool is_at_ub = (m_basis_index[nbi] == NONBASIC_UB);
            if(m_is_minimization ^ (negative == is_at_ub)) {
                result.emplace_back(i_in_nb);
            }
        }
        return result;
    }

    template<typename LPBackend, typename ObjectiveCoefficientType>
    std::vector<std::pair<std::size_t, ObjectiveCoefficientType>>
    compute_index_in_nonbasis_of_dual_infeasibilities_symbolically(
        const ImpreciseModel<LPBackend, ObjectiveCoefficientType> &inexact) {
        return compute_index_in_nonbasis_of_dual_infeasibilities_symbolically(
            inexact,
            [](const NumberType &n, const ObjectiveCoefficientType &o) -> ObjectiveCoefficientType { return n * o; });
    }

    template<typename InexactModel, typename ExactDualInfeasibilities>
    void primal_simplex_iteration(const InexactModel &inexact_model, ExactDualInfeasibilities &dual_infeasibilities) {
        auto [entering_var, entering_in_nonbasis] = p_choose_entering_first_index(dual_infeasibilities);
        bool moving_down = (m_basis_index[entering_var] == NONBASIC_UB);
        p_compute_primal_step_direction(entering_var, entering_in_nonbasis, moving_down);
        auto [leaving_var, leaving_in_basis, hit_ub, step_length] = p_compute_leaving_var();
        p_compute_dual_step_direction(leaving_var, leaving_in_basis);
        p_dual_step_length_to_table(entering_in_nonbasis);
        p_update_dual_values(inexact_model, entering_in_nonbasis, hit_ub);
        p_update_basic_values(step_length, entering_var, leaving_in_basis, moving_down);
        p_update_indices(leaving_var, leaving_in_basis, entering_var, entering_in_nonbasis, hit_ub);
        p_primal_dir_and_leaving_to_eta(leaving_in_basis);
    }

    void filter_certain_dual_infeasibilities(std::vector<std::size_t> &indices_in_nonbasis_out) {
        indices_in_nonbasis_out.clear();
        for(std::size_t i_in_nb : m_potential_dual_infeasibilities) {
            auto iv = m_nonbasic_objective_intervals[i_in_nb];
            if(iv.inf() >= 0 || iv.sup() <= 0) {
                indices_in_nonbasis_out.push_back(i_in_nb);
            }
        }
    }

  private:
    /**
     * Used to re-use the LU decomposition for more
     * than one primal simplex iteration.
     */
    struct EtaEntry {
        EtaEntry(std::size_t leaving_in_basis, std::vector<NumberType> &primal_dir)
            : index_in_basis(leaving_in_basis), x_i(std::move(primal_dir[leaving_in_basis])) {
            for(std::size_t i = 0; i < leaving_in_basis; ++i) {
                if(!primal_dir[i])
                    continue;
                eta_values.push_back(Entry{i, std::move(primal_dir[i])});
            }
            eta_values.push_back(Entry{leaving_in_basis, x_i - 1});
            for(std::size_t i = leaving_in_basis + 1, k = primal_dir.size(); i < k; ++i) {
                if(!primal_dir[i])
                    continue;
                eta_values.push_back(Entry{i, std::move(primal_dir[i])});
            }
        }

        std::size_t index_in_basis; //< the index i of the column that was replaced
        NumberType x_i;             //< the value in x_B of the column that was replaced
        ColumnType eta_values;      //< the vector (\Delta x_B - e_i)
    };

    void p_dual_step_length_to_table(std::size_t entering_in_nonbasis) {
        m_symbolic_table->read_symbolic_value(m_nonbasic_objective_coefficients[entering_in_nonbasis]);
        m_symbolic_table->divide_coefficients(m_dual_dir[entering_in_nonbasis]);
    }

    template<typename InexactModel>
    void p_update_dual_values(const InexactModel &inexact, std::size_t entering_in_nonbasis, bool basic_moving_up) {
        for(std::size_t i = 0, k = m_nonbasis.size(); i < k; ++i) {
            m_symbolic_table->subtract_scaled_from(m_nonbasic_objective_coefficients[i], m_dual_dir[i]);
        }
        if(!m_nonbasic_objective_coefficients[entering_in_nonbasis].empty()) {
            throw std::logic_error("Coefficient of entering variable in dual is not zero!");
        }
        m_symbolic_table->output_coefficients_to_vector(m_nonbasic_objective_coefficients[entering_in_nonbasis], true);
        p_update_dual_intervals(inexact);
    }

    template<typename InexactModel> void p_update_dual_intervals(const InexactModel &inexact) {
        CGAL::Protect_FPU_rounding protect;
        for(std::size_t i = 0, k = m_nonbasis.size(); i < k; ++i) {
            auto &interval = m_nonbasic_objective_intervals[i];
            interval = substitute_symbols(
                m_nonbasic_objective_coefficients[i], inexact.objective_approx,
                [&](const NumberType &n, const Interval &o) -> Interval { return mwt::exact_to_interval(n) * o; });
        }
    }

    void p_update_basic_values(const NumberType &step_length, std::size_t entering_var, std::size_t leaving_in_basis,
                               bool moving_down) {
        if(step_length) {
            for(std::size_t i = 0; i < m_num_constraints; ++i) {
                m_basic_values[i] -= m_primal_dir[i] * step_length;
            }
        }
        if(moving_down) {
            m_basic_values[leaving_in_basis] = p_ub_by_index(entering_var) - step_length;
        } else {
            m_basic_values[leaving_in_basis] = p_lb_by_index(entering_var) + step_length;
        }
    }

    void p_update_indices(std::size_t leaving_var, std::size_t leaving_in_basis, std::size_t entering_var,
                          std::size_t entering_in_nonbasis, bool hit_ub) {
        m_basis_index[leaving_var] = hit_ub ? NONBASIC_UB : NONBASIC_LB;
        m_basis_index[entering_var] = leaving_in_basis;
        m_nonbasic_index[leaving_var] = entering_in_nonbasis;
        m_nonbasic_index[entering_var] = std::numeric_limits<std::size_t>::max();
        m_basis[leaving_in_basis] = entering_var;
        m_nonbasis[entering_in_nonbasis] = leaving_var;
    }

    template<typename ExactDualInfeasibilities>
    std::pair<std::size_t, std::size_t> p_choose_entering_first_index(const ExactDualInfeasibilities &di) {
        std::size_t min_nbi = std::numeric_limits<std::size_t>::max();
        std::size_t index_in_nonbasis = std::numeric_limits<std::size_t>::max();
        for(const auto &e : di) {
            std::size_t nbi;
            std::size_t i_in_nb;
            if constexpr(std::is_convertible_v<decltype(e), std::size_t>) {
                i_in_nb = e;
                nbi = m_nonbasis[i_in_nb];
            } else {
                i_in_nb = e.first;
                nbi = m_nonbasis[i_in_nb];
            }
            if(nbi < min_nbi) {
                min_nbi = nbi;
                index_in_nonbasis = i_in_nb;
            }
        }
        return {min_nbi, index_in_nonbasis};
    }

    void p_compute_primal_step_direction(std::size_t entering_col, std::size_t entering_in_nb, bool downwards) {
        m_tmp_rhs.assign(m_num_constraints, 0);
        for(const auto &entry : p_column_by_index(entering_col)) {
            m_tmp_rhs[entry.index] = entry.value;
        }
        m_lu->solve_linear_system(m_tmp_rhs, m_primal_dir, &m_tmp_buffer);
        for(const auto &e : m_eta_entries) {
            p_apply_eta(e, m_primal_dir);
        }
        if(downwards) {
            for(auto &v : m_primal_dir) {
                v.negate();
            }
        }
    }

    std::tuple<std::size_t, std::size_t, bool, NumberType> p_compute_leaving_var() {
        NumberType step_length(0);
        NumberType this_length(0);
        std::size_t leaving_var = std::numeric_limits<std::size_t>::max();
        std::size_t leaving_in_basis = std::numeric_limits<std::size_t>::max();
        bool hit_ub = false;
        for(std::size_t index_in_basis = 0; index_in_basis < m_num_constraints; ++index_in_basis) {
            const auto &dir_entry = m_primal_dir[index_in_basis];
            CGAL::Sign sgn = dir_entry.sign();
            if(sgn == CGAL::ZERO)
                continue;
            std::size_t bindex = m_basis[index_in_basis];
            auto [lb, ub] = p_bounds_by_index(bindex);
            if(sgn == CGAL::POSITIVE) {
                if(!lb)
                    continue;
                this_length = m_basic_values[index_in_basis];
                this_length -= *lb;
                this_length /= dir_entry;
                if(leaving_in_basis > m_num_constraints || this_length <= step_length) {
                    if(this_length == step_length && leaving_var < bindex)
                        continue;
                    step_length = std::move(this_length);
                    leaving_var = bindex;
                    leaving_in_basis = index_in_basis;
                    hit_ub = false;
                }
            } else {
                if(!ub)
                    continue;
                this_length = m_basic_values[index_in_basis];
                this_length -= *ub;
                this_length /= dir_entry;
                if(leaving_in_basis > m_num_constraints || this_length <= step_length) {
                    if(this_length == step_length && leaving_var < bindex)
                        continue;
                    step_length = std::move(this_length);
                    leaving_var = bindex;
                    leaving_in_basis = index_in_basis;
                    hit_ub = true;
                }
            }
        }
        return {leaving_var, leaving_in_basis, hit_ub, std::move(step_length)};
    }

    const NumberType &p_lb_by_index(std::size_t index) const {
        if(index < m_original_vars) {
            return m_var_lb[index].value();
        } else {
            return m_lhs[index - m_original_vars].value();
        }
    }

    const NumberType &p_ub_by_index(std::size_t index) const {
        if(index < m_original_vars) {
            return m_var_ub[index].value();
        } else {
            return m_rhs[index - m_original_vars].value();
        }
    }

    std::pair<const NumberType *, const NumberType *> p_bounds_by_index(std::size_t index) {
        auto opt_ptr = [](const BoundType &b) -> const NumberType * { return b ? &*b : nullptr; };
        if(index < m_original_vars) {
            const auto &lb = m_var_lb[index];
            const auto &ub = m_var_ub[index];
            return {opt_ptr(lb), opt_ptr(ub)};
        } else {
            const auto &lb = m_lhs[index - m_original_vars];
            const auto &ub = m_rhs[index - m_original_vars];
            return {opt_ptr(lb), opt_ptr(ub)};
        }
    }

    void p_compute_dual_step_direction(std::size_t leaving_var, std::size_t leaving_in_basis) {
        m_tmp_rhs.assign(m_num_constraints, 0);
        m_tmp_rhs[leaving_in_basis] = 1;
        for(auto eta_i = m_eta_entries.rbegin(); eta_i != m_eta_entries.rend(); ++eta_i) {
            p_apply_eta_transposed(*eta_i, m_tmp_rhs);
        }
        m_lu->solve_transposed_linear_system(m_tmp_rhs, m_dual_dir_no_N, &m_tmp_buffer);
        m_dual_dir.resize(m_nonbasis.size());
        for(std::size_t index_in_nonbasis = 0, k = m_nonbasis.size(); index_in_nonbasis < k; ++index_in_nonbasis) {
            std::size_t nbi = m_nonbasis[index_in_nonbasis];
            if(nbi < m_original_vars) {
                m_dual_dir[index_in_nonbasis] = 0;
                for(const auto &entry : m_original_matrix[nbi]) {
                    m_dual_dir[index_in_nonbasis] += entry.value * m_dual_dir_no_N[entry.index];
                }
            } else {
                m_dual_dir[index_in_nonbasis] = -m_dual_dir_no_N[nbi - m_original_vars];
            }
        }
    }

    void p_apply_eta(const EtaEntry &eta, std::vector<NumberType> &vec) const {
        std::size_t index = eta.index_in_basis;
        if(!vec[index])
            return;
        NumberType factor = vec[index];
        factor /= eta.x_i;
        for(const auto &entry : eta.eta_values) {
            vec[entry.index] -= entry.value * factor;
        }
    }

    void p_apply_eta_transposed(const EtaEntry &eta, std::vector<NumberType> &vec) const {
        std::size_t index = eta.index_in_basis;
        NumberType scalprod(0);
        for(const auto &entry : eta.eta_values) {
            scalprod += entry.value * vec[entry.index];
        }
        scalprod /= eta.x_i;
        vec[index] -= scalprod;
    }

    void p_primal_dir_and_leaving_to_eta(std::size_t leaving_in_basis) {
        if(m_eta_entries.size() >= 1000) {
            m_eta_entries.clear();
            p_recompute_lu();
            return;
        }
        m_eta_entries.emplace_back(leaving_in_basis, m_primal_dir);
    }

    void p_compute_cb_interval(const std::vector<Interval> &objective) {
        m_cb_interval.clear();
        m_cb_interval.reserve(m_basis.size());
        for(std::size_t vb : m_basis) {
            if(vb >= m_original_vars) {
                m_cb_interval.emplace_back(0.0);
            } else {
                m_cb_interval.emplace_back(objective[vb]);
            }
        }
    }

    void p_compute_cb_symbolic() {
        m_cb_symbolic.resize(m_basis.size());
        std::size_t index = 0;
        for(std::size_t vb : m_basis) {
            if(vb >= m_original_vars) {
                m_cb_symbolic[index++] = std::nullopt;
            } else {
                m_cb_symbolic[index++] = vb;
            }
        }
    }

    void p_ensure_have_table() {
        if(!m_symbolic_table) {
            m_symbolic_table.emplace(m_original_vars, m_num_constraints);
        }
    }

    void p_recompute_lu() {
        if(!m_lu) {
            m_lu.emplace(m_num_constraints, typename RationalSparseLU<NumberType>::ColumnMajor{},
                         p_basic_column_range());
        } else {
            m_lu->reset_columns(m_num_constraints, p_basic_column_range());
        }
        m_lu->compute_decomposition();
    }

    auto p_basic_column_range() const {
        auto transform = [this](std::size_t index) -> const ColumnType & { return p_column_by_index(index); };
        auto tbegin = boost::make_transform_iterator(m_basis.begin(), transform);
        auto tend = boost::make_transform_iterator(m_basis.end(), transform);
        return boost::make_iterator_range(tbegin, tend);
    }

    const ColumnType &p_column_by_index(std::size_t index) const {
        if(index < m_original_vars) {
            return m_original_matrix[index];
        } else {
            m_tmp_column.clear();
            m_tmp_column.push_back(Entry{index - m_original_vars, -1});
            return m_tmp_column;
        }
    }

    /**
     * LP information.
     */
    bool m_is_minimization;
    std::size_t m_original_vars{0};            //< the original number of variables
    std::size_t m_num_constraints{0};          //< the number of constraints
    std::vector<ColumnType> m_original_matrix; //< the matrix A of our LP
    std::vector<BoundType> m_lhs;              //< the LHS vector a of our LP (a <= Ax <= b)
    std::vector<BoundType> m_rhs;              //< the RHS vector b of our LP (a <= Ax <= b)
    std::vector<BoundType> m_var_lb;           //< the lower bound vector l of the variables (l <= x <= u)
    std::vector<BoundType> m_var_ub;           //< the upper bound vector u of the variables (l <= x <= u)
    std::vector<bool> m_is_equality;           //< whether the i'th constraint is an equality

    /**
     * Current solution information.
     */
    constexpr static std::size_t NONBASIC_LB = std::numeric_limits<std::size_t>::max();
    constexpr static std::size_t NONBASIC_UB = std::numeric_limits<std::size_t>::max() - 1;

    std::vector<std::size_t> m_basis_index; //< the index in the basis of each variable, or NONBASIC_LB/UB
    std::vector<std::size_t> m_basis;       //< the column indices of the basis columns
    std::vector<NumberType> m_basic_values; //< the values of the basic variables

    std::vector<std::size_t> m_nonbasic_index; //< the index in the non-basic variables
    std::vector<std::size_t> m_nonbasis;       //< the column indices of the non-basic columns

    /**
     * Interval objective related variables/information.
     */
    std::vector<Interval> m_cb_interval;     //< the interval coefficients of the basic variables in the objective
    std::vector<Interval> m_interval_buffer; //< buffer for temporary interval values
    std::vector<Interval> m_nonbasic_objective_intervals; //< intervals containing the nonbasic objective coefficients
    std::vector<std::size_t>
        m_potential_dual_infeasibilities; //< the index in the nonbasis of potential dual infeasibilities

    /**
     * Symbolic objective related variables/information.
     */
    std::vector<std::vector<SymEntry>> m_symbolic_buffer; //< buffer to solve symbolic systems into
    std::vector<std::vector<SymEntry>>
        m_nonbasic_objective_coefficients; //< the symbolic coefficients of the non-basic variables in the objective
    std::vector<std::optional<std::size_t>> m_cb_symbolic; //< the symbolic coefficients of the original basic variables
    std::optional<SymTable> m_symbolic_table;              //< the symbolic table we use for computations

    /**
     * Buffers and decompositions used for solving linear systems.
     */
    std::optional<RationalSparseLU<NumberType>> m_lu; //< if we have it, the LU decomposition of the basis matrix
    mutable ColumnType m_tmp_column;                  //< a temporary column buffer
    std::vector<NumberType> m_tmp_rhs;                //< a temporary RHS buffer
    std::vector<NumberType> m_primal_dir;             //< a buffer for the primal step direction
    std::vector<NumberType> m_dual_dir;               //< a buffer for the dual step direction
    std::vector<NumberType> m_tmp_buffer;             //< temporary buffer for forward/backward substitution
    std::vector<NumberType> m_dual_dir_no_N;          //< a buffer for the dual step direction without the N^T part
    std::vector<EtaEntry> m_eta_entries;              //< set of eta entries used for re-using the LU decomposition

    /**
     * Create the variables in the imprecise model.
     */
    template<typename LPBackend, typename ObjectiveCoefficientType>
    void p_imprecise_create_vars(ImpreciseModel<LPBackend, ObjectiveCoefficientType> &imodel) {
        CGAL::Protect_FPU_rounding protect;
        imodel.vars.reserve(m_original_vars);
        for(VarIndex i = 0; i < m_original_vars; ++i) {
            auto lb = LPBackend::infinite_lb();
            auto ub = LPBackend::infinite_ub();
            if(m_var_lb[i]) {
                lb = mwt::to_interval(*m_var_lb[i]).inf();
            }
            if(m_var_ub[i]) {
                ub = mwt::to_interval(*m_var_ub[i]).sup();
            }
            imodel.vars.push_back(LPBackend::create_continuous(*imodel.model, lb, ub));
        }
    }

    /**
     * Create the constraints in the imprecise model.
     */
    template<typename LPBackend, typename ObjectiveCoefficientType>
    void p_imprecise_create_constraints(ImpreciseModel<LPBackend, ObjectiveCoefficientType> &imodel) {
        CGAL::Protect_FPU_rounding protect;
        imodel.constrs.reserve(m_num_constraints);
        auto expr = LPBackend::empty_expression();
        std::vector<std::vector<std::pair<double, VarIndex>>> rows;
        rows.resize(m_num_constraints);
        for(VarIndex v = 0; v < m_original_vars; ++v) {
            const auto &col = m_original_matrix[v];
            for(const auto &entry : col) {
                rows[entry.index].emplace_back(mwt::to_interval(entry.value).inf(), v);
            }
        }
        for(ConstrIndex i = 0; i < m_num_constraints; ++i) {
            const auto &row = rows[i];
            LPBackend::clear_expression(expr);
            for(const auto &entry : row) {
                LPBackend::add_to_expression(expr, entry.first, imodel.vars[entry.second]);
            }
            if(!m_lhs[i] && !m_rhs[i]) {
                throw std::logic_error("Bidirectionally unbounded superfluous constraint");
            }
            if(!m_lhs[i]) {
                imodel.constrs.push_back(
                    LPBackend::add_less_equal(*imodel.model, expr, mwt::to_interval(*m_rhs[i]).sup()));
                continue;
            } else if(!m_rhs[i]) {
                imodel.constrs.push_back(
                    LPBackend::add_greater_equal(*imodel.model, expr, mwt::to_interval(*m_lhs[i]).inf()));
                continue;
            }
            if(m_is_equality[i]) {
                auto value = mwt::to_interval(*m_lhs[i]).inf();
                imodel.constrs.push_back(LPBackend::add_equal(*imodel.model, expr, value));
            } else {
                auto lb = mwt::to_interval(*m_lhs[i]).inf();
                auto ub = mwt::to_interval(*m_rhs[i]).sup();
                imodel.constrs.push_back(LPBackend::add_range(*imodel.model, expr, lb, ub));
            }
        }
    }

    /**
     * Set the objective in the imprecise model.
     */
    template<typename LPBackend, typename ObjectiveCoefficientType>
    void p_imprecise_set_objective(ImpreciseModel<LPBackend, ObjectiveCoefficientType> &imodel) {
        auto expr = LPBackend::empty_expression();
        for(VarIndex i = 0; i < m_original_vars; ++i) {
            LPBackend::add_to_expression(expr, imodel.objective_approx[i].inf(), imodel.vars[i]);
        }
        if(m_is_minimization) {
            LPBackend::minimize(*imodel.model, expr);
        } else {
            LPBackend::maximize(*imodel.model, expr);
        }
    }
};

} // namespace mwt

#endif
