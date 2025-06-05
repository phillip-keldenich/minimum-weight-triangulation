#ifndef CGAL_MWT_EXACT_LP_FACE_TRIANGULATOR_H_INCLUDED_
#define CGAL_MWT_EXACT_LP_FACE_TRIANGULATOR_H_INCLUDED_

#include "Dynamic_program_utils.h"
#include "LP_backends.h"
#include "Rational_or_int.h"
#include "Rational_sparse_matrix.h"
#include "Scaled_weight_sign.h"
#include "Triangle_based_linear_model.h"
#include "time_util.h"
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel_with_sqrt.h>
#include <nlohmann/json.hpp>

namespace mwt {

template<typename Skeleton_, typename LPBackend_> class Exact_LP_face_triangulator {
  public:
    using LPBackend = LPBackend_;
    using Skeleton = Skeleton_;
    using Analyzer = Face_analyzer<Skeleton>;
    using Face = mwt::Face<Skeleton>;
    using HalfedgeHandle = const typename Skeleton::Halfedge *;
    using Interval = CGAL::Interval_nt_advanced;
    using Rational = CGAL::Exact_rational;
    using RationalOrInt = mwt::RationalOrInt<Rational>;
    using SqrtNT = CGAL::Exact_predicates_exact_constructions_kernel_with_sqrt::FT;
    using ExactLP = RationalSparseLPWithSymbolicObjective<RationalOrInt>;
    using SymEntry = typename ExactLP::SymEntry;
    using ObjectiveCoefficientType = CGAL::Exact_predicates_exact_constructions_kernel_with_sqrt::FT;
    using InexactLP = typename ExactLP::template ImpreciseModel<LPBackend, ObjectiveCoefficientType>;
    using LazyExactRational = CGAL::Lazy_exact_nt<Rational>;
    using ScaledWS = Scaled_weight_sign<typename Skeleton::Halfedge, Rational, LazyExactRational>;

    explicit Exact_LP_face_triangulator(FaceTriangulatorOptions options) : m_options(options) {}

    void operator()(Face &&face) {
        m_face_analysis_time = measure_time([&]() { m_analyzer.analyze_face(face); });
        m_num_triangle_vars = m_analyzer.empty_triangles().size();
        m_rank_reduction_time = measure_time([&]() { p_build_rank_reduced_constraints(); });
        p_create_imprecise_model();
        m_imprecise_solve_time = measure_time([&]() { m_lp->solve_imprecise_model(*m_inexact); });
        if(!p_check_imprecise_is_integral()) {
            throw std::runtime_error("Inexact solution is not close to integral; not supported in exact mode!");
        }
        if(!p_exact_check_primal_feasible()) {
            throw std::runtime_error("Basis from LP solver is primally infeasible -- this is most likely a bug!");
        }
        if(!p_check_dual_feasibility_with_intervals()) {
            m_symbolic_objective_coefficient_initialization_time =
                measure_time([&]() { m_lp->compute_symbolic_objective_coefficients(); });
            m_exact_simplex_time = measure_time([&]() { p_simplex_till_exact_optimality(); });
        }
        if(!p_check_exact_is_integral()) {
            throw std::runtime_error("Exact solution is not integral -- this is a very interesting instance!");
        }
        p_solution_value_to_status();
    }

    void combine_stats(nlohmann::json &overall_stats) {
        nlohmann::json &stat_list = overall_stats["nonsimple_faces"];
        nlohmann::json &stats = stat_list.emplace_back();
        stats["exact_simplex_iterations"] = m_exact_simplex_iterations;
        stats["exact_dual_infeasibility_computations"] = m_exact_dual_infeasibility_computations;
        stats["exact_simplex_time"] = m_exact_simplex_time;
        stats["rank_reduction_time"] = m_rank_reduction_time;
        stats["rows_before_reduction"] = m_rows_before_reduction;
        stats["rows_after_reduction"] = m_rows_after_reduction;
        stats["num_triangle_vars"] = m_num_triangle_vars;
        stats["face_analysis_time"] = m_face_analysis_time;
        stats["exact_coefficient_initialization_time"] = m_exact_coefficient_initialization_time;
        stats["imprecise_model_creation_time"] = m_imprecise_model_creation_time;
        stats["imprecise_lp_solve_time"] = m_imprecise_solve_time;
        stats["imprecise_solve_stats"] = LPBackend::get_solution_process_stats(*m_inexact->model);
        stats["exact_startup_time"] = m_exact_starting_time;
        stats["initial_interval_dual_feasibility_check_time"] = m_initial_interval_dual_feasibility_check_time;
        stats["intervals_show_inexact_lp_was_optimal"] = m_intervals_show_inexact_lp_was_optimal;
        stats["symbolic_objective_coefficient_initialization_time"] =
            m_symbolic_objective_coefficient_initialization_time;
    }

  private:
    void p_simplex_till_exact_optimality() {
        std::vector<std::size_t> definite_dual_infeasibilities;

        for(;; m_lp->compute_index_in_nonbasis_of_potential_dual_infeasibilities()) {
            // if interval arithmetic has shown that there are no dual
            // infeasibilities, we are optimal and can stop
            const auto &pdis = m_lp->index_in_nonbasis_of_potential_dual_infeasibilities();
            if(pdis.empty()) {
                return;
            }

            // if interval arithmetic has shown definite dual infeasibilities,
            // we do not have to resort to CGAL's exact type with square root
            // to continue primal simplex pivots
            m_lp->filter_certain_dual_infeasibilities(definite_dual_infeasibilities);
            if(!definite_dual_infeasibilities.empty()) {
                ++m_exact_simplex_iterations;
                m_lp->primal_simplex_iteration(*m_inexact, definite_dual_infeasibilities);
                continue;
            }

            // we have to resort to exact real arithmetic to continue
            ++m_exact_dual_infeasibility_computations;
            auto compute_sign = [this](const std::vector<SymEntry> &symbolic) {
                return p_compute_sign_of_symbolic_expression(symbolic);
            };
            auto exact_infeasibilities =
                m_lp->compute_index_in_nonbasis_of_dual_infeasibilities_with_custom_sign(*m_inexact, compute_sign);
            if(exact_infeasibilities.empty()) {
                // computing exactly, we find the solution to be optimal
                return;
            }
            ++m_exact_simplex_iterations;
            m_lp->primal_simplex_iteration(*m_inexact, exact_infeasibilities);
        }
    }

    CGAL::Sign p_compute_sign_of_symbolic_expression(const std::vector<SymEntry> &symbolic) {
        assert(CGAL::FPU_get_cw() == CGAL_FE_UPWARD);
        const auto &triangles = m_analyzer.empty_triangles();
        m_scaled_weight_sign.clear();
        for(const auto &entry : symbolic) {
            const auto &triangle = triangles[entry.symbol];
            for(auto *halfedge_ptr : triangle.edges) {
                if(halfedge_ptr->status == LMTStatus::Possible) {
                    m_scaled_weight_sign.add_halfedge_with_coefficient(entry.value, halfedge_ptr);
                }
            }
        }
        return m_scaled_weight_sign.compute_sign();
    }

    void p_solution_value_to_status() {
        p_mark_chosen_as_certain();
        detail::mark_remaining_as_impossible(m_analyzer.face().inner_edges);
    }

    void p_mark_chosen_as_certain() {
        const auto &triangles = m_analyzer.empty_triangles();
        for(std::size_t i = 0, nt = triangles.size(); i < nt; ++i) {
            if(!!m_lp->solution_value(i)) {
                for(auto *h : triangles[i].edges) {
                    if(h->status == LMTStatus::Possible) {
                        h->status = LMTStatus::Certain;
                        h->twin()->status = LMTStatus::Certain;
                    }
                }
            }
        }
    }

    bool p_check_dual_feasibility_with_intervals() {
        m_initial_interval_dual_feasibility_check_time = measure_time([&]() {
            m_lp->compute_interval_objective_coefficients(*m_inexact);
            m_intervals_show_inexact_lp_was_optimal =
                m_lp->compute_index_in_nonbasis_of_potential_dual_infeasibilities().empty();
        });
        return m_intervals_show_inexact_lp_was_optimal;
    }

    bool p_exact_check_primal_feasible() {
        m_exact_starting_time = measure_time([&]() {
            m_lp->imprecise_load_basis(*m_inexact);
            m_lp->compute_basic_values();
        });
        return m_lp->is_primal_feasible();
    }

    bool p_check_exact_is_integral() {
        const auto &basic_values = m_lp->basic_values();
        return std::all_of(basic_values.begin(), basic_values.end(), [](const auto &val) { return val.is_integer(); });
    }

    bool p_check_imprecise_is_integral() {
        const auto &values = m_inexact->get_variable_values();
        return std::all_of(values.begin(), values.end(),
                           [](double v) { return (v >= -0.05 && v <= 0.05) || (v >= 0.95 && v <= 1.05); });
    }

    void p_create_imprecise_model() {
        CGAL::Protect_FPU_rounding protect_rounding;
        std::vector<ObjectiveCoefficientType> objective_coefficients;
        std::vector<Interval> objective_approx;
        m_exact_coefficient_initialization_time = measure_time([&]() {
            objective_coefficients.reserve(m_analyzer.empty_triangles().size());
            for(const auto &triangle : m_analyzer.empty_triangles()) {
                objective_coefficients.push_back(triangle.template possible_weight<ObjectiveCoefficientType>());
                objective_approx.push_back(triangle.template possible_weight<Interval>());
            }
        });
        m_imprecise_model_creation_time = measure_time([&]() {
            m_inexact.emplace(m_lp->to_imprecise_model<LPBackend>(std::move(objective_coefficients),
                                                                  std::move(objective_approx), m_options.lp_verbose));
        });
    }

    void p_build_rank_reduced_constraints() {
        const auto &triangles = m_analyzer.empty_triangles();
        RationalSparseMatrix<RationalOrInt> rank_reduction(triangles.size());
        for(auto *container :
            {&m_analyzer.face().boundary, &m_analyzer.face().hole_boundaries, &m_analyzer.face().inner_edges}) {
            for(HalfedgeHandle h : *container) {
                rank_reduction.start_new_row();
                create_halfedge_constraint<int>(
                    m_analyzer, h,
                    [&](int coeff, std::size_t var_index) { rank_reduction.new_row_add_entry(var_index, coeff); },
                    [&](int rhs) { rank_reduction.finish_row(rhs); });
            }
        }
        m_rows_before_reduction = rank_reduction.num_rows();
        m_rows_after_reduction = rank_reduction.eliminate();
        rank_reduction.remove_zero_rows();
        m_rows_after_reduction = rank_reduction.num_rows();
        m_lp.emplace(true);
        m_lp->add_variables(triangles.size(), 0, {});
        m_lp->add_equalities(std::move(rank_reduction));
    }

    Analyzer m_analyzer;
    std::optional<ExactLP> m_lp;
    std::optional<InexactLP> m_inexact;
    nlohmann::json m_run_stats;

    FaceTriangulatorOptions m_options;
    std::size_t m_exact_simplex_iterations = 0;
    std::size_t m_exact_dual_infeasibility_computations = 0;
    double m_face_analysis_time;
    double m_rank_reduction_time;
    double m_exact_coefficient_initialization_time;
    double m_imprecise_model_creation_time;
    double m_imprecise_solve_time;
    double m_exact_starting_time;
    double m_initial_interval_dual_feasibility_check_time;
    double m_symbolic_objective_coefficient_initialization_time = 0.0;
    double m_exact_simplex_time = 0.0;
    std::size_t m_rows_before_reduction;
    std::size_t m_rows_after_reduction;
    std::size_t m_num_triangle_vars;
    bool m_intervals_show_inexact_lp_was_optimal = false;
    ScaledWS m_scaled_weight_sign;
};

} // namespace mwt

#endif
