#include <CGAL_MWT/Gurobi_LP_backend.h>
#include <CGAL_MWT/Null_LP_backend.h>
#include <CGAL_MWT/Rational_sparse_matrix.h>
#include <doctest/doctest.h>

TEST_CASE("[RationalSparseLP] simple LP test - values after basis loading") {
    using InnerNT = CGAL::Exact_rational;
    using NT = mwt::RationalOrInt<InnerNT>;
    using LP = mwt::RationalSparseLPWithSymbolicObjective<NT>;
    using Entry = LP::Entry;
    using LENT = CGAL::Exact_predicates_exact_constructions_kernel::FT;

    SUBCASE("known LP, only <= constraints") {
        std::vector<LENT> objective_coefficients{5, 4, 3};
        LP lp(false);
        auto x1 = lp.add_variable(0, std::nullopt);
        auto x2 = lp.add_variable(0, std::nullopt);
        auto x3 = lp.add_variable(0, std::nullopt);
        auto w1 = lp.add_constraint(std::nullopt, 5, std::initializer_list<Entry>{{x1, 2}, {x2, 3}, {x3, 1}});
        auto w2 = lp.add_constraint(std::nullopt, 11, std::initializer_list<Entry>{{x1, 4}, {x2, 1}, {x3, 2}});
        auto w3 = lp.add_constraint(std::nullopt, 8, std::initializer_list<Entry>{{x1, 3}, {x2, 4}, {x3, 2}});
        auto inexact = lp.to_imprecise_model<mwt::GurobiLPBackend>(objective_coefficients);
        lp.solve_imprecise_model(inexact);
        lp.imprecise_load_basis(inexact);
        lp.compute_basic_values();
        const auto &basis = lp.basis_column_indices();
        CHECK(basis.size() == 3);
        CHECK(basis[0] == x1);
        CHECK(basis[1] == x3);
        CHECK(basis[2] == 3 + w2);
        const auto &basic_values = lp.basic_values();
        CHECK(basic_values.size() == 3);
        CHECK(basic_values[0] == 2);
        CHECK(basic_values[1] == 1);
        CHECK(basic_values[2] == 10);
        lp.compute_interval_objective_coefficients(inexact);
        const auto &pdi = lp.compute_index_in_nonbasis_of_potential_dual_infeasibilities();
        CHECK(pdi.empty());
        lp.compute_symbolic_objective_coefficients();
        const auto &obj_coeffs = lp.nonbasic_objective_coefficients();
        const auto &obj_intervals = lp.nonbasic_objective_intervals();
        CHECK(obj_coeffs.size() == 3);
        CHECK(obj_intervals.size() == 3);
        CHECK(lp.nonbasis_column_indices().size() == 3);
        CHECK(lp.nonbasis_column_indices() == std::vector<std::size_t>{x2, 3 + w1, 3 + w3});
        auto coefficient_times_value = [](const NT &coeff, const LENT &value) -> LENT {
            return (coeff * CGAL::exact(value)).cast_to<LENT>();
        };
        CHECK(substitute_symbols(obj_coeffs[0], objective_coefficients, coefficient_times_value) == -3);
        CHECK(substitute_symbols(obj_coeffs[1], objective_coefficients, coefficient_times_value) == 1);
        CHECK(substitute_symbols(obj_coeffs[2], objective_coefficients, coefficient_times_value) == 1);
        CHECK(std::round(obj_intervals[0].inf()) == -3);
        CHECK(std::round(obj_intervals[1].inf()) == 1);
        CHECK(std::round(obj_intervals[2].inf()) == 1);
        CHECK(obj_intervals[0].inf() <= -3);
        CHECK(obj_intervals[0].sup() >= -3);
        CHECK(obj_intervals[1].inf() <= 1);
        CHECK(obj_intervals[1].sup() >= 1);
        CHECK(obj_intervals[2].inf() <= 1);
        CHECK(obj_intervals[2].sup() >= 1);
    }

    SUBCASE("known LP, >= and <= constraints") {
        // max x1 + 2x2 - x3
        // w1 = x1 + x2 >= 7
        // w2 = x1 + x3 >= 1
        // w3 = x2 - x3 <= 4
        // 0 <= x1 <= 5, 0 <= x2 <= 4, -4 <= x3 <= 3
        // optimal basis: x3, w1, w2
        // x1 = 5, x2 = 4, x3 = 0
        // w1 = 9, w2 = 5, w3 = 4
        std::vector<LENT> objective_coefficients2{1, 2, -1};
        LP lp2(false);
        auto x1 = lp2.add_variable(0, 5);
        auto x2 = lp2.add_variable(0, 4);
        auto x3 = lp2.add_variable(-4, 3);
        auto w1 = lp2.add_constraint(7, std::nullopt, std::initializer_list<Entry>{{x1, 1}, {x2, 1}});
        auto w2 = lp2.add_constraint(1, std::nullopt, std::initializer_list<Entry>{{x1, 1}, {x3, 1}});
        auto w3 = lp2.add_constraint(std::nullopt, 4, std::initializer_list<Entry>{{x2, 1}, {x3, -1}});
        auto inexact2 = lp2.to_imprecise_model<mwt::GurobiLPBackend>(objective_coefficients2);
        lp2.solve_imprecise_model(inexact2);
        lp2.imprecise_load_basis(inexact2);
        lp2.compute_basic_values();
        const auto &basis2 = lp2.basis_column_indices();
        CHECK(basis2.size() == 3);
        CHECK(basis2[0] == x3);
        CHECK(basis2[1] == 3 + w1);
        CHECK(basis2[2] == 3 + w2);
        const auto &basic_values2 = lp2.basic_values();
        CHECK(basic_values2.size() == 3);
        CHECK(basic_values2[0] == 0);
        CHECK(basic_values2[1] == 9);
        CHECK(basic_values2[2] == 5);
        lp2.compute_interval_objective_coefficients(inexact2);
        CHECK(lp2.compute_index_in_nonbasis_of_potential_dual_infeasibilities().empty());
    }

    SUBCASE("known LP, >= and == constraints") {
        // max x1 - x2 - x3
        // w1 = x1 + x2 == 4
        // w2 = x1 + x3 >= 3
        // w3 = x2 + x3 == 5
        // x1 >= 0, x2 >= 0, x3 >= 0
        // opt: x1 = 4, x2 = 0, x3 = 5
        std::vector<LENT> objective_coefficients3{1, -1, -1};
        LP lp3(false);
        auto x1 = lp3.add_variable(0, std::nullopt);
        auto x2 = lp3.add_variable(0, std::nullopt);
        auto x3 = lp3.add_variable(0, std::nullopt);
        auto w1 = lp3.add_constraint(4, 4, std::initializer_list<Entry>{{x1, 1}, {x2, 1}});
        auto w2 = lp3.add_constraint(3, std::nullopt, std::initializer_list<Entry>{{x1, 1}, {x3, 1}});
        auto w3 = lp3.add_constraint(5, 5, std::initializer_list<Entry>{{x2, 1}, {x3, 1}});
        auto inexact3 = lp3.to_imprecise_model<mwt::GurobiLPBackend>(objective_coefficients3);
        lp3.solve_imprecise_model(inexact3);
        auto inexact_values = inexact3.get_variable_values();
        CHECK(std::abs(inexact_values[x1] - 4.0) < 0.0001);
        CHECK(std::abs(inexact_values[x2]) < 0.0001);
        CHECK(std::abs(inexact_values[x3] - 5.0) < 0.0001);
        lp3.imprecise_load_basis(inexact3);
        lp3.compute_basic_values();
        const auto &basis3 = lp3.basis_column_indices();
        CHECK(basis3.size() == 3);
        CHECK(basis3[0] == x1);
        CHECK(basis3[1] == x3);
        CHECK(basis3[2] == 3 + w2);
        const auto &basic_values3 = lp3.basic_values();
        CHECK(basic_values3.size() == 3);
        CHECK(basic_values3[0] == 4);
        CHECK(basic_values3[1] == 5);
        CHECK(basic_values3[2] == 9);
        CHECK(lp3.solution_value(x1) == 4);
        CHECK(lp3.solution_value(x2) == 0);
        CHECK(lp3.solution_value(x3) == 5);
        CHECK(lp3.solution_value(3 + w1) == 4);
        CHECK(lp3.solution_value(3 + w2) == 9);
        CHECK(lp3.solution_value(3 + w3) == 5);
        lp3.compute_interval_objective_coefficients(inexact3);
        const auto &pdi = lp3.compute_index_in_nonbasis_of_potential_dual_infeasibilities();
        CHECK(pdi.empty());
    }
}

TEST_CASE("[RationalSparseLP] primal general simplex test") {
    using InnerNT = CGAL::Exact_rational;
    using NT = mwt::RationalOrInt<InnerNT>;
    using LP = mwt::RationalSparseLPWithSymbolicObjective<NT>;
    using Entry = LP::Entry;

    SUBCASE("known general LP, with steps") {
        LP lp(false);
        auto x1 = lp.add_variable(-2, std::nullopt);
        auto x2 = lp.add_variable(0, 6);
        auto w1 = 2 + lp.add_constraint(1, 5, std::initializer_list<Entry>{{x1, -1}, {x2, 1}});
        auto w2 = 2 + lp.add_constraint(2, 10, std::initializer_list<Entry>{{x1, -3}, {x2, 2}});
        auto w3 = 2 + lp.add_constraint({}, 0, std::initializer_list<Entry>{{x1, 2}, {x2, -1}});
        std::vector<mwt::BasicStatus> var_basis{mwt::BasicStatus::LOWER_BOUND, mwt::BasicStatus::LOWER_BOUND};
        std::vector<mwt::BasicStatus> con_basis{mwt::BasicStatus::BASIC, mwt::BasicStatus::BASIC,
                                                mwt::BasicStatus::BASIC};
        std::vector<NT> objective_coefficients{3, -1};
        auto inexact = lp.to_imprecise_model<mwt::NullLPBackend>(objective_coefficients);
        lp.load_basis(var_basis, con_basis);
        lp.compute_basic_values();
        const auto &basis = lp.basis_column_indices();
        CHECK(basis == std::vector<std::size_t>{w1, w2, w3});
        const auto &nonbasis = lp.nonbasis_column_indices();
        CHECK(nonbasis == std::vector<std::size_t>{x1, x2});
        CHECK(lp.is_primal_feasible());
        const auto &basic_values = lp.basic_values();
        CHECK(basic_values == std::vector<NT>{2, 6, -4});
        lp.compute_interval_objective_coefficients(inexact);
        auto pdis = lp.compute_index_in_nonbasis_of_potential_dual_infeasibilities();
        CHECK(pdis.size() == 1);
        CHECK(pdis[0] == 0);
        auto intervals = lp.nonbasic_objective_intervals();
        CHECK(intervals[0].inf() <= 3);
        CHECK(intervals[0].sup() >= 3);
        CHECK(intervals[1].inf() <= -1);
        CHECK(intervals[1].sup() >= -1);
        lp.compute_symbolic_objective_coefficients();
        const auto &obj_coeffs = lp.nonbasic_objective_coefficients();
        CHECK(obj_coeffs.size() == 2);
        CHECK(substitute_symbols(obj_coeffs[0], objective_coefficients) == 3);
        CHECK(substitute_symbols(obj_coeffs[1], objective_coefficients) == -1);
        lp.primal_simplex_iteration(inexact, pdis);
        CHECK(lp.is_primal_feasible());
        CHECK(lp.basic_values() == std::vector<NT>{-1, 3, -2});
        CHECK(lp.basis_column_indices() == std::vector<std::size_t>{x1, w2, w3});
        CHECK(lp.nonbasis_column_indices() == std::vector<std::size_t>{w1, x2});
        const auto &obj_coeffs2 = lp.nonbasic_objective_coefficients();
        CHECK(&obj_coeffs == &obj_coeffs2);
        CHECK(substitute_symbols(obj_coeffs[0], objective_coefficients) == -3);
        CHECK(substitute_symbols(obj_coeffs[1], objective_coefficients) == 2);
        CHECK(lp.nonbasic_objective_intervals()[0].inf() <= -3);
        CHECK(lp.nonbasic_objective_intervals()[0].sup() >= -3);
        CHECK(lp.nonbasic_objective_intervals()[1].inf() <= 2);
        CHECK(lp.nonbasic_objective_intervals()[1].sup() >= 2);
        pdis = lp.compute_index_in_nonbasis_of_potential_dual_infeasibilities();
        CHECK(pdis.size() == 1);
        CHECK(pdis[0] == 1);
        CHECK(lp.nonbasis_column_indices()[pdis[0]] == x2);
        lp.primal_simplex_iteration(inexact, pdis);
        CHECK(lp.is_primal_feasible());
        CHECK(lp.basis_column_indices() == std::vector<std::size_t>{x1, x2, w3});
        CHECK(lp.nonbasis_column_indices() == std::vector<std::size_t>{w1, w2});
        CHECK(lp.basic_values() == std::vector<NT>{0, 1, -1});
        CHECK(lp.solution_value(w1) == 1);
        CHECK(lp.solution_value(w2) == 2);
        CHECK(lp.nonbasic_objective_intervals()[0].inf() <= 3);
        CHECK(lp.nonbasic_objective_intervals()[0].sup() >= 3);
        CHECK(lp.nonbasic_objective_intervals()[1].inf() <= -2);
        CHECK(lp.nonbasic_objective_intervals()[1].sup() >= -2);
        pdis = lp.compute_index_in_nonbasis_of_potential_dual_infeasibilities();
        CHECK(pdis.size() == 1);
        CHECK(pdis[0] == 0);
        lp.primal_simplex_iteration(inexact, pdis);
        CHECK(lp.is_primal_feasible());
        CHECK(lp.basis_column_indices() == std::vector<std::size_t>{x1, x2, w1});
        CHECK(lp.nonbasis_column_indices() == std::vector<std::size_t>{w3, w2});
        CHECK(lp.basic_values() == std::vector<NT>{2, 4, 2});
        CHECK(lp.solution_value(w3) == 0);
        CHECK(lp.solution_value(w2) == 2);
        CHECK(lp.nonbasic_objective_intervals()[0].inf() <= 3);
        CHECK(lp.nonbasic_objective_intervals()[0].sup() >= 3);
        CHECK(lp.nonbasic_objective_intervals()[1].inf() <= 1);
        CHECK(lp.nonbasic_objective_intervals()[1].sup() >= 1);
        pdis = lp.compute_index_in_nonbasis_of_potential_dual_infeasibilities();
        CHECK(pdis.size() == 1);
        CHECK(pdis[0] == 1);
        lp.primal_simplex_iteration(inexact, pdis);
        CHECK(lp.is_primal_feasible());
        CHECK(lp.basis_column_indices() == std::vector<std::size_t>{x1, w2, w1});
        CHECK(lp.nonbasis_column_indices() == std::vector<std::size_t>{w3, x2});
        CHECK(lp.solution_value(w3) == 0);
        CHECK(lp.solution_value(x2) == 6);
        CHECK(lp.nonbasic_objective_intervals()[0].inf() <= 1.5);
        CHECK(lp.nonbasic_objective_intervals()[0].sup() >= 1.5);
        CHECK(lp.nonbasic_objective_intervals()[1].inf() <= 0.5);
        CHECK(lp.nonbasic_objective_intervals()[1].sup() >= 0.5);
        CHECK(lp.basic_values() == std::vector<NT>{3, 3, 3});
        CHECK(lp.compute_index_in_nonbasis_of_potential_dual_infeasibilities().empty());
    }
}
