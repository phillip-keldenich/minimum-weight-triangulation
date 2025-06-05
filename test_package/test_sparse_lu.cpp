#include <CGAL/Exact_rational.h>
#include <CGAL_MWT/Rational_sparse_matrix.h>
#include <doctest/doctest.h>

TEST_CASE("[RationalSparseLU] simple SparseLU test") {
    using InnerNT = CGAL::Exact_rational;
    using NT = mwt::RationalOrInt<InnerNT>;
    using SparseLU = mwt::RationalSparseLU<NT>;
    using Entry = SparseLU::Entry;

    SUBCASE("known linear system, normal and symbolic") {
        std::vector<std::vector<Entry>> sparse_cols{{{0, NT(2)}, {1, 3}, {2, -1}},
                                                    {{1, 1}, {3, -1}},
                                                    {{0, 4}, {2, -1}, {4, 1}},
                                                    {{1, 1}},
                                                    {{0, -2}, {2, -2}, {3, -6}, {4, 4}}};
        std::vector<NT> rhs{7, -2, 0, 3, 0};
        std::vector<NT> result;
        SparseLU solver(5, SparseLU::ColumnMajor{}, sparse_cols);
        solver.compute_decomposition();
        solver.solve_linear_system(rhs, result);
        CHECK(result.size() == 5);
        CHECK(result[0] == -1);
        CHECK(result[1] == 0);
        CHECK(result[2] == 2);
        CHECK(result[3] == 1);
        CHECK(result[4] == NT(-1, 2));

        using SymEntry = SparseLU::SymbolicEntry;
        mwt::SymbolicCalcTable<NT> table(3, 5);
        std::vector<std::vector<SymEntry>> result_sym;
        solver.solve_linear_system_with_symbols(table, result_sym, {0, 1, {}, 2, {}});
        CHECK(result_sym.size() == 5);

        std::vector<NT> symbol_values{7, -2, 3};
        CHECK(substitute_symbols(result_sym[0], symbol_values) == -1);
        CHECK(substitute_symbols(result_sym[1], symbol_values) == 0);
        CHECK(substitute_symbols(result_sym[2], symbol_values) == 2);
        CHECK(substitute_symbols(result_sym[3], symbol_values) == 1);
        CHECK(substitute_symbols(result_sym[4], symbol_values) == NT(-1, 2));
    }

    SUBCASE("known linear system, transposed, normal and symbolic") {
        std::vector<std::vector<Entry>> sparse_rows{{{0, 2}, {2, 4}, {4, -2}},
                                                    {{0, 3}, {1, 1}, {3, 1}},
                                                    {{0, -1}, {2, -1}, {4, -2}},
                                                    {{1, -1}, {4, -6}},
                                                    {{2, 1}, {4, 4}}};
        std::vector<std::vector<Entry>> copy = sparse_rows;
        std::vector<NT> rhs{7, -2, 0, 3, 0};
        std::vector<NT> rhscopy = rhs;
        std::vector<NT> result, result2;
        SparseLU solver(5, SparseLU::RowMajor{}, sparse_rows);
        solver.compute_decomposition();
        solver.solve_linear_system(rhs, result);
        CHECK(result.size() == 5);
        CHECK(result[0] == -1);
        CHECK(result[1] == 0);
        CHECK(result[2] == 2);
        CHECK(result[3] == 1);
        CHECK(result[4] == NT(-1, 2));
        CHECK(copy == sparse_rows);
        CHECK(rhscopy == rhs);

        SparseLU solver_trans(5, SparseLU::ColumnMajor{}, sparse_rows);
        solver_trans.compute_decomposition();
        solver_trans.solve_transposed_linear_system(rhs, result2);
        CHECK(result2.size() == 5);
        CHECK(result2[0] == -1);
        CHECK(result2[1] == 0);
        CHECK(result2[2] == 2);
        CHECK(result2[3] == 1);
        CHECK(result2[4] == NT(-1, 2));

        using SymEntry = SparseLU::SymbolicEntry;
        mwt::SymbolicCalcTable<NT> table(3, 5);
        std::vector<std::vector<SymEntry>> result_sym;
        solver_trans.solve_transposed_linear_system_with_symbols(table, result_sym, {0, 1, {}, 2, {}});
        CHECK(result_sym.size() == 5);
        std::vector<NT> symbol_values{7, -2, 3};
        CHECK(substitute_symbols(result_sym[0], symbol_values) == -1);
        CHECK(substitute_symbols(result_sym[1], symbol_values) == 0);
        CHECK(substitute_symbols(result_sym[2], symbol_values) == 2);
        CHECK(substitute_symbols(result_sym[3], symbol_values) == 1);
        CHECK(substitute_symbols(result_sym[4], symbol_values) == NT(-1, 2));
    }
}
