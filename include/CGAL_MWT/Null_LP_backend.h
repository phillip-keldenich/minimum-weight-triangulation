#ifndef CGAL_MWT_NULL_LP_BACKEND_H_INCLUDED_
#define CGAL_MWT_NULL_LP_BACKEND_H_INCLUDED_

#include "LP_backends.h"
#include <exception>
#include <limits>
#include <nlohmann/json.hpp>
#include <stdexcept>
#include <vector>

namespace mwt {

/**
 * A null LP backend.
 * It cannot solve anything and can be used when
 * no LP solver is actually needed, and is also
 * useful for testing.
 * Its expressions and variables can be created
 * and manipulated, but all these operations are
 * effectively no-ops.
 * Calling any solving method will throw an exception.
 */
struct NullLPBackend {
    struct Empty {};

    using Model = Empty;
    using Var = Empty;
    using Constr = Empty;
    using LinExpr = Empty;

    /**
     * Return a value that encodes the absence of a bound.
     */
    static double infinite_lb() noexcept { return -std::numeric_limits<double>::infinity(); }
    static double infinite_ub() noexcept { return std::numeric_limits<double>::infinity(); }

    static LinExpr empty_expression() { return {}; }

    static void clear_expression(LinExpr & /*expr*/) {}

    static void add_to_expression(LinExpr & /*expr*/, double /*value*/) {}

    static void add_to_expression(LinExpr & /*expr*/, double /*coefficient*/, const Var & /*var*/) {}

    static void model_construction_done(Model & /*model*/) {}

    static Model create_model(bool /*silent*/) { return {}; }

    static Var create_continuous(Model & /*model*/, double /*lb*/, double /*ub*/) { return {}; }

    static Var create_boolean(Model & /*model*/) { return {}; }

    static Var create_integer(Model & /*model*/, double /*lb*/, double /*ub*/) { return {}; }

    static void make_var_boolean(Model & /*model*/, Var & /*var*/) { return; }

    static Constr add_greater_equal(Model & /*model*/, const LinExpr & /*expr*/, double /*rhs*/) { return {}; }

    static Constr add_equal(Model & /*model*/, const LinExpr & /*expr*/, double /*rhs*/) { return {}; }

    static Constr add_less_equal(Model & /*model*/, const LinExpr & /*expr*/, double /*rhs*/) { return {}; }

    static Constr add_range(Model & /*model*/, const LinExpr & /*expr*/, double /*lhs*/, double /*rhs*/) { return {}; }

    /**
     * Set the objective to minimize the given expression.
     */
    static void minimize(Model & /*model*/, const LinExpr & /*expr*/) {}

    /**
     * Set the objective to maximize the given expression.
     */
    static void maximize(Model & /*model*/, const LinExpr & /*expr*/) {}

    /**
     * Optimize the given model and return
     * the resulting status.
     */
    static LPStatus optimize(Model &model) { throw std::logic_error("NullLPBackend cannot optimize!"); }

    static double objective_value(const Model &model) {
        throw std::logic_error("NullLPBackend cannot get objective value!");
    }

    /**
     * Get solution values of a set of variables into a preallocated vector.
     */
    static void get_solution(Model & /*model*/, const std::vector<Var> & /*vars*/, std::vector<double> & /*values*/) {
        throw std::logic_error("NullLPBackend cannot get solution values!");
    }

    /**
     * Get solution values of a set of variables in a new vector.
     */
    static std::vector<double> get_solution(Model & /*model*/, const std::vector<Var> & /*vars*/) {
        throw std::logic_error("NullLPBackend cannot get solution values!");
    }

    /**
     * Get dual values of a set of constraints into a preallocated vector.
     */
    static void get_dual_values(Model & /*model*/, const std::vector<Constr> & /*constrs*/,
                                std::vector<double> & /*values*/) {
        throw std::logic_error("NullLPBackend cannot get dual values!");
    }

    /**
     * Get dual values of a set of constraints into a new vector.
     */
    static std::vector<double> get_dual_values(Model & /*model*/, const std::vector<Constr> & /*constrs*/) {
        throw std::logic_error("NullLPBackend cannot get dual values!");
    }

    /**
     * Get the basis status for variables or constraints.
     * BasisObject must be Var or Constr, otherwise compilation will fail.
     * BasicStatusType must be BasicStatus or ExtendedBasicStatus, otherwise compilation will fail.
     */
    template<typename BasisObject, typename BasicStatusType>
    static void get_basis_status(Model & /*model*/, const std::vector<BasisObject> & /*objs*/,
                                 std::vector<BasicStatusType> & /*statuses*/) {
        throw std::logic_error("NullLPBackend cannot get basis status!");
    }

    /**
     * Get the basis status for variables or constraints.
     * BasisObject must be Var or Constr, otherwise compilation will fail.
     * If any object is SUPERBASIC in Gurobi (only non-linear models),
     * raise an error.
     */
    template<typename BasisObject>
    static std::vector<BasicStatus> get_basis_status(Model & /*model*/, const std::vector<BasisObject> & /*objs*/) {
        throw std::logic_error("NullLPBackend cannot get basis status!");
    }

    /**
     * Get the basis status for variables or constraints.
     * BasisObject must be Var or Constr, otherwise compilation will fail.
     * Does not raise errors on superbasic objects.
     */
    template<typename BasisObject>
    static std::vector<ExtendedBasicStatus> get_extended_basis_status(Model & /*model*/,
                                                                      const std::vector<BasisObject> & /*objs*/) {
        throw std::logic_error("NullLPBackend cannot get basis status!");
    }

    static nlohmann::json get_solution_process_stats(const Model & /*model*/) { return {}; }

    /**
     * Check whether this backend is available at runtime.
     * Involves checking whether the Gurobi library is available
     * and does not raise errors, e.g., due to license issues.
     */
    static bool test_runtime_availability(bool warn_on_error) { return false; }
};

} // namespace mwt

#endif
