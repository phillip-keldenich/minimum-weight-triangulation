#ifndef CGAL_MWT_GUROBI_LP_BACKEND_H_INCLUDED_
#define CGAL_MWT_GUROBI_LP_BACKEND_H_INCLUDED_

#include "LP_backends.h"

// TODO remove/make optional
#define CGAL_MWT_HAVE_GUROBI 1

#if defined(CGAL_MWT_HAVE_GUROBI)

#include <algorithm>
#include <gurobi_c++.h>
#include <map>
#include <memory>
#include <nlohmann/json.hpp>
#include <vector>

namespace mwt {

struct GurobiLPBackend {
    using Model = GRBModel;
    using Var = GRBVar;
    using Constr = GRBConstr;
    using LinExpr = GRBLinExpr;

    /**
     * Return a value that encodes the absence of a bound.
     */
    static double infinite_lb() noexcept { return -GRB_INFINITY; }
    static double infinite_ub() noexcept { return GRB_INFINITY; }

    static LinExpr empty_expression() { return 0; }

    static void clear_expression(LinExpr &expr) { expr.clear(); }

    static void add_to_expression(LinExpr &expr, double value) { expr += value; }

    static void add_to_expression(LinExpr &expr, double coefficient, const Var &var) { expr += coefficient * var; }

    static void model_construction_done(Model &model) { model.update(); }

    static Model create_model(bool silent) { return Model(p_get_env(silent)); }

    static Var create_continuous(Model &model, double lb, double ub) {
        return model.addVar(lb, ub, 0.0, GRB_CONTINUOUS);
    }

    static Var create_boolean(Model &model) { return model.addVar(0.0, 1.0, 0.0, GRB_BINARY); }

    static Var create_integer(Model &model, double lb, double ub) { return model.addVar(lb, ub, 0.0, GRB_INTEGER); }

    static void make_var_boolean(Model & /*model*/, Var &var) { var.set(GRB_CharAttr_VType, GRB_BINARY); }

    static Constr add_greater_equal(Model &model, const LinExpr &expr, double rhs) {
        return model.addConstr(expr >= rhs);
    }

    static Constr add_equal(Model &model, const LinExpr &expr, double rhs) { return model.addConstr(expr == rhs); }

    static Constr add_less_equal(Model &model, const LinExpr &expr, double rhs) { return model.addConstr(expr <= rhs); }

    static Constr add_range(Model &model, const LinExpr &expr, double lhs, double rhs) {
        return model.addRange(expr, lhs, rhs);
    }

    /**
     * Set the objective to minimize the given expression.
     */
    static void minimize(Model &model, const LinExpr &expr) { model.setObjective(expr, GRB_MINIMIZE); }

    /**
     * Set the objective to maximize the given expression.
     */
    static void maximize(Model &model, const LinExpr &expr) { model.setObjective(expr, GRB_MAXIMIZE); }

    /**
     * Optimize the given model and return
     * the resulting status.
     */
    static LPStatus optimize(Model &model) {
        model.optimize();
        auto status = model.get(GRB_IntAttr_Status);
        switch(status) {
        case GRB_OPTIMAL:
            return LPStatus::OPTIMAL;
        case GRB_INFEASIBLE:
            return LPStatus::INFEASIBLE;
        case GRB_UNBOUNDED:
            return LPStatus::UNBOUNDED;

        default:
        case GRB_LOADED:
        case GRB_INF_OR_UNBD:
        case GRB_SUBOPTIMAL:
        case GRB_INPROGRESS:
        case GRB_CUTOFF:
        case GRB_NUMERIC:
        case GRB_USER_OBJ_LIMIT: {
            std::map<int, std::string> lookup{{GRB_LOADED, "GRB_LOADED"},         {GRB_INF_OR_UNBD, "GRB_INF_OR_UNBD"},
                                              {GRB_SUBOPTIMAL, "GRB_SUBOPTIMAL"}, {GRB_INPROGRESS, "GRB_INPROGRESS"},
                                              {GRB_CUTOFF, "GRB_CUTOFF"},         {GRB_NUMERIC, "GRB_NUMERIC"}};
            if(lookup.count(status)) {
                throw std::runtime_error("Unexpected Gurobi status " + lookup[status]);
            } else {
                throw std::runtime_error("Unknown and unexpected Gurobi status " + std::to_string(status));
            }
        }

        case GRB_NODE_LIMIT:
        case GRB_TIME_LIMIT:
        case GRB_SOLUTION_LIMIT:
        case GRB_INTERRUPTED:
        case GRB_ITERATION_LIMIT:
        case GRB_WORK_LIMIT:
        case GRB_MEM_LIMIT: {
            return LPStatus::TIMEOUT;
        }
        }
    }

    static double objective_value(const Model &model) { return model.get(GRB_DoubleAttr_ObjVal); }

    /**
     * Get solution values of a set of variables into a preallocated vector.
     */
    static void get_solution(Model &model, const std::vector<Var> &vars, std::vector<double> &values) {
        if(vars.size() > std::size_t(std::numeric_limits<int>::max())) {
            throw std::runtime_error("Too many variables to get solution for");
        }
        std::unique_ptr<double[]> data;
        data.reset(model.get(GRB_DoubleAttr_X, vars.data(), static_cast<int>(vars.size())));
        values.assign(data.get(), data.get() + vars.size());
    }

    /**
     * Get solution values of a set of variables in a new vector.
     */
    static std::vector<double> get_solution(Model &model, const std::vector<Var> &vars) {
        std::vector<double> values;
        get_solution(model, vars, values);
        return values;
    }

    /**
     * Get dual values of a set of constraints into a preallocated vector.
     */
    static void get_dual_values(Model &model, const std::vector<Constr> &constrs, std::vector<double> &values) {
        if(constrs.size() > std::size_t(std::numeric_limits<int>::max())) {
            throw std::runtime_error("Too many constraints to get dual values for");
        }
        std::unique_ptr<double[]> data;
        data.reset(model.get(GRB_DoubleAttr_Pi, constrs.data(), static_cast<int>(constrs.size())));
        values.assign(data.get(), data.get() + constrs.size());
    }

    /**
     * Get dual values of a set of constraints into a new vector.
     */
    static std::vector<double> get_dual_values(Model &model, const std::vector<Constr> &constrs) {
        std::vector<double> values;
        get_dual_values(model, constrs, values);
        return values;
    }

    /**
     * Get the basis status for variables or constraints.
     * BasisObject must be Var or Constr, otherwise compilation will fail.
     * BasicStatusType must be BasicStatus or ExtendedBasicStatus, otherwise compilation will fail.
     */
    template<typename BasisObject, typename BasicStatusType>
    static void get_basis_status(Model &model, const std::vector<BasisObject> &objs,
                                 std::vector<BasicStatusType> &statuses) {
        auto data = p_get_bstatus(model, objs);
        statuses.clear();
        statuses.reserve(objs.size());
        std::transform(data.get(), data.get() + objs.size(), std::back_inserter(statuses),
                       [](int status) { return p_convert_status<BasicStatusType>(status); });
    }

    /**
     * Get the basis status for variables or constraints.
     * BasisObject must be Var or Constr, otherwise compilation will fail.
     * If any object is SUPERBASIC in Gurobi (only non-linear models),
     * raise an error.
     */
    template<typename BasisObject>
    static std::vector<BasicStatus> get_basis_status(Model &model, const std::vector<BasisObject> &objs) {
        std::vector<BasicStatus> statuses;
        get_basis_status(model, objs, statuses);
        return statuses;
    }

    /**
     * Get the basis status for variables or constraints.
     * BasisObject must be Var or Constr, otherwise compilation will fail.
     * Does not raise errors on superbasic objects.
     */
    template<typename BasisObject>
    static std::vector<ExtendedBasicStatus> get_extended_basis_status(Model &model,
                                                                      const std::vector<BasisObject> &objs) {
        std::vector<ExtendedBasicStatus> statuses;
        get_basis_status(model, objs, statuses);
        return statuses;
    }

    /**
     * Check whether this backend is available at runtime.
     * Involves checking whether the Gurobi library is available
     * and does not raise errors, e.g., due to license issues.
     */
    static bool test_runtime_availability(bool warn_on_error) {
        try {
            GRBEnv env = p_create_env(true);
            GRBModel model(env);
            GRBVar var = model.addVar(0, 1, 2, GRB_CONTINUOUS);
            model.optimize();
            if(model.get(GRB_IntAttr_Status) != GRB_OPTIMAL) {
                if(warn_on_error) {
                    std::cerr << "Gurobi test optimization did not return optimal status; "
                                 "something seems to be wrong with Gurobi!"
                              << std::endl;
                }
                return false;
            }
            return true;
        } catch(const GRBException &exception) {
            if(warn_on_error) {
                std::cerr << "Gurobi seems to be unavailable: " << exception.getMessage() << std::endl;
            }
            return false;
        }
    }

    static nlohmann::json get_solution_process_stats(const Model &model) {
        return nlohmann::json{{"obj_value", model.get(GRB_DoubleAttr_ObjVal)},
                              {"runtime", model.get(GRB_DoubleAttr_Runtime)},
                              {"simplex_iterations", model.get(GRB_DoubleAttr_IterCount)},
                              {"barrier_iterations", model.get(GRB_IntAttr_BarIterCount)},
                              {"condition_estimate", model.get(GRB_DoubleAttr_Kappa)},
                              {"num_nonzeros", model.get(GRB_IntAttr_NumNZs)},
                              {"num_constraints", model.get(GRB_IntAttr_NumConstrs)},
                              {"num_variables", model.get(GRB_IntAttr_NumVars)}};
    }

  private:
    template<typename BasisObject>
    static std::unique_ptr<int[]> p_get_bstatus(Model &model, const std::vector<BasisObject> &vars) {
        if(vars.size() > std::size_t(std::numeric_limits<int>::max())) {
            throw std::runtime_error("Too many objects to get basis status for");
        }

        std::unique_ptr<int[]> result;
        if constexpr(std::is_same_v<BasisObject, Var>) {
            result.reset(model.get(GRB_IntAttr_VBasis, vars.data(), static_cast<int>(vars.size())));
        } else {
            static_assert(std::is_same_v<BasisObject, Constr>, "BasisObject must be Var or Constr");
            result.reset(model.get(GRB_IntAttr_CBasis, vars.data(), static_cast<int>(vars.size())));
        }
        return result;
    }

    template<typename StatusType> static StatusType p_convert_status(int status) {
        switch(status) {
        case GRB_NONBASIC_LOWER:
            return StatusType::LOWER_BOUND;
        case GRB_NONBASIC_UPPER:
            return StatusType::UPPER_BOUND;
        case GRB_BASIC:
            return StatusType::BASIC;
        case GRB_SUPERBASIC:
            if constexpr(std::is_same_v<StatusType, ExtendedBasicStatus>) {
                return StatusType::SUPERBASIC;
            } else {
                throw std::runtime_error("Unexpected Gurobi basis status SUPERBASIC");
            }
        default:
            throw std::runtime_error("Unknown Gurobi basis status " + std::to_string(status));
        }
    }

    static GRBEnv &p_get_env(bool quiet) {
        thread_local GRBEnv env = p_create_env(quiet);
        return env;
    }

    static GRBEnv p_create_env(bool quiet) {
        GRBEnv env(true);
        env.set(GRB_IntParam_OutputFlag, quiet ? 0 : 1);
        env.start();
        return env;
    }
};

} // namespace mwt

#endif

#endif
