#ifndef CGAL_MWT_LP_FACE_TRIANGULATOR_H_INCLUDED_
#define CGAL_MWT_LP_FACE_TRIANGULATOR_H_INCLUDED_

#include "Gurobi_LP_backend.h"
#include "Triangle_based_linear_model.h"
#include <CGAL/Exact_rational.h>
#include <CGAL/Triangle_2.h>
#include <boost/container/small_vector.hpp>
#include <boost/container/static_vector.hpp>
#include <nlohmann/json.hpp>
#include <unordered_map>
#include <unordered_set>
#include <variant>

namespace mwt {

namespace detail {

struct FractionalTriangleInfo {
    std::size_t triangle_index;
    double solution_value;
};

template<typename Point_2_, typename Triangle_> class Cover_cut_generator {
  public:
    using Triangle = Triangle_;
    using Point_2 = Point_2_;
    using TriangleIndices = boost::container::small_vector<std::size_t, 5>;
    using TriangleMap = std::unordered_map<std::size_t, TriangleIndices>;

    Cover_cut_generator(const Point_2 *center, const FractionalTriangleInfo *ft_info,
                        const std::vector<Triangle> *all_triangles, std::size_t num_triangles,
                        const std::vector<double> *triangle_values)
        : m_center_vertex(center), m_ft_info(ft_info), m_all_triangles(*all_triangles), m_num_triangles(num_triangles),
          m_triangle_values(*triangle_values) {
        p_triangle_conflict_graph();
        p_make_cyclic_order_map();
    }

    nlohmann::json json_description() {
        nlohmann::json desc;
        auto &all_triangles = desc["all_triangles"];
        for(std::size_t i = 0, num_all_triangles = m_all_triangles.size(); i < num_all_triangles; ++i) {
            const auto &t = m_all_triangles[i];
            nlohmann::json points;
            for(int j = 0; j < 3; ++j) {
                points.push_back(nlohmann::json{t.edges[j]->source().x(), t.edges[j]->source().y()});
            }
            all_triangles.emplace_back(std::move(points));
        }
        desc["center_point"] = nlohmann::json{m_center_vertex->x(), m_center_vertex->y()};
        auto &fractional_triangles = desc["fractional_triangles"];
        for(std::size_t i = 0; i < m_num_triangles; ++i) {
            fractional_triangles.push_back(nlohmann::json{{"triangle_index", m_ft_info[i].triangle_index},
                                                          {"solution_value", m_ft_info[i].solution_value}});
        }
        return desc;
    }

    std::size_t min_forward_clique_size() {
        std::size_t min_clique_size = std::numeric_limits<std::size_t>::max();
        for(const auto &e : m_triangle_events) {
            if(e.is_begin) {
                std::size_t clique_size = p_make_forward_clique(e.triangle_index).size();
                if(clique_size < min_clique_size) {
                    min_clique_size = clique_size;
                }
            }
        }
        return min_clique_size;
    }

    double max_forward_clique_weight() const noexcept { return m_max_forward_clique_weight; }

    const std::vector<std::size_t> &max_weight_forward_clique() const noexcept { return m_max_weight_forward_clique; }

    bool graph_is_odd_cycle() {
        if(m_num_edges % 2 == 0) {
            return false;
        }
        if(!std::all_of(m_conflict_neighbors.begin(), m_conflict_neighbors.end(),
                        [](const auto &v) { return v.second.size() == 2; })) {
            return false;
        }
        std::size_t first_triangle = m_ft_info[0].triangle_index;
        std::unordered_set<std::size_t> visited;
        visited.reserve(m_num_edges);
        visited.insert(first_triangle);
        std::size_t current_triangle = first_triangle;
        for(;;) {
            const auto &neighbors = m_conflict_neighbors[current_triangle];
            if(!visited.insert(neighbors[0]).second) {
                if(!visited.insert(neighbors[1]).second) {
                    break;
                }
                current_triangle = neighbors[1];
            } else {
                current_triangle = neighbors[0];
            }
        }
        return visited.size() == m_num_edges;
    }

    std::vector<std::size_t> graph_find_odd_cycle() {
        std::unordered_map<std::size_t, ParityAndParent> parity_map;
        std::vector<std::size_t> queue;
        for(std::size_t toffs = 0; toffs < m_num_triangles; ++toffs) {
            std::size_t current_root = m_ft_info[toffs].triangle_index;
            if(parity_map.count(current_root))
                continue;
            parity_map[current_root] = ParityAndParent{std::numeric_limits<std::size_t>::max(), 0};
            queue.assign(1, current_root);
            std::size_t queue_pos = 0;
            std::size_t v_candidate = std::numeric_limits<std::size_t>::max();
            std::size_t w_candidate = std::numeric_limits<std::size_t>::max();
            while(queue_pos < queue.size()) {
                std::size_t current_triangle = queue[queue_pos++];
                std::uint8_t current_parity = parity_map[current_triangle].parity;
                for(std::size_t neighbor : m_conflict_neighbors[current_triangle]) {
                    auto [iter, inserted] = parity_map.try_emplace(
                        neighbor, ParityAndParent{current_triangle, std::uint8_t(1 - current_parity)});
                    if(inserted) {
                        queue.push_back(neighbor);
                        continue;
                    }
                    std::uint8_t neighbor_parity = iter->second.parity;
                    if(neighbor_parity == current_parity) {
                        // r --> current -- neighbor --> r is a (possibly non-simple) odd cycle
                        // current and neighbor must be on the same BFS tree level
                        v_candidate = current_triangle;
                        w_candidate = neighbor;
                    }
                }
            }
            if(v_candidate != std::numeric_limits<std::size_t>::max()) {
                return p_make_simple_cycle(current_root, v_candidate, w_candidate, parity_map);
            }
        }
        return {};
    }

  private:
    bool p_is_forward_triangle(std::size_t current_root, std::size_t crb, std::size_t cre, std::size_t triangle) const {
        std::size_t triangle_beg = m_cyclic_order_map.at(triangle);
        if(crb < cre) {
            return triangle_beg >= crb && triangle_beg <= cre;
        } else {
            return triangle_beg >= crb || triangle_beg <= cre;
        }
    }

    bool p_ends_after_root(std::size_t current_root, std::size_t crb, std::size_t cre, std::size_t triangle) const {
        std::size_t triangle_end = m_cyclic_order_ends.at(triangle);
        if(crb < cre) {
            return triangle_end >= cre || triangle_end <= crb;
        } else {
            return triangle_end >= cre && triangle_end <= crb;
        }
    }

    std::vector<std::size_t> p_make_forward_clique(std::size_t start) {
        std::vector<std::size_t> clique;
        double clique_weight = 0.0;
        clique.push_back(start);
        clique_weight += m_triangle_values[start];
        const std::size_t crb = m_cyclic_order_map.at(start);
        const std::size_t cre = m_cyclic_order_ends.at(start);
        for(std::size_t neigh : m_conflict_neighbors.at(start)) {
            if(!p_is_forward_triangle(start, crb, cre, neigh))
                continue;
            if(!p_ends_after_root(start, crb, cre, neigh))
                continue;
            clique.push_back(neigh);
            clique_weight += m_triangle_values[neigh];
        }
        if(clique_weight > m_max_forward_clique_weight) {
            m_max_forward_clique_weight = clique_weight;
            m_max_weight_forward_clique = clique;
        }
        return clique;
    }

    struct ParityAndParent {
        std::size_t parent;
        std::uint8_t parity;
    };

    void p_make_cyclic_order_map() {
        m_cyclic_order_map.clear();
        m_cyclic_order_ends.clear();
        m_cyclic_order_map.reserve(m_num_triangles);
        m_cyclic_order_ends.reserve(m_num_triangles);
        std::size_t o = 0;
        for(const auto &e : m_triangle_events) {
            if(e.is_begin) {
                m_cyclic_order_map[e.triangle_index] = o++;
            } else {
                m_cyclic_order_ends[e.triangle_index] = o;
            }
        }
    }

    std::vector<std::size_t> p_make_simple_cycle(std::size_t bfs_root, std::size_t v, std::size_t w,
                                                 std::unordered_map<std::size_t, ParityAndParent> &parents) {
        std::vector<std::size_t> cycle_result;
        // walk v -> root in parents, erasing entries
        std::size_t current = v;
        while(current != bfs_root) {
            cycle_result.push_back(current);
            current = parents.extract(current).mapped().parent;
        }
        cycle_result.push_back(bfs_root);
        parents.erase(bfs_root);
        std::size_t first_path_length = cycle_result.size();
        // walk from w towards root until hitting a missing element (common ancestor)
        current = w;
        std::size_t common_ancestor;
        for(;;) {
            cycle_result.push_back(current);
            auto it = parents.find(current);
            if(it == parents.end()) {
                common_ancestor = current;
                break;
            }
            current = it->second.parent;
        }
        std::reverse(cycle_result.begin() + first_path_length, cycle_result.end());
        // cycle_result now: v -> [common_ancestor -> root] | common_ancestor -> w
        // remove the part in [...]
        auto it = cycle_result.begin() + first_path_length - 1;
        assert(std::find(cycle_result.begin(), it + 1, common_ancestor) != it + 1);
        while(*it != common_ancestor) {
            --it;
        }
        it = std::move(cycle_result.begin() + first_path_length, cycle_result.end(), it);
        cycle_result.erase(it, cycle_result.end());
        if(!(cycle_result.size() % 2)) {
            throw std::logic_error("Odd cycle construction yielded even length!");
        }
        return cycle_result;
    }

    struct TriangleEvent {
        std::size_t triangle_index;
        const Point_2 *other_point;
        bool is_begin;
    };

    void p_add_to_active_set(std::size_t triangle_index) {
        for(std::size_t i : m_active_triangles) {
            p_add_conflict_edge(i, triangle_index);
        }
        m_active_triangles.insert(triangle_index);
    }

    void p_create_triangle_events(const Point_2 *op1, const Point_2 *op2, std::size_t triangle_index) {
        if(CGAL::orientation(*m_center_vertex, *op1, *op2) == CGAL::LEFT_TURN) {
            std::swap(op1, op2);
        }
        bool spans_begin = (op1->x() < m_center_vertex->x()) && (op2->x() >= m_center_vertex->x());
        if(spans_begin) {
            m_active_triangles.insert(triangle_index);
        }
        m_triangle_events.push_back({triangle_index, op1, true});
        m_triangle_events.push_back({triangle_index, op2, false});
    }

    void p_create_triangle_events() {
        m_triangle_events.clear();
        m_triangle_events.reserve(2 * m_num_triangles);
        for(std::size_t i = 0; i < m_num_triangles; ++i) {
            std::size_t triangle_index = m_ft_info[i].triangle_index;
            const Triangle &triangle = m_all_triangles[triangle_index];
            const Point_2 *other_points[2];
            for(int j = 0, i = 0; j < 3; ++j) {
                const Point_2 *other_point = triangle.edges[j]->source_handle();
                if(other_point == m_center_vertex)
                    continue;
                other_points[i++] = other_point;
            }
            p_create_triangle_events(other_points[0], other_points[1], triangle_index);
        }
    }

    void p_sort_triangle_events() {
        const auto center_x = m_center_vertex->x();
        const auto center_y = m_center_vertex->y();
        auto x_partition = [&](const TriangleEvent &e) { return e.other_point->x() > center_x; };
        auto right_y_partition = [&](const TriangleEvent &e) { return e.other_point->y() >= center_y; };
        auto left_y_partition = [&](const TriangleEvent &e) { return e.other_point->y() < center_y; };
        auto left_begin = std::partition(m_triangle_events.begin(), m_triangle_events.end(), x_partition);
        auto right_bot_begin = std::partition(m_triangle_events.begin(), left_begin, right_y_partition);
        auto left_top_begin = std::partition(left_begin, m_triangle_events.end(), left_y_partition);
        auto compare_events = [this](const TriangleEvent &e1, const TriangleEvent &e2) {
            if(e1.other_point == e2.other_point) {
                return e1.is_begin < e2.is_begin;
            }
            return CGAL::right_turn(*m_center_vertex, *e1.other_point, *e2.other_point);
        };
        std::sort(m_triangle_events.begin(), right_bot_begin, compare_events);
        std::sort(right_bot_begin, left_begin, compare_events);
        std::sort(left_begin, left_top_begin, compare_events);
        std::sort(left_top_begin, m_triangle_events.end(), compare_events);
    }

    void p_triangle_conflict_graph() {
        p_create_triangle_events();
        p_sort_triangle_events();
        for(const TriangleEvent &e : m_triangle_events) {
            if(e.is_begin) {
                p_add_to_active_set(e.triangle_index);
            } else {
                m_active_triangles.erase(e.triangle_index);
            }
        }
    }

    void p_add_conflict_edge(std::size_t i, std::size_t j) {
        m_conflict_neighbors[i].push_back(j);
        m_conflict_neighbors[j].push_back(i);
        ++m_num_edges;
    }

    const Point_2 *m_center_vertex;
    const FractionalTriangleInfo *m_ft_info;
    const std::vector<Triangle> &m_all_triangles;
    std::size_t m_num_triangles;
    const std::vector<double> &m_triangle_values;
    TriangleMap m_conflict_neighbors;
    std::unordered_map<std::size_t, std::size_t> m_cyclic_order_map;
    std::unordered_map<std::size_t, std::size_t> m_cyclic_order_ends;
    std::vector<TriangleEvent> m_triangle_events;
    std::unordered_set<std::size_t> m_active_triangles;
    std::size_t m_num_edges = 0;
    double m_max_forward_clique_weight = 0.0;
    std::vector<std::size_t> m_max_weight_forward_clique;
};

template<typename Skeleton_, typename LPBackend_> class Inexact_nonsimple_face_triangulator_with_max_gap {
  public:
    using LPBackend = LPBackend_;
    using Model = typename LPBackend::Model;
    using Var = typename LPBackend::Var;
    using Constr = typename LPBackend::Constr;
    using LinExpr = typename LPBackend::LinExpr;
    using Skeleton = Skeleton_;
    using Analyzer = Face_analyzer<Skeleton>;
    using Face = mwt::Face<Skeleton>;
    using HalfedgeHandle = const typename Skeleton::Halfedge *;
    using Interval = CGAL::Interval_nt_advanced;
    using HVec = boost::container::static_vector<std::size_t, 3>;
    using Point_2 = typename Skeleton::Point_2;
    using Triangle = typename Analyzer::Triangle;

    using FractionalTriangles = std::vector<FractionalTriangleInfo>;
    using FractionalTrianglesMap = std::unordered_map<const Point_2 *, FractionalTriangles>;

    explicit Inexact_nonsimple_face_triangulator_with_max_gap(const Analyzer *analyzer, bool silent, bool cuts)
        : m_analyzer(analyzer), m_model(LPBackend::create_model(silent)), m_use_cuts(cuts) {}

    /**
     * Error raised when the solution is not integral.
     */
    struct NonIntegral : public std::runtime_error {
        NonIntegral() : std::runtime_error("Gurobi_LP_nonsimple_face_triangulator: non-integral solution") {}
    };

    /**
     * Build the model.
     */
    void build_model() {
        p_build_model();
        m_model.update();
    }

    /**
     * Solve the linear model.
     * If the relaxed solution is not integral, raises a NonIntegral error.
     * If the relaxed solution is integral, updates the status of the edges
     * in the face to LMTStatus::Certain and LMTStatus::Impossible;
     * does not remove/unlink edges, as this would be impossible to do
     * in a parallelized manner.
     */
    void solve() {
        m_var_values.clear();
        m_dual_values.clear();
        if(LPBackend::optimize(m_model) != LPStatus::OPTIMAL) {
            throw std::runtime_error("Inexact_nonsimple_face_triangulator_with_max_gap: optimization failed");
        }
        p_buffer_solution();
        if(p_check_integral()) {
            // solution is integral;
            // round & verify that it is truly feasible
            p_check_rounded_solution();
            m_compute_gap_time = measure_time([&]() { p_compute_maximum_gap(); });
        } else {
            // store the dual values for later; when
            // we move to MIP, they will no longer be available
            p_dual_values();
            if(m_use_cuts) {
                std::size_t iter = 0;
                auto search_cuts = [&]() {
                    bool res;
                    m_cut_search_time += measure_time([&]() { res = p_search_odd_wheel_cuts(); });
                    return res;
                };
                m_cut_phase_time = measure_time([&]() {
                    while(search_cuts()) {
                        if(LPBackend::optimize(m_model) != LPStatus::OPTIMAL) {
                            throw std::runtime_error(
                                "Inexact_nonsimple_face_triangulator_with_max_gap: post-cut optimization failed");
                        }
                        p_buffer_solution();
                        if(p_check_integral()) {
                            m_integral_by_cuts = true;
                            p_check_rounded_solution();
                            m_compute_gap_time = measure_time([&]() { p_compute_maximum_gap(); });
                            p_integral_solution_to_status();
                            break;
                        }
                        if(++iter > 10)
                            break;
                    }
                });
            }
            if(!m_integral_by_cuts) {
                m_mip_phase_time = measure_time([&]() { p_solve_as_mip(); });
            }
        }
        p_integral_solution_to_status();
    }

    /**
     * Get the maximum gap between our solution
     * and a truly (mathematically) optimal solution.
     */
    double get_maximum_gap() const noexcept {
        CGAL_assertion(!m_dual_values.empty());
        return m_maximum_gap;
    }

    double get_gap_time() const noexcept { return m_compute_gap_time; }

    bool get_integral_by_cuts() const noexcept { return m_integral_by_cuts; }

    double get_cut_search_time() const noexcept { return m_cut_search_time; }

    double get_cut_phase_time() const noexcept { return m_cut_phase_time; }

    double get_mip_phase_time() const noexcept { return m_mip_phase_time; }

  private:
    nlohmann::json p_solution_description() {
        nlohmann::json desc;
        auto &triangles = desc["triangles"];
        for(std::size_t i = 0; i < m_triangle_vars.size(); ++i) {
            const auto &t = m_analyzer->empty_triangles()[i];
            triangles.push_back(nlohmann::json{
                {"points",
                 {{CGAL::to_double(t.edges[0]->source().x()), CGAL::to_double(t.edges[0]->source().y())},
                  {CGAL::to_double(t.edges[1]->source().x()), CGAL::to_double(t.edges[1]->source().y())},
                  {CGAL::to_double(t.edges[2]->source().x()), CGAL::to_double(t.edges[2]->source().y())}}},
                {"relaxed_value", m_var_values[i]}});
        }
        return desc;
    }

    bool p_compute_forward_clique_based_cut() {
        std::size_t max_frac_count = 0;
        auto max_frac_count_iter = m_incident_ft.end();
        for(auto it = m_incident_ft.begin(), e = m_incident_ft.end(); it != e; ++it) {
            std::size_t fcount = it->second.size();
            if(fcount > max_frac_count) {
                max_frac_count = fcount;
                max_frac_count_iter = it;
            }
        }
        if(max_frac_count < 3) {
            return false;
        }
        detail::Cover_cut_generator<Point_2, Triangle> generator{
            max_frac_count_iter->first, max_frac_count_iter->second.data(), &m_analyzer->empty_triangles(),
            max_frac_count, &m_var_values};
        std::size_t fcs = generator.min_forward_clique_size();
        std::size_t max_num_selected = max_frac_count / fcs;
        double selected_weight = 0.0;
        for(const auto &info_entry : max_frac_count_iter->second) {
            selected_weight += info_entry.solution_value;
        }
        if(selected_weight - 0.05 > max_num_selected) {
            auto expr = LPBackend::empty_expression();
            for(const auto &info_entry : max_frac_count_iter->second) {
                LPBackend::add_to_expression(expr, 1.0, m_triangle_vars[info_entry.triangle_index]);
            }
            LPBackend::add_less_equal(m_model, expr, max_num_selected);
            return true;
        }
        double max_weight = generator.max_forward_clique_weight();
        if(max_weight >= 1.05) {
            auto expr = LPBackend::empty_expression();
            for(std::size_t ti : generator.max_weight_forward_clique()) {
                LPBackend::add_to_expression(expr, 1.0, m_triangle_vars[ti]);
            }
            LPBackend::add_less_equal(m_model, expr, 1);
            return true;
        }
        return false;
    }

    bool p_search_odd_wheel_cuts() {
        bool found = false;
        p_compute_vertices_with_incident_fractional_triangles();
        if(p_compute_forward_clique_based_cut()) {
            return true;
        }
        for(const auto &e : m_incident_ft) {
            if(e.second.size() < 3)
                continue;
            detail::Cover_cut_generator<Point_2, Triangle> generator{
                e.first, e.second.data(), &m_analyzer->empty_triangles(), e.second.size(), &m_var_values};
            if(generator.graph_is_odd_cycle()) {
                double sum_weight = 0;
                for(const auto &info_entry : e.second) {
                    sum_weight += info_entry.solution_value;
                }
                if(sum_weight > e.second.size() / 2 + 0.05) {
                    found = true;
                    auto expr = LPBackend::empty_expression();
                    for(const auto &tinfo : e.second) {
                        LPBackend::add_to_expression(expr, 1.0, m_triangle_vars[tinfo.triangle_index]);
                    }
                    LPBackend::add_less_equal(m_model, expr, e.second.size() / 2);
                }
            }
            if(!found) {
                std::vector<std::size_t> cycle = generator.graph_find_odd_cycle();
                if(!cycle.empty()) {
                    double sum_weight = std::reduce(cycle.begin(), cycle.end(), 0.0,
                                                    [&](double acc, std::size_t ti) { return acc + m_var_values[ti]; });
                    if(sum_weight > cycle.size() / 2 + 0.05) {
                        found = true;
                        auto expr = LPBackend::empty_expression();
                        for(std::size_t ti : cycle) {
                            LPBackend::add_to_expression(expr, 1.0, m_triangle_vars[ti]);
                        }
                        LPBackend::add_less_equal(m_model, expr, cycle.size() / 2);
                    }
                }
            }
        }
        return found;
    }

    void p_add_fractional_triangle(const Point_2 *point, std::size_t triangle_index, double frac_value) {
        auto [iter, _] = m_incident_ft.try_emplace(point);
        iter->second.push_back({triangle_index, frac_value});
    }

    void p_compute_vertices_with_incident_fractional_triangles() {
        m_incident_ft.clear();
        const auto &triangles = m_analyzer->empty_triangles();
        const std::size_t nt = triangles.size();
        for(std::size_t ti = 0; ti < nt; ++ti) {
            double v = m_var_values[ti];
            if(p_is_integral(v))
                continue;
            const auto &t = triangles[ti];
            for(HalfedgeHandle h : t.edges) {
                const Point_2 *source = h->source_handle();
                p_add_fractional_triangle(source, ti, v);
            }
        }
    }

    /**
     * Switch to MIP and solve the problem
     * despite its non-integral relaxation.
     * The numerical gap will be measured to the
     * true dual feasible solution of the root relaxation,
     * so it could be a bit larger.
     */
    void p_solve_as_mip() {
        for(Var &v : m_triangle_vars) {
            LPBackend::make_var_boolean(m_model, v);
        }
        auto status = LPBackend::optimize(m_model);
        if(status != LPStatus::OPTIMAL) {
            throw std::runtime_error("MIP solver failed to find optimal solution!");
        }
        p_buffer_solution();
        if(!p_check_integral()) {
            throw std::runtime_error("MIP solver failed to produce a solution that looks integral!");
        }
        p_check_rounded_solution();
        p_compute_maximum_gap();
    }

    /**
     * Build the LP model.
     */
    void p_build_model() {
        p_add_triangle_vars();
        for(auto *container :
            {&m_analyzer->face().boundary, &m_analyzer->face().hole_boundaries, &m_analyzer->face().inner_edges}) {
            for(HalfedgeHandle h : *container) {
                p_add_halfedge_constraint(h);
                m_halfedges.push_back(h);
            }
        }
    }

    /**
     * Create a variable for each empty triangle.
     */
    void p_add_triangle_vars() {
        CGAL::Protect_FPU_rounding protect;
        const auto &triangles = m_analyzer->empty_triangles();
        m_triangle_vars.reserve(triangles.size());
        m_triangle_weights.reserve(triangles.size());
        m_triangle_lhe_indices.assign(triangles.size(), HVec());
        m_triangle_rhe_indices.assign(triangles.size(), HVec());
        auto obj = LPBackend::empty_expression();
        for(std::size_t i = 0, n = triangles.size(); i < n; ++i) {
            const auto &t = triangles[i];
            m_triangle_weights.push_back(t.template possible_weight<Interval>());
            m_triangle_vars.push_back(LPBackend::create_continuous(m_model, 0.0, LPBackend::infinite_ub()));
            LPBackend::add_to_expression(obj, m_triangle_weights.back().inf(), m_triangle_vars.back());
        }
        LPBackend::minimize(m_model, obj);
    }

    /**
     * Add a constraint for a halfedge, depending on its status.
     */
    void p_add_halfedge_constraint(HalfedgeHandle h) {
        std::size_t he_index = m_constraints.size();
        LPBackend::clear_expression(m_tmp);
        auto add_entry = [&](double coefficient, std::size_t index) {
            LPBackend::add_to_expression(m_tmp, coefficient, m_triangle_vars[index]);
            if(coefficient < 0) {
                m_triangle_rhe_indices[index].push_back(he_index);
            } else {
                m_triangle_lhe_indices[index].push_back(he_index);
            }
        };
        auto finalize = [&](double rhs) { m_constraints.push_back(LPBackend::add_equal(m_model, m_tmp, rhs)); };
        create_halfedge_constraint<double>(*m_analyzer, h, add_entry, finalize);
    }

    /**
     * Get the dual values.
     */
    const std::vector<double> &p_dual_values() {
        if(m_dual_values.empty()) {
            LPBackend::get_dual_values(m_model, m_constraints, m_dual_values);
        }
        return m_dual_values;
    }

    /**
     * Fill the vector of variable values.
     */
    void p_buffer_solution() { LPBackend::get_solution(m_model, m_triangle_vars, m_var_values); }

    /**
     * Check if the solution is (close enough to) integral.
     * Also extracts the rounded solution to m_rounded.
     */
    bool p_check_integral() {
        const std::size_t n = m_triangle_vars.size();
        m_rounded.clear();
        m_rounded.reserve(n);
        for(std::size_t i = 0; i < n; ++i) {
            double v = m_var_values[i];
            if(v < 0.0)
                v = 0.0;
            if(v > 1.0)
                v = 1.0;
            m_var_values[i] = v;
            if(!p_is_integral(v)) {
                return false;
            }
            m_rounded.push_back(std::uint8_t(v > 0.5));
        }
        return true;
    }

    static bool p_is_integral(double v) noexcept { return v <= 1.0e-5 || v >= 1.0 - 1.0e-5; }

    /**
     * Check that the solution rounded to integer values
     * is actually primally feasible in a mathematically rigorous manner.
     */
    void p_check_rounded_solution() {
        for(auto *container :
            {&m_analyzer->face().boundary, &m_analyzer->face().hole_boundaries, &m_analyzer->face().inner_edges}) {
            for(HalfedgeHandle h : *container) {
                p_check_halfedge_constraint(h);
            }
        }
    }

    /**
     * Check a constraint against the rounded solution.
     */
    void p_check_halfedge_constraint(HalfedgeHandle h) {
        if(h->status == LMTStatus::Possible) {
            p_check_possible_halfedge_constraint(h);
        } else {
            p_check_certain_halfedge_constraint(h);
        }
    }

    /**
     * Check a constraint for a 'Possible' halfedge against the rounded solution.
     */
    void p_check_possible_halfedge_constraint(HalfedgeHandle h) {
        if(!h->is_primary())
            return; // we only need to check the primary halfedge

        HalfedgeHandle t = h->twin();
        const auto &analyzer_info = m_analyzer->edge_info();
        const auto &p_info = analyzer_info[h];
        const auto &t_info = analyzer_info[t];
        auto is_chosen = [&](std::size_t tind) { return m_rounded[tind] == 1; };
        auto p_chosen = std::count_if(p_info.left_triangles.begin(), p_info.left_triangles.end(), is_chosen);
        auto t_chosen = std::count_if(t_info.left_triangles.begin(), t_info.left_triangles.end(), is_chosen);
        if(p_chosen != t_chosen) {
            throw std::logic_error(
                "The rounded solution violates a constraint for a possible halfedge (sides unequal)!");
        }
        if(p_chosen > 1) {
            throw std::logic_error(
                "The rounded solution violates a constraint for a possible halfedge (too many triangles chosen)!");
        }
    }

    /**
     * Check a constraint for a 'Certain' halfedge
     * against the rounded solution.
     */
    void p_check_certain_halfedge_constraint(HalfedgeHandle h) {
        const auto &analyzer_info = m_analyzer->edge_info();
        const auto &h_info = analyzer_info[h];
        auto is_chosen = [&](std::size_t tind) { return m_rounded[tind] == 1; };
        auto begin = h_info.left_triangles.begin();
        auto end = h_info.left_triangles.end();
        auto first = std::find_if(begin, end, is_chosen);
        if(first == end || std::any_of(first + 1, end, is_chosen)) {
            throw std::logic_error("The rounded solution violates a constraint for a certain halfedge!");
        }
    }

    /**
     * After solving and obtaining an integral solution,
     * get the maximum gap between the solution's objective value
     * and the mathematically true value of an optimal solution,
     * including all possible rounding errors.
     * This is done by creating a truly dual feasible solution in interval
     * arithmetic by modifying the values of the dual variables
     * corresponding to the <= 1 constraints.
     */
    void p_compute_maximum_gap() {
        CGAL::Protect_FPU_rounding protect;
        Interval sol_weight = p_solution_weight();
        const auto &dual = p_dual_values();
        Interval dual_sum = p_dual_certain_sum();
        for(std::size_t ti = 0, n = m_triangle_vars.size(); ti < n; ++ti) {
            Interval lhs = 0;
            for(std::size_t hi : m_triangle_lhe_indices[ti]) {
                lhs += dual[hi];
            }
            for(std::size_t hi : m_triangle_rhe_indices[ti]) {
                lhs -= dual[hi];
            }
            Interval rhs = m_triangle_weights[ti];
            rhs -= lhs;
            if(rhs.inf() < 0.0) {
                // the dual constraint is possibly violated;
                // this can be either due to imprecision or because
                // Gurobi's solution actually uses a dual value >0
                // for the dual variable of ti's <=1 constraint
                dual_sum += rhs.inf(); // reduce the dual objective value by the maximum amount of violation
            }
        }
        // now, dual sum is the weight of a truly dually feasible solution!
        // the 0.5 is because our objective is twice the actual objective.
        Interval gap = 0.5 * (sol_weight - dual_sum);
        CGAL_assertion(gap.sup() >= 0.0);
        m_maximum_gap = gap.sup();
    }

    /**
     * Compute 2 * the weight of all possible
     * edges in the (rounded, integral, checked) solution.
     */
    Interval p_solution_weight() {
        CGAL::Protect_FPU_rounding protect;
        Interval solution_weight = 0;
        for(std::size_t i = 0, n = m_triangle_vars.size(); i < n; ++i) {
            if(m_rounded[i]) {
                solution_weight += m_triangle_weights[i];
            }
        }
        return solution_weight;
    }

    /**
     * Compute the sum of the dual values of all
     * certain halfedges; this is not the dual objective value
     * since some dual constraints only hold due to variable
     * bound <=1-dual variables.
     */
    Interval p_dual_certain_sum() {
        Interval dual_sum = 0;
        const auto &dual_values = p_dual_values();
        for(std::size_t i = 0, n = m_constraints.size(); i < n; ++i) {
            if(m_halfedges[i]->status == LMTStatus::Certain || m_halfedges[i]->status == LMTStatus::CH) {
                dual_sum += dual_values[i];
            }
        }
        return dual_sum;
    }

    /**
     * Use an integral solution to update
     * the edge status of edges in the solution
     * to 'Certain', followed by flagging
     * all remaining 'Possible' edges as 'Impossible'.
     */
    void p_integral_solution_to_status() {
        p_mark_chosen_as_certain();
        mark_remaining_as_impossible(m_analyzer->face().inner_edges);
    }

    /**
     * Mark all edges of chosen triangles as certain.
     */
    void p_mark_chosen_as_certain() {
        CGAL_assertion(m_rounded.size() == m_triangle_vars.size());
        const auto &triangles = m_analyzer->empty_triangles();
        for(std::size_t i = 0, nt = m_triangle_vars.size(); i < nt; ++i) {
            if(m_rounded[i] != 1)
                continue;
            for(auto *h : triangles[i].edges) {
                if(h->status == LMTStatus::Possible) {
                    h->status = LMTStatus::Certain;
                    h->twin()->status = LMTStatus::Certain;
                }
            }
        }
    }

    Model m_model;                            //< Inexact LP model
    const Analyzer *m_analyzer;               //< Face analyzer for the face to triangulate
    std::vector<HVec> m_triangle_lhe_indices; //< Indices of halfedges having each triangle on the left
    std::vector<HVec> m_triangle_rhe_indices; //< Indices of halfedges having each triangle on the right
    std::vector<Var> m_triangle_vars;         //< variables for the triangles
    std::vector<Interval> m_triangle_weights; //< weights for the triangles
    std::vector<Constr> m_constraints;        //< all constraints
    std::vector<HalfedgeHandle> m_halfedges;  //< all halfedges, parallel to m_constraints
    std::vector<std::uint8_t> m_rounded;      //< buffer for rounded values
    LinExpr m_tmp;                            //< buffer for expressions
    std::vector<double> m_var_values;         //< buffer for variable values
    std::vector<double> m_dual_values;        //< buffer for dual values
    FractionalTrianglesMap m_incident_ft;     //< fractional triangles incident to each vertex
    double m_maximum_gap;                     //< maximum gap between the solution and the true optimal solution
    double m_compute_gap_time;                //< time spent computing the gap
    double m_cut_search_time = 0.0;           //< time spent searching for cuts
    double m_cut_phase_time = 0.0;            //< total time spent in the cut phase
    double m_mip_phase_time = 0.0;
    bool m_integral_by_cuts = false; //< whether the integral solution was found by cuts
    bool m_use_cuts;                 //< use custom cuts before resorting to MIP
};

} // namespace detail

/**
 * Face triangulator that can handle most non-simple faces (so long as the
 * LP relaxation of the problem modeled using triangle variables and edge constraints
 * is integral, which is the case for most).
 * Otherwise, it will throw an exception of type NonIntegral; then, the face will have
 * be triangulated using a different method (usually, a MIP).
 */
template<typename Skeleton_, typename LPBackend_> class Inexact_LP_face_triangulator_with_max_gap {
  public:
    static constexpr bool accepts_silent_flag = true;

    explicit Inexact_LP_face_triangulator_with_max_gap(FaceTriangulatorOptions options)
        : m_lp_quiet(!options.lp_verbose), m_use_cuts(options.use_cuts) {}

    using Skeleton = Skeleton_;
    using LPBackend = LPBackend_;
    using Face = mwt::Face<Skeleton>;

    void operator()(Face &&face) {
        m_face_analysis_time = measure_time([&]() { m_analyzer.analyze_face(face); });
        detail::Inexact_nonsimple_face_triangulator_with_max_gap<Skeleton, LPBackend> impl(&m_analyzer, m_lp_quiet,
                                                                                           m_use_cuts);
        m_model_build_type = measure_time([&]() { impl.build_model(); });
        m_solve_time = measure_time([&]() { impl.solve(); });
        m_solve_time -= impl.get_gap_time();
        m_gap_time = impl.get_gap_time();
        m_gap = impl.get_maximum_gap();
        m_cut_search_time = impl.get_cut_search_time();
        m_cut_phase_time = impl.get_cut_phase_time();
        m_mip_phase_time = impl.get_mip_phase_time();
        m_integral_by_cuts = impl.get_integral_by_cuts();
    }

    double get_maximum_gap() const noexcept { return m_gap; }

    void combine_stats(nlohmann::json &overall_stats) {
        nlohmann::json &stat_list = overall_stats["nonsimple_faces"];
        nlohmann::json &stats = stat_list.emplace_back();
        stats["gap"] = m_gap;
        stats["face_analysis_time"] = m_face_analysis_time;
        stats["model_build_time"] = m_model_build_type;
        stats["solve_time"] = m_solve_time;
        stats["gap_time"] = m_gap_time;
        stats["use_cuts"] = m_use_cuts;
        stats["cut_search_time"] = m_cut_search_time;
        stats["cut_phase_time"] = m_cut_phase_time;
        stats["mip_phase_time"] = m_mip_phase_time;
        stats["integral_by_cuts"] = m_integral_by_cuts;
    }

  private:
    Face_analyzer<Skeleton> m_analyzer;
    double m_gap;
    double m_face_analysis_time;
    double m_model_build_type;
    double m_solve_time;
    double m_gap_time;
    double m_cut_search_time;
    double m_cut_phase_time;
    double m_mip_phase_time;
    bool m_lp_quiet;
    bool m_use_cuts;
    bool m_integral_by_cuts;
};

} // namespace mwt

#endif
