#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL_MWT/Face_analyzer.h>
#include <CGAL_MWT/Inexact_LP_face_triangulator_with_max_gap.h>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <doctest/doctest.h>
#include <nlohmann/json.hpp>
#include <sstream>

extern std::stringstream debug_cover_cut_generator_data();

nlohmann::json debug_cover_cut_json() {
    std::stringstream raw = debug_cover_cut_generator_data();
    boost::iostreams::filtering_istream in;
    in.push(boost::iostreams::bzip2_decompressor());
    in.push(raw);

    nlohmann::json j;
    in >> j;
    return j;
}

template<typename Point_2_> class FakeHalfedge {
  public:
    using Point_2 = Point_2_;

    FakeHalfedge(const Point_2 *src, const Point_2 *tgt) : m_source(src), m_target(tgt) {}

    Point_2 source() const noexcept { return *m_source; }

    const Point_2 *source_handle() const noexcept { return m_source; }

    Point_2 target() const noexcept { return *m_target; }

    const Point_2 *target_handle() const noexcept { return m_target; }

    bool operator<(const FakeHalfedge &other) const noexcept {
        return std::tie(m_source, m_target) < std::tie(other.m_source, other.m_target);
    }

    bool operator==(const FakeHalfedge &other) const noexcept {
        return std::tie(m_source, m_target) == std::tie(other.m_source, other.m_target);
    }

    bool operator!=(const FakeHalfedge &other) const noexcept { return !(*this == other); }

  private:
    const Point_2 *m_source;
    const Point_2 *m_target;
};

TEST_CASE("[CoverCutGenerator] Cover Cut Generator Regression #1") {
    // ---------- loading data and creating mock structures ----------
    auto jsdata = debug_cover_cut_json();
    using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
    using Point_2 = Kernel::Point_2;
    using Halfedge = FakeHalfedge<Point_2>;
    using Triangle = mwt::Face_triangle<Halfedge>;

    std::vector<Point_2> points;
    std::vector<Halfedge> halfedges;
    std::vector<Triangle> triangles;

    auto to_point = [](const nlohmann::json &j) { return Point_2(j[0].get<double>(), j[1].get<double>()); };

    for(const auto &t : jsdata["all_triangles"]) {
        points.push_back(to_point(t[0]));
        points.push_back(to_point(t[1]));
        points.push_back(to_point(t[2]));
    }

    std::sort(points.begin(), points.end());
    points.erase(std::unique(points.begin(), points.end()), points.end());

    auto point_to_handle = [&](const Point_2 &p) {
        auto it = std::lower_bound(points.begin(), points.end(), p);
        if(it == points.end() || *it != p) {
            throw std::runtime_error("Point not found");
        }
        return &(*it);
    };

    for(const auto &t : jsdata["all_triangles"]) {
        halfedges.emplace_back(point_to_handle(to_point(t[0])), point_to_handle(to_point(t[1])));
        halfedges.emplace_back(point_to_handle(to_point(t[1])), point_to_handle(to_point(t[0])));
        halfedges.emplace_back(point_to_handle(to_point(t[1])), point_to_handle(to_point(t[2])));
        halfedges.emplace_back(point_to_handle(to_point(t[2])), point_to_handle(to_point(t[1])));
        halfedges.emplace_back(point_to_handle(to_point(t[0])), point_to_handle(to_point(t[2])));
        halfedges.emplace_back(point_to_handle(to_point(t[2])), point_to_handle(to_point(t[0])));
    }

    std::sort(halfedges.begin(), halfedges.end());
    halfedges.erase(std::unique(halfedges.begin(), halfedges.end()), halfedges.end());

    auto halfedge_to_handle = [&](const Halfedge &h) {
        auto it = std::lower_bound(halfedges.begin(), halfedges.end(), h);
        if(it == halfedges.end() || *it != h) {
            throw std::runtime_error("Halfedge not found");
        }
        return &(*it);
    };

    std::vector<double> all_triangle_values;
    for(const auto &t : jsdata["all_triangles"]) {
        triangles.emplace_back();
        triangles.back().edges[0] =
            halfedge_to_handle(Halfedge(point_to_handle(to_point(t[0])), point_to_handle(to_point(t[1]))));
        triangles.back().edges[1] =
            halfedge_to_handle(Halfedge(point_to_handle(to_point(t[1])), point_to_handle(to_point(t[2]))));
        triangles.back().edges[2] =
            halfedge_to_handle(Halfedge(point_to_handle(to_point(t[2])), point_to_handle(to_point(t[0]))));
        all_triangle_values.push_back(0.0);
    }

    std::vector<mwt::detail::FractionalTriangleInfo> fractional_triangles;
    for(const auto &ft : jsdata["fractional_triangles"]) {
        std::size_t tindex = ft["triangle_index"].get<std::size_t>();
        double value = ft["solution_value"].get<double>();
        fractional_triangles.push_back({tindex, value});
        all_triangle_values[tindex] = value;
    }

    // ---------- actual tests ----------
    CHECK(fractional_triangles.size() == 4);
    const Point_2 *center = point_to_handle(to_point(jsdata["center_point"]));
    mwt::detail::Cover_cut_generator<Point_2, Triangle> generator{center, fractional_triangles.data(), &triangles,
                                                                  fractional_triangles.size(), &all_triangle_values};
    CHECK(!generator.graph_is_odd_cycle());
    auto odd_cycle = generator.graph_find_odd_cycle();
    CHECK(odd_cycle.empty());
    CHECK(generator.min_forward_clique_size() == 1);
}
