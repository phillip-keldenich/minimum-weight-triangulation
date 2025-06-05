#ifndef CGAL_MWT_VALIDATE_H_INCLUDED_
#define CGAL_MWT_VALIDATE_H_INCLUDED_

#include "Crossing_free.h"
#include "Total_weight.h"
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel_with_sqrt.h>
#include <CGAL/FPU.h>
#include <CGAL/Surface_sweep_2_algorithms.h>
#include <CGAL/Uncertain.h>
#include <CGAL/convex_hull_2.h>
#include <boost/config.hpp>
#include <chrono>
#include <climits>
#include <random>

namespace mwt {

/**
 * Quick (and optionally strong, but slow)
 * verifier for MWT computed by the MWT solver process itself.
 */
template<typename Skeleton> class MWTVerifier {
  public:
    using Kernel = typename Skeleton::Kernel;
    using Delaunay = CGAL::Delaunay_triangulation_2<Kernel>;
    using Circulator = typename Delaunay::Vertex_circulator;

    explicit MWTVerifier(const Skeleton *skeleton) : skeleton(skeleton) {}

    /**
     * Get the time it took to compute the
     * Delaunay triangulation.
     */
    double get_delaunay_time() {
        p_delaunay();
        return delaunay_time;
    }

    /**
     * Get the time it took to check for intersections.
     */
    double get_intersection_time() {
        if(intersection_check_time < 0.0) {
            p_check_no_intersections();
        }
        return intersection_check_time;
    }

    /**
     * Do a 'quick' check of the triangulation.
     * This still involves computing a Delaunay triangulation.
     * Checks the number of edges, both in the convex hull
     * and the total number of edges; does not check intersection-freeness.
     */
    void quick_check() {
        p_check_edge_count();
        p_check_total_weight_double_interval();
    }

    /**
     * Do a 'strong' check of the triangulation; basically,
     * everything the quick check does, plus intersection-freeness.
     */
    void strong_check() {
        p_check_edge_count();
        p_check_no_intersections();
        p_check_total_weight_double_interval();
    }

  private:
    /**
     * Check that the total weight of the edges in the MWT
     * is (possibly) at most the total weight of the edges
     * in the Delaunay triangulation (using the given interval type).
     */
    template<typename NT> bool p_check_weights_with_number_type() {
        NT total_mwt_weight = total_weight<NT>(*skeleton);
        NT total_del_weight = total_weight<NT>(p_delaunay());
        CGAL::Uncertain uncertain = (total_mwt_weight <= total_del_weight);
        if(CGAL::certainly_not(uncertain)) {
            throw std::logic_error("Delaunay triangulation definitely has less weight than computed MWT!");
        }
        return CGAL::is_certain(uncertain);
    }

    /**
     * Check that the total weight of the edges in the MWT
     * is (possibly) at most the total weight of the edges
     * in the Delaunay triangulation (in double interval arithmetic).
     * Returns true if the check was definitive.
     */
    bool p_check_total_weight_double_interval() {
        CGAL::Protect_FPU_rounding setrounding;
        return p_check_weights_with_number_type<CGAL::Interval_nt_advanced>();
    }

    /**
     * Check that there are no intersections.
     * This is slow (usually much slower than MWT computation)
     * because CGAL's method requires exact constructions.
     */
    void p_check_no_intersections() {
        auto before = std::chrono::steady_clock::now();
        std::vector<CGAL::Segment_2<Kernel>> segments;
        for(auto he = skeleton->halfedges_begin(); he != skeleton->halfedges_end(); ++he) {
            if(!he->is_primary()) {
                continue;
            }
            if(he->status == LMTStatus::Impossible) {
                continue;
            }
            segments.emplace_back(he->source(), he->target());
        }
        CrossingFreenessVerifier<Kernel> verifier{std::move(segments)};
        if(verifier()) {
            auto s = verifier.get_intersecting_segments();
            std::stringstream msg;
            msg << "Validation error: Segments " << std::setprecision(19) << s.first << " and " << s.second
                << " intersect!";
            throw std::runtime_error(msg.str());
        }
        auto after = std::chrono::steady_clock::now();
        intersection_check_time = std::chrono::duration_cast<std::chrono::duration<double>>(after - before).count();
    }

    /**
     * Check that the number of edges in the triangulation
     * matches the number of edges in the Delaunay triangulation
     * (CGAL's best out-of-the-box method to compute the convex hull including
     * collinear points is to compute the Delaunay triangulation).
     * Also check that the total weight of the Delaunay edges
     * is (possibly, i.e., in interval arithmetic)
     * at most the total weight of our edges.
     */
    void p_check_edge_count() {
        p_count_statuses();
        p_delaunay_count_chull();
        p_delaunay_count_edges();
        if(mwt_num_chull != del_num_chull) {
            throw std::logic_error("Number of convex hull edges in MWT (" + std::to_string(mwt_num_chull) +
                                   ") does not match number of convex hull edges in Delaunay triangulation (" +
                                   std::to_string(del_num_chull) + ")!");
        }
        if(mwt_num_edges != del_num_edges) {
            throw std::logic_error("Number of edges in MWT (" + std::to_string(mwt_num_edges) +
                                   ") does not match number of edges in Delaunay triangulation (" +
                                   std::to_string(del_num_edges) + ")!");
        }
    }

    /**
     * Return a reference to a Delaunay triangulation of the points.
     * If it has not been computed yet, compute it.
     */
    Delaunay &p_delaunay() {
        if(!delaunay) {
            auto before = std::chrono::steady_clock::now();
            const auto &tree = skeleton->get_tree();
            delaunay.emplace(tree.points_begin(), tree.points_end());
            delaunay->finite_edges_begin();
            auto after = std::chrono::steady_clock::now();
            delaunay_time = std::chrono::duration_cast<std::chrono::duration<double>>(after - before).count();
        }
        return *delaunay;
    }

    /**
     * Examine and verify the edge status flags in the MWT.
     */
    void p_count_statuses() {
        mwt_num_chull = mwt_num_edges = 0;
        for(auto he = skeleton->halfedges_begin(), hend = skeleton->halfedges_end(); he != hend; ++he) {
            if(he->status != he->twin()->status) {
                throw std::logic_error("Halfedge and its twin have different status!");
            }
            if(!he->is_primary()) {
                continue;
            }
            switch(he->status) {
            default:
                throw std::logic_error("Unknown halfedge status!");
            case LMTStatus::Impossible:
                continue;
            case LMTStatus::Possible:
                throw std::logic_error("Possible halfedge should not be left after completing triangulation!");
            case LMTStatus::Certain:
                ++mwt_num_edges;
                break;
            case LMTStatus::CH:
                ++mwt_num_edges;
                ++mwt_num_chull;
                break;
            }
        }
    }

    /**
     * Use the Delaunay triangulation to find the number of edges
     * in the convex hull of the points.
     */
    void p_delaunay_count_chull() {
        del_num_chull = 0;
        Circulator vc = p_delaunay().incident_vertices(p_delaunay().infinite_vertex()), done(vc);
        do {
            ++del_num_chull;
            ++vc;
        } while(vc != done);
    }

    /**
     * Count the number of finite edges in the Delaunay triangulation.
     */
    void p_delaunay_count_edges() {
        del_num_edges = std::distance(p_delaunay().finite_edges_begin(), p_delaunay().finite_edges_end());
    }

    const Skeleton *skeleton;
    std::optional<Delaunay> delaunay;
    double delaunay_time = -1.0;
    double intersection_check_time = -1.0;
    std::size_t mwt_num_chull = 0;
    std::size_t mwt_num_edges = 0;
    std::size_t del_num_chull = 0;
    std::size_t del_num_edges = 0;
};

/**
 * MWT verifier for an MWT computed by another process
 * and loaded from a file.
 */
template<typename Kernel> class LoadedMWTVerifier {
  public:
    using Point_2 = CGAL::Point_2<Kernel>;
    using Segment_2 = CGAL::Segment_2<Kernel>;
    using Delaunay = CGAL::Delaunay_triangulation_2<Kernel>;

    LoadedMWTVerifier(std::vector<Point_2> all_points, std::vector<std::array<std::uint64_t, 2>> all_edges)
        : m_points(std::move(all_points)), m_edges(std::move(all_edges)) {}

    void check_edge_count() {
        std::size_t num_edges_loaded = m_edges.size();
        Delaunay delaunay{m_points.begin(), m_points.end()};
        std::size_t num_edges_delaunay = std::distance(delaunay.finite_edges_begin(), delaunay.finite_edges_end());
        if(num_edges_delaunay != num_edges_loaded) {
            throw std::runtime_error("Validation failed - number of edges in loaded MWT (" +
                                     std::to_string(num_edges_loaded) +
                                     ") does not match number of edges in Delaunay triangulation (" +
                                     std::to_string(num_edges_delaunay) + ")!");
        }
    }

    void check_no_intersections_cgal_exact() {
        using EPoint = CGAL::Point_2<CGAL::Epeck>;
        using ESegment = CGAL::Segment_2<CGAL::Epeck>;
        auto ep = [&](std::size_t i1) { return EPoint{m_points[i1].x(), m_points[i1].y()}; };
        std::vector<ESegment> segments;
        for(const auto &edge : m_edges) {
            segments.emplace_back(ep(edge[0]), ep(edge[1]));
        }
        if(CGAL::do_curves_intersect(segments.begin(), segments.end())) {
            throw std::runtime_error("Validation failed - there are some intersecting edges!");
        }
    }

    void check_no_intersections() {
        std::vector<CGAL::Segment_2<Kernel>> segments;
        for(const auto &edge : m_edges) {
            segments.emplace_back(m_points[edge[0]], m_points[edge[1]]);
        }
        CrossingFreenessVerifier<Kernel> verifier{std::move(segments)};
        if(verifier()) {
            auto [seg1, seg2] = verifier.get_intersecting_segments();
            std::stringstream errmsg;
            errmsg << "Validation failed - edges (" << seg1.source() << " -> " << seg1.target() << ") and ("
                   << seg2.source() << " -> " << seg2.target() << ") intersect!";
            throw std::runtime_error(errmsg.str());
        }
    }

    const std::vector<Point_2> &get_points() const noexcept { return m_points; }

    const std::vector<std::array<std::uint64_t, 2>> &get_edges() const noexcept { return m_edges; }

  private:
    std::vector<Point_2> m_points;
    std::vector<std::array<std::uint64_t, 2>> m_edges;
};

} // namespace mwt

#endif
