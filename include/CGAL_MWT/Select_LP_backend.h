#ifndef CGAL_MWT_SELECT_LP_BACKEND_H_INCLUDED_
#define CGAL_MWT_SELECT_LP_BACKEND_H_INCLUDED_

#include "Exact_LP_face_triangulator.h"
#include "Face_triangulator.h"
#include "Gurobi_LP_backend.h"
#include "LP_backends.h"

namespace mwt {

template<typename Skeleton_> struct BackendsExhaustedNonsimpleFaceTriangulator {
    using Kernel = typename Skeleton_::Kernel;
    using FT = typename Kernel::FT;
    using Skeleton = Skeleton_;
    using Face = mwt::Face<Skeleton>;
    using Halfedge = typename Skeleton::Halfedge;
    using Handle = Halfedge *;
    using Interval = CGAL::Interval_nt_advanced;
    using Iterator = typename std::vector<Handle>::iterator;

    BackendsExhaustedNonsimpleFaceTriangulator(FaceTriangulatorOptions) {}

    void operator()(Face &&face) {
        throw std::runtime_error("No LP backend is available at runtime to triangulate a nonsimple face");
    }

    double get_maximum_gap() const noexcept { return std::numeric_limits<double>::infinity(); }
};

template<typename FaceCollector>
nlohmann::json triangulate_nonsimple_faces_with_lp(typename FaceCollector::Skeleton &skeleton,
                                                   FaceTriangulatorOptions options) {
    using Skeleton = typename FaceCollector::Skeleton;
    using SimpleTriangulator = Exact_simple_face_triangulator<Skeleton>;

    /**
     * Try using Gurobi first, if we are compiled with it.
     */
#if defined(CGAL_MWT_HAVE_GUROBI)
    if(is_runtime_available<GurobiLPBackend>()) {
        try {
            using GurobiBackendNonsimpleTriangulator = Exact_LP_face_triangulator<Skeleton, GurobiLPBackend>;
            return triangulate_faces<FaceCollector, SimpleTriangulator, GurobiBackendNonsimpleTriangulator>(skeleton,
                                                                                                            options);
        } catch(const GRBException &e) {
            throw std::runtime_error("Gurobi exception (GRBException): " + e.getMessage() + " (" +
                                     std::to_string(e.getErrorCode()) + ")");
        }
    }
#endif

    using BackendsExhausted = BackendsExhaustedNonsimpleFaceTriangulator<Skeleton>;
    return triangulate_faces<FaceCollector, SimpleTriangulator, BackendsExhausted>(skeleton, options);
}

template<typename FaceCollector>
nlohmann::json triangulate_nonsimple_faces_with_gap(typename FaceCollector::Skeleton &skeleton, double *max_gap,
                                                    FaceTriangulatorOptions options) {
    using Skeleton = typename FaceCollector::Skeleton;
    using SimpleTriangulator = Exact_simple_face_triangulator<Skeleton>;

#if defined(CGAL_MWT_HAVE_GUROBI)
    if(is_runtime_available<GurobiLPBackend>()) {
        using GurobiBackendNonsimpleTriangulator = Inexact_LP_face_triangulator_with_max_gap<Skeleton, GurobiLPBackend>;
        return triangulate_faces_with_gap<FaceCollector, SimpleTriangulator, GurobiBackendNonsimpleTriangulator>(
            skeleton, max_gap, options);
    }
#endif

    using BackendsExhausted = BackendsExhaustedNonsimpleFaceTriangulator<Skeleton>;
    return triangulate_faces_with_gap<FaceCollector, SimpleTriangulator, BackendsExhausted>(skeleton, max_gap, options);
}

} // namespace mwt

#endif
