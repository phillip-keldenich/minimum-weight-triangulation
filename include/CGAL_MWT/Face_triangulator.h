#ifndef CGAL_MWT_FACE_TRIANGULATOR_H_INCLUDED_
#define CGAL_MWT_FACE_TRIANGULATOR_H_INCLUDED_

#include "Exact_simple_face_triangulator.h"
#include "Face_collector.h"
#include "Inexact_LP_face_triangulator_with_max_gap.h"
#include <CGAL/FPU.h>
#include <nlohmann/json.hpp>
#include <tbb/parallel_pipeline.h>
#include <tbb/tbb.h>

namespace mwt {

namespace triangulator_has_stats_impl {

template<typename FaceTriangulator,
         typename = decltype(std::declval<FaceTriangulator &>().combine_stats(std::declval<nlohmann::json &>()))>
std::true_type has_stats_impl(const FaceTriangulator *);
std::false_type has_stats_impl(...);

} // namespace triangulator_has_stats_impl

template<typename FaceTriangulator>
constexpr bool triangulator_has_stats_v =
    decltype(triangulator_has_stats_impl::has_stats_impl(std::declval<const FaceTriangulator *>()))::value;

template<typename FaceCollector> auto make_tbb_face_producer(FaceCollector &collector) {
    using Face = typename FaceCollector::Face;
    return tbb::make_filter<void, Face>(tbb::filter_mode::serial_out_of_order, [&](tbb::flow_control &fc) -> Face {
        if(collector.next()) {
            return collector.get_current_face_copy();
        } else {
            fc.stop();
            return Face{};
        }
    });
}

template<typename FaceCollector, typename SimpleFaceTriangulator, typename NonsimpleFaceTriangulator>
nlohmann::json triangulate_faces(typename FaceCollector::Skeleton &skeleton, FaceTriangulatorOptions options) {
    using Face = typename FaceCollector::Face;
    FaceCollector collector{skeleton};
    std::unique_ptr<std::mutex> stat_lock;
    nlohmann::json stats;
    if constexpr(triangulator_has_stats_v<SimpleFaceTriangulator> ||
                 triangulator_has_stats_v<NonsimpleFaceTriangulator>) {
        stat_lock = std::make_unique<std::mutex>();
    }

    auto handle_face = [&](Face face) {
        CGAL::Protect_FPU_rounding protect_rounding;
        if(face.is_simple()) {
            SimpleFaceTriangulator triangulator(options);
            triangulator(std::move(face));
            if constexpr(triangulator_has_stats_v<SimpleFaceTriangulator>) {
                std::unique_lock l{*stat_lock};
                triangulator.combine_stats(stats);
            }
        } else {
            NonsimpleFaceTriangulator triangulator(options);
            triangulator(std::move(face));
            if constexpr(triangulator_has_stats_v<NonsimpleFaceTriangulator>) {
                std::unique_lock l{*stat_lock};
                triangulator.combine_stats(stats);
            }
        }
    };
    auto producer = make_tbb_face_producer(collector);
    auto consumer = tbb::make_filter<Face, void>(tbb::filter_mode::parallel, handle_face);
    tbb::parallel_pipeline(tbb::this_task_arena::max_concurrency(), producer & consumer);
    return stats;
}

template<typename FaceCollector, typename SimpleTriangulator, typename NonsimpleTriangulator>
nlohmann::json triangulate_faces_with_gap(typename FaceCollector::Skeleton &skeleton, double *max_gap,
                                          FaceTriangulatorOptions options) {
    using Skeleton = typename FaceCollector::Skeleton;
    using Face = typename FaceCollector::Face;
    using Interval = CGAL::Interval_nt_advanced;

    FaceCollector collector{skeleton};
    std::mutex gap_lock;
    nlohmann::json stats;

    auto producer = make_tbb_face_producer(collector);
    auto consumer = tbb::make_filter<Face, void>(tbb::filter_mode::parallel, [&](Face face) {
        if(face.is_simple()) {
            CGAL::Protect_FPU_rounding protect_thread;
            SimpleTriangulator triangulator(options);
            triangulator(std::move(face));
            if constexpr(triangulator_has_stats_v<SimpleTriangulator>) {
                std::unique_lock l{gap_lock};
                triangulator.combine_stats(stats);
            }
        } else {
            NonsimpleTriangulator triangulator(options);
            triangulator(std::move(face));
            double this_gap = triangulator.get_maximum_gap();
            if constexpr(triangulator_has_stats_v<NonsimpleTriangulator>) {
                std::unique_lock l{gap_lock};
                triangulator.combine_stats(stats);
                if(this_gap > 0) {
                    CGAL::Protect_FPU_rounding protect_thread;
                    Interval v{*max_gap};
                    v += this_gap;
                    *max_gap = v.sup();
                }
            } else {
                if(this_gap > 0) {
                    CGAL::Protect_FPU_rounding protect_thread;
                    std::unique_lock l{gap_lock};
                    Interval v{*max_gap};
                    v += this_gap;
                    *max_gap = v.sup();
                }
            }
        }
    });

    tbb::parallel_pipeline(tbb::this_task_arena::max_concurrency(), producer & consumer);
    return stats;
}

} // namespace mwt

#endif
