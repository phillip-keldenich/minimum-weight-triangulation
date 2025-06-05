add_library(cgal_mwt STATIC version.cpp)
target_link_libraries(cgal_mwt PUBLIC CGAL::CGAL)
target_link_libraries(cgal_mwt PUBLIC onetbb::onetbb)
target_link_libraries(cgal_mwt PUBLIC gurobi::gurobi)
target_link_libraries(cgal_mwt PUBLIC Boost::boost Boost::iostreams)

target_compile_features(cgal_mwt PUBLIC cxx_std_17)

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  target_compile_options(cgal_mwt PUBLIC "-frounding-math")
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
  target_compile_options(cgal_mwt PUBLIC "/fp:strict" "/fp:except-")
endif()

target_include_directories(cgal_mwt PUBLIC 
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/../include> # include the public header directory in the source tree when building the lib
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)    # include the include directory when using the installed lib

target_compile_definitions(cgal_mwt PRIVATE "CGAL_MWT_VERSION=\"${CGAL_MWT_FILE_VERSION}\"")

add_executable(solve_mwt solve_mwt.cpp)
target_link_libraries(solve_mwt PRIVATE cgal_mwt Boost::program_options
                      nlohmann_json::nlohmann_json)
target_compile_features(solve_mwt PUBLIC cxx_std_17)

add_executable(verify_triangulation verify_triangulation.cpp)
target_link_libraries(verify_triangulation PRIVATE cgal_mwt
                      nlohmann_json::nlohmann_json)
target_compile_features(verify_triangulation PUBLIC cxx_std_17)

add_executable(compare_triangulations compare_triangulations.cpp)
target_link_libraries(compare_triangulations PRIVATE cgal_mwt
                      nlohmann_json::nlohmann_json)
target_compile_features(compare_triangulations PUBLIC cxx_std_17)

set(PUBLIC_HEADERS "Are_all_collinear.h"
"CGAL_Rational_aux.h"
"Compare_weights.h"
"Crossing_free.h"
"Dead_sector.h"
"Diamond_filter.h"
"Directional_filter.h"
"Dynamic_program_utils.h"
"Exact_LP_face_triangulator.h"
"Exact_diamond_filtered_search_driver.h"
"Exact_simple_face_triangulator.h"
"Face_analyzer.h"
"Face_collector.h"
"Face_triangulator.h"
"Filtered_incremental_search.h"
"Generate_random_integral.h"
"Gurobi_LP_backend.h"
"Inexact_LP_face_triangulator_with_max_gap.h"
"Interval_set.h"
"Is_interval.h"
"LMT_halfedge.h"
"LMT_skeleton.h"
"LP_backends.h"
"Mutable_index_priority_queue.h"
"Mwt_traits.h"
"Neighbor_set.h"
"Null_LP_backend.h"
"Quadtree_aux.h"
"Rational_or_int.h"
"Rational_sparse_LU.h"
"Rational_sparse_gaussian.h"
"Rational_sparse_matrix.h"
"Read_instance.h"
"Scaled_weight_sign.h"
"Search_aux.h"
"Select_LP_backend.h"
"Static_quadtree.h"
"Total_weight.h"
"Triangle_based_linear_model.h"
"Unique_hash_map.h"
"Validate.h"
"output.h"
"time_util.h"
"version.h")

list(TRANSFORM PUBLIC_HEADERS PREPEND "${CMAKE_CURRENT_LIST_DIR}/../include/CGAL_MWT/")
set_target_properties(cgal_mwt PROPERTIES PUBLIC_HEADER "${PUBLIC_HEADERS}")

install(TARGETS cgal_mwt EXPORT cgal_mwt_targets PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/CGAL_MWT")
