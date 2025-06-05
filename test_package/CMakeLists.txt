cmake_minimum_required(VERSION 3.19)

add_executable(test_cgal_mwt test_main.cpp test_interval_set.cpp test_static_quadtree.cpp
                             test_mwt_traits.cpp test_directional_filter.cpp test_incremental_search.cpp
                             test_filtered_search.cpp test_lmt_skeleton.cpp test_face_collector.cpp
                             test_interval_face_triangulator.cpp test_rational_or_int.cpp
                             test_number_types.cpp test_sparse_lu.cpp 
                             test_exact_lp.cpp test_rational_sparse_matrix.cpp
                             test_mutable_index_pq.cpp debug_cover_cut_generator_data.cpp
                             test_cover_cut_generator.cpp test_intersection_check.cpp
                             test_read_instance.cpp test_scaled_weight_sign.cpp)

get_directory_property(HAS_PARENT PARENT_DIRECTORY)
if(HAS_PARENT)
  # we are not the top-level file (i.e., add_subdirectory'd from main project)
  target_link_libraries(test_cgal_mwt PRIVATE cgal_mwt nlohmann_json::nlohmann_json)
else()
  # we are the top-level file (i.e., called during conan test_package)
  project(test_cgal_mwt CXX)
  find_package(doctest REQUIRED)
  find_package(cgal REQUIRED)
  find_package(TBB REQUIRED)
  find_package(cgal_mwt REQUIRED)
  find_package(nlohmann_json REQUIRED)
  target_link_libraries(test_cgal_mwt PRIVATE cgal_mwt::cgal_mwt CGAL::CGAL onetbb::onetbb nlohmann_json::nlohmann_json)
  target_compile_features(test_cgal_mwt PRIVATE cxx_std_17)
  include(CTest)
endif()
target_link_libraries(test_cgal_mwt PRIVATE doctest::doctest)
add_test(NAME run_test_cgal_mwt WORKING_DIRECTORY "${CMAKE_BINARY_DIR}" COMMAND test_cgal_mwt)
