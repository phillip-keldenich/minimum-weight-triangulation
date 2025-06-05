cmake_minimum_required(VERSION 3.19)

# try to auto-detect the conan toolchain file, if present
if(NOT DEFINED CMAKE_TOOLCHAIN_FILE)
  if(EXISTS "${CMAKE_BINARY_DIR}/conan_toolchain.cmake")
    set(CMAKE_TOOLCHAIN_FILE "${CMAKE_BINARY_DIR}/conan_toolchain.cmake")
  endif()
  if(EXISTS "${CMAKE_BINARY_DIR}/generators/conan_toolchain.cmake")
    set(CMAKE_TOOLCHAIN_FILE "${CMAKE_BINARY_DIR}/generators/conan_toolchain.cmake")
  endif()
endif()

file(READ "${CMAKE_CURRENT_LIST_DIR}/version.txt" CGAL_MWT_FILE_VERSION)
string(STRIP "${CGAL_MWT_FILE_VERSION}" CGAL_MWT_FILE_VERSION)
project(cgal_mwt LANGUAGES CXX VERSION ${CGAL_MWT_FILE_VERSION})

option(CGAL_MWT_ENABLE_TESTING "Enable testing the code using doctest" OFF)

include(GNUInstallDirs)
find_package(CGAL REQUIRED)
find_package(TBB REQUIRED)
find_package(gurobi REQUIRED)
find_package(Boost REQUIRED COMPONENTS program_options iostreams)
find_package(nlohmann_json REQUIRED)

add_subdirectory("src")

if(CGAL_MWT_ENABLE_TESTING)
  find_package(doctest REQUIRED)
  include(CTest)
  add_subdirectory("test_package")
endif()
