cmake_minimum_required(VERSION 3.15)
project(TDK LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Eigen3 3.3 REQUIRED NO_MODULE)
find_path(GMSH_INCLUDE_DIR gmsh.h)
find_library(GMSH_LIBRARY gmsh)

if(NOT GMSH_INCLUDE_DIR OR NOT GMSH_LIBRARY)
  message(FATAL_ERROR "Gmsh library not found")
endif()

add_executable(tdk src/main.cpp)
target_include_directories(tdk PRIVATE ${GMSH_INCLUDE_DIR})
target_link_libraries(tdk PRIVATE Eigen3::Eigen ${GMSH_LIBRARY})
