cmake_minimum_required(VERSION 3.15)
project(TDK VERSION 0.1 LANGUAGES CXX)

# Set default build type to Debug if not specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE "Debug" CACHE STRING "Choose the type of build." FORCE)
endif()

# Require modern C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Export compile_commands.json for better IDE integration
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Set library paths for Homebrew installations
set(EIGEN3_INCLUDE_DIR "/opt/homebrew/opt/eigen/include/eigen3")
set(GMSH_INCLUDE_DIR "/opt/homebrew/opt/gmsh/include")
set(GMSH_LIB_DIR "/opt/homebrew/opt/gmsh/lib")

# Check if directories exist
if(NOT EXISTS ${EIGEN3_INCLUDE_DIR})
    message(FATAL_ERROR "Eigen3 not found at ${EIGEN3_INCLUDE_DIR}")
endif()

if(NOT EXISTS ${GMSH_INCLUDE_DIR})
    message(FATAL_ERROR "Gmsh not found at ${GMSH_INCLUDE_DIR}")
endif()

message(STATUS "Found Eigen3 at: ${EIGEN3_INCLUDE_DIR}")
message(STATUS "Found Gmsh at: ${GMSH_INCLUDE_DIR}")

# Add the main executable
add_executable(tdk
    src/main.cpp
)

# Set compile options for debugging information
target_compile_options(tdk PRIVATE
    $<$<CONFIG:Debug>:-g>
)

target_include_directories(tdk
    PRIVATE
        ${PROJECT_SOURCE_DIR}/src
        ${EIGEN3_INCLUDE_DIR}
        ${GMSH_INCLUDE_DIR}
)

target_link_libraries(tdk 
    PRIVATE
        ${GMSH_LIB_DIR}/libgmsh.dylib
)

install(TARGETS tdk DESTINATION bin)