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
)

install(TARGETS tdk DESTINATION bin)