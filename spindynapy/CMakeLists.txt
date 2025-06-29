# ===================================================================
#  Project Configuration for Spin Dynamics C++ Backend
# ===================================================================

cmake_minimum_required(VERSION 3.5)
project(core)

# -------------------------------------------------------------------
#  Language Standard and General Compiler Settings
# -------------------------------------------------------------------

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Output directories
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}")

# Use cached packages if available (CPM)
set(CPM_USE_LOCAL_PACKAGES ON CACHE BOOL "Use local packages if available")
set(CPM_SOURCE_CACHE "${CMAKE_SOURCE_DIR}/build-deps/" CACHE PATH "Path to CPM source cache")

# Generate compile_commands.json (used by clangd, IDEs, tooling)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# -------------------------------------------------------------------
#  Compiler Flags (Optimization, Debug, Diagnostics)
# -------------------------------------------------------------------

# Release build flags
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -march=native -flto")
#   -O3            : Maximum compiler optimization
#   -march=native  : Optimize for the host CPU
#   -flto          : Enable link-time optimization

# Debug build flags
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -O0 -fsanitize=address -static-libasan -fsanitize=undefined")
#   -g                        : Generate debug information
#   -O0                       : Disable optimization for easier debugging
#   -fsanitize=address        : Detect memory issues (ASan)
#   -fsanitize=undefined      : Detect undefined behavior (UBSan)
#   -static-libasan           : Static link for ASan (e.g., Python C extensions)

# Additional diagnostics and static analysis flags
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -Wextra -Wpedantic -Wshadow -Wconversion -Wduplicated-cond -Wduplicated-branches -Wlogical-op -Wuseless-cast")
#   -Wall                     : Basic warnings
#   -Wextra                   : Extra warnings
#   -Wpedantic                : Enforce strict compliance with the standard
#   -Wshadow                  : Warn on variable shadowing (e.g., lambdas, loops)
#   -Wconversion              : Warn on implicit conversions
#   -Wduplicated-cond         : Duplicate conditions in if/else
#   -Wduplicated-branches     : Duplicate code in if/else
#   -Wlogical-op              : Suspicious logical expressions
#   -Wuseless-cast            : Redundant type casts (relevant for pybind11)

# Global warnings for all build types
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")
#   -Wall                     : Basic warnings
#   -Wextra                   : Extra warnings

# -------------------------------------------------------------------
#  Dependencies
# -------------------------------------------------------------------

# Python 3.12 headers (mandatory) <!!!>
find_package(Python 3.12 REQUIRED COMPONENTS Interpreter Development)

# OpenMP (multi-threading support)
find_package(OpenMP REQUIRED)
if(OpenMP_CXX_FOUND)
    message(STATUS "OpenMP found: ${OpenMP_CXX_LIBRARIES}")
else()
    message(FATAL_ERROR "OpenMP not found. Please install OpenMP (e.g., libomp-dev on Ubuntu).")
endif()

# CPM - dependency manager
include(${CMAKE_SOURCE_DIR}/cpm/cmake/CPM.cmake)

# pybind11 (Python bindings)
CPMAddPackage(
    NAME pybind11
    GITHUB_REPOSITORY pybind/pybind11
    VERSION 2.13.0
)
include_directories(${pybind11_SOURCE_DIR})

# Eigen (linear algebra)
CPMAddPackage(
    NAME eigen
    GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
    GIT_TAG 3.4.0
)
include_directories(${eigen_SOURCE_DIR})

# fmt (string formatting)
CPMAddPackage(
    NAME fmt
    GITHUB_REPOSITORY fmtlib/fmt
    GIT_TAG 11.2.0
)
include_directories(${fmt_SOURCE_DIR})

# -------------------------------------------------------------------
#  Source Files
# -------------------------------------------------------------------

file(GLOB HEADERS 
    "core_src/*.hpp"
    "core_src/**/*.hpp"
)

file(GLOB SOURCES 
    "core_src/*.cpp"
    "core_src/**/*.cpp"
)

# -------------------------------------------------------------------
#  pybind11 Module Definition and Properties
# -------------------------------------------------------------------

pybind11_add_module(
    ${PROJECT_NAME}
    ${HEADERS}
    ${SOURCES}
)

target_include_directories(${PROJECT_NAME} PUBLIC
    ${eigen_SOURCE_DIR}
    ${pybind11_SOURCE_DIR}
)

target_link_libraries(${PROJECT_NAME} PUBLIC eigen pybind11::module OpenMP::OpenMP_CXX fmt::fmt-header-only)

# Enable link-time optimization for release builds
set_property(TARGET ${PROJECT_NAME} PROPERTY INTERPROCEDURAL_OPTIMIZATION True)

# -------------------------------------------------------------------
#  Tooling and Post-Build Steps
# -------------------------------------------------------------------

# Copy compile_commands.json to source root for clangd, IDE support (./.clandg)
add_custom_command(
    TARGET ${PROJECT_NAME} POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy
        "${CMAKE_BINARY_DIR}/compile_commands.json"
        "${CMAKE_SOURCE_DIR}/compile_commands.json"
    COMMENT "Copying compile_commands.json to source directory"
)
