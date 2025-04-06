# Сборка CPP-backend библиотеки

cmake_minimum_required(VERSION 3.5)

project(spindynapy_core)

# Преднастройки проекта
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/lib")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/lib")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/bin")

# python3-dev (!!!)
# Явно указываем Python 3.12
find_package(Python 3.12 REQUIRED COMPONENTS Interpreter Development)

# Зависимости
include(${CMAKE_SOURCE_DIR}/cpm/cmake/CPM.cmake)
CPMAddPackage(
    NAME pybind11  # pybind для работы с CPython API
    GITHUB_REPOSITORY pybind/pybind11
    VERSION 2.13.0
)

# Создаём модуль с pybind11_add_module
pybind11_add_module(
    spindynapy_core
    core/bindings.cpp
    core/hello.cpp
)
