# Сборка CPP-backend библиотеки

cmake_minimum_required(VERSION 3.5)

project(core)

# Преднастройки проекта
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}")

# python3-dev headers (!!!)
# Явно указываем Python 3.12
find_package(Python 3.12 REQUIRED COMPONENTS Interpreter Development)

# Зависимости
include(${CMAKE_SOURCE_DIR}/cpm/cmake/CPM.cmake)
CPMAddPackage(
    NAME pybind11  # pybind для работы с CPython API
    GITHUB_REPOSITORY pybind/pybind11
    VERSION 2.13.0
)

# Собираем все .hpp файлы
file(GLOB HEADERS 
    "core_src/*.hpp"
    "core_src/**/*.hpp"
)

# Собираем все .cpp файлы
file(GLOB SOURCES 
    "core_src/*.cpp"
    "core_src/**/*.cpp"
)

# Создаём модуль с pybind11_add_module
pybind11_add_module(
    ${PROJECT_NAME}
    ${HEADERS}
    ${SOURCES}
)
