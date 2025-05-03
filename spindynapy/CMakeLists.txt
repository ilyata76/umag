# Сборка CPP-backend библиотеки

cmake_minimum_required(VERSION 3.5)

project(core)

# Преднастройки проекта
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}")
set(CPM_USE_LOCAL_PACKAGES ON CACHE BOOL "Use local packages if available")
set(CPM_SOURCE_CACHE "${CMAKE_SOURCE_DIR}/build-deps/" CACHE PATH "Path to CPM source cache")
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Флаги для Release сборки
# -O3 - максимальная оптимизация
# -march=native - оптимизация под архитектуру текущего процессора (может снизить переносимость бинарника)
# -flto - Link-Time Optimization (оптимизация на этапе линковки)
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -march=native -flto")

# Флаги для Debug сборки
# -g - добавить отладочные символы (обязательно для GDB/LLDB)
# -O0 - отключить оптимизацию (упрощает отладку, т.к. код не переупорядочивается)
# -fsanitize=address - включить AddressSanitizer (ASan) для поиска ошибок памяти
# -fsanitize=undefined - включить UndefinedBehaviorSanitizer (UBSan) для поиска неопределенного поведения
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -O0 -fsanitize=address -static-libasan -fsanitize=undefined")

# Флаги для Debug сборки для cland-сервера
# -Wall - Базовые предупреждения
# -Wextra - Дополнительные полезные предупреждения
# -Wpedantic - Строгое соблюдение стандарта
# -Wshadow - Тени переменных (особенно актуально в lambdas и циклах)
# -Wconversion - неявные преобразования
# -Wduplicated-cond - Повторяющиеся условия в if/else
# -Wduplicated-branches - Повторяющиеся блоки if/else
# -Wlogical-op - Подозрительные логические операции
# -Wuseless-cast - Бесполезные касты (особенно при привязке к pybind)
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -Wextra -Wpedantic -Wshadow -Wconversion -Wduplicated-cond -Wduplicated-branches -Wlogical-op -Wuseless-cast")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")

# python3-dev headers (!!!)
# Явно указываем Python 3.12
find_package(Python 3.12 REQUIRED COMPONENTS Interpreter Development)

# Зависимости

# OpenMP должен поставляться через API компилятора
find_package(OpenMP REQUIRED)
if(OpenMP_CXX_FOUND)
    message(STATUS "OpenMP found: ${OpenMP_CXX_LIBRARIES}")
else()
    message(FATAL_ERROR "OpenMP not found. Please install OpenMP (e.g., libomp-dev on Ubuntu).")
endif()

include(${CMAKE_SOURCE_DIR}/cpm/cmake/CPM.cmake)

CPMAddPackage(
    NAME pybind11  # pybind для работы с CPython API
    GITHUB_REPOSITORY pybind/pybind11
    VERSION 2.13.0
)
include_directories(${pybind11_SOURCE_DIR})

CPMAddPackage(
    NAME eigen  # eigen для работы с линейной алгеброй
    GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
    GIT_TAG 3.4.0
)
include_directories(${eigen_SOURCE_DIR})

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

target_include_directories(${PROJECT_NAME} PUBLIC
    ${eigen_SOURCE_DIR}
    ${pybind11_SOURCE_DIR}
)

target_link_libraries(${PROJECT_NAME} PUBLIC eigen pybind11::module OpenMP::OpenMP_CXX)

# Включаем LTO для Release сборок через свойство цели (замена -flto)
set_property(TARGET ${PROJECT_NAME} PROPERTY INTERPROCEDURAL_OPTIMIZATION True)

# Копирует compile_commands.json в корень проекта после конфигурации
add_custom_command(
    TARGET ${PROJECT_NAME} POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy
        "${CMAKE_BINARY_DIR}/compile_commands.json"
        "${CMAKE_SOURCE_DIR}/compile_commands.json"
    COMMENT "Copying compile_commands.json to source directory"
)

