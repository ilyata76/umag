/*
 * Точка входа в библиотеку, сборка с PyBind11 всех модулей, подмодулей,
 * функций и классов.
 */

#include "constants.hpp"
#include "geometries.hpp"
#include "interactions.hpp"
#include "logger.hpp"
#include "printer.hpp"
#include "registries.hpp"
#include "simulation.hpp"
#include "solvers.hpp"
#include "types.hpp"

#include <pybind11/functional.h>
#include <pybind11/iostream.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

PYBIND11_MODULE(core, module) {
    py::add_ostream_redirect(module, "ostream_redirect");
    pyBindConstants(module);
    pyBindTypes(module);
    pyBindRegistries(module);
    pyBindGeometries(module);
    pyBindInteractions(module);
    pyBindSolvers(module);
    pyBindSimulation(module);
    pyBindPrinter(module);
    pyBindLogger(module); // регистрация логгера
};
