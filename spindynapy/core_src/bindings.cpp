/*
 * Точка входа в библиотеку, сборка с PyBind11 всех модулей, подмодулей,
 * функций и классов.
 */

#include "constants.hpp"
#include "geometries.hpp"
#include "interactions.hpp"
#include "registries.hpp"
#include "simulation.hpp"
#include "solvers.hpp"
#include "types.hpp"

#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

PYBIND11_MODULE(core, module) {
    pyBindConstants(module);
    pyBindTypes(module);
    pyBindRegistries(module);
    pyBindGeometries(module);
    pyBindInteractions(module);
    pyBindSolvers(module);
    pyBindSimulation(module);
};
