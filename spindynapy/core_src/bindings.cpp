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

namespace py = pybind11;

PYBIND11_MODULE(core, module) {

    // types
    py::module_ types_module = module.def_submodule("types");
    pyBindTypes(types_module);

    // geometries
    py::module_ geometries_module = module.def_submodule("geometries");
    pyBindGeometries(geometries_module);

    // solvers
    py::module_ solvers_module = module.def_submodule("solvers");
    pyBindSolvers(solvers_module);

    // constants
    py::module_ constants_module = module.def_submodule("constants");
    pyBindConstants(constants_module);

    // interactions
    py::module_ interaction_module = module.def_submodule("interactions");
    pyBindInteractions(interaction_module);

    // registries
    py::module_ registries_module = module.def_submodule("registries");
    pyBindRegistries(registries_module);

    // simulation
    py::module_ simulation_module = module.def_submodule("simulation");
    pyBindSimulation(simulation_module);

}; // ! PYBIND11_MODULE

// TODO: комментарии нормальные ко всем модулям, функциям и т.д.
