/**
 * @file    bindings.cpp
 * @brief   PyBind11 entry point for the SpinDynaPy C++ core library.
 *
 * This translation unit glues together every C++ subsystem (constants,
 *   types, registries, geometries, interactions, solvers, simulation, …) and
 *   publishes them under the top‑level Python module name **``spindynapy.core``**.
 * All values are exposed in SI units and documented via docstrings so that
 *   they are introspectable from Python (`help(core.constants)`).
 *
 * Typical usage from Python:
 * ```python
 *   import spindynapy.core as core
 *   geom      = core.geometries.cartesian.Geometry(...)
 *   materials = core.registries.MaterialRegistry(...)
 *   solver    = core.solvers.cartesian.LLGSolver()
 *   sim       = core.simulation.cartesian.Simulation(geom, solver, ...)
 *   sim.simulate_many_steps(100)
 * ```
 *
 * Sub‑module layout created here (low‑level → high‑level):
 *   - **constants**    – fundamental physical constants (see constants.hpp).
 *   - **types**        – abstract base types (Moments, Coordinates, …).
 *   - **registries**   – thread‑safe object stores & caches.
 *   - **geometries**   – spatial distributions of magnetic moments.
 *   - **interactions** – exchange, demagnetisation, Zeeman, ….
 *   - **solvers**      – numerical integrators and field update strategies (Euler, Heun, …).
 *   - **simulation**   – orchestration layer that brings everything together.
 *   - **printer**      – text/visualisation helpers.
 *   - **logger**       – logging utilities.
 *
 * Keep the binding order logical: low‑level modules should be registered
 *   before high‑level ones to avoid missing type conversions.
 *
 * @copyright 2025 SpinDynaPy
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

//! CPython entry‑point
PYBIND11_MODULE(core, module) {

    module.doc() = "High-performance backend for atomistic spin-dynamics simulations.  All heavy\n"
                   "numerics are executed in C++20 with optional parallelism and exposed to\n"
                   "Python via PyBind11.";

    // ---------------------------------------------------------------------
    //  I/O redirection (std::cout/std::cerr → sys.stdout/sys.stderr)
    // ---------------------------------------------------------------------
    py::add_ostream_redirect(module, "ostream_redirect");

    // ---------------------------------------------------------------------
    //  Register sub‑modules
    // ---------------------------------------------------------------------
    pyBindConstants(module);    // fundamental physical constants
    pyBindTypes(module);        // basic abstract data types
    pyBindRegistries(module);   // registries / caches
    pyBindGeometries(module);   // geometries & spatial containers
    pyBindInteractions(module); // physical interactions
    pyBindSolvers(module);      // numerical integrators
    pyBindSimulation(module);   // simulation orchestrator
    pyBindPrinter(module);      // text & file output helpers
    pyBindLogger(module);       // logging utilities (depends on others)
};
