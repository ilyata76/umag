#ifndef __SIMULATION_HPP__
#define __SIMULATION_HPP__

/**
 * @file   simulation.hpp
 * @brief  Simulation interface and core integration logic for spin dynamics systems.
 *
 * This header defines the core simulation controller that integrates:
 *   - spatial geometry and magnetic moment layout,
 *   - solver strategies for time evolution,
 *   - interaction registry (exchange, anisotropy, dipolar, etc.),
 *   - material registry (atomic parameters),
 *   - time integration loop (with macrocell refresh, caching, etc.).
 *
 * Functional overview:
 *
 * - `Simulation<CoordSystem>`:
 *   - Owns the full simulation state: geometry, solver, interactions, materials.
 *   - Advances the system state using a time integrator (`ISolver`).
 *   - Manages step-by-step snapshots (`SimulationStepData`).
 *   - Provides thread-safe execution with optional caching and macrocell updates.
 *
 * - `SimulationStepData<CoordSystem>`:
 *   - Records solver buffer, geometry snapshot, magnetization, and energy at each saved step.
 *   - Used for trajectory postprocessing, diagnostics, and time-series export.
 *
 * Python bindings:
 * - Submodule: `simulation.cartesian`
 * - Macroc:
 *   - `SIMULATION_TEMPLATE_BINDINGS(cls)`
 *   - `SIMULATIONSTEPDATA_TEMPLATE_BINDINGS(cls)`
 * - Types exposed:
 *   - `Simulation` (core loop),
 *   - `SimulationStepData` (individual step snapshot).
 *
 * Interface:
 * - `Simulation<CoordSystem>`     – time integration and state evolution interface,
 * - `SimulationStepData<CoordSystem>` – snapshot of state at a simulation step.
 *
 * Concrete:
 * - `cartesian::Simulation`       – simulation in Cartesian coordinate system,
 * - `cartesian::SimulationStepData` – corresponding Cartesian snapshot.
 *
 * @note All units are SI unless explicitly stated.
 *       Time is in [s], magnetic field in [T], energy in [J], position in [m].
 *
 * @note Geometry, solver, and interaction registry must be preconfigured before simulation starts.
 *
 * @copyright 2025 SpinDynaPy
 */

#include "geometries.hpp"
#include "interactions.hpp"
#include "logger.hpp"
#include "registries.hpp"
#include "solvers.hpp"
#include "types.hpp"

#include <algorithm>
#include <ctime>
#include <memory>
#include <pybind11/chrono.h>
#include <pybind11/detail/common.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/iostream.h>
#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>
#include <pybind11/stl.h>
#include <stdexcept>
#include <string>

namespace py = pybind11;

namespace PYTHON_API spindynapy {

// ===========================================================================
//  Step data
// ===========================================================================

/**
 * @struct SimulationStepData
 * @brief Snapshot of simulation state at a given time step.
 *
 * Represents a saved snapshot of the spin system at a discrete time step during simulation.
 * @details Inherits from `SolverData` and extends it with temporal metadata and a full geometry snapshot.
 *
 * This structure contains:
 *   - Simulation time and step number,
 *   - Full clone of geometry at this time step,
 *   - Solver-specific data (effective fields, energies, etc.) for saved step.
 *
 * Used for:
 *   - Trajectory logging,
 *   - Postprocessing (e.g., time-series export),
 *   - State resumption,
 *   - Visualization.
 *
 * @tparam CoordSystem The coordinate system in which the simulation is run.
 *
 * @note All fields are stored in SI units.
 */
template <CoordSystemConcept CoordSystem> struct SimulationStepData : public SolverData<CoordSystem> {
  public:
    /**
     * @brief Simulation time at this step [s] from zero.
     */
    double time;

    /**
     * @brief Step index (discrete time index, starting from zero).
     */
    uint step;

    /**
     * @brief Snapshot of the geometry (moments and layout) at this step.
     *
     * This is a deep clone (without cache) of the original geometry, used for replay and analysis.
     */
    std::unique_ptr<IGeometry<CoordSystem>> geometry;

    /**
     * @brief Default constructor.
     *
     * Initializes an empty (null) snapshot. Used for initialization and
     *      subsequent population.
     * @warning Fields must be manually assigned before use.
     */
    SimulationStepData() {};

    /**
     * @brief Fully construct a snapshot from simulation data.
     *
     * @param time          Simulation time at this step [s].
     * @param step          Step index.
     * @param geometry      Cloned geometry at this step. Ensure that Macrocells are updated if used.
     * @param solver_data   Copy of effective fields and energies.
     */
    SimulationStepData(
        double time,
        uint step,
        std::unique_ptr<IGeometry<CoordSystem>> geometry,
        SolverData<CoordSystem> solver_data
    )
        : SolverData<CoordSystem>(solver_data), time(time), step(step), geometry(std::move(geometry)) {};

    /**
     * @brief Compute the mean magnetization vector ⟨𝐌⟩ of the system.
     *
     * Calculates the vector average of all spin directions across the geometry.
     * @details Formula: ⟨𝐌⟩ = (1/N) ∑ᵢ 𝐒ᵢ
     *
     * @returns Average magnetization unit-vector.
     */
    PYTHON_API Magnetization getMeanMagnetization() const {
        size_t moments_size = this->geometry->size();
        Eigen::Vector3d mean_magnetization = Eigen::Vector3d::Zero();
        if (moments_size != 0) {
            for (size_t i = 0; i < moments_size; ++i) {
                mean_magnetization += this->geometry->operator[](i).getDirection().asVector();
            }
            mean_magnetization /= double(moments_size);
        }
        return mean_magnetization;
    }

    /**
     * @brief Compute the norm of mean magnetization ∥⟨𝐌⟩∥.
     *
     * Indicates how aligned the system is.
     * Value ranges from 0 (disordered) to 1 (fully aligned).
     *
     * @returns Norm (length) of the average magnetization vector.
     */
    PYTHON_API double getMeanMagnetizationNorm() const { return this->getMeanMagnetization().norm(); }

    /**
     * @brief Compute the total energy of the system [J].
     *
     * Sums up all scalar energy contributions from all interactions:
     *   E_total = ∑ᵢ Eᵢ
     *
     * @returns Scalar total energy [J].
     */
    PYTHON_API double getEnergy() const {
        return std::accumulate(this->energies.begin(), this->energies.end(), 0.0);
    }

    /**
     * @brief Get energy breakdown per interaction [J].
     *
     * Returns a map from interaction index to its total contribution.
     *   Eᵢ = ∑ⱼ Eᵢⱼ for interaction i and moment j.
     *
     * @returns Map from interaction ID to total energy [J].
     */
    PYTHON_API std::unordered_map<regnum, double> getEnergyByInteraction() const {
        std::unordered_map<regnum, double> energy_by_interaction;
        for (const auto &[reg, energies] : this->interaction_energies) {
            energy_by_interaction[reg] = std::accumulate(energies.begin(), energies.end(), 0.0);
        }
        return energy_by_interaction;
    }
};

namespace cartesian {

/**
 * @struct AbstractSimulationStepData
 * @brief Snapshot of simulation state at a given time step.
 *
 * Represents a saved snapshot of the spin system at a discrete time step during simulation.
 * @details Inherits from `SolverData` and extends it with temporal metadata and a full geometry snapshot.
 *
 * This structure contains:
 *   - Simulation time and step number,
 *   - Full clone of geometry at this time step,
 *   - Solver-specific data (effective fields, energies, etc.) for saved step.
 *
 * Used for:
 *   - Trajectory logging,
 *   - Postprocessing (e.g., time-series export),
 *   - State resumption,
 *   - Visualization.
 *
 * @note All fields are stored in SI units.
 */
using AbstractSimulationStepData = SimulationStepData<NamespaceCoordSystem>;

}; // namespace cartesian

// ===========================================================================
//  Simulator
// ===========================================================================

/**
 * @class Simulation
 * @brief Simulation orchestrator for time evolution of spin systems.
 *
 * Central entry point for performing spin dynamics simulations.
 *
 * Encapsulates the evolving state of the system (geometry, solver, interactions)
 *   and provides methods to step the simulation forward in time.
 *
 * The class manages:
 *   - The evolving spin geometry (positions, directions, materials),
 *   - A solver instance responsible for time integration,
 *   - Registries for materials and magnetic interactions,
 *   - Discrete time tracking (`step`, `time`, `dt`),
 *   - Buffered results for post-analysis or output (via `SimulationStepData`).
 *
 * This class implements a minimal but thread-safe interface for:
 *   - Initial system preparation (geometry and interactions),
 *   - Advancing the system by one or more steps,
 *   - Saving snapshots of the current system state,
 *   - Accessing simulation history or clearing memory.
 *
 * @tparam CoordSystem The coordinate system used throughout the simulation (e.g., Cartesian).
 *
 * @note Units are SI unless otherwise noted.
 * @note Geometry, solver, and registries are passed as shared pointers to allow reuse and sharing.
 */
template <CoordSystemConcept CoordSystem> class Simulation {
  protected:
    /**
     * @brief Geometry of the spin system (shared and mutable).
     *
     * Contains moment data and provides spatial access, macrocell control, etc.
     * Owned externally and passed by shared pointer.
     */
    std::shared_ptr<IGeometry<CoordSystem>> _geometry;

    /**
     * @brief Solver responsible for numerical time integration.
     *
     * Implements specific integration method (e.g., Euler, Heun, etc.).
     * @details Contains methods to update moment states and compute effective fields inside.
     */
    std::shared_ptr<ISolver<CoordSystem>> _solver;

    /**
     * @brief Registry of materials used in the system.
     *
     * Provides magnetic constants and identifiers for all atoms/spins.
     */
    std::shared_ptr<MaterialRegistry> _material_registry;

    /**
     * @brief Registry of interactions (exchange, dipole, etc.).
     *
     * Determines which interaction contributions are computed at each step.
     */
    std::shared_ptr<InteractionRegistry<CoordSystem>> _interaction_registry;

    /**
     * @brief Data buffer from the latest solver step.
     *
     * Contains effective fields and energy contributions at current state/step.
     */
    SolverData<CoordSystem> _step_solver_data;

    /**
     * @brief Accumulated trajectory of simulation states.
     *
     * Each entry contains a cloned geometry and solver data.
     */
    std::vector<SimulationStepData<CoordSystem>> steps;

    /**
     * @brief Number of simulation steps performed so far.
     */
    uint _step = 0;

    /**
     * @brief Current simulation time [s].
     */
    double _current_time = 0.0;

    /**
     * @brief Time step for numerical integration [s].
     *
     * Constant across all simulation steps.
     */
    double _dt;

    /**
     * @brief Whether the system has been prepared for interaction evaluation.
     *
     * Set to true after calling `prepareSystem()`.
     */
    bool _system_prepared = false;

    /**
     * @brief Internal preparation of system components before simulation.
     *
     * Prepares geometry layout and notifies all interactions to initialize.
     * Called automatically on first step if not prepared yet.
     *
     * @returns void – mutates internal state.
     */
    void prepareSystem() {
        SCOPED_LOG_TIMER_PRINT("Preparing system for simulation");
        {
            SCOPED_LOG_TIMER_DEBUG("├─ Preparing geometry");
            this->_geometry->prepare();
        }
        {
            SCOPED_LOG_TIMER_DEBUG("├─ Preparing interactions");
            for (auto &[_, interaction] : *_interaction_registry) {
                SCOPED_LOG_TIMER_DEBUG("│  ├─ Preparing interaction: " + interaction->getName());
                interaction->prepare(*this->_geometry, *this->_material_registry);
            }
        }
        this->_system_prepared = true;
    }

    /**
     * @brief Save a snapshot of the current state to internal buffer.
     *
     * Pushes a new `SimulationStepData` containing time, step, geometry, and solver data.
     *
     * @returns void – appends to `steps`.
     */
    void saveStep() {
        this->steps.push_back(SimulationStepData<CoordSystem>(
            this->_current_time, this->_step, this->_geometry->clone(false, true), this->_step_solver_data
        ));
    }

  public:
    /**
     * @brief Construct a simulation instance.
     *
     * @param geometry             Geometry (initial spin state and layout).
     * @param solver               Time integrator (Euler, Heun, etc.) with Field update methods.
     * @param material_registry    Registry of materials used by spins.
     * @param interaction_registry Registry of magnetic interactions.
     * @param dt                   Time step for evolution [s].
     *
     * @throws std::invalid_argument if any component is invalid or empty.
     *
     * @note All components are shared by pointer and reused.
     * @note Initial snapshot is saved automatically.
     */
    Simulation(
        std::shared_ptr<IGeometry<CoordSystem>> geometry,
        std::shared_ptr<ISolver<CoordSystem>> solver,
        std::shared_ptr<MaterialRegistry> material_registry,
        std::shared_ptr<InteractionRegistry<CoordSystem>> interaction_registry,
        double dt = 1e-13
    )
        : _geometry(geometry),
          _solver(solver),
          _material_registry(material_registry),
          _interaction_registry(interaction_registry),
          _dt(dt) {
        // проверки
        if (!_geometry)
            throw std::invalid_argument("Геометрия не может быть None");
        if (!_solver)
            throw std::invalid_argument("Решатель не может быть None");
        if (!_interaction_registry || _material_registry->isEmpty())
            throw std::invalid_argument("Регистр взаимодействий не может быть None");
        if (!_material_registry || _material_registry->isEmpty())
            throw std::invalid_argument("Регистр материалов не может быть None");
        if (_dt <= 0)
            throw std::invalid_argument("Time step dt must be positive.");

        this->_step_solver_data = SolverData<CoordSystem>();
        // Инициализируем буфер эффективных полей нужного размера
        this->_step_solver_data.clear(this->_geometry->size(), *this->_interaction_registry);
        this->_step_solver_data.correct(this->_geometry->size(), *this->_interaction_registry);
        this->_step = 0; // нулевой шаг - начальная конфигурация
        // сохранить начальную конфигурацию
        this->saveStep();
    };

    /**
     * @brief Human-readable description.
     *
     * @returns String describing parameters.
     */
    std::string __str__() const { return _geometry->__str__(); };

    /**
     * @brief Perform one simulation step (advance time by dt).
     *
     * Updates geometry via solver, optionally updates macrocells,
     *   and optionally saves current step to history.
     *
     * @param save_step          Whether to save the current state snapshot.
     * @param update_macrocells  Whether to recompute macrocell averages.
     *
     * @throws std::exception on solver or geometry failure.
     *
     * @returns void – mutates system state.
     */
    void simulateOneStep(bool save_step = false, bool update_macrocells = true) {
        this->_step += 1;
        try {
            SCOPED_LOG_TIMER_PRINT(" ============ | Simulation step " + std::to_string(this->_step));
            if (!this->_system_prepared)
                this->prepareSystem();
            if (update_macrocells) {
                this->_geometry->updateMacrocells();
            }
            this->_current_time += this->_dt;
            this->_step_solver_data = this->_solver->updateMoments(
                *this->_geometry, *this->_interaction_registry, *this->_material_registry, this->_dt
            );
            if (save_step) {
                this->saveStep();
            }
        } catch (const std::exception &e) {
            LOG_MSG_PRINT("Error during simulation step: " + std::string(e.what()));
            throw e;
        }
    };

    /**
     * @brief Perform multiple simulation steps in a loop.
     *
     * Allows controlling the frequency of snapshot saving and macrocell updates.
     *
     * @param sim_steps                    Number of time steps to simulate.
     * @param save_every_step              Save snapshot every N steps (default: 1).
     * @param update_macrocells_every_step Update macrocell averages every N steps (default: 1).
     *
     * @returns void – mutates system state.
     */
    void simulateManySteps(uint sim_steps, uint save_every_step = 1, uint update_macrocells_every_step = 1) {
        for (uint i = 0; i < sim_steps; ++i) {
            simulateOneStep(i % save_every_step == 0, i % update_macrocells_every_step == 0);
        }
    };

    /**
     * @brief Access stored simulation steps (trajectory).
     *
     * @returns Reference to internal buffer of snapshots.
     */
    std::vector<SimulationStepData<CoordSystem>> &getSteps() { return this->steps; }

    /**
     * @brief Clear all stored simulation steps (snapshots).
     *
     * @note Useful to free memory in long simulations.
     *
     * @returns void – mutates internal buffer.
     */
    void clearSteps() { return this->steps.clear(); }

    // TODO: pop_back last step...
    // TODO: getStepData(uint step)
    // TODO: оптимизация через batching-вызовы взаимодействий, решателя и пр.
};

namespace cartesian {

/**
 * @class AbstractSimulation
 * @brief Simulation orchestrator for time evolution of spin systems.
 *
 * Central entry point for performing spin dynamics simulations.
 *
 * Encapsulates the evolving state of the system (geometry, solver, interactions)
 *   and provides methods to step the simulation forward in time.
 *
 * The class manages:
 *   - The evolving spin geometry (positions, directions, materials),
 *   - A solver instance responsible for time integration,
 *   - Registries for materials and magnetic interactions,
 *   - Discrete time tracking (`step`, `time`, `dt`),
 *   - Buffered results for post-analysis or output (via `SimulationStepData`).
 *
 * This class implements a minimal but thread-safe interface for:
 *   - Initial system preparation (geometry and interactions),
 *   - Advancing the system by one or more steps,
 *   - Saving snapshots of the current system state,
 *   - Accessing simulation history or clearing memory.
 *
 * @note Units are SI unless otherwise noted.
 * @note Geometry, solver, and registries are passed as shared pointers to allow reuse and sharing.
 */
using AbstractSimulation = Simulation<NamespaceCoordSystem>;

}; // namespace cartesian

}; // namespace PYTHON_API spindynapy

// ===========================================================================
//  Python bindings
// ===========================================================================

#define SIMULATION_STEPDATA_TEMPLATE_BINDINGS(cls)                                                           \
    .def_readonly("time", &cls::time, py::doc("@brief Simulation time at this step [s] from zero."))         \
        .def_readonly(                                                                                       \
            "step", &cls::step, py::doc("@brief Step index (discrete time index, starting from zero).")      \
        )                                                                                                    \
        .def_property_readonly(                                                                              \
            "geometry",                                                                                      \
            [](const cls &self) { return self.geometry.get(); },                                             \
            py::return_value_policy::reference_internal,                                                     \
            py::doc("@brief Snapshot of the geometry (moments and layout) at this step.\n"                   \
                    "\n"                                                                                     \
                    "This is a deep clone (without cache) of the original geometry, used for replay and "    \
                    "analysis.")                                                                             \
        )                                                                                                    \
        .def(                                                                                                \
            "get_mean_magnetization",                                                                        \
            &cls::getMeanMagnetization,                                                                      \
            py::doc("Get the mean magnetization of the system at this step.\n"                               \
                    "\n"                                                                                     \
                    "Calculates the vector average of all spin directions across the geometry.\n"            \
                    "Formula: ⟨𝐌⟩ = (1/N) ∑ᵢ 𝐒ᵢ\n"                                           \
                    "@returns Average magnetization unit-vector.")                                           \
        )                                                                                                    \
        .def(                                                                                                \
            "get_mean_magnetization_norm",                                                                   \
            &cls::getMeanMagnetizationNorm,                                                                  \
            py::doc("Get the norm of the mean magnetization of the system at this step.\n"                   \
                    "\n"                                                                                     \
                    "Indicates how aligned the system is.\n"                                                 \
                    "Value ranges from 0 (disordered) to 1 (fully aligned).\n"                               \
                    "@returns Norm (length) of the average magnetization vector.")                           \
        )                                                                                                    \
        .def(                                                                                                \
            "get_energy",                                                                                    \
            &cls::getEnergy,                                                                                 \
            py::doc("Get the total energy of the system at this step.\n"                                     \
                    "\n"                                                                                     \
                    "Sums up all scalar energy contributions from all interactions:\n"                       \
                    "E_total = ∑ᵢ Eᵢ\n"                                                                \
                    "@returns Scalar total energy [J].")                                                     \
        )                                                                                                    \
        .def(                                                                                                \
            "get_energy_by_interaction",                                                                     \
            &cls::getEnergyByInteraction,                                                                    \
            py::doc("Get energy breakdown per interaction at this step.\n"                                   \
                    "\n"                                                                                     \
                    "Returns a map from interaction index to its total contribution:\n"                      \
                    "Eᵢ = ∑ⱼ Eᵢⱼ for interaction i and moment j.\n"                                \
                    "@returns Map from interaction ID to total energy [J].")                                 \
        )

#define SIMULATION_TEMPLATE_BINDINGS(cls)                                                                    \
    .def(                                                                                                    \
        "__str__",                                                                                           \
        &cls::__str__,                                                                                       \
        py::doc("@brief Return a string representation of the simulation parameters.\n"                      \
                "@returns Human-readable parameter summary.")                                                \
    )                                                                                                        \
        .def(                                                                                                \
            "simulate_one_step",                                                                             \
            &cls::simulateOneStep,                                                                           \
            py::arg("save_step") = false,                                                                    \
            py::arg("update_macrocells") = true,                                                             \
            py::call_guard<py::gil_scoped_release>(),                                                        \
            py::doc("@brief Perform one simulation step.\n"                                                  \
                    "\n"                                                                                     \
                    "Advances the simulation by one time step (dt),\n"                                       \
                    "  updating geometry and solver state.\n"                                                \
                    "\n"                                                                                     \
                    "@param save_step          Whether to save the current state snapshot.\n"                \
                    "@param update_macrocells  Whether to recompute macrocell averages.\n")                  \
        )                                                                                                    \
        .def(                                                                                                \
            "simulate_many_steps",                                                                           \
            &cls::simulateManySteps,                                                                         \
            py::arg("steps"),                                                                                \
            py::arg("save_every_step") = 1,                                                                  \
            py::arg("update_macrocells_every_step") = 1,                                                     \
            py::call_guard<py::gil_scoped_release>(),                                                        \
            py::doc("@brief Perform multiple simulation steps in a loop.\n"                                  \
                    "\n"                                                                                     \
                    "Allows controlling the frequency of snapshot saving and macrocell updates.\n"           \
                    "\n"                                                                                     \
                    "@param steps                   Number of time steps to simulate.\n"                     \
                    "@param save_every_step         Save snapshot every N steps (default: 1).\n"             \
                    "@param update_macrocells_every_step Update macrocell averages every N steps (default: " \
                    "1).\n")                                                                                 \
        )                                                                                                    \
        .def(                                                                                                \
            "get_steps",                                                                                     \
            &cls::getSteps,                                                                                  \
            py::return_value_policy::reference,                                                              \
            py::doc("@brief Access stored simulation steps (trajectory).\n"                                  \
                    "\n"                                                                                     \
                    "@returns Reference to internal buffer of snapshots.")                                   \
        )                                                                                                    \
        .def(                                                                                                \
            "clear_steps",                                                                                   \
            &cls::clearSteps,                                                                                \
            py::doc("@brief Clear all stored simulation steps (snapshots).\n"                                \
                    "\n"                                                                                     \
                    "@note Useful to free memory in long simulations.\n"                                     \
                    "@returns void – mutates internal buffer.")                                              \
        )

/**
 * @brief Bind the simulation utilities to a Python sub-module.
 *
 * @param module Parent PyBind11 module (usually the core extension module).
 * @returns void – extends the parent module.
 */
inline void pyBindSimulation(py::module_ &module) {
    using namespace spindynapy;

    // -------- | SIMULATION | --------
    py::module_ simulation_module = module.def_submodule("simulation");

    simulation_module.doc() =
        "@brief  Simulation interface and core integration logic for spin dynamics systems.\n"
        "\n"
        "This header defines the core simulation controller that integrates:\n"
        "  - spatial geometry and magnetic moment layout,\n"
        "  - solver strategies for time evolution,\n"
        "  - interaction registry (exchange, anisotropy, dipolar, etc.),\n"
        "  - material registry (atomic parameters),\n"
        "  - time integration loop (with macrocell refresh, caching, etc.).\n"
        "\n"
        "Functional overview:\n"
        "\n"
        "- `Simulation<CoordSystem>`:\n"
        "  - Owns the full simulation state: geometry, solver, interactions, materials.\n"
        "  - Advances the system state using a time integrator (`ISolver`).\n"
        "  - Manages step-by-step snapshots (`SimulationStepData`).\n"
        "  - Provides thread-safe execution with optional caching and macrocell updates.\n"
        "\n"
        "- `SimulationStepData<CoordSystem>`:\n"
        "  - Records solver buffer, geometry snapshot, magnetization, and energy at each saved step.\n"
        "  - Used for trajectory postprocessing, diagnostics, and time-series export.\n"
        "\n"
        "Python bindings:\n"
        "- Submodule: `simulation.cartesian`\n"
        "- Types exposed:\n"
        "  - `Simulation` (core loop),\n"
        "  - `SimulationStepData` (individual step snapshot).\n"
        "\n"
        "Interface:\n"
        "- `Simulation<CoordSystem>`     – time integration and state evolution interface,\n"
        "- `SimulationStepData<CoordSystem>` – snapshot of state at a simulation step.\n"
        "\n"
        "Concrete:\n"
        "- `cartesian::Simulation`       – simulation in Cartesian coordinate system,\n"
        "- `cartesian::SimulationStepData` – corresponding Cartesian snapshot.\n"
        "\n"
        "@note All units are SI unless explicitly stated.\n"
        "      Time is in [s], magnetic field in [T], energy in [J], position in [m].\n"
        "\n"
        "@note Geometry, solver, and interaction registry must be preconfigured before simulation starts.\n";

    // -------- | CARTESIAN SIMULATION | --------
    py::module_ cartesian = simulation_module.def_submodule("cartesian");
    {
        using cartesian::AbstractGeometry;
        using cartesian::AbstractInteractionRegistry;
        using cartesian::AbstractSimulation;
        using cartesian::AbstractSimulationStepData;
        using cartesian::AbstractSolver;
        using cartesian::AbstractSolverData;

        py::class_<AbstractSimulationStepData, AbstractSolverData>(cartesian, "SimulationStepData")
            SIMULATION_STEPDATA_TEMPLATE_BINDINGS(AbstractSimulationStepData)
                .doc() =
            "@struct AbstractSimulationStepData\n"
            "@brief Snapshot of simulation state at a given time step.\n"
            "\n"
            "Represents a saved snapshot of the spin system at a discrete time step during simulation.\n"
            "@details Inherits from `SolverData` and extends it with temporal metadata and a full geometry "
            "snapshot.\n"
            "\n"
            "This structure contains:\n"
            "  - Simulation time and step number,\n"
            "  - Full clone of geometry at this time step,\n"
            "  - Solver-specific data (effective fields, energies, etc.) for saved step.\n"
            "\n"
            "Used for:\n"
            "  - Trajectory logging,\n"
            "  - Postprocessing (e.g., time-series export),\n"
            "  - State resumption,\n"
            "  - Visualization.\n"
            "\n"
            "@note All fields are stored in SI units.\n";

        py::class_<AbstractSimulation>(cartesian, "Simulation")
            .def(
                py::init<
                    std::shared_ptr<AbstractGeometry>,
                    std::shared_ptr<AbstractSolver>,
                    std::shared_ptr<MaterialRegistry>,
                    std::shared_ptr<AbstractInteractionRegistry>,
                    double>(),
                py::arg("geometry"),
                py::arg("solver"),
                py::arg("material_registry"),
                py::arg("interaction_registry"),
                py::arg("dt") = 1e-13
            ) SIMULATION_TEMPLATE_BINDINGS(AbstractSimulation)
            .doc() = "@class AbstractSimulation\n"
                     "@brief Simulation orchestrator for time evolution of spin systems.\n"
                     "\n"
                     "Central entry point for performing spin dynamics simulations.\n"
                     "\n"
                     "Encapsulates the evolving state of the system (geometry, solver, interactions)\n"
                     "  and provides methods to step the simulation forward in time.\n"
                     "\n"
                     "The class manages:\n"
                     "  - The evolving spin geometry (positions, directions, materials),\n"
                     "  - A solver instance responsible for time integration,\n"
                     "  - Registries for materials and magnetic interactions,\n"
                     "  - Discrete time tracking (`step`, `time`, `dt`),\n"
                     "  - Buffered results for post-analysis or output (via `SimulationStepData`).\n"
                     "\n"
                     "This class implements a minimal but thread-safe interface for:\n"
                     "  - Initial system preparation (geometry and interactions),\n"
                     "  - Advancing the system by one or more steps,\n"
                     "  - Saving snapshots of the current system state,\n"
                     "  - Accessing simulation history or clearing memory.\n"
                     "\n"
                     "@note Units are SI unless otherwise noted.\n"
                     "@note Geometry, solver, and registries are passed as shared pointers to allow reuse "
                     "and sharing.\n";
    }
}

#endif // ! __SIMULATION_HPP__
