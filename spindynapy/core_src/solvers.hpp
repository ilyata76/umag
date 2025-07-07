#ifndef __SOLVERS_HPP__
#define __SOLVERS_HPP__

/**
 * @file   solvers.hpp
 * @brief  Explicit time-integration solvers for atomistic spin dynamics (LLG-based).
 *
 * This header defines the solver interface (`ISolver`) and concrete integrators
 *   (`LLGSolver`, etc.).
 * It also introduces `IFieldUpdater` for computing effective fields and concrete updaters.
 *
 * Conceptual overview:
 * - A **Solver** implements a time-stepping algorithm (Euler, Heun, etc.) that updates the
 *   directions of magnetic moments according to computed effective fields.
 * - A **FieldUpdater** computes effective fields and energy densities for all spins
 *   given the current geometry, interaction set, and material parameters.
 * - The `SolverData` structure acts as a per-step container for intermediate results
 *   (net field, total and per-interaction energy), decoupling computation from output.
 *
 * All solvers operate on geometries satisfying the `IGeometry` interface and are
 *   parametrised by a coordinate system (`CoordSystemConcept`) to allow future extension.
 *
 * Exposed entities:
 * - `SolverStrategy` – enum describing integration schemes (Euler, Heun…)
 * - `SolverData<CoordSystem>` – storage class for fields and energies
 * - `Solver<CoordSystem>` – base class for solvers
 * - `FieldUpdater<CoordSystem>` – base class for field evaluators
 *
 * Concrete:
 * - `LLGSolver` – default LLG-based integrator using Euler or Heun method
 * - `OMPFieldUpdater` – OpenMP-based field updater implementation
 *
 * - Python binding macros: `SOLVER_TEMPLATE_BINDINGS`, `SOLVERDATA_TEMPLATE_BINDINGS`,
 *   `FIELDUPDATER_TEMPLATE_BINDINGS`
 * - Python submodule binder: `pyBindSolvers()`
 *
 * @note All units are SI. Solver classes mutate geometry in-place.
 *       The update steps are not atomic – parallelisation responsibility lies within the updater.
 *
 * @copyright 2025 SpinDynaPy
 */

#include "constants.hpp"
#include "geometries.hpp"
#include "interactions.hpp"
#include "logger.hpp"
#include "registries.hpp"
#include "types.hpp"

#include <Eigen/Dense>
#include <cassert>
#include <cmath>
#include <memory>
#ifdef _OPENMP // если OpenMP доступен (флаг компиляции, заголовки)
#include <omp.h>
#endif
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <stdexcept>

namespace py = pybind11;

namespace PYTHON_API spindynapy {

/**
 * @enum   SolverStrategy
 * @brief  Numerical integration scheme used by a solver.
 *
 * The supported explicit integrators are ordered by increasing accuracy (and computational cost):
 *   - EULER    First-order Euler step.
 *   - MIDPOINT Classical mid-point (Euler predictor, no extra field eval).
 *   - HEUN     Second-order Heun predictor–corrector (two field evals).
 */
enum class PYTHON_API SolverStrategy {
    EULER = 0,    ///< First-order Euler method.
    MIDPOINT = 1, ///< Mid-point method ( (F+S)/2 ).
    HEUN = 2      ///< Heun (predictor-corrector) method.
};

// ==========================================================================
//  SolverData
// ==========================================================================

/**
 * @struct SolverData
 * @brief  Per-step buffer for effective fields and energy terms.
 *
 * @tparam CoordSystem  Any type satisfying @ref CoordSystemConcept (e.g. Cartesian).
 *
 * A solver typically performs two phases per time-step:
 * 1. Field update – compute H_eff and interaction-resolved contributions.
 * 2. Moment update – integrate the LLG equation using those fields.
 *
 * `SolverData` stores all intermediate results so they can be inspected, logged, or fed into a
 *   post-processing pipeline.
 *
 * Layout
 * - `effective_fields[i]`         → net H_eff acting on spin *i*.
 * - `energies[i]`                 → net energy density for spin *i*.
 * - `interaction_effective_fields[reg][i]` → per-interaction H_eff.
 * - `interaction_energies[reg][i]`         → per-interaction energy term.
 */
template <CoordSystemConcept CoordSystem> struct PYTHON_API SolverData {
  public:
    /** @brief Net effective field (sum) on every spin (indexed by geometry order). */
    PYTHON_API EffectiveFieldVector effective_fields;
    /** @brief Net energy (sum) for every spin (J). */
    PYTHON_API std::vector<double> energies;
    /**
     * @brief Map: interaction registry id → vector of field contributions on every *i* spin.
     * Each entry mirrors `effective_fields` but for a single interaction only.
     */
    PYTHON_API std::unordered_map<regnum, EffectiveFieldVector> interaction_effective_fields;
    /**
     * @brief Map: interaction registry id → vector of energy contributions on every *i* spin.
     * Units are joules; length matches `energies`.
     */
    PYTHON_API std::unordered_map<regnum, std::vector<double>> interaction_energies;

    /**
     * @brief  Default constructor – creates *empty* buffers.
     *
     * @note The struct is intentionally initialised with **zero-sized containers**.
     *  The owning solver must call `clear()`/`correct()` prior to
     *  use so that the internal buffers match the geometry size.
     *
     * @returns Newly constructed *SolverData* with zero-initialised members.
     */
    SolverData() {};

    /**
     * @brief  Fully parameterised constructor (shallow copy of buffers).
     *
     * All buffers are taken *by value* to avoid dangling references
     *
     * @param effective_fields             Effective fields per spin.
     * @param energies                     Energies per spin.
     * @param interaction_effective_fields Detailed field map per interaction.
     * @param interaction_energies         Detailed energy map per interaction.
     *
     * @returns New *SolverData* populated with provided buffers.
     */
    SolverData(
        const EffectiveFieldVector &effective_fields,
        std::vector<double> energies,
        const std::unordered_map<regnum, EffectiveFieldVector> &interaction_effective_fields,
        const std::unordered_map<regnum, std::vector<double>> &interaction_energies
    )
        : effective_fields(effective_fields),
          energies(energies),
          interaction_effective_fields(interaction_effective_fields),
          interaction_energies(interaction_energies) {}

    /**
     * @brief  Reset all buffers to zero-filled state Sized for a geometry.
     *
     * @param moments_size         Number of spins in the geometry.
     * @param interaction_registry Registry of active interactions (used to size the per-interaction maps).
     *
     * @returns void – *mutates the internal state* (buffers are resized and zero-initialised).
     */
    PYTHON_API void clear(size_t moments_size, InteractionRegistry<CoordSystem> &interaction_registry) {
        this->effective_fields = EffectiveFieldVector(moments_size, EffectiveField::Zero());
        this->energies = std::vector<double>(moments_size, 0.0);
        this->interaction_effective_fields.clear();
        this->interaction_energies.clear();
        for (auto &[interaction_regnum, interaction] : interaction_registry) {
            this->interaction_effective_fields[interaction_regnum] =
                EffectiveFieldVector(moments_size, EffectiveField::Zero());
            this->interaction_energies[interaction_regnum] = std::vector<double>(moments_size, 0.0);
        }
    };

    /**
     * @brief  Adjust buffer sizes when geometry size changes.
     *
     * If the geometry grows/shrinks (e.g., after dynamic loading or a mask
     * update) this method resizes all internal containers accordingly while
     * preserving existing data.
     *
     * @param moments_size         Number of spins in the geometry.
     * @param interaction_registry Registry of active interactions (used to size the per-interaction maps).
     *
     * @returns void – *mutates the internal state* by re-allocating buffers if needed.
     */
    PYTHON_API void correct(size_t moments_size, InteractionRegistry<CoordSystem> &interaction_registry) {
        if (moments_size == this->effective_fields.size())
            return;
        this->effective_fields.clear();
        this->effective_fields.resize(moments_size);
        this->energies.clear();
        this->energies.resize(moments_size);
        this->interaction_effective_fields.clear();
        this->interaction_energies.clear();
        for (auto &[interaction_regnum, interaction] : interaction_registry) {
            this->interaction_effective_fields[interaction_regnum].resize(moments_size);
            this->interaction_energies[interaction_regnum].resize(moments_size);
        }
    }
};

namespace cartesian {

/**
 * @struct AbstractSolverData
 * @brief  Per-step buffer for effective fields and energy terms.
 *
 * Cartesian coord system.
 *
 * A solver typically performs two phases per time-step:
 * 1. Field update – compute H_eff and interaction-resolved contributions.
 * 2. Moment update – integrate the LLG equation using those fields.
 *
 * `SolverData` stores all intermediate results so they can be inspected, logged, or fed into a
 *   post-processing pipeline.
 *
 * Layout
 * - `effective_fields[i]`         → net H_eff acting on spin *i*.
 * - `energies[i]`                 → net energy density for spin *i*.
 * - `interaction_effective_fields[reg][i]` → per-interaction H_eff.
 * - `interaction_energies[reg][i]`         → per-interaction energy term.
 */
using AbstractSolverData = PYTHON_API SolverData<NamespaceCoordSystem>;

} // namespace cartesian

// ==========================================================================
//  Field updater
// ==========================================================================

/**
 * @class  IFieldUpdater
 * @brief  Strategy interface that calculates **effective fields & energies**.
 *
 * @tparam CoordSystem Coordinate-system tag.
 *
 * A field-updater encapsulates the sum of all interactions for a given
 *   geometry.  Different implementations may employ OpenMP, CUDA, TBB…; the
 *   solver is agnostic.
 */
template <CoordSystemConcept CoordSystem> class PYTHON_API IFieldUpdater {
  protected:
    /**
     * @brief Protected default constructor.
     */
    IFieldUpdater() = default;

  public:
    /**
     * @brief Virtual destructor (default).
     */
    virtual ~IFieldUpdater() = default;

    /**
     * @brief Compute effective fields & energies for *current* geometry state.
     *
     * @param geometry             Geometry for which to compute fields.
     * @param interaction_registry Active interactions.
     * @param material_registry    Materials referenced by `geometry`.
     *
     * @returns `SolverData` containing freshly calculated fields and energies.
     */
    PYTHON_API virtual SolverData<CoordSystem> calculateFields(
        IGeometry<CoordSystem> &geometry,
        InteractionRegistry<CoordSystem> &interaction_registry,
        MaterialRegistry &material_registry
    ) = 0;
};

namespace cartesian {

/**
 * @class  AbstractFieldUpdater
 * @brief  Strategy interface that calculates **effective fields & energies**.
 *
 * Cartesian coord system.
 *
 * A field-updater encapsulates the sum of all interactions for a given
 *   geometry.  Different implementations may employ OpenMP, CUDA, TBB…; the
 *   solver is agnostic.
 */
using AbstractFieldUpdater = PYTHON_API IFieldUpdater<NamespaceCoordSystem>;

/**
 * @class  OMPFieldUpdater
 * @brief  OpenMP-parallel field updater for Cartesian geometries.
 *
 * A field-updater encapsulates the sum of all interactions for a given geometry.
 * The class parallelises *per-spin* interaction contributions with OMP pragma.
 */
class PYTHON_API OMPFieldUpdater : public AbstractFieldUpdater {
  protected:
    AbstractSolverData _step_data; ///< Internal buffer reused every call (minimises re-allocations).

    /**
     * @brief Calculate field & energy contribution from a single interaction.
     *
     * Helper invoked inside the outer loop over `interaction_registry`.
     * Thread-safe: writes into local vectors then enters a critical section for
     *   the final reduction.
     *
     * @param geometry           Active geometry.
     * @param interaction_regnum Registry id of the interaction.
     * @param interaction        Pointer to the interaction implementation.
     * @param material_registry  Material look-up table.
     *
     * @returns void – mutates internal data.
     */
    void calculateFieldContribution(
        AbstractGeometry &geometry,
        regnum interaction_regnum,
        AbstractInteraction *interaction,
        MaterialRegistry &material_registry
    ) {
        SCOPED_LOG_TIMER_DEBUG(
            "│  ├─ Calculate effective fields and energies for interaction " + interaction->getName()
        );
        // --- предсоздать буфферы ---

        size_t moments_size = geometry.size(); // Размер системы
        EffectiveFieldVector contribution_vector(moments_size
        ); // Вектор для полевого вклада текущего interaction
        std::vector<double> energy(moments_size); // Вектор для энергетич. вклада текущего interaction

        // --- посчитать энергии и эфф. поля в CPU-параллели ---

        // clang-format off
        #pragma omp parallel for schedule(dynamic)
        // clang-format on
        for (size_t i = 0; i < moments_size; ++i) {
            contribution_vector[i] = interaction->calculateFieldContribution(i, geometry, material_registry);
            energy[i] = interaction->calculateEnergy(geometry[i], contribution_vector[i]);
        }

        // --- проверить на правильность результатов ---

        for (size_t i = 0; i < moments_size; ++i) {
            if (contribution_vector[i].hasNaN()) {
                throw std::runtime_error("Invalid contribution_vector data (NaN detected)");
            }
        }

        // --- записать в текущее шаговое состояние (гарантировать последовательность записи) ---

        // clang-format off
        #pragma omp critical
        // clang-format on
        for (size_t i = 0; i < moments_size; ++i) {
            this->_step_data.effective_fields[i] += contribution_vector[i];
            this->_step_data.interaction_effective_fields[interaction_regnum][i] = contribution_vector[i];
            // и энергии
            this->_step_data.energies[i] += energy[i];
            this->_step_data.interaction_energies[interaction_regnum][i] = energy[i];
        }
    }

  public:
    /** Default constructor – no state. */
    PYTHON_API OMPFieldUpdater() {};

    /**
     * @brief Compute effective fields & energies for *current* geometry state.
     *
     * @param geometry             Geometry for which to compute fields.
     * @param interaction_registry Active interactions.
     * @param material_registry    Materials referenced by `geometry`.
     *
     * @returns `SolverData` containing freshly calculated fields and energies.
     */
    PYTHON_API virtual AbstractSolverData calculateFields(
        AbstractGeometry &geometry,
        AbstractInteractionRegistry &interaction_registry,
        MaterialRegistry &material_registry
    ) override {
        SCOPED_LOG_TIMER_DEBUG("├─ Calculate Effective Fields and Energies");
        size_t moments_size = geometry.size();                      // Размер системы
        this->_step_data.clear(moments_size, interaction_registry); // Занулить буфер
        this->_step_data.correct(moments_size, interaction_registry); // Не потерять входящие изменения

        // обсчёт полей и энергий от зарегистрированных взаимодействий с записью во внутреннее состояние
        for (auto &[interaction_regnum, interaction] : interaction_registry) {
            this->calculateFieldContribution(
                geometry, interaction_regnum, interaction.get(), material_registry
            );
        }

        return this->_step_data;
    }
};

}; // namespace cartesian

// ==========================================================================
//  Solvers
// ==========================================================================

/**
 * @class  ISolver
 * @brief  Abstract interface for spin dynamics solvers in a given coordinate system.
 *
 * @tparam CoordSystem  Coordinate-system tag that satisfies @ref CoordSystemConcept.
 *
 * Conceptual overview:
 *  A solver integrates the system's magnetic moments over time by computing and evolving
 *    the spin directions according to a chosen time-integration scheme (Euler,
 *    Heun, etc.). This interface decouples the integration policy from the
 *    geometry and interaction definitions.
 *
 * The solver modifies the geometry in-place and returns a full snapshot of the
 *   effective fields and energy contributions used in the step. This enables
 *   inspection, logging, or post-processing without recomputing.
 *
 * @note Concrete implementations must implement the method `updateMoments()`.
 *       They are responsible for allocating intermediate data and managing
 *       consistency between field updates and geometry states.
 */
template <CoordSystemConcept CoordSystem> class PYTHON_API ISolver {
  protected:
    /**
     * @brief Protected default constructor to prevent direct instantiation.
     */
    ISolver() = default;

  public:
    /**
     * @brief Virtual destructor (default).
     */
    virtual ~ISolver() = default;

    /**
     * @brief Advance the geometry by one time-step.
     *
     * The implementation must: (i) compute effective fields via its
     * `IFieldUpdater`, (ii) integrate equations or others according to its
     * `SolverStrategy`, and (iii) write updated spin directions back into
     * `geometry`.
     *
     * @param geometry             Target geometry (spins mutated in-place).
     * @param interaction_registry Registry of active interactions.
     * @param material_registry    Registry of materials (look-ups only).
     * @param dt                   Time-step (s).
     *
     * @returns A fully populated `SolverData` snapshot corresponding to the
     *          fields/energies *used* during this step.
     */
    PYTHON_API virtual SolverData<CoordSystem> updateMoments(
        IGeometry<CoordSystem> &geometry,
        InteractionRegistry<CoordSystem> &interaction_registry,
        MaterialRegistry &material_registry,
        double dt
    ) = 0;
};

namespace cartesian {

/**
 * @class  AbstractSolver
 * @brief  Abstract interface for spin dynamics solvers in a given coordinate system.
 *
 * Cartesian coord system.
 *
 * Conceptual overview:
 *  A solver integrates the system's magnetic moments over time by computing and evolving
 *    the spin directions according to a chosen time-integration scheme (Euler,
 *    Heun, etc.). This interface decouples the integration policy from the
 *    geometry and interaction definitions.
 *
 * The solver modifies the geometry in-place and returns a full snapshot of the
 *   effective fields and energy contributions used in the step. This enables
 *   inspection, logging, or post-processing without recomputing.
 *
 * @note Concrete implementations must implement the method `updateMoments()`.
 *       They are responsible for allocating intermediate data and managing
 *       consistency between field updates and geometry states.
 */
using AbstractSolver = PYTHON_API ISolver<NamespaceCoordSystem>;

/**
 * @class LLGSolver
 * @brief Solver for time evolution of magnetic systems via the Landau-Lifshitz-Gilbert (LLG) equation.
 *
 * The LLGSolver class implements the core of explicit time-integration for atomistic spin dynamics
 *   using the Landau-Lifshitz-Gilbert formalism. It combines a numerical integrator (Euler, Heun)
 *   with a field update strategy to simulate the evolution of spin systems in time.
 *
 *  Formula: dS/dt = −γ / (1 + α²) [ S × H_eff + α S × (S × H_eff) ]
 *
 * This solver mutates the geometry in place, updating the directions of magnetic moments using
 *   effective fields derived from registered interactions. It is responsible for both computing
 *   the spin derivatives (`dS/dt`) and applying them to the geometry using a selected time-stepping method.
 *
 * Internally, the solver delegates effective field computation to an Field Updater and
 *   integrates the LLG equation using either Strategy methods.
 *
 * @note All physical units are assumed to be in SI. Spin vectors are expected to be normalized.
 */
class PYTHON_API LLGSolver : public AbstractSolver {
  protected:
    SolverStrategy _strategy; ///< Numerical integration scheme used by this solver (EULER, HEUN, etc.).
    std::unique_ptr<AbstractFieldUpdater>
        _field_updater; ///< Field updater strategy used to compute effective fields.

    /**
     * @brief Compute the time derivative of a spin (dS/dt) using the LLG equation.
     *
     * This method implements the core of the LLG right-hand side:
     *     dS/dt = −γ / (1 + α²) [ S × H_eff + α S × (S × H_eff) ]
     *   where γ is the gyromagnetic ratio and α is the damping constant.
     *
     * @param material         Material parameters (gyromagnetic ratio, damping).
     * @param moment_vector    Current spin vector (unit length).
     * @param effective_field  Local effective magnetic field (SI units).
     *
     * @returns The instantaneous spin change vector dS/dt (unitless).
     */
    Eigen::Vector3d calculateLLGMomentChange(
        Material &material, Eigen::Vector3d moment_vector, Eigen::Vector3d effective_field
    ) {
        auto prefix_term = material.getGyromagneticRatio() / (1.0 + pow(material.getDampingConstant(), 2));

        Eigen::Vector3d moment_change =
            -(prefix_term * moment_vector.cross(effective_field) +
              prefix_term * material.getDampingConstant() *
                  moment_vector.cross(moment_vector.cross(effective_field)));

        return moment_change;
    }

    /**
     * @brief Compute LLG moment changes for the entire geometry.
     *
     * Iterates over each spin in the geometry, using its material properties and the local
     *   effective field (from `SolverData`) to compute the LLG derivative vector.
     *
     * This method implements the core of the LLG right-hand side:
     *     dS/dt = −γ / (1 + α²) [ S × H_eff + α S × (S × H_eff) ]
     *   where γ is the gyromagnetic ratio and α is the damping constant.
     *
     * @param geometry  Magnetic geometry containing spin moments.
     * @param data      Snapshot of effective fields and energies.
     *
     * @returns Vector of dS/dt changes for each spin.
     */
    virtual std::vector<Eigen::Vector3d> calculateLLGMomentsChange(
        IGeometry<NamespaceCoordSystem> &geometry, AbstractSolverData &data
    ) {
        SCOPED_LOG_TIMER_DEBUG("├─ Calculate LLG moments change (all spins)");
        auto moments_size = geometry.size();
        std::vector<Eigen::Vector3d> result;
        result.resize(moments_size, Eigen::Vector3d::Zero());
        for (size_t i = 0; i < moments_size; ++i) {
            auto &moment = geometry[i];
            auto &base_moment_vector = moment.getDirection().asVector();

            result[i] = this->calculateLLGMomentChange(
                moment.getMaterial(), base_moment_vector, data.effective_fields[i]
            );
        }
        return result;
    }

    /**
     * @brief Apply spin updates to geometry.
     *
     * Uses finite difference: S(t+dt) = S(t) + dt * dS/dt
     *
     * @param geometry        The system geometry to mutate.
     * @param moment_changes  Computed dS/dt values.
     * @param dt              Time-step (s).
     *
     * @returns Mutates the geometry in place.
     */
    virtual void updateMomentsInGeometry(
        IGeometry<NamespaceCoordSystem> &geometry, std::vector<Eigen::Vector3d> &moment_changes, double dt
    ) {
        SCOPED_LOG_TIMER_DEBUG("├─ Update moments in geometry (all spins)");
        auto moments_size = geometry.size();
        if (moment_changes.size() != moments_size) {
            throw std::invalid_argument("Invalid moment changes size");
        }
        for (size_t i = 0; i < moments_size; ++i) {
            auto &moment = geometry[i];
            Eigen::Vector3d new_moment = moment.getDirection().asVector() + moment_changes[i] * dt;
            moment.setDirection(new_moment);
        }
    }

    /**
     * @brief Perform a single integration step using the Euler method.
     *
     * Effective fields are computed once, then the LLG equation is integrated explicitly
     *   using the computed dS/dt to advance the spin directions.
     *
     * @param geometry             System geometry to mutate.
     * @param interaction_registry Interaction registry providing contributions.
     * @param material_registry    Material lookup registry.
     * @param dt                   Time-step (s).
     *
     * @returns `SolverData` corresponding to the computed effective fields.
     */
    virtual AbstractSolverData updateMomentsViaEuler(
        IGeometry<NamespaceCoordSystem> &geometry,
        InteractionRegistry<NamespaceCoordSystem> &interaction_registry,
        MaterialRegistry &material_registry,
        double dt
    ) {
        SCOPED_LOG_TIMER_DEBUG("EULER STEP");
        auto solver_data =
            this->_field_updater->calculateFields(geometry, interaction_registry, material_registry);
        auto moments_changes = this->calculateLLGMomentsChange(geometry, solver_data);
        this->updateMomentsInGeometry(geometry, moments_changes, dt);
        return solver_data;
    }

    /**
     * @brief Perform a single integration step using the Heun (predictor-corrector) method.
     *
     * Two field evaluations are used: one on the base geometry and one on the predicted geometry
     *   to increase accuracy. The spin change is then averaged and applied to the original state.
     *
     * @param geometry             System geometry to mutate.
     * @param interaction_registry Interaction registry providing contributions.
     * @param material_registry    Material lookup registry.
     * @param dt                   Time-step (s).
     *
     * @returns `SolverData` from the predictor pass.
     */
    virtual AbstractSolverData updateMomentsViaHeun(
        IGeometry<NamespaceCoordSystem> &geometry,
        InteractionRegistry<NamespaceCoordSystem> &interaction_registry,
        MaterialRegistry &material_registry,
        double dt
    ) {
        SCOPED_LOG_TIMER_DEBUG("HEUN STEP");
        // рассчитать поля и изменения для предиктора
        auto predictor_data =
            this->_field_updater->calculateFields(geometry, interaction_registry, material_registry);
        auto predictor_geometry = geometry.clone(true);
        auto predictor_moments_changes = this->calculateLLGMomentsChange(geometry, predictor_data);

        // рассчитать новое направление спина (predictor_geometry изменилась)
        this->updateMomentsInGeometry(*predictor_geometry, predictor_moments_changes, dt);
        // рассчитать новые поля для корректора на основе изменённой геометрии предиктором
        auto corrector_data = this->_field_updater->calculateFields(
            *predictor_geometry, interaction_registry, material_registry
        );
        auto corrector_moments_changes = this->calculateLLGMomentsChange(*predictor_geometry, corrector_data);

        // усреднить изменения dS
        std::vector<Eigen::Vector3d> avg_moments_changes(geometry.size());
        for (size_t i = 0; i < geometry.size(); ++i) {
            avg_moments_changes[i] = (predictor_moments_changes[i] + corrector_moments_changes[i]) * 0.5;
        }

        // посчитать изменение на обобщённом dS и начальной геометрии
        this->updateMomentsInGeometry(geometry, avg_moments_changes, dt);

        // вернуть данные, которые привели к текущей ориентации (предиктор) TODO: усреднить пред+корр
        return predictor_data;
    }

  public:
    /**
     * @brief Construct a new LLGSolver.
     *
     * Initializes the solver with a given integration strategy and a default OpenMP-based field updater.
     *
     * @param strategy  Integration method to use (defaults to Euler).
     */
    PYTHON_API LLGSolver(SolverStrategy strategy = SolverStrategy::EULER)
        : _strategy(strategy), _field_updater(std::make_unique<OMPFieldUpdater>()) {};

    /**
     * @brief Advance the system by one time-step.
     *
     * Delegates field computation to the field updater, then integrates using
     *   the configured strategy (Euler or Heun). The spin moments are updated in place.
     * @note Mutates the geometry in-place.
     *
     * @param geometry             Magnetic geometry to update.
     * @param interaction_registry Registry of interactions (used for field computation).
     * @param material_registry    Registry of material properties (read-only).
     * @param dt                   Simulation time-step (s).
     *
     * @returns Snapshot of the solver state (fields and energies) that were used.
     */
    PYTHON_API virtual AbstractSolverData updateMoments(
        IGeometry<NamespaceCoordSystem> &geometry,
        InteractionRegistry<NamespaceCoordSystem> &interaction_registry,
        MaterialRegistry &material_registry,
        double dt
    ) override {
        //
        if (this->_strategy == SolverStrategy::EULER) {
            return this->updateMomentsViaEuler(geometry, interaction_registry, material_registry, dt);
        } else if (this->_strategy == SolverStrategy::HEUN) {
            return this->updateMomentsViaHeun(geometry, interaction_registry, material_registry, dt);
        };

        throw std::invalid_argument("Invalid solver strategy enum value / doesnt support");
    };
};

}; // namespace cartesian

}; // namespace PYTHON_API spindynapy

// ==========================================================================
//  Python bindings
// ==========================================================================

/**
 * @def   SOLVER_TEMPLATE_BINDINGS
 * @brief PyBind11 boilerplate for exposing Solver methods.
 */
#define SOLVER_TEMPLATE_BINDINGS(cls)                                                                        \
    .def(                                                                                                    \
        "update_moments",                                                                                    \
        &cls::updateMoments,                                                                                 \
        py::arg("geometry"),                                                                                 \
        py::arg("interaction_registry"),                                                                     \
        py::arg("material_registry"),                                                                        \
        py::arg("dt"),                                                                                       \
        py::doc("@brief Advance the geometry by one time-step.\n"                                            \
                "\n"                                                                                         \
                "The implementation must: (i) compute effective fields via its\n"                            \
                "`IFieldUpdater`, (ii) integrate equations or others according to its\n"                     \
                "`SolverStrategy`, and (iii) write updated spin directions back into\n"                      \
                "`geometry`.\n"                                                                              \
                "\n"                                                                                         \
                "@param geometry             Target geometry (spins mutated in-place).\n"                    \
                "@param interaction_registry Registry of active interactions.\n"                             \
                "@param material_registry    Registry of materials (look-ups only).\n"                       \
                "@param dt                   Time-step (s).\n"                                               \
                "\n"                                                                                         \
                "@returns A fully populated `SolverData` snapshot corresponding to the\n"                    \
                "        fields/energies *used* during this step.\n")                                        \
    )

/**
 * @def   SOLVERDATA_TEMPLATE_BINDINGS
 * @brief PyBind11 boilerplate for exposing the members of @ref SolverData.
 */
#define SOLVERDATA_TEMPLATE_BINDINGS(cls)                                                                    \
    .def_readwrite("effective_fields", &cls::effective_fields)                                               \
        .def_readwrite("energies", &cls::energies)                                                           \
        .def_readwrite("interaction_effective_fields", &cls::interaction_effective_fields)                   \
        .def_readwrite("interaction_energies", &cls::interaction_energies)                                   \
        .def(                                                                                                \
            "clear",                                                                                         \
            &cls::clear,                                                                                     \
            py::arg("moments_size"),                                                                         \
            py::arg("interaction_registry"),                                                                 \
            py::doc(                                                                                         \
                "@brief  Reset all buffers to zero-filled state Sized for a geometry.\n"                     \
                "\n"                                                                                         \
                "@param moments_size         Number of spins in the geometry.\n"                             \
                "@param interaction_registry Registry of active interactions (used to size the "             \
                "per-interaction maps).\n"                                                                   \
                "\n"                                                                                         \
                "@returns void – *mutates the internal state* (buffers are resized and zero-initialised)." \
            )                                                                                                \
        )                                                                                                    \
        .def(                                                                                                \
            "correct",                                                                                       \
            &cls::correct,                                                                                   \
            py::arg("moments_size"),                                                                         \
            py::arg("interaction_registry"),                                                                 \
            py::doc("@brief  Adjust buffer sizes when geometry size changes.\n"                              \
                    "\n"                                                                                     \
                    "If the geometry grows/shrinks (e.g., after dynamic loading or a mask\n"                 \
                    "update) this method resizes all internal containers accordingly while\n"                \
                    "preserving existing data.\n"                                                            \
                    "\n"                                                                                     \
                    "@param moments_size         Number of spins in the geometry.\n"                         \
                    "@param interaction_registry Registry of active interactions (used to size the "         \
                    "per-interaction maps).\n"                                                               \
                    "\n"                                                                                     \
                    "@returns void – *mutates the internal state* by re-allocating buffers if needed.")      \
        )

/**
 * @def   FIELDUPDATER_TEMPLATE_BINDINGS
 * @brief PyBind11 boilerplate for exposing Updater methods.
 */
#define FIELDUPDATER_TEMPLATE_BINDINGS(cls)                                                                  \
    .def(                                                                                                    \
        "calculate_fields",                                                                                  \
        &cls::calculateFields,                                                                               \
        py::arg("geometry"),                                                                                 \
        py::arg("interaction_registry"),                                                                     \
        py::arg("material_registry"),                                                                        \
        py::doc("@brief Compute effective fields & energies for *current* geometry state.\n"                 \
                "\n"                                                                                         \
                "@param geometry             Geometry for which to compute fields.\n"                        \
                "@param interaction_registry Active interactions.\n"                                         \
                "@param material_registry    Materials referenced by `geometry`.\n"                          \
                "\n"                                                                                         \
                "@returns `SolverData` containing freshly calculated fields and energies.")                  \
    )

inline void pyBindSolvers(py::module_ &module) {
    using namespace spindynapy;

    // ---------------------------------------------------------------------
    //  Create submodule
    // ---------------------------------------------------------------------

    py::module_ solvers_module = module.def_submodule("solvers");

    solvers_module.doc() =
        "@brief  Explicit time-integration solvers for atomistic spin dynamics (LLG-based).\n"
        "\n"
        "This header defines the solver interface (`ISolver`) and concrete integrators\n"
        "  (`LLGSolver`, etc.).\n"
        "It also introduces `IFieldUpdater` for computing effective fields and concrete updaters.\n"
        "\n"
        "Conceptual overview:\n"
        "- A **Solver** implements a time-stepping algorithm (Euler, Heun, etc.) that updates the\n"
        "  directions of magnetic moments according to computed effective fields.\n"
        "- A **FieldUpdater** computes effective fields and energy densities for all spins\n"
        "  given the current geometry, interaction set, and material parameters.\n"
        "- The `SolverData` structure acts as a per-step container for intermediate results\n"
        "  (net field, total and per-interaction energy), decoupling computation from output.\n"
        "All solvers operate on geometries satisfying the `IGeometry` interface and are\n"
        "  parametrised by a coordinate system (`CoordSystemConcept`) to allow future extension.\n"
        "\n"
        "Exposed entities:\n"
        "- `SolverStrategy` – enum describing integration schemes (Euler, Heun…)\n"
        "- `SolverData<CoordSystem>` – storage class for fields and energies\n"
        "- `Solver<CoordSystem>` – base class for solvers\n"
        "- `FieldUpdater<CoordSystem>` – base class for field evaluators\n"
        "Concrete:\n"
        "- `LLGSolver` – default LLG-based integrator using Euler or Heun method\n"
        "- `OMPFieldUpdater` – OpenMP-based field updater implementation\n"
        "\n"
        "- Python binding macros: `SOLVER_TEMPLATE_BINDINGS`, `SOLVERDATA_TEMPLATE_BINDINGS`, "
        "`FIELDUPDATER_TEMPLATE_BINDINGS`\n"
        "- Python submodule binder: `pyBindSolvers()`\n"
        "\n"
        "@note All units are SI. Solver classes mutate geometry in-place.\n"
        "      The update steps are not atomic – parallelisation responsibility lies within the updater.\n";

    // ---------------------------------------------------------------------
    //  Strategy for solvers
    // ---------------------------------------------------------------------

    py::enum_<SolverStrategy>(solvers_module, "SolverStrategy")
        .value("EULER", SolverStrategy::EULER, "First-order Euler method.")
        .value("HEUN", SolverStrategy::HEUN, "Heun (predictor-corrector) method.");

    // ---------------------------------------------------------------------
    //  Cartesian submodule bindings
    // ---------------------------------------------------------------------

    py::module_ cartesian = solvers_module.def_submodule("cartesian");
    cartesian.doc() = solvers_module.doc();

    {
        using cartesian::AbstractFieldUpdater;
        using cartesian::AbstractSolver;
        using cartesian::AbstractSolverData;
        using cartesian::LLGSolver;
        using cartesian::OMPFieldUpdater;

        py::class_<AbstractSolverData>(cartesian, "SolverData")
            SOLVERDATA_TEMPLATE_BINDINGS(AbstractSolverData)
                .doc() =
            "@struct AbstractSolverData\n"
            "@brief  Per-step buffer for effective fields and energy terms.\n"
            "\n"
            "Cartesian coord system.\n"
            "A solver typically performs two phases per time-step:\n"
            "1. Field update – compute H_eff and interaction-resolved contributions.\n"
            "2. Moment update – integrate the LLG equation using those fields.\n"
            "\n"
            "`SolverData` stores all intermediate results so they can be inspected, logged, or fed into a\n"
            "  post-processing pipeline.\n"
            "\n"
            "Layout\n"
            "- `effective_fields[i]`         → net H_eff acting on spin *i*.\n"
            "- `energies[i]`                 → net energy density for spin *i*.\n"
            "- `interaction_effective_fields[reg][i]` → per-interaction H_eff.\n"
            "- `interaction_energies[reg][i]`         → per-interaction energy term.";

        py::class_<AbstractFieldUpdater, std::shared_ptr<AbstractFieldUpdater>>(
            cartesian, "AbstractFieldUpdater"
        ) FIELDUPDATER_TEMPLATE_BINDINGS(AbstractFieldUpdater)
            .doc() = "Strategy interface that calculates **effective fields & energies**.\n"
                     "\n"
                     "Cartesian coord system.\n"
                     "A field-updater encapsulates the sum of all interactions for a given\n"
                     "  geometry.  Different implementations may employ OpenMP, CUDA, TBB…; the\n"
                     "  solver is agnostic.";

        py::class_<OMPFieldUpdater, AbstractFieldUpdater, std::shared_ptr<OMPFieldUpdater>>(
            cartesian, "OMPFieldUpdater"
        )
            .def(py::init<>())
            .doc() = "OpenMP-parallel field updater for Cartesian geometries.\n"
                     "\n"
                     "A field-updater encapsulates the sum of all interactions for a given geometry.\n"
                     "The class parallelises *per-spin* interaction contributions with OMP pragma.\n";

        py::class_<AbstractSolver, std::shared_ptr<AbstractSolver>>(cartesian, "AbstractSolver")
            SOLVER_TEMPLATE_BINDINGS(AbstractSolver)
                .doc() =
            "Abstract interface for spin dynamics solvers in a given coordinate system.\n"
            "\n"
            "Cartesian coord system.\n"
            "Conceptual overview:\n"
            "  A solver integrates the system's magnetic moments over time by computing and evolving\n"
            "    the spin directions according to a chosen time-integration scheme (Euler,\n"
            "    Heun, etc.). This interface decouples the integration policy from the\n"
            "    geometry and interaction definitions.\n"
            "\n"
            "The solver modifies the geometry in-place and returns a full snapshot of the\n"
            "  effective fields and energy contributions used in the step. This enables\n"
            "  inspection, logging, or post-processing without recomputing.\n"
            "\n"
            "@note Concrete implementations must implement the method `updateMoments()`.\n"
            "      They are responsible for allocating intermediate data and managing\n"
            "      consistency between field updates and geometry states.";

        py::class_<LLGSolver, AbstractSolver, std::shared_ptr<LLGSolver>>(cartesian, "LLGSolver")
            .def(py::init<SolverStrategy>(), py::arg("strategy") = SolverStrategy::EULER)
            .doc() =
            "Solver for time evolution of magnetic systems via the Landau-Lifshitz-Gilbert (LLG) equation.\n"
            "\n"
            "The LLGSolver class implements the core of explicit time-integration for atomistic spin "
            "dynamics\n"
            "  using the Landau-Lifshitz-Gilbert formalism. It combines a numerical integrator (Euler, "
            "Heun)\n"
            "  with a field update strategy to simulate the evolution of spin systems in time.\n"
            "\n"
            " Formula: dS/dt = -γ / (1 + α²) [ S × H_eff + α S × (S × H_eff) ]\n"
            "\n"
            "This solver mutates the geometry in place, updating the directions of magnetic moments using\n"
            "  effective fields derived from registered interactions. It is responsible for both computing\n"
            "  the spin derivatives (`dS/dt`) and applying them to the geometry using a selected "
            "time-stepping method.\n"
            "\n"
            "Internally, the solver delegates effective field computation to an Field Updater and\n"
            "  integrates the LLG equation using either Strategy methods.\n"
            "\n"
            "@note All physical units are assumed to be in SI. Spin vectors are expected to be normalized.";
    }
}

#endif // ! __SOLVERS_HPP__
