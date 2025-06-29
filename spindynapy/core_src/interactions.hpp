#ifndef __INTERACTIONS_HPP__
#define __INTERACTIONS_HPP__

/**
 * @file   interactions.hpp
 * @brief  Interaction interfaces and implementations for spin dynamics simulations.
 *
 * This header defines the abstract interface Interaction and
 *   concrete realisations of magnetic interactions for the Cartesian coordinate system.
 *
 * Functional overview:
 *
 * - An **Interaction** represents any physical mechanism that contributes to the
 *   effective magnetic field on each spin (and optionally, energy terms).
 *
 * - Each interaction defines:
 *   - `prepare(...)`                – optional pre-pass for neighbour caching or validation,
 *   - `calculateFieldContribution(...)` – computes 𝐁ᵢ interaction on a given spin,
 *   - `calculateEnergy(...)`       – returns energy contribution for a given spin and field,
 *   - `getName()` / `__str__()`    – for introspection, serialization and diagnostics.
 *
 * - Typical interactions include:
 *   - exchange (Heisenberg-type),
 *   - external field,
 *   - anisotropy,
 *   - dipole-dipole,
 *   - demagnetization field,
 *   - thermal fluctuation (stochastic noise).
 *
 * Interface:
 * - `IInteraction<CoordSystem>`       – abstract interface for interaction definitions,
 * - `InteractionRegistry<CoordSystem>`– registry for managing interaction instances.
 *
 * Concrete:
 * - `ExchangeInteraction`            – isotropic Heisenberg exchange (cutoff-based),
 * - `ExternalInteraction`           – constant applied magnetic field,
 * - `AnisotropyInteraction`         – uniaxial anisotropy (K·(S·n)² type),
 * - `DipoleDipoleInteraction`       – pairwise dipolar field (cutoff or macrocells),
 * - `DemagnetizationInteraction`    – shape-field approximation with self-term subtraction,
 * - `ThermalInteraction`            – stochastic thermal noise per LLG-FDT model.
 *
 * Python bindings:
 * - Macros: `INTERACTION_TEMPLATE_BINDINGS`
 * - Submodule binder: `pyBindInteractions(py::module_ &)`
 *
 * @note All units are SI unless explicitly stated. Effective fields are returned in [T],
 *       energies in [J]. Positions and directions are dimensionless in interface.
 *
 * @note All interactions are read-only after construction and can be safely shared
 *       across threads via the `InteractionRegistry`.
 *
 * @copyright 2025 SpinDynaPy
 */

#include "constants.hpp"
#include "geometries.hpp"
#include "registries.hpp"
#include "types.hpp"

#include <cmath>
#include <format>
#include <memory>
#include <pybind11/pybind11.h>
#include <random>
#include <stdexcept>
#include <string>

namespace py = pybind11;

namespace PYTHON_API spindynapy {

/**
 * @class  IInteraction
 * @brief  Abstract interface for physical interactions in a spin system.
 *
 * @tparam CoordSystem  Coordinate system type satisfying `CoordSystemConcept`.
 *
 * This interface defines the contract for implementing physical interactions
 *   between magnetic moments in a given geometry. It abstracts away the details
 *   of how effective fields and energy contributions are computed based on the
 *   spatial and material context of the system.
 *
 * Core responsibilities:
 *   - optionally precompute geometry-dependent data (e.g. neighbour lists),
 *   - compute the effective field contributed by the interaction at a given moment,
 *   - compute the interaction energy for a moment in a given field.
 *
 * @note The `prepare()` step is optional and may be skipped depending on simulation context.
 *       Implementations should not assume it will always be called.
 *
 * Derived classes must be registered via `InteractionRegistry`.
 *
 * @see InteractionRegistry
 */
template <CoordSystemConcept CoordSystem> class PYTHON_API IInteraction {
  protected:
    /**
     * @brief Protected default constructor (interface only).
     * @details Prevents direct instantiation.
     */
    IInteraction() noexcept = default;

  public:
    /**
     * @brief Virtual destructor.
     */
    virtual ~IInteraction() noexcept = default;

    /**
     * @brief Optional interaction precomputation step.
     *
     * This hook is invoked once before the simulation begins, allowing
     *   the interaction to cache neighbour data, allocate buffers, or
     *   validate geometry/material parameters.
     *
     * @param geometry          Geometry of the system.
     * @param material_registry Reference to material registry.
     *
     * @returns void - mutates internal or external (e.g geometry) state.
     */
    PYTHON_API virtual void prepare(
        [[maybe_unused]] IGeometry<CoordSystem> &geometry,
        [[maybe_unused]] MaterialRegistry &material_registry
    ) = 0;

    /**
     * @brief Compute the effective field induced by this interaction.
     *
     * For the specified moment, calculates the vector contribution to
     *   the total effective field due to this interaction only.
     *
     * @param moment_index       Index of the moment within the geometry.
     * @param geometry           Geometry of the system.
     * @param material_registry  Reference to material registry.
     *
     * @returns Effective field vector (SI units, typically [T]).
     */
    PYTHON_API virtual EffectiveField calculateFieldContribution(
        size_t moment_index, IGeometry<CoordSystem> &geometry, MaterialRegistry &material_registry
    ) const = 0;

    /**
     * @brief Compute interaction energy for a given moment and effective field.
     *
     * This evaluates the contribution of this interaction to the total system
     *   energy at the level of a single moment.
     *
     * @param moment Moment object to evaluate.
     * @param field  Effective field applied to the moment.
     *
     * @returns Energy contribution (joules).
     */
    PYTHON_API virtual double calculateEnergy(CoordSystem::Moment &moment, EffectiveField &field) const = 0;

    /**
     * @brief Get the canonical name of the interaction.
     *
     * Useful for serialization, registry display, and diagnostics.
     *
     * @returns Short interaction name (e.g., "EXCHANGE").
     */
    PYTHON_API virtual std::string getName() const = 0;

    /**
     * @brief Return a string representation of the interaction parameters.
     * @returns Human-readable parameter summary.
     */
    PYTHON_API virtual std::string __str__() const = 0;
};

/**
 * @typedef InteractionRegistry
 * @brief  Registry of physical interactions for a given coordinate system.
 *
 * The `InteractionRegistry` stores instances of `Interaction` and provides
 *   global access to the active set of interactions in the system.
 *
 * Registries are read-only after construction and therefore safe for concurrent access.
 *
 * Intended usage:
 *   - Maintain a globally available set of interactions used in the simulation.
 *   - Allow simulation loops to iterate over all interactions when calculating
 *     effective fields and energy contributions.
 *
 * @note Interaction instances are stored as shared pointers, allowing polymorphic usage
 *       and dynamic dispatch at runtime.
 *
 * @see IInteraction
 * @see Registry
 */
template <CoordSystemConcept CoordSystem>
using InteractionRegistry = PYTHON_API Registry<IInteraction<CoordSystem>>;

// ============================================================================
//  Cartesian interactions
// ============================================================================

namespace PYTHON_API cartesian {

/**
 * @typedef  AbstractInteraction
 * @brief  Abstract interface for physical interactions in a spin system.
 *
 * This interface defines the contract for implementing physical interactions
 *   between magnetic moments in a given geometry. It abstracts away the details
 *   of how effective fields and energy contributions are computed based on the
 *   spatial and material context of the system.
 *
 * Core responsibilities:
 *   - optionally precompute geometry-dependent data (e.g. neighbour lists),
 *   - compute the effective field contributed by the interaction at a given moment,
 *   - compute the interaction energy for a moment in a given field.
 *
 * @note The `prepare()` step is optional and may be skipped depending on simulation context.
 *       Implementations should not assume it will always be called.
 *
 * Derived classes must be registered via `InteractionRegistry`.
 *
 * @see InteractionRegistry
 */
using AbstractInteraction = PYTHON_API IInteraction<NamespaceCoordSystem>;

/**
 * @class  ExchangeInteraction
 * @brief  Heisenberg exchange interaction between magnetic moments.
 *
 * This class implements a short-range isotropic exchange interaction.
 *   The exchange field on a moment is computed as the weighted sum of its neighbors'
 *   directions within a specified cutoff radius.
 *
 * This follows the classical Heisenberg model with pairwise interactions:
 *
 *     B_exch(i) = ∑_j ( J_ij · S_j )
 *
 * where J_ij is selected based on whether moments i and j belong to the same material.
 *
 * The energy contribution per moment is computed as:
 *
 *     E = −½ · μ_s · S · B_exch
 *
 * The symmetry factor ½ ensures that each moment pair contributes only once.
 */
class PYTHON_API ExchangeInteraction : public AbstractInteraction {
  protected:
    /**
     * @brief Cutoff radius (in meters) for neighbor interaction.
     *
     * All neighbors within this radius are considered during exchange calculation.
     */
    double _cutoff_radius;

  public:
    /**
     * @brief Construct an exchange interaction with given cutoff.
     *
     * @param cutoff_radius Cutoff radius in meters (SI units) for only neighbor interactions.
     */
    PYTHON_API ExchangeInteraction(double cutoff_radius) : _cutoff_radius(cutoff_radius) {};

    /**
     * @brief Optional interaction precomputation step.
     *
     * This hook is invoked once before the simulation begins, allowing
     *   the interaction to cache neighbour data, allocate buffers, or
     *   validate geometry/material parameters.
     *
     * @param geometry          Geometry of the system.
     * @param material_registry Reference to material registry.
     *
     * @returns void - mutates internal or external (e.g geometry) state.
     */
    PYTHON_API virtual void prepare(
        [[maybe_unused]] IGeometry<NamespaceCoordSystem> &geometry,
        [[maybe_unused]] MaterialRegistry &material_registry
    ) override {
        size_t moments_size = geometry.size();
        for (size_t i = 0; i < moments_size; ++i) {
            // TODO: сделать это на стороне prepare() геометрии
            geometry.getNeighbors(i, this->_cutoff_radius);
        }
        return;
    };

    /**
     * @brief Compute interaction energy for a given moment and effective field.
     *
     * This evaluates the contribution of this interaction to the total system
     *   energy at the level of a single moment.
     *
     * @param moment Moment object to evaluate.
     * @param field  Effective field applied to the moment.
     *
     * @returns Energy contribution (joules).
     */
    PYTHON_API double calculateEnergy(NamespaceCoordSystem::Moment &moment, EffectiveField &field)
        const override {
        const auto &current_material = moment.getMaterial();
        const auto &direction = moment.getDirection().asVector();
        return -0.5 * current_material.atomic_magnetic_saturation_absolute *
               direction.dot(field); // TODO сделать через флаг pairwise на вызывающей стороне
    }

    /**
     * @brief Compute the effective field induced by this interaction.
     *
     * For the specified moment, calculates the vector contribution to
     *   the total effective field due to this interaction only.
     *
     * @param moment_index       Index of the moment within the geometry.
     * @param geometry           Geometry of the system.
     * @param material_registry  Reference to material registry.
     *
     * @returns Effective field vector (SI units, typically [T]).
     */
    PYTHON_API virtual EffectiveField
    calculateFieldContribution(size_t moment_index, IGeometry<NamespaceCoordSystem> &geometry, MaterialRegistry &)
        const override {
        EffectiveField exchange_field = EffectiveField::Zero();
        auto &current_material = geometry[moment_index].getMaterial();

        auto neighbor_indices = geometry.getNeighbors(moment_index, this->_cutoff_radius);
        for (size_t neighbor_index : neighbor_indices) {
            exchange_field +=
                (geometry[neighbor_index].getMaterial() == current_material
                     ? current_material.exchange_monomaterial_prefix
                     : current_material
                           .exchange_monomaterial_prefix /* тут можно будет интерфейс использовать */
                ) *
                geometry[neighbor_index].getDirection().asVector();
        }
        return exchange_field;
    }

    /**
     * @brief Get the canonical name of the interaction.
     *
     * Useful for serialization, registry display, and diagnostics.
     *
     * @returns Short interaction name (e.g., "EXCHANGE").
     */
    PYTHON_API virtual std::string getName() const override { return "EXCHANGE"; }

    /**
     * @brief Return a string representation of the interaction parameters.
     * @returns Human-readable parameter summary.
     */
    PYTHON_API virtual std::string __str__() const override {
        return std::format("ExchangeInteraction(r={})", _cutoff_radius);
    };
};

/**
 * @class  ExternalInteraction
 * @brief  Uniform magnetic field interaction (external Zeeman term).
 *
 * This interaction applies a constant external magnetic field to all moments
 *   in the system. The same field vector is used for every moment, regardless
 *   of material or position.
 *
 * The contribution to the effective field is:
 *
 *     B_ext(i) = B_external
 *
 * where B_ext is the externally applied field (in T).
 *
 * The energy of a moment in the external field is:
 *
 *     E = −μ_s · S · B_ext
 *
 * @note This interaction is stateless and does not depend on geometry.
 *
 * @see IInteraction
 */
class PYTHON_API ExternalInteraction : public AbstractInteraction {
  protected:
    /**
     * @brief Constant external magnetic field vector (in tesla).
     *
     * This field is applied uniformly to all moments in the system.
     */
    Eigen::Vector3d external_field;

  public:
    /**
     * @brief Construct the interaction with a 3D external field vector.
     *
     * @param external_field  Constant field vector to apply (in tesla).
     */
    PYTHON_API ExternalInteraction(const Eigen::Vector3d &external_field) : external_field(external_field) {};

    /**
     * @brief Construct the interaction with individual field components.
     *
     * @param sx  Field component along x-axis (in T).
     * @param sy  Field component along y-axis (in T).
     * @param sz  Field component along z-axis (in T).
     */
    PYTHON_API ExternalInteraction(double sx, double sy, double sz) : external_field(sx, sy, sz) {};

    /**
     * @brief Optional interaction precomputation step.
     *
     * This hook is invoked once before the simulation begins, allowing
     *   the interaction to cache neighbour data, allocate buffers, or
     *   validate geometry/material parameters.
     *
     * @param geometry          Geometry of the system.
     * @param material_registry Reference to material registry.
     *
     * @returns void - mutates internal or external (e.g geometry) state.
     */
    PYTHON_API virtual void prepare(
        [[maybe_unused]] IGeometry<NamespaceCoordSystem> &geometry,
        [[maybe_unused]] MaterialRegistry &material_registry
    ) override {
        return;
    };

    /**
     * @brief Compute interaction energy for a given moment and effective field.
     *
     * This evaluates the contribution of this interaction to the total system
     *   energy at the level of a single moment.
     *
     * @param moment Moment object to evaluate.
     * @param field  Effective field applied to the moment.
     *
     * @returns Energy contribution (joules).
     */
    PYTHON_API double calculateEnergy(NamespaceCoordSystem::Moment &moment, EffectiveField &field)
        const override {
        const auto &current_material = moment.getMaterial();
        const auto &direction = moment.getDirection().asVector();
        const auto current_moment_norm =
            current_material.atomic_magnetic_saturation_magnetization * constants::BOHR_MAGNETON;
        return -current_moment_norm * direction.dot(field);
    }

    /**
     * @brief Compute the effective field induced by this interaction.
     *
     * For the specified moment, calculates the vector contribution to
     *   the total effective field due to this interaction only.
     *
     * @param moment_index       Index of the moment within the geometry.
     * @param geometry           Geometry of the system.
     * @param material_registry  Reference to material registry.
     *
     * @returns Effective field vector (SI units, typically [T]).
     */
    PYTHON_API virtual EffectiveField
    calculateFieldContribution(size_t, IGeometry<NamespaceCoordSystem> &, MaterialRegistry &) const override {
        return this->external_field;
    }

    /**
     * @brief Get the canonical name of the interaction.
     *
     * Useful for serialization, registry display, and diagnostics.
     *
     * @returns Short interaction name (e.g., "EXCHANGE").
     */
    PYTHON_API virtual std::string getName() const override { return "EXTERNAL"; }

    /**
     * @brief Return a string representation of the interaction parameters.
     * @returns Human-readable parameter summary.
     */
    PYTHON_API virtual std::string __str__() const override {
        return std::format(
            "ExternalInteraction(field=({}, {}, {}))",
            this->external_field.x(),
            this->external_field.y(),
            this->external_field.z()
        );
    };
};

/**
 * @class  AnisotropyInteraction
 * @brief  Magnetocrystalline anisotropy interaction.
 *
 * This interaction models the internal energy and effective field resulting
 *   from material-specific anisotropy. It supports uniaxial anisotropy and is
 *   extensible to others (e.g., cubic).
 *
 * The effective field for uniaxial anisotropy is computed as:
 *
 *     B_aniso = (2 · k_u / μ_s) · (S · n) · n
 *
 * where:
 *   - k_u is the uniaxial anisotropy constant [J/atom],
 *   - n is the anisotropy axis (unit vector).
 *
 * The energy associated with this field is (uniaxial case):
 *
 *     E = −½ · μ_s · S · B_aniso
 *
 * The ½ factor arises due to differentiation of the quadratic energy form.
 *
 * @note This interaction expects material objects to hold a valid `Anisotropy` instance.
 *       If no anisotropy is defined, this interaction contributes nothing.
 *
 * @throws std::invalid_argument if anisotropy type is unsupported.
 *
 * @see UniaxialAnisotropy
 */
class PYTHON_API AnisotropyInteraction : public AbstractInteraction {
  public:
    /**
     * @brief Construct an anisotropy interaction (default, stateless).
     */
    PYTHON_API AnisotropyInteraction() {};

    /**
     * @brief Optional interaction precomputation step.
     *
     * This hook is invoked once before the simulation begins, allowing
     *   the interaction to cache neighbour data, allocate buffers, or
     *   validate geometry/material parameters.
     *
     * @param geometry          Geometry of the system.
     * @param material_registry Reference to material registry.
     *
     * @returns void - mutates internal or external (e.g geometry) state.
     */
    PYTHON_API virtual void prepare(
        [[maybe_unused]] IGeometry<NamespaceCoordSystem> &geometry,
        [[maybe_unused]] MaterialRegistry &material_registry
    ) override {
        return;
    };

    /**
     * @brief Compute interaction energy for a given moment and effective field.
     *
     * This evaluates the contribution of this interaction to the total system
     *   energy at the level of a single moment.
     *
     * @param moment Moment object to evaluate.
     * @param field  Effective field applied to the moment.
     *
     * @returns Energy contribution (joules).
     */
    PYTHON_API double calculateEnergy(NamespaceCoordSystem::Moment &moment, EffectiveField &field)
        const override {
        const auto &current_material = moment.getMaterial();
        const auto &direction = moment.getDirection().asVector();
        const auto current_moment_norm =
            current_material.atomic_magnetic_saturation_magnetization * constants::BOHR_MAGNETON;

        // из-за взятия производной от квадрата => / 2
        return -current_moment_norm / 2 * direction.dot(field);
    }

    /**
     * @brief Compute the effective field induced by this interaction.
     *
     * For the specified moment, calculates the vector contribution to
     *   the total effective field due to this interaction only.
     *
     * @param moment_index       Index of the moment within the geometry.
     * @param geometry           Geometry of the system.
     * @param material_registry  Reference to material registry.
     *
     * @returns Effective field vector (SI units, typically [T]).
     */
    PYTHON_API virtual EffectiveField
    calculateFieldContribution(size_t moment_index, IGeometry<NamespaceCoordSystem> &geometry, MaterialRegistry &)
        const override {
        //
        auto &moment = geometry[moment_index];
        auto &material = moment.getMaterial();
        auto anisotropy = material.anisotropy;
        auto atomic_magnetic_moments_norm =
            material.atomic_magnetic_saturation_magnetization * constants::BOHR_MAGNETON;

        if (!anisotropy) {
            return EffectiveField::Zero(); // Нет анизотропии
        }

        if (auto uniaxial = dynamic_cast<UniaxialAnisotropy *>(anisotropy.get())) {
            return 2 * uniaxial->constant / atomic_magnetic_moments_norm *
                   moment.getDirection().asVector().dot(uniaxial->axis) * uniaxial->axis;
        }

        throw std::invalid_argument("Unsupported anisotropy type");
    }

    /**
     * @brief Get the canonical name of the interaction.
     *
     * Useful for serialization, registry display, and diagnostics.
     *
     * @returns Short interaction name (e.g., "EXCHANGE").
     */
    PYTHON_API virtual std::string getName() const override { return "ANISOTROPY"; }

    /**
     * @brief Return a string representation of the interaction parameters.
     * @returns Human-readable parameter summary.
     */
    PYTHON_API virtual std::string __str__() const override { return "AnisotropyInteraction()"; };
};

using DistanceNormCache = Cache<std::pair<size_t, size_t>, double>;
using DistancePow3Cache = Cache<std::pair<size_t, size_t>, double>;
using DistancePow5Cache = Cache<std::pair<size_t, size_t>, double>;

/**
 * @class  DipoleDipoleInteraction
 * @brief  Long-range magnetostatic (dipole–dipole) interaction.
 *
 * This interaction accounts for the dipolar magnetic coupling between all
 *   magnetic moments. It models the classical field of a magnetic dipole
 *   acting on another dipole, and is responsible for magnetostatic effects
 *   such as demagnetizing fields and long-range anisotropies.
 *
 * The effective field at position **rᵢ** due to a moment **mⱼ** located at **rⱼ**
 *   is computed using the formula:
 *
 *     B_dd(rᵢ) = (μ₀ / 4π) ∑ⱼ [ (3(mⱼ·rᵢⱼ)rᵢⱼ / |rᵢⱼ|⁵) − (mⱼ / |rᵢⱼ|³) ]
 *
 * where:
 *   - rᵢⱼ = rⱼ − rᵢ is the displacement vector from i to j,
 *   - μ₀ is the magnetic permeability of vacuum.
 *
 * Energy of interaction is given by:
 *
 *     E = −½ · μ_s · S · B_dd
 *
 * Strategy options:
 * - `"cutoff"`: sum over neighbors within a spherical radius (explicit pairwise summation);
 * - `"macrocells"`: group distant moments into representative macro-averaged cells
 *                   to accelerate long-range computations (using a macrocell grid + amplitude).
 *
 * @note This implementation does not use fast multipole methods or FFT;
 *       complexity is depending on cutoff and macrocell density.
 *
 * @see DemagnetizationInteraction
 */
class PYTHON_API DipoleDipoleInteraction : public AbstractInteraction {
  protected:
    /**
     * @brief Evaluation strategy used for neighbor selection.
     *
     * Controls how the dipole-dipole interaction is computed:
     * - `"cutoff"` — direct summation over neighbors within radius;
     * - `"macrocells"` — average over macrocell representatives.
     */
    std::string _strategy;

    /**
     * @brief Cutoff radius in meters.
     *
     * Used for both `"cutoff"` and `"macrocells"` strategies to determine interaction scope.
     */
    double _cutoff_radius;

    /**
     * @brief Compute total dipolar field from a set of source moments.
     *
     * This helper method applies the classical dipole field formula to compute
     *   the net contribution of all `calculation_moments` on the given `current_moment`.
     *
     * @param current_moment       Moment receiving the field (target).
     * @param current_material     Material of the target moment (used for μ₀).
     * @param calculation_moments  Moments contributing to the dipolar field.
     *
     * @returns Total dipolar field vector acting on the target moment.
     */
    EffectiveField calculate(
        Moment &current_moment, Material &, MomentsContainer<NamespaceCoordSystem> calculation_moments
    ) const {
        EffectiveField dipole_field = EffectiveField::Zero();

        for (auto moment : calculation_moments) {
            auto distance_vector =
                moment->getCoordinates().asVector() - current_moment.getCoordinates().asVector();

            auto neighbor_atomistic_moment = moment->amplitude *
                                             moment->getMaterial().atomic_magnetic_saturation_magnetization *
                                             constants::BOHR_MAGNETON;
            auto distance_norm = distance_vector.norm();

            if (distance_norm < 1e-30) {
                // предотвратить деление на ноль (превращение в NaN)
                continue;
            }

            dipole_field += neighbor_atomistic_moment *
                            ((3 * moment->getDirection().asVector().dot(distance_vector) * distance_vector) /
                                 pow(distance_norm, 5) -
                             moment->getDirection().asVector() / pow(distance_norm, 3));
        }

        dipole_field = (constants::VACUUM_MAGNETIC_PERMEABILITY / 4 / constants::NUMBER_PI) * dipole_field;

        return dipole_field;
    }

  public:
    /**
     * @brief Construct dipole-dipole interaction with a given cutoff and strategy.
     *
     * @param cutoff_radius Cutoff distance in meters.
     * @param strategy      One of: `"cutoff"` or `"macrocells"`.
     *
     * @throws std::invalid_argument if strategy is not recognized.
     */
    PYTHON_API DipoleDipoleInteraction(double cutoff_radius, std::string strategy = "cutoff")
        : _strategy(strategy), _cutoff_radius(cutoff_radius) {
        if (strategy != "cutoff" && strategy != "macrocells") {
            throw std::invalid_argument("Invalid strategy string");
        }
    };

    /**
     * @brief Optional interaction precomputation step.
     *
     * This hook is invoked once before the simulation begins, allowing
     *   the interaction to cache neighbour data, allocate buffers, or
     *   validate geometry/material parameters.
     *
     * @param geometry          Geometry of the system.
     * @param material_registry Reference to material registry.
     *
     * @returns void - mutates internal or external (e.g geometry) state.
     */
    PYTHON_API virtual void prepare(
        [[maybe_unused]] IGeometry<NamespaceCoordSystem> &geometry,
        [[maybe_unused]] MaterialRegistry &material_registry
    ) override {
        size_t moments_size = geometry.size();
        for (size_t i = 0; i < moments_size; ++i) {
            // TODO: сделать это на стороне prepare() геометрии
            geometry.getNeighbors(i, this->_cutoff_radius);
        }
        return;
    };

    /**
     * @brief Compute interaction energy for a given moment and effective field.
     *
     * This evaluates the contribution of this interaction to the total system
     *   energy at the level of a single moment.
     *
     * @param moment Moment object to evaluate.
     * @param field  Effective field applied to the moment.
     *
     * @returns Energy contribution (joules).
     */
    PYTHON_API double calculateEnergy(NamespaceCoordSystem::Moment &moment, EffectiveField &field)
        const override {
        const auto &current_material = moment.getMaterial();
        const auto &direction = moment.getDirection().asVector();
        const auto current_moment_norm =
            current_material.atomic_magnetic_saturation_magnetization * constants::BOHR_MAGNETON;
        return -0.5 * current_moment_norm * moment.amplitude *
               direction.dot(field); // TODO сделать через флаг pairwise на вызывающей стороне
    }

    /**
     * @brief Compute the effective field induced by this interaction.
     *
     * For the specified moment, calculates the vector contribution to
     *   the total effective field due to this interaction only.
     *
     * @param moment_index       Index of the moment within the geometry.
     * @param geometry           Geometry of the system.
     * @param material_registry  Reference to material registry.
     *
     * @returns Effective field vector (SI units, typically [T]).
     */
    PYTHON_API virtual EffectiveField
    calculateFieldContribution(size_t moment_index, IGeometry<NamespaceCoordSystem> &geometry, MaterialRegistry &)
        const override {
        auto &current_moment = geometry[moment_index];
        auto &current_material = current_moment.getMaterial();

        if (this->_strategy == "cutoff") {
            auto neighbor_indices = geometry.getNeighbors(moment_index, this->_cutoff_radius);
            auto calculation_moments = geometry.getFromIndexes(neighbor_indices);
            return this->calculate(current_moment, current_material, calculation_moments);
        } else if (this->_strategy == "macrocells") {
            // TODO как учесть внутри? получать список внутри ячейки находящихся
            auto calculation_moments = geometry.getMomentsFromMacrocells(moment_index, this->_cutoff_radius);
            return this->calculate(current_moment, current_material, calculation_moments);
        }
        // Если не удалось найти подходящую стратегию
        throw std::invalid_argument("Invalid strategy string");
    }

    /**
     * @brief Get the canonical name of the interaction.
     *
     * Useful for serialization, registry display, and diagnostics.
     *
     * @returns Short interaction name (e.g., "EXCHANGE").
     */
    PYTHON_API virtual std::string getName() const override { return "DIPOLE-DIPOLE"; }

    /**
     * @brief Return a string representation of the interaction parameters.
     * @returns Human-readable parameter summary.
     */
    PYTHON_API virtual std::string __str__() const override {
        return std::format(
            "DipoleDipoleInteraction(cutoff_radius={}, strategy={})", this->_cutoff_radius, this->_strategy
        );
    };
};

/**
 * @class  DemagnetizationInteraction
 * @brief  Effective demagnetizing field interaction based on dipole-dipole approximation.
 *
 * This interaction models the **macroscopic demagnetizing field** as arising from
 *   long-range dipole-dipole effects that cannot be fully resolved within practical
 *   cutoffs or coarse-grained representations (not all cases).
 *
 * Unlike `DipoleDipoleInteraction`, this class explicitly introduces a **self-term correction**
 *   to account for **shape-based demagnetization**. Its role is to compensate for the dipoles
 *   that are *not* explicitly included in the local sum (either due to cutoff or coarse macrocelling).
 *
 * ❱❱ Total field:
 *
 *     𝐁_total = 𝐁_ext + 𝐁_demag
 *
 * ❱❱ Dipole field:
 *
 *     𝐁_dip,i = (μ₀ / 4π) ∑ⱼ≠ᵢ [ 3(mⱼ·rᵢⱼ)rᵢⱼ / |rᵢⱼ|⁵ − mⱼ / |rᵢⱼ|³ ]
 *
 * ❱❱ Demagnetizing correction (macrocell case, shape-factor):
 *
 *     𝐁_self = - (μ₀ / 3) · (𝐌 / V)
 *
 * where:
 *  - 𝐌 = average magnetic moment vector of macrocell [A·m²],
 *  - V = effective macrocell volume ≈ N · a³ (a = atomic cell size),
 *  - 𝐁_demag = − 𝐁_self + 𝐁_dip
 *
 * ❱❱ In shape-based demagnetization, we define:
 *
 *     ℋ_demag = {
 *         −ℋ_dip           // fully computed case (honest sum over j ≠ i)
 *         −ℋ_self + ℋ_dip  // subtract only shape term from partial sum
 *     }
 *
 * ❱❱ The **self-term correction** is:
 *
 *     𝐁_self = μ₀ / 3 · (𝐌 / V)     (macrocell form)
 *     𝐁_self = μ₀ / 3 · (μ_s · 𝐒 / V₀) (cutoff, single-atom effective volume)
 *
 * where V₀ ≈ atomic volume for single spin.
 *
 * This shape self-field implements the **−𝑁 · 𝐌** contribution where 𝑁 is the demagnetizing tensor
 *   for a uniform ellipsoid. Here, we approximate it using the `1/3` isotropic term.
 *
 * @note This class explicitly **subtracts short-range dipole fields** already computed as
 *       `DipoleDipoleInteraction`, and adds the missing **demagnetizing self-term**
 *       as a correction due to NOT CALCULATED OTHERS DIPOLES. <!!!>
 *
 * @note In practice, this interaction improves convergence of magnetostatic fields
 *       in heterogeneous geometries or when using coarse macrocells.
 *
 * @see DipoleDipoleInteraction
 */
class DemagnetizationInteraction : public DipoleDipoleInteraction {
  public:
    /**
     * @brief Construct a demagnetizing-field interaction with strategy and cutoff radius.
     *
     * This behaves like a dipole-dipole interaction but subtracts the average shape self-term
     *   based on the selected strategy.
     *
     * @param cutoff_radius  Radius used for neighbor/macrocell selection (in meters).
     * @param strategy       Either `"cutoff"` (local spins) or `"macrocells"` (coarse averaging).
     *
     * @throws std::invalid_argument if strategy is not recognized.
     */
    PYTHON_API DemagnetizationInteraction(double cutoff_radius, std::string strategy = "macrocells")
        : DipoleDipoleInteraction(cutoff_radius, strategy) {};

    /**
     * @brief Compute the effective field induced by this interaction.
     *
     * For the specified moment, calculates the vector contribution to
     *   the total effective field due to this interaction only.
     *
     * @param moment_index       Index of the moment within the geometry.
     * @param geometry           Geometry of the system.
     * @param material_registry  Reference to material registry.
     *
     * @returns Effective field vector (SI units, typically [T]).
     */
    PYTHON_API virtual EffectiveField
    calculateFieldContribution(size_t moment_index, IGeometry<NamespaceCoordSystem> &geometry, MaterialRegistry &)
        const override {
        auto &current_moment = geometry[moment_index];
        auto &current_material = current_moment.getMaterial();

        if (this->_strategy == "cutoff") {
            auto neighbor_indices = geometry.getNeighbors(moment_index, this->_cutoff_radius);
            auto calculation_moments = geometry.getFromIndexes(neighbor_indices);
            auto self_term = (constants::VACUUM_MAGNETIC_PERMEABILITY / 3.0) *
                             (current_material.atomic_magnetic_saturation_magnetization *
                              constants::BOHR_MAGNETON * current_moment.getDirection().asVector()) /
                             current_material.atom_cell_size;
            return this->calculate(current_moment, current_material, calculation_moments) - self_term;
        } else if (this->_strategy == "macrocells") {

            auto calculation_moments = geometry.getMomentsFromMacrocells(moment_index, this->_cutoff_radius);

            auto macrocell = geometry.getMacrocells()[geometry.getMacrocellIndexBySpin(moment_index)];
            auto macrocell_moment = macrocell.avg_moment;

            auto moment_term = macrocell_moment->amplitude * macrocell_moment->getDirection().asVector() *
                               macrocell_moment->getMaterial().atomic_magnetic_saturation_magnetization *
                               constants::BOHR_MAGNETON;

            auto self_term =
                (constants::VACUUM_MAGNETIC_PERMEABILITY / 3.0) *
                (moment_term / (macrocell.moment_indices.size() * current_material.atom_cell_size));

            return this->calculate(current_moment, current_material, calculation_moments) - self_term;
        }
        // Если не удалось найти подходящую стратегию
        throw std::invalid_argument("Invalid strategy string");
    }

    /**
     * @brief Get the canonical name of the interaction.
     *
     * Useful for serialization, registry display, and diagnostics.
     *
     * @returns Short interaction name (e.g., "EXCHANGE").
     */
    PYTHON_API virtual std::string getName() const override { return "DEMAGNETIZATION"; }

    /**
     * @brief Return a string representation of the interaction parameters.
     * @returns Human-readable parameter summary.
     */
    PYTHON_API virtual std::string __str__() const override {
        return std::format(
            "DemagnetizationInteraction(cutoff_radius={}, strategy={})", this->_cutoff_radius, this->_strategy
        );
    };
};

/**
 * @class  ThermalInteraction
 * @brief  Stochastic thermal field contribution for finite-temperature LLG dynamics.
 *
 * This interaction injects random magnetic field fluctuations into the effective field
 *   according to the **fluctuation–dissipation theorem** at a given temperature.
 *
 * This is required to simulate **finite-temperature effects** such as thermal agitation,
 *   relaxation, and switching due to noise.
 *
 * ❱❱ Model:
 *
 * The thermal field 𝐁ₜₕ is defined by:
 *
 *     𝐁ₜₕ = √[ 2α k_B T / (γ μ_s Δt) ] · ξ,     ξ ~ N(0,1)
 *
 * where:
 * - α — Gilbert damping constant (unitless),
 * - k_B — Boltzmann constant [J/K],
 * - T — temperature [K],
 * - γ — gyromagnetic ratio [rad·s⁻¹·T⁻¹],
 * - μ_s — magnetic moment of spin [A·m²],
 * - Δt — time step [s],
 * - ξ — vector of 3 normally distributed random numbers N(0,1).
 *
 * This expression ensures:
 * - ⟨𝐁ₜₕ(t)⟩ = 0 (mean zero),
 * - ⟨𝐁ₜₕ(t) · 𝐁ₜₕ(t')⟩ ∝ δ(t − t') (delta-correlated).
 *
 * @note This interaction does **not** contribute to total energy (it is purely statistical).
 *
 * @note Internally uses `std::mt19937` and `std::normal_distribution<double>` for sampling.
 *       RNG is seeded with a user-provided seed or system random device (if zero).
 */
class PYTHON_API ThermalInteraction : public AbstractInteraction {
  private:
    /**
     * @brief Temperature of the thermal bath [K].
     */
    double _temperature;

    /**
     * @brief Time step used in the numerical solver [s].
     */
    double _dt;

    /**
     * @brief Random number generator (Mersenne Twister engine).
     *
     * @details Used for generating Gaussian-distributed thermal noise.
     */
    mutable std::mt19937 _rng;

    /**
     * @brief Standard normal distribution N(0,1) for sampling noise.
     *
     * @details Used to generate thermal field components independently.
     *          Bound to the `_rng` generator.
     */
    mutable std::normal_distribution<double> _norm {0.0, 1.0};

  public:
    /**
     * @brief Construct a thermal stochastic interaction.
     *
     * @param temperature  Temperature in kelvin [K].
     * @param dt           Time step used in solver [s].
     * @param seed         Optional RNG seed. If 0, uses `std::random_device`.
     *
     * @returns Constructed interaction instance.
     */
    PYTHON_API ThermalInteraction(double temperature, double dt, uint32_t seed = 0)
        : _temperature(temperature), _dt(dt) {
        if (seed <= 0) {
            _rng.seed(std::random_device {}());
        } else {
            _rng.seed(seed);
        }
    }

    /**
     * @brief Optional interaction precomputation step.
     *
     * This hook is invoked once before the simulation begins, allowing
     *   the interaction to cache neighbour data, allocate buffers, or
     *   validate geometry/material parameters.
     *
     * @param geometry          Geometry of the system.
     * @param material_registry Reference to material registry.
     *
     * @returns void - mutates internal or external (e.g geometry) state.
     */
    PYTHON_API virtual void prepare(
        [[maybe_unused]] IGeometry<NamespaceCoordSystem> &geometry,
        [[maybe_unused]] MaterialRegistry &material_registry
    ) override {
        return;
    };

    /**
     * @brief Compute the effective field induced by this interaction.
     *
     * For the specified moment, calculates the vector contribution to
     *   the total effective field due to this interaction only.
     *
     * @param moment_index       Index of the moment within the geometry.
     * @param geometry           Geometry of the system.
     * @param material_registry  Reference to material registry.
     *
     * @returns Effective field vector (SI units, typically [T]).
     */
    PYTHON_API EffectiveField calculateFieldContribution(
        size_t moment_index, IGeometry<NamespaceCoordSystem> &geometry, MaterialRegistry & /*unused*/
    ) const override {
        auto &moment = geometry[moment_index];
        auto &material = moment.getMaterial();
        // учтём, что amplitude может быть >1
        const double mu_s =
            moment.amplitude * material.atomic_magnetic_saturation_magnetization * constants::BOHR_MAGNETON;

        const double pref = std::sqrt(
            2.0 * material.damping_constant * constants::BOLTZMANN_CONSTANT * _temperature /
            (material.gyromagnetic_ratio * mu_s * _dt)
        );

        /* ---------- потокобезопасный RNG (по копии seed-а) ---------- */
        thread_local static std::mt19937 local_rng(std::random_device {}()); // TODO переделать
        thread_local static std::normal_distribution<double> local_norm(0.0, 1.0);
        /* ------------------------------------------------------------ */
        auto x = local_norm(local_rng);
        auto y = local_norm(local_rng);
        auto z = local_norm(local_rng);

        return {pref * x, pref * y, pref * z}; // Возвращаем эффективное поле
    }

    /**
     * @brief Compute interaction energy for a given moment and effective field.
     *
     * This evaluates the contribution of this interaction to the total system
     *   energy at the level of a single moment.
     *
     * @param moment Moment object to evaluate.
     * @param field  Effective field applied to the moment.
     *
     * @returns Energy contribution (joules).
     */
    PYTHON_API double calculateEnergy(NamespaceCoordSystem::Moment &, EffectiveField &) const override {
        return 0.0; // Тепловое поле не имеет собственной энергии в макроскопическом выводе
    }

    /**
     * @brief Get the canonical name of the interaction.
     *
     * Useful for serialization, registry display, and diagnostics.
     *
     * @returns Short interaction name (e.g., "EXCHANGE").
     */
    PYTHON_API std::string getName() const override { return "THERMAL"; }

    /**
     * @brief Return a string representation of the interaction parameters.
     * @returns Human-readable parameter summary.
     */
    PYTHON_API std::string __str__() const override {
        return std::format("ThermalInteraction(T={} K, dt={} s)", _temperature, _dt);
    }
};

/**
 * @typedef AbstractInteractionRegistry
 * @brief  Registry of physical interactions for a given coordinate system.
 *
 * The `InteractionRegistry` stores instances of `Interaction` and provides
 *   global access to the active set of interactions in the system.
 *
 * Registries are read-only after construction and therefore safe for concurrent access.
 *
 * Intended usage:
 *   - Maintain a globally available set of interactions used in the simulation.
 *   - Allow simulation loops to iterate over all interactions when calculating
 *     effective fields and energy contributions.
 *
 * @note Interaction instances are stored as shared pointers, allowing polymorphic usage
 *       and dynamic dispatch at runtime.
 *
 * @see IInteraction
 * @see Registry
 */
using AbstractInteractionRegistry = PYTHON_API InteractionRegistry<NamespaceCoordSystem>;

}; // namespace PYTHON_API cartesian

}; // namespace PYTHON_API spindynapy

// ===========================================================================
//  Python bindings
// ===========================================================================

#define INTERACTION_TEMPLATE_BINDINGS(cls)                                                                   \
    .def(                                                                                                    \
        "__str__",                                                                                           \
        &cls::__str__,                                                                                       \
        py::doc("@brief Return a string representation of the interaction parameters.\n"                     \
                "@returns Human-readable parameter summary.")                                                \
    )                                                                                                        \
        .def(                                                                                                \
            "calculate_field_contribution",                                                                  \
            &cls::calculateFieldContribution,                                                                \
            py::arg("moment_index"),                                                                         \
            py::arg("geometry"),                                                                             \
            py::arg("material_registry"),                                                                    \
            py::doc("@brief Compute the effective field induced by this interaction.\n"                      \
                    "\n"                                                                                     \
                    "For the specified moment, calculates the vector contribution to\n"                      \
                    "  the total effective field due to this interaction only.\n"                            \
                    "\n"                                                                                     \
                    "@param moment_index       Index of the moment within the geometry.\n"                   \
                    "@param geometry           Geometry of the system.\n"                                    \
                    "@param material_registry  Reference to material registry.\n"                            \
                    "\n"                                                                                     \
                    "@returns Effective field vector (SI units, typically [T]).\n")                          \
        )                                                                                                    \
        .def(                                                                                                \
            "calculate_energy",                                                                              \
            &cls::calculateEnergy,                                                                           \
            py::arg("moment"),                                                                               \
            py::arg("field"),                                                                                \
            py::doc("@brief Compute interaction energy for a given moment and effective field.\n"            \
                    "\n"                                                                                     \
                    "This evaluates the contribution of this interaction to the total system\n"              \
                    "  energy at the level of a single moment.\n"                                            \
                    "\n"                                                                                     \
                    "@param moment Moment object to evaluate.\n"                                             \
                    "@param field  Effective field applied to the moment.\n"                                 \
                    "\n"                                                                                     \
                    "@returns Energy contribution (joules).")                                                \
        )                                                                                                    \
        .def(                                                                                                \
            "prepare",                                                                                       \
            &cls::prepare,                                                                                   \
            py::arg("geometry"),                                                                             \
            py::arg("material_registry"),                                                                    \
            py::doc("@brief Optional interaction precomputation step.\n"                                     \
                    "\n"                                                                                     \
                    "This hook is invoked once before the simulation begins, allowing\n"                     \
                    "  the interaction to cache neighbour data, allocate buffers, or\n"                      \
                    "  validate geometry/material parameters.\n"                                             \
                    "\n"                                                                                     \
                    "@param geometry          Geometry of the system.\n"                                     \
                    "@param material_registry Reference to material registry.\n"                             \
                    "\n"                                                                                     \
                    "@returns void - mutates internal or external (e.g geometry) state.")                    \
        )                                                                                                    \
        .def(                                                                                                \
            "get_name",                                                                                      \
            &cls::getName,                                                                                   \
            py::doc("@brief Get the canonical name of the interaction.\n"                                    \
                    "\n"                                                                                     \
                    "Useful for serialization, registry display, and diagnostics.\n"                         \
                    "\n"                                                                                     \
                    "@returns Short interaction name (e.g., \"EXCHANGE\").")                                 \
        )

/**
 * @brief Bind the interactions utilities to a Python sub‑module.
 *
 * @param module Parent PyBind11 module (usually the core extension module).
 * @returns void – extends the parent module.
 */
inline void pyBindInteractions(py::module_ &module) {
    using namespace spindynapy;

    // -------- | INTERACTIONS | --------
    py::module_ interaction_module = module.def_submodule("interactions");

    interaction_module.doc() =
        "@brief  Interaction interfaces and implementations for spin dynamics simulations.\n"
        "\n"
        "This header defines the abstract interface Interaction and\n"
        "  concrete realisations of magnetic interactions for the Cartesian coordinate system.\n"
        "\n"
        "Functional overview:\n"
        "\n"
        "- An **Interaction** represents any physical mechanism that contributes to the\n"
        "  effective magnetic field on each spin (and optionally, energy terms).\n"
        "\n"
        "- Each interaction defines:\n"
        "  - `prepare(...)`                – optional pre-pass for neighbour caching or validation,\n"
        "  - `calculateFieldContribution(...)` – computes 𝐁ᵢ interaction on a given spin,\n"
        "  - `calculateEnergy(...)`       – returns energy contribution for a given spin and field,\n"
        "  - `getName()` / `__str__()`    – for introspection, serialization and diagnostics.\n"
        "\n"
        "- Typical interactions include:\n"
        "  - exchange (Heisenberg-type),\n"
        "  - external field,\n"
        "  - anisotropy,\n"
        "  - dipole-dipole,\n"
        "  - demagnetization field,\n"
        "  - thermal fluctuation (stochastic noise).\n"
        "\n"
        "Interface:\n"
        "- `IInteraction<CoordSystem>`       – abstract interface for interaction definitions,\n"
        "- `InteractionRegistry<CoordSystem>`– registry for managing interaction instances.\n"
        "\n"
        "Concrete:\n"
        "- `ExchangeInteraction`            – isotropic Heisenberg exchange (cutoff-based),\n"
        "- `ExternalInteraction`           – constant applied magnetic field,\n"
        "- `AnisotropyInteraction`         – uniaxial anisotropy (K·(S·n)² type),\n"
        "- `DipoleDipoleInteraction`       – pairwise dipolar field (cutoff or macrocells),\n"
        "- `DemagnetizationInteraction`    – shape-field approximation with self-term subtraction,\n"
        "- `ThermalInteraction`            – stochastic thermal noise per LLG-FDT model.\n"
        "\n"
        "Python bindings:\n"
        "- Macros: `INTERACTION_TEMPLATE_BINDINGS`\n"
        "- Submodule binder: `pyBindInteractions(py::module_ &)`\n"
        "\n"
        "@note All units are SI unless explicitly stated. Effective fields are returned in [T],\n"
        "      energies in [J]. Positions and directions are dimensionless in interface.\n"
        "\n"
        "@note All interactions are read-only after construction and can be safely shared\n"
        "      across threads via the `InteractionRegistry`.";

    // -------- | CARTESIAN INTERACTIONS | --------
    py::module_ cartesian = interaction_module.def_submodule("cartesian");

    cartesian.doc() = interaction_module.doc();

    {
        using cartesian::AbstractInteraction;
        using cartesian::AbstractInteractionRegistry;
        using cartesian::AnisotropyInteraction;
        using cartesian::DemagnetizationInteraction;
        using cartesian::DipoleDipoleInteraction;
        using cartesian::ExchangeInteraction;
        using cartesian::ExternalInteraction;
        using cartesian::ThermalInteraction;

        py::class_<AbstractInteraction, std::shared_ptr<AbstractInteraction>>(
            cartesian, "AbstractInteraction"
        ) INTERACTION_TEMPLATE_BINDINGS(AbstractInteraction)
            .doc() = "@class  ExchangeInteraction\n"
                     "@brief  Heisenberg exchange interaction between magnetic moments.\n"
                     "\n"
                     "This class implements a short-range isotropic exchange interaction.\n"
                     "  The exchange field on a moment is computed as the weighted sum of its neighbors'\n"
                     "  directions within a specified cutoff radius.\n"
                     "\n"
                     "This follows the classical Heisenberg model with pairwise interactions:\n"
                     "\n"
                     "    B_exch(i) = ∑_j ( J_ij · S_j )\n"
                     "\n"
                     "where J_ij is selected based on whether moments i and j belong to the same material.\n"
                     "\n"
                     "The energy contribution per moment is computed as:\n"
                     "\n"
                     "    E = −½ · μ_s · S · B_exch\n"
                     "\n"
                     "The symmetry factor ½ ensures that each moment pair contributes only once.";

        py::class_<ExchangeInteraction, AbstractInteraction, std::shared_ptr<ExchangeInteraction>>(
            cartesian, "ExchangeInteraction"
        )
            .def(py::init<double>(), py::arg("cutoff_radius"))
            .doc() = "@class  ExchangeInteraction\n"
                     "@brief  Heisenberg exchange interaction between magnetic moments.\n"
                     "\n"
                     "This class implements a short-range isotropic exchange interaction.\n"
                     "  The exchange field on a moment is computed as the weighted sum of its neighbors'\n"
                     "  directions within a specified cutoff radius.\n"
                     "\n"
                     "This follows the classical Heisenberg model with pairwise interactions:\n"
                     "\n"
                     "    B_exch(i) = ∑_j ( J_ij · S_j )\n"
                     "\n"
                     "where J_ij is selected based on whether moments i and j belong to the same material.\n"
                     "\n"
                     "The energy contribution per moment is computed as:\n"
                     "\n"
                     "    E = −½ · μ_s · S · B_exch\n"
                     "\n"
                     "The symmetry factor ½ ensures that each moment pair contributes only once.";

        py::class_<ExternalInteraction, AbstractInteraction, std::shared_ptr<ExternalInteraction>>(
            cartesian, "ExternalInteraction"
        )
            .def(py::init<const Eigen::Vector3d &>(), py::arg("external_field"))
            .def(py::init<double, double, double>(), py::arg("sx"), py::arg("sy"), py::arg("sz"))
            .doc() = "@class  ExternalInteraction\n"
                     "@brief  Uniform magnetic field interaction (external Zeeman term).\n"
                     "\n"
                     "This interaction applies a constant external magnetic field to all moments\n"
                     "  in the system. The same field vector is used for every moment, regardless\n"
                     "  of material or position.\n"
                     "\n"
                     "The contribution to the effective field is:\n"
                     "\n"
                     "    B_ext(i) = B_external\n"
                     "\n"
                     "where B_ext is the externally applied field (in T).\n"
                     "\n"
                     "The energy of a moment in the external field is:\n"
                     "\n"
                     "    E = −μ_s · S · B_ext\n"
                     "\n"
                     "@note This interaction is stateless and does not depend on geometry.";

        py::class_<AnisotropyInteraction, AbstractInteraction, std::shared_ptr<AnisotropyInteraction>>(
            cartesian, "AnisotropyInteraction"
        )
            .def(py::init())
            .doc() =
            "@class  AnisotropyInteraction\n"
            "@brief  Magnetocrystalline anisotropy interaction.\n"
            "\n"
            "This interaction models the internal energy and effective field resulting\n"
            "  from material-specific anisotropy. It supports uniaxial anisotropy and is\n"
            "  extensible to others (e.g., cubic).\n"
            "\n"
            "The effective field for uniaxial anisotropy is computed as:\n"
            "\n"
            "    B_aniso = (2 · k_u / μ_s) · (S · n) · n\n"
            "\n"
            "where:\n"
            "  - k_u is the uniaxial anisotropy constant [J/atom],\n"
            "  - n is the anisotropy axis (unit vector).\n"
            "\n"
            "The energy associated with this field is (uniaxial case):\n"
            "\n"
            "    E = −½ · μ_s · S · B_aniso\n"
            "\n"
            "The ½ factor arises due to differentiation of the quadratic energy form.\n"
            "\n"
            "@note This interaction expects material objects to hold a valid `Anisotropy` instance.\n"
            "      If no anisotropy is defined, this interaction contributes nothing.\n"
            "\n"
            "@throws std::invalid_argument if anisotropy type is unsupported.\n"
            "\n"
            "@see UniaxialAnisotropy\n";

        py::class_<DipoleDipoleInteraction, AbstractInteraction, std::shared_ptr<DipoleDipoleInteraction>>(
            cartesian, "DipoleDipoleInteraction"
        )
            .def(py::init<double, std::string>(), py::arg("cutoff_radius"), py::arg("strategy") = "cutoff")
            .doc() =
            "@class  DipoleDipoleInteraction\n"
            "@brief  Long-range magnetostatic (dipole–dipole) interaction.\n"
            "\n"
            "This interaction accounts for the dipolar magnetic coupling between all\n"
            "  magnetic moments. It models the classical field of a magnetic dipole\n"
            "  acting on another dipole, and is responsible for magnetostatic effects\n"
            "  such as demagnetizing fields and long-range anisotropies.\n"
            "\n"
            "The effective field at position **rᵢ** due to a moment **mⱼ** located at **rⱼ**\n"
            "  is computed using the formula:\n"
            "\n"
            "    B_dd(rᵢ) = (μ₀ / 4π) ∑ⱼ [ (3(mⱼ·rᵢⱼ)rᵢⱼ / |rᵢⱼ|⁵) − (mⱼ / |rᵢⱼ|³) ]\n"
            "\n"
            "where:\n"
            "  - rᵢⱼ = rⱼ − rᵢ is the displacement vector from i to j,\n"
            "  - μ₀ is the magnetic permeability of vacuum.\n"
            "\n"
            "Energy of interaction is given by:\n"
            "\n"
            "    E = −½ · μ_s · S · B_dd\n"
            "\n"
            "Strategy options:\n"
            "- `\"cutoff\"`: sum over neighbors within a spherical radius (explicit pairwise summation);\n"
            "- `\"macrocells\"`: group distant moments into representative macro-averaged cells\n"
            "                  to accelerate long-range computations (using a macrocell grid + amplitude).\n"
            "\n"
            "@note This implementation does not use fast multipole methods or FFT;\n"
            "      complexity is depending on cutoff and macrocell density.";

        py::class_<
            DemagnetizationInteraction,
            DipoleDipoleInteraction,
            std::shared_ptr<DemagnetizationInteraction>>(cartesian, "DemagnetizationInteraction")
            .def(py::init<double, std::string>(), py::arg("cutoff_radius"), py::arg("strategy") = "cutoff")
            .doc() =
            "@class  DemagnetizationInteraction\n"
            "@brief  Effective demagnetizing field interaction based on dipole-dipole approximation.\n"
            "\n"
            "This interaction models the **macroscopic demagnetizing field** as arising from\n"
            "  long-range dipole-dipole effects that cannot be fully resolved within practical\n"
            "  cutoffs or coarse-grained representations (not all cases).\n"
            "\n"
            "Unlike `DipoleDipoleInteraction`, this class explicitly introduces a **self-term correction**\n"
            "  to account for **shape-based demagnetization**. Its role is to compensate for the dipoles\n"
            "  that are *not* explicitly included in the local sum (either due to cutoff or coarse "
            "macrocelling).\n"
            "\n"
            "❱❱ Total field:\n"
            "\n"
            "    𝐁_total = 𝐁_ext + 𝐁_demag\n"
            "\n"
            "❱❱ Dipole field:\n"
            "\n"
            "    𝐁_dip,i = (μ₀ / 4π) ∑ⱼ≠ᵢ [ 3(mⱼ·rᵢⱼ)rᵢⱼ / |rᵢⱼ|⁵ − mⱼ / |rᵢⱼ|³ ]\n"
            "\n"
            "❱❱ Demagnetizing correction (macrocell case, shape-factor):\n"
            "\n"
            "    𝐁_self = - (μ₀ / 3) · (𝐌 / V)\n"
            "\n"
            "where:\n"
            " - 𝐌 = average magnetic moment vector of macrocell [A·m²],\n"
            " - V = effective macrocell volume ≈ N · a³ (a = atomic cell size),\n"
            " - 𝐁_demag = − 𝐁_self + 𝐁_dip\n"
            "\n"
            "❱❱ In shape-based demagnetization, we define:\n"
            "\n"
            "    ℋ_demag = {\n"
            "        −ℋ_dip           // fully computed case (honest sum over j ≠ i)\n"
            "        −ℋ_self + ℋ_dip  // subtract only shape term from partial sum\n"
            "    }\n"
            "\n"
            "❱❱ The **self-term correction** is:\n"
            "\n"
            "    𝐁_self = μ₀ / 3 · (𝐌 / V)     (macrocell form)\n"
            "    𝐁_self = μ₀ / 3 · (μ_s · 𝐒 / V₀) (cutoff, single-atom effective volume)\n"
            "\n"
            "where V₀ ≈ atomic volume for single spin.\n"
            "\n"
            "This shape self-field implements the **−𝑁 · 𝐌** contribution where 𝑁 is the demagnetizing "
            "tensor\n"
            "  for a uniform ellipsoid. Here, we approximate it using the `1/3` isotropic term.\n"
            "\n"
            "@note This class explicitly **subtracts short-range dipole fields** already computed as\n"
            "      `DipoleDipoleInteraction`, and adds the missing **demagnetizing self-term**\n"
            "      as a correction due to NOT CALCULATED OTHERS DIPOLES. <!!!>\n"
            "\n"
            "@note In practice, this interaction improves convergence of magnetostatic fields\n"
            "      in heterogeneous geometries or when using coarse macrocells.";

        py::class_<ThermalInteraction, AbstractInteraction, std::shared_ptr<ThermalInteraction>>(
            cartesian, "ThermalInteraction"
        )
            .def(
                py::init<double, double, uint32_t>(),
                py::arg("temperature"),
                py::arg("dt"),
                py::arg("seed") = 0
            )
            .doc() =
            "@class  ThermalInteraction\n"
            "@brief  Stochastic thermal field contribution for finite-temperature LLG dynamics.\n"
            "\n"
            "This interaction injects random magnetic field fluctuations into the effective field\n"
            "  according to the **fluctuation–dissipation theorem** at a given temperature.\n"
            "\n"
            "This is required to simulate **finite-temperature effects** such as thermal agitation,\n"
            "  relaxation, and switching due to noise.\n"
            "\n"
            "❱❱ Model:\n"
            "The thermal field 𝐁ₜₕ is defined by:\n"
            "\n"
            "    𝐁ₜₕ = √[ 2α k_B T / (γ μ_s Δt) ] · ξ,     ξ ~ N(0,1)\n"
            "\n"
            "where:\n"
            "- α — Gilbert damping constant (unitless),\n"
            "- k_B — Boltzmann constant [J/K],\n"
            "- T — temperature [K],\n"
            "- γ — gyromagnetic ratio [rad·s⁻¹·T⁻¹],\n"
            "- μ_s — magnetic moment of spin [A·m²],\n"
            "- Δt — time step [s],\n"
            "- ξ — vector of 3 normally distributed random numbers N(0,1).\n"
            "\n"
            "This expression ensures:\n"
            "- ⟨𝐁ₜₕ(t)⟩ = 0 (mean zero),\n"
            "- ⟨𝐁ₜₕ(t) · 𝐁ₜₕ(t')⟩ ∝ δ(t − t') (delta-correlated).\n"
            "\n"
            "@note This interaction does **not** contribute to total energy (it is purely statistical).\n"
            "@note Internally uses `std::mt19937` and `std::normal_distribution<double>` for sampling.\n"
            "      RNG is seeded with a user-provided seed or system random device (if zero)."

            ;

        py::class_<AbstractInteractionRegistry, std::shared_ptr<AbstractInteractionRegistry>>(
            cartesian, "InteractionRegistry"
        )
            .def(py::init<>())
            .def(py::init<RegistryContainer<AbstractInteraction>>(), py::arg("container"))
                REGISTRY_TEMPLATE_BINDINGS(AbstractInteractionRegistry)
            .doc() = "@typedef AbstractInteractionRegistry\n"
                     "@brief  Registry of physical interactions for a given coordinate system.\n"
                     "\n"
                     "The `InteractionRegistry` stores instances of `Interaction` and provides\n"
                     "  global access to the active set of interactions in the system.\n"
                     "\n"
                     "Registries are read-only after construction and therefore safe for concurrent access.\n"
                     "\n"
                     "Intended usage:\n"
                     "  - Maintain a globally available set of interactions used in the simulation.\n"
                     "  - Allow simulation loops to iterate over all interactions when calculating\n"
                     "    effective fields and energy contributions.\n"
                     "\n"
                     "@note Interaction instances are stored as shared pointers, allowing polymorphic usage\n"
                     "      and dynamic dispatch at runtime.";
    }
}

#endif // ! __INTERACTIONS_HPP__
