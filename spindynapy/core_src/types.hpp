#ifndef __SPINDYNAPY_TYPES_HPP__
#define __SPINDYNAPY_TYPES_HPP__

/**
 * @file   types.hpp
 * @brief  Fundamental data types used by the SpinDynaPy atomistic core.
 *
 * The header provides lightweight abstractions that represent materials,
 * magnetic moments, coordinate systems and helper concepts.  All physical
 * quantities are expressed in **SI units** unless stated otherwise.
 *
 * Exposed entities
 * - Scalar aliases: `regnum`, `EffectiveField`, `EffectiveFieldVector`.
 * - `Anisotropy` / `UniaxialAnisotropy`.
 * - `Material` – immutable magnetic-property bundle.
 * - Namespaces `cartesian` / `spherical` with `Coordinates`, `Direction`,
 *     `Moment` tags or concrete classes.
 * - Tag structs `CartesianCoordSystem`, `SphericalCoordSystem`.
 * - `CoordSystemConcept` – compile-time interface check.
 * - Helper macro `BIND_STR_REPR`.
 * - Function `pyBindTypes()` – publishes the above to Python.
 *
 * @note  This file deliberately avoids heavy algorithms, just bricks.
 *
 * @copyright 2025 SpinDynaPy
 */

#include "constants.hpp"

#include <format>
#include <memory>
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>
#include <stdexcept>
#include <string>

namespace py = pybind11;

namespace spindynapy {

// ==========================================================================
//  Scalar aliases
// ==========================================================================

/**
 * @typedef regnum
 * @brief   Registry identifier used as a primary key throughout the code-base.
 */
typedef short regnum;

/** @brief  Vector representing the net effective magnetic field */
using EffectiveField = Eigen::Vector3d;

/** @brief  Contiguous array of effective fields (one per spin). */
using EffectiveFieldVector = std::vector<EffectiveField>;

// ==========================================================================
//  Anisotropy hierarchy
// ==========================================================================

/**
 * @class  Anisotropy
 * @brief  Abstract interface for magnetic anisotropy terms.
 *
 * Concrete implementations compute energy contributions based on the
 * orientation of a spin relative to crystal axes.
 *
 * @note   Default constructor is `protected`, preventing direct instantiation.
 */
class Anisotropy {
  protected:
    /** @brief Protected default constructor – class is abstract. */
    Anisotropy() = default;

  public:
    /** @brief destructor (default) */
    virtual ~Anisotropy() = default;
    /**
     * @brief Human-readable description.
     *
     * Must be overridden by every concrete anisotropy.
     *
     * @returns String describing parameters.
     */
    virtual std::string __str__() const { throw std::logic_error("Method __str__ Not implemented"); };
};

/**
 * @class  UniaxialAnisotropy
 * @brief  Uniaxial (universal easy axis) magnetic anisotropy.
 *
 * @param axis     Unit axis specifying the easy axis direction.
 * @param constant Per-atom anisotropy constant (J).
 */
class UniaxialAnisotropy : public Anisotropy {
  public:
    Eigen::Vector3d axis; ///< Normalised easy axis.
    double constant;      ///< Anisotropy constant per atom (J).

    /** @brief Construct with axis `(x,y,z)` (m) and constant `k` = `K` per atom (J). */
    UniaxialAnisotropy(const Eigen::Vector3d &axis, double constant)
        : axis(axis.normalized()), constant(constant) {}

    /**
     * @brief Human-readable description.
     *
     * @returns String describing parameters.
     */
    std::string __str__() const override {
        return std::format(
            "UniaxialAnisotropy(axis=({}, {}, {}), constant={:.2e})", axis.x(), axis.y(), axis.z(), constant
        );
    }
};

// ==========================================================================
//  Material
// ==========================================================================

/**
 * @class  Material
 * @brief  Immutable bundle of intrinsic magnetic properties.
 *
 * A `Material` instance is stored once in a registry and shared by reference
 * throughout the simulation to avoid duplication.
 *
 * @param exchange_constant_J Exchange constant (J) – strength of exchange interaction.
 * @param atomic_magnetic_saturation_magnetization Magnetic saturation (Bohr magnetons) – maximum magnetic
 *   moment per atom.
 * @param unit_cell_size Size of the unit cell (m) – defines the periodicity of the material.
 * @param atom_cell_size Size of the atomic cell (m) – fraction of the unit cell occupied by a single atom.
 * @param gyromagnetic_ratio Gyromagnetic ratio (rad/s/T) – relates magnetic moment to angular momentum.
 * @param damping_constant Damping constant (dimensionless) – describes energy dissipation in the system.
 * @param anisotropy Anisotropy instance (optional, Anisotropy class) – describes the material's magnetic
 *   anisotropy.
 * @param atomic_magnetic_saturation_absolute Absolute magnetic saturation (A·m²) – derived from
 * @param exchange_monomaterial_prefix Exchange monomaterial prefix (T) – derived from the exchange constant
 *      for a single material (no interface).
 */
class Material {
  public:
    regnum material_number;                          ///< Registry index.
    double exchange_constant_J;                      ///< Heisenberg exchange integral (*J*).
    double atomic_magnetic_saturation_magnetization; ///< Saturation μ (Bohr magnetons).
    std::shared_ptr<Anisotropy> anisotropy;          ///< Optional intrinsic anisotropy.
    double gyromagnetic_ratio;                       ///< γ Gyromagnetic ratio (rad/s/T)
    double damping_constant;                         ///< Gilbert damping α.
    double unit_cell_size;                           ///< Edge length of unit cell (m).
    double atom_cell_size;                           ///< Edge length of atomic cell (m).

    //

    double atomic_magnetic_saturation_absolute; ///< μ_s = M_s·μ_B      [A·m²]
    double exchange_monomaterial_prefix;        ///< C   = J / μ_s      [T]

    /**
     * @brief Construct a fully-specified material.
     *
     * @param material_number  Registry id.
     * @param exchange_constant_J  Heisenberg exchange integral (*J*).
     * @param atomic_magnetic_saturation_magnetization  Saturation moment (μ_B).
     * @param unit_cell_size   Unit cell size (m).
     * @param atom_cell_size   Atomic cell size (m).
     * @param gyromagnetic_ratio  γ (rad s⁻¹ T⁻¹). Defaults to free-spin value.
     * @param damping_constant     Gilbert damping α. Defaults to 0.1.
     * @param anisotropy           Optional intrinsic anisotropy (may be `nullptr`).
     */
    Material(
        regnum material_number,
        double exchange_constant_J,
        double atomic_magnetic_saturation_magnetization,
        double unit_cell_size,
        double atom_cell_size,
        double gyromagnetic_ratio = constants::FREE_SPIN_GYROMAGNETIC_RATIO, // ~=~ Co Ni Fe
        double damping_constant = 0.1,
        std::shared_ptr<Anisotropy> anisotropy = nullptr
    )
        : material_number(material_number),
          exchange_constant_J(exchange_constant_J),
          atomic_magnetic_saturation_magnetization(atomic_magnetic_saturation_magnetization),
          anisotropy(anisotropy),
          gyromagnetic_ratio(gyromagnetic_ratio),
          damping_constant(damping_constant),
          unit_cell_size(unit_cell_size),
          atom_cell_size(atom_cell_size) {
        this->atomic_magnetic_saturation_absolute =
            atomic_magnetic_saturation_magnetization * constants::BOHR_MAGNETON;
        this->exchange_monomaterial_prefix = exchange_constant_J / this->atomic_magnetic_saturation_absolute;
    };

    /**
     * @brief Registry id accessor.
     *
     * @returns `regnum` unique number.
     */
    regnum getNumber() const { return this->material_number; };

    /**
     * @brief Human-readable description.
     *
     * @returns String describing parameters.
     */
    std::string __str__() const { return std::to_string(this->material_number); };

    /**
     * @brief Python-like string representation.
     *
     * @returns Formatted string with all properties.
     */
    std::string __repr__() const {
        return std::format(
            "Material(material_number={}, exchange_constant_J={}, "
            "atomic_magnetic_saturation_magnetization={})",
            this->material_number,
            this->exchange_constant_J,
            this->atomic_magnetic_saturation_magnetization
        );
    }

    /**
     * @brief Equality operator.
     *
     * @param other Material to compare with.
     * @details Compares only the `material_number` field.
     *
     * @returns `true` if materials match.
     */
    bool operator==(const Material &other) const { return material_number == other.material_number; }
};

/**
 * @class  MaterialInterface
 * @brief  Interface between two materials.
 */
class MaterialInterface {};

// ==========================================================================
//  Region
// ==========================================================================

/**
 * @class  Region
 * @brief  Generic spatial region used for constraints / masks / etc.
 */
class Region {
  protected:
    double _accuracy; ///< Tolerance parameter

  public:
    Region(double accuracy) : _accuracy(accuracy) {};
};

// ==========================================================================
//  Cartesian coordinate system
// ==========================================================================

namespace cartesian {

/**
 * @class  Coordinates
 * @brief  Immutable 3D coordinates in Cartesian space.
 *
 * Represents a point in 3D space with methods for distance calculations.
 */
class Coordinates {
  protected:
    Eigen::Vector3d _coords; ///< Internal storage (m).

  public:
    /** @brief Construct from components *(x,y,z)* (m). */
    Coordinates(double x, double y, double z) noexcept : _coords(x, y, z) {};
    /** @brief Construct from existing vector (m). */
    Coordinates(const Eigen::Vector3d &coords) noexcept : _coords(coords) {};

    /**
     * @brief Human-readable description.
     *
     * @returns String describing parameters.
     */
    std::string __str__() const {
        return std::format("{:.3f}\t{:.3f}\t{:.3f}", _coords(0) * 1e10, _coords(1) * 1e10, _coords(2) * 1e10);
    };

    /**
     * @brief Python-like string representation.
     *
     * @returns Formatted string with all properties.
     */
    std::string __repr__() const {
        return std::format("Coordinates(x={:.3f}, y={:.3f}, z={:.3f})", _coords(0), _coords(1), _coords(2));
    };

    /**
     * @brief Euclidean distance to another point.
     *
     * @param other Reference point.
     *
     * @returns Distance (m).
     */
    double getDistanceFrom(const Coordinates &other) const { return (this->_coords - other._coords).norm(); };

    /**
     * @brief Squared distance (avoids `sqrt` if only comparison is needed).
     *
     * @param other Reference point.
     *
     * @returns Distance² (m²).
     */
    double getDistanceSquareFrom(const Coordinates &other) const {
        return (this->_coords - other._coords).squaredNorm();
    };

    /**
     * @brief Overwrite the coordinates with explicit components.
     *
     * After the call the internal vector is `(x, y, z)` **metres**.
     *
     * @param x New *x*-component (m).
     * @param y New *y*-component (m).
     * @param z New *z*-component (m).
     *
     * @returns void - mutates the internal state.
     */
    void setCoordinates(double x, double y, double z) {
        _coords[0] = x;
        _coords[1] = y;
        _coords[2] = z;
        return;
    }

    /**
     * @brief Overwrite the coordinates with a full vector.
     *
     * The argument is copied as-is.
     *
     * @param coords 3-vector in meters.
     *
     * @returns void - mutates the internal state.
     */
    void setCoordinates(Eigen::Vector3d coords) {
        _coords[0] = coords[0];
        _coords[1] = coords[1];
        _coords[2] = coords[2];
        return;
    }

    /**
     * @brief Read-only component access.
     *
     * @param idx Index 0→*x*, 1→*y*, 2→*z*.
     *
     * @returns Requested component (m).
     */
    double operator[](size_t index) const { return _coords[index]; }

    /**
     * @brief Direct, mutable access to the underlying vector.
     *
     * This is a low-level operation that allows direct manipulation of the coordinates.
     * @note Use with care – callers must keep the quantity in metres.
     *
     * @returns Reference to the internal `Eigen::Vector3d`.
     */
    Eigen::Vector3d &asVector() { return this->_coords; }
};

/**
 * @class  Direction
 * @brief  Immutable 3D direction vector in Cartesian space.
 *
 * Represents a unit vector in 3D space with methods for normalization and access.
 */
class Direction {
  protected:
    Eigen::Vector3d _vector; ///< Internal storage (unit vector).

  public:
    /** @brief Construct from components; auto-normalises. */
    Direction(double x, double y, double z) noexcept : _vector(x, y, z) { _vector.normalize(); };
    /** @brief Construct from existing vector; auto-normalises. */
    Direction(const Eigen::Vector3d &vector) noexcept : _vector(vector) { _vector.normalize(); };

    /**
     * @brief Human-readable description.
     *
     * @returns String describing parameters.
     */
    std::string __str__() const {
        return std::format("{:.3f}\t{:.3f}\t{:.3f}", _vector(0), _vector(1), _vector(2));
    };

    /**
     * @brief Python-like string representation.
     *
     * @returns Formatted string with all properties.
     */
    std::string __repr__() const {
        return std::format("Direction(x={:.3f}, y={:.3f}, z={:.3f})", _vector(0), _vector(1), _vector(2));
    };

    /**
     * @brief Replace the direction using explicit components.
     *
     * @details The vector is immediately normalised to unit length.
     *
     * @param x New *x*-component.
     * @param y New *y*-component.
     * @param z New *z*-component.
     *
     * @returns void.
     */
    void setDirection(double x, double y, double z) {
        _vector[0] = x;
        _vector[1] = y;
        _vector[2] = z;
        _vector.normalize();
        return;
    }

    /**
     * @brief Replace the direction using a full vector.
     *
     * @details The vector is immediately normalised to unit length.
     *
     * @param vector New direction vector.
     *
     * @returns void.
     */
    void setDirection(Eigen::Vector3d vector) {
        _vector = vector;
        _vector.normalize();
        return;
    }

    /**
     * @brief Read-only component access.
     *
     * @param index Index 0→*x*, 1→*y*, 2→*z*.
     *
     * @returns Requested component.
     */
    double operator[](size_t index) const { return _vector[index]; }

    /**
     * @brief Direct, mutable access to the underlying vector.
     *
     * This is a low-level operation that allows direct manipulation of the direction.
     * @note Use with care – callers must keep the quantity as a unit vector.
     *
     * @returns Reference to the internal `Eigen::Vector3d`.
     */
    Eigen::Vector3d &asVector() { return this->_vector; }

    /**
     * @brief Re-normalise if external code modified the vector length.
     *
     * @returns void - mutates the internal state.
     */
    void normalize() { return _vector.normalize(); }
};

/**
 * @class  Moment
 * @brief  Immutable magnetic moment in Cartesian space.
 *
 * Represents a magnetic moment with its position, direction, and material properties.
 */
class Moment {
  protected:
    std::unique_ptr<Coordinates> _coordinates; ///< Position.
    std::unique_ptr<Direction> _direction;     ///< Unit spin vector.
    std::shared_ptr<Material> _material;       ///< Intrinsic material properties.

  public:
    double amplitude = 1.0; ///< Scale factor (e.g. macro-spin grouping). Represents correlating between
                            ///< the moment's directions in group.

  public:
    /** @brief Construct from lvalue references. */
    Moment(
        const Coordinates &coordinates, const Direction &direction, std::shared_ptr<Material> material
    ) noexcept
        : _coordinates(std::make_unique<Coordinates>(coordinates)),
          _direction(std::make_unique<Direction>(direction)),
          _material(material) {};

    /** @brief Construct from rvalue references (move). */
    Moment(Coordinates &&coordinates, Direction &&direction, std::shared_ptr<Material> material) noexcept
        : _coordinates(std::make_unique<Coordinates>(std::move(coordinates))),
          _direction(std::make_unique<Direction>(std::move(direction))),
          _material(std::move(material)) {};

    /**
     * @brief Access the spin direction.
     * @returns Modifiable reference to the underlying @ref Direction.
     */
    Direction &getDirection() { return *_direction; };

    /**
     * @brief Overwrite the spin direction from components.
     * @param sx,sy,sz Components of the new vector (will be normalised).
     * @returns void - mutates the internal state.
     */
    void setDirection(double sx, double sy, double sz) { this->_direction->setDirection(sx, sy, sz); }

    /**
     * @brief Overwrite the spin direction from an existing vector.
     * @param new_direction_vector 3-vector (will be normalised internally).
     * @returns void - mutates the internal state.
     */
    void setDirection(Eigen::Vector3d new_direction_vector) {
        this->_direction->setDirection(new_direction_vector);
    }

    /**
     * @brief Update the spatial coordinates of the moment.
     * @param new_coordinates_vector 3-vector in metres.
     * @returns void - mutates the internal state.
     */
    void setCoordinates(Eigen::Vector3d new_coordinates_vector) {
        this->_coordinates->setCoordinates(new_coordinates_vector);
    }

    /**
     * @brief Set the material pointer.
     * @param new_material Shared pointer to a @ref Material.
     * @returns void - mutates the internal state.
     */
    void setMaterial(std::shared_ptr<Material> new_material) { this->_material = new_material; }

    /**
     * @brief Access the coordinates.
     * @returns Modifiable reference to the underlying @ref Coordinates.
     */
    Coordinates &getCoordinates() { return *_coordinates; };

    /**
     * @brief Access the associated material.
     * @returns Modifiable reference to the underlying @ref Material.
     */
    Material &getMaterial() { return *_material; };

    /**
     * @brief Shared-ptr accessor.
     * @returns Shared pointer to the material.
     */
    std::shared_ptr<Material> &getMaterialSharedPtr() { return this->_material; };

    /**
     * @brief Deep copy.
     * @returns New `Moment` allocated on the heap.
     */
    std::shared_ptr<Moment> clone() const {
        return std::make_shared<Moment>(*_coordinates, *_direction, _material);
    }

    /**
     * @brief Human-readable description.
     *
     * @returns String describing parameters.
     */
    std::string __str__() const {
        return _material->__str__() + "\t0\t" + _coordinates->__str__() + "\t" + _direction->__str__();
    };

    /**
     * @brief Python-like string representation.
     *
     * @returns Formatted string with all properties.
     */
    std::string __repr__() const {
        return std::format(
            "Moment(coordinates={}, direction={}, material={})",
            _coordinates->__repr__(),
            _direction->__repr__(),
            _material->__repr__()
        );
    };
};

} // namespace cartesian

// ==========================================================================
//  Spherical coordinate system (marker only for now)
// ==========================================================================

namespace spherical {

/**
 * @struct Coordinates
 * @brief  Marker for future spherical coordinates *(r,θ,φ)*.
 */
class Coordinates {};

/**
 * @struct Direction
 * @brief  Marker for future spherical direction vector.
 */
class Direction {};

/**
 * @struct Moment
 * @brief  Marker for future spherical moment implementation.
 */
class Moment {};

} // namespace spherical

// ==========================================================================
//  Coordinate-system tags and concept
// ==========================================================================

/**
 * @concept CoordSystemConcept
 * @brief   Compile-time check for coordinate-system compatibility.
 *
 * Requirements:
 * - Nested aliases `Moment`, `Coordinates`, `Direction`
 * - Constexpr field `name` convertible to `const char*`
 */
template <typename CoordSystem>
concept CoordSystemConcept = requires {
    typename CoordSystem::Moment;
    typename CoordSystem::Coordinates;
    typename CoordSystem::Direction;
    { CoordSystem {}.name } -> std::convertible_to<const char *>;
};

/**
 * @struct  CoordSystem
 * @brief   Base class for coordinate systems.
 *
 * Provides a common interface for coordinate systems.
 */
struct CoordSystem {};

/**
 * @struct  SphericalCoordSystem
 * @brief   Aggregates spherical marker types.
 */
struct SphericalCoordSystem : public CoordSystem {
    using Moment = spherical::Moment;
    using Coordinates = spherical::Coordinates;
    using Direction = spherical::Direction;
    const char *name = "SphericalCoordSystem";
};

/**
 * @struct  CartesianCoordSystem
 * @brief   Aggregates cartesian concrete types.
 */
struct CartesianCoordSystem : public CoordSystem {
    using Moment = cartesian::Moment;
    using Coordinates = cartesian::Coordinates;
    using Direction = cartesian::Direction;
    const char *name = "CoordSystem";
};

// --------------------------------------------------------------------------
namespace cartesian {
using NamespaceCoordSystem = CartesianCoordSystem;
}
namespace spherical {
using NamespaceCoordSystem = SphericalCoordSystem;
}
// --------------------------------------------------------------------------

}; // namespace spindynapy

// ==========================================================================
//  Python bindings
// ==========================================================================

/**
 * @def   BIND_STR_REPR
 * @brief Convenience macro binding both `__str__` and `__repr__`.
 * @param cls  C++ class exposing member functions `__str__` and `__repr__` to Python.
 */
#define BIND_STR_REPR(cls) .def("__str__", &cls::__str__).def("__repr__", &cls::__repr__)

/**
 * @brief Bind basic C++ data types to a Python sub-module.
 *
 * @param module Parent PyBind11 module (usually the core extension module).
 * @returns void – extends the parent module hierarchy.
 */
inline void pyBindTypes(py::module_ &module) {
    using namespace spindynapy;

    // ---------------------------------------------------------------------
    //  Create submodule
    // ---------------------------------------------------------------------

    py::module_ types_module = module.def_submodule("types");
    types_module.doc() = "@brief  Fundamental data types used by the SpinDynaPy atomistic core.\n"
                         "\n"
                         "The header provides lightweight abstractions that represent materials,\n"
                         "magnetic moments, coordinate systems and helper concepts.  All physical\n"
                         "quantities are expressed in **SI units** unless stated otherwise.\n"
                         "\n"
                         "Exposed entities\n"
                         "- Scalar aliases: `regnum`, `EffectiveField`, `EffectiveFieldVector`.\n"
                         "- `Anisotropy` / `UniaxialAnisotropy`.\n"
                         "- `Material` – immutable magnetic-property bundle.\n"
                         "- Namespaces `cartesian` / `spherical` with `Coordinates`, `Direction`,\n"
                         "    `Moment` tags or concrete classes.\n"
                         "- Tag structs `CartesianCoordSystem`, `SphericalCoordSystem`.\n"
                         "- `CoordSystemConcept` – compile-time interface check.\n"
                         "- Helper macro `BIND_STR_REPR`.\n"
                         "- Function `pyBindTypes()` – publishes the above to Python.\n"
                         "\n"
                         "@note  This file deliberately avoids heavy algorithms, just bricks.";

    // ---------------------------------------------------------------------
    //  Anisotropy bindings
    // ---------------------------------------------------------------------

    py::class_<Anisotropy, std::shared_ptr<Anisotropy>>(types_module, "Anisotropy")
        .def("__str__", &Anisotropy::__str__)
        .doc() = "@class  Anisotropy\n"
                 "@brief  Abstract interface for magnetic anisotropy terms.\n"
                 "\n"
                 "Concrete implementations compute energy contributions based on the\n"
                 "orientation of a spin relative to crystal axes.\n"
                 "\n"
                 "@note   Default constructor is `protected`, preventing direct instantiation.";

    py::class_<UniaxialAnisotropy, Anisotropy, std::shared_ptr<UniaxialAnisotropy>>(
        types_module, "UniaxialAnisotropy"
    )
        .def(py::init<const Eigen::Vector3d &, double>(), py::arg("axis"), py::arg("constant"))
        .def_readwrite("axis", &UniaxialAnisotropy::axis)
        .def_readwrite("constant", &UniaxialAnisotropy::constant)
        .doc() = "@class  UniaxialAnisotropy\n"
                 "@brief  Uniaxial (universal easy axis) magnetic anisotropy.\n"
                 "\n"
                 "@param axis     Unit axis specifying the easy axis direction.\n"
                 "@param constant Per-atom anisotropy constant (J).\n"
                 "\n"
                 "This class implements a uniaxial anisotropy model, where the energy depends on the angle\n"
                 "between the spin direction and the specified easy axis.";

    // ---------------------------------------------------------------------
    //  Material bindings
    // ---------------------------------------------------------------------

    py::class_<Material, std::shared_ptr<Material>>(types_module, "Material") BIND_STR_REPR(Material)
        .def(
            py::init<regnum, double, double, double, double, double, double, std::shared_ptr<Anisotropy>>(),
            py::arg("material_number"),
            py::arg("exchange_constant_J"),
            py::arg("atomic_magnetic_saturation_magnetization"),
            py::arg("unit_cell_size"),
            py::arg("atom_cell_size"),
            py::arg("gyromagnetic_ratio") = constants::FREE_SPIN_GYROMAGNETIC_RATIO,
            py::arg("damping_constant") = 0.1,
            py::arg("anisotropy").none(true) = pybind11::none()
        )
        .def(
            "get_number",
            &Material::getNumber,
            py::doc("@brief Registry id accessor.\n"
                    "\n"
                    "@returns `regnum` unique number.")
        )
        .def_readwrite("material_number", &Material::material_number)
        .def_readwrite("exchange_constant_J", &Material::exchange_constant_J)
        .def_readwrite(
            "atomic_magnetic_saturation_magnetization", &Material::atomic_magnetic_saturation_magnetization
        )
        .def_readwrite("anisotropy", &Material::anisotropy)
        .def_readwrite("gyromagnetic_ratio", &Material::gyromagnetic_ratio)
        .def_readwrite("damping_constant", &Material::damping_constant)
        .def_readwrite("unit_cell_size", &Material::unit_cell_size)
        .def_readwrite("atom_cell_size", &Material::atom_cell_size)
        .def_readwrite("atomic_magnetic_saturation_absolute", &Material::atomic_magnetic_saturation_absolute)
        .def_readwrite("exchange_monomaterial_prefix", &Material::exchange_monomaterial_prefix)
        .doc() =
        "@class  Material\n"
        "@brief  Immutable bundle of intrinsic magnetic properties.\n"
        "\n"
        "A `Material` instance is stored once in a registry and shared by reference\n"
        "throughout the simulation to avoid duplication.\n"
        "\n"
        "@param exchange_constant_J Exchange constant (J) – strength of exchange interaction.\n"
        "@param atomic_magnetic_saturation_magnetization Magnetic saturation (Bohr magnetons) – maximum "
        "magnetic\n"
        "   moment per atom.\n"
        "@param unit_cell_size Size of the unit cell (m) – defines the periodicity of the material.\n"
        "@param atom_cell_size Size of the atomic cell (m) – fraction of the unit cell occupied by a single "
        "atom.\n"
        "@param gyromagnetic_ratio Gyromagnetic ratio (rad/s/T) – relates magnetic moment to angular "
        "momentum.\n"
        "@param damping_constant Damping constant (dimensionless) – describes energy dissipation in the "
        "system.\n"
        "@param anisotropy Anisotropy instance (optional, Anisotropy class) – describes the material's "
        "magnetic\n"
        "   anisotropy.";

    // ---------------------------------------------------------------------
    //  Region bindings
    // ---------------------------------------------------------------------

    py::class_<Region>(module, "Region").def(py::init<double>(), py::arg("accuracy")).doc() =
        "@class  Region\n"
        "@brief  Generic spatial region used for constraints / masks / etc.\n"
        "\n"
        "This class represents a spatial region with an accuracy parameter that can be used for various "
        "purposes "
        "  such as defining boundaries, constraints, or regions of interest in the simulation.\n"
        "\n"
        "@param accuracy Tolerance parameter.";

    // ---------------------------------------------------------------------
    //  CoordSystems bindings
    // ---------------------------------------------------------------------

    py::class_<CartesianCoordSystem, std::shared_ptr<CartesianCoordSystem>>(module, "CartesianCoordSystem")
        .def(py::init<>())
        .def_readonly("name", &CartesianCoordSystem::name)
        .doc() =
        "@struct  CartesianCoordSystem\n"
        "@brief   Aggregates cartesian concrete types.\n"
        "\n"
        "This struct serves as a tag for the Cartesian coordinate system, aggregating the concrete types "
        "`Moment`, `Coordinates`, and `Direction` from the `cartesian` namespace.";

    py::class_<SphericalCoordSystem, std::shared_ptr<SphericalCoordSystem>>(module, "SphericalCoordSystem")
        .def(py::init<>())
        .def_readonly("name", &SphericalCoordSystem::name)
        .doc() =
        "@struct  SphericalCoordSystem\n"
        "@brief   Aggregates spherical marker types.\n"
        "\n"
        "This struct serves as a tag for the Spherical coordinate system, aggregating the marker types "
        "`Moment`, `Coordinates`, and `Direction` from the `spherical` namespace.";

    // ---------------------------------------------------------------------
    //  Cartesian submodule bindings
    // ---------------------------------------------------------------------

    py::module_ cartesian = types_module.def_submodule("cartesian");

    cartesian.doc() =
        "@brief  Cartesian coordinate system types.\n"
        "\n"
        "This submodule contains concrete implementations of the Cartesian coordinate system, "
        "  including `Coordinates`, `Direction`, and `Moment` classes.\n"
        "\n"
        "These types are used to represent points, directions, and magnetic moments in a 3D Cartesian "
        "  space, providing methods for distance calculations and vector operations.";

    using cartesian::Coordinates;
    using cartesian::Direction;
    using cartesian::Moment;

    py::class_<Coordinates>(cartesian, "Coordinates")
        .def(py::init<double, double, double>(), py::arg("x"), py::arg("y"), py::arg("z"))
        .def(py::init<const Eigen::Vector3d &>(), py::arg("coords"))
        .def(
            "get_distance_from",
            &Coordinates::getDistanceFrom,
            py::arg("other"),
            py::doc("@brief Calculate Euclidean distance to another point.\n"
                    "\n"
                    "@param other Reference point.\n"
                    "\n"
                    "@returns Distance (m).")
        )
        .def(
            "get_distance_square_from",
            &Coordinates::getDistanceSquareFrom,
            py::arg("other"),
            py::doc("@brief Calculate squared distance to another point (avoids `sqrt` if only comparison is "
                    "needed).\n"
                    "\n"
                    "@param other Reference point.\n"
                    "\n"
                    "@returns Distance² (m²).")
        )
        .def(
            "__getitem__",
            &Coordinates::operator[],
            py::arg("index"),
            py::doc("@brief Read-only component access.\n"
                    "\n"
                    "@param index Index 0→*x*, 1→*y*, 2→*z*.\n"
                    "\n"
                    "@returns Requested component (m).")
        )
        .def_property_readonly("x", [](const Coordinates &v) { return v[0]; })
        .def_property_readonly("y", [](const Coordinates &v) { return v[1]; })
        .def_property_readonly("z", [](const Coordinates &v) { return v[2]; }) BIND_STR_REPR(Coordinates)
        .doc() = "@class  Coordinates\n"
                 "@brief  Immutable 3D coordinates in Cartesian space.\n"
                 "\n"
                 "Represents a point in 3D space with methods for distance calculations.\n"
                 "\n"
                 "@param x, y, z Coordinates in metres (m).\n"
                 "\n"
                 "Provides methods to calculate distances and access individual components.";

    py::class_<Direction>(cartesian, "Direction")
        .def(py::init<double, double, double>(), py::arg("x"), py::arg("y"), py::arg("z"))
        .def(py::init<const Eigen::Vector3d &>(), py::arg("vector")) BIND_STR_REPR(Direction)
        .def(
            "__getitem__",
            &Direction::operator[],
            py::arg("index"),
            py::doc("Read-only component access.\n"
                    "\n"
                    "@param index Index 0→*x*, 1→*y*, 2→*z*.\n"
                    "\n"
                    "@returns Requested component.")
        )
        .def_property_readonly("x", [](const Direction &v) { return v[0]; })
        .def_property_readonly("y", [](const Direction &v) { return v[1]; })
        .def_property_readonly("z", [](const Direction &v) { return v[2]; })
        .doc() = "@class  Direction\n"
                 "@brief  Immutable 3D direction vector in Cartesian space.\n"
                 "\n"
                 "Represents a unit vector in 3D space with methods for normalization and access.\n"
                 "\n"
                 "@param x, y, z Direction components (will be normalised to unit length).\n"
                 "\n"
                 "Provides methods to set the direction and access individual components.";

    py::class_<Moment, std::shared_ptr<Moment>>(cartesian, "Moment")
        .def(
            py::init([](const Coordinates &coords, const Direction &dir, Material &mat) {
                return new Moment(coords, dir, std::make_shared<Material>(mat));
            }),
            py::arg("coordinates"),
            py::arg("direction"),
            py::arg("material")
        ) BIND_STR_REPR(Moment)
        .def(
            "get_material",
            &Moment::getMaterial,
            py::return_value_policy::reference,
            py::doc("@brief Access the associated material.\n"
                    "\n"
                    "@returns Modifiable reference to the underlying @ref Material.")
        )
        .def(
            "get_coordinates",
            &Moment::getCoordinates,
            py::return_value_policy::reference,
            py::doc("@brief Access the coordinates.\n"
                    "\n"
                    "@returns Modifiable reference to the underlying @ref Coordinates.")
        )
        .def(
            "get_direction",
            &Moment::getDirection,
            py::return_value_policy::reference,
            py::doc("@brief Access the spin direction.\n"
                    "\n"
                    "@returns Modifiable reference to the underlying @ref Direction.")
        )
        .def_readonly("amplitude", &Moment::amplitude)
        .doc() = "@class  Moment\n"
                 "@brief  Immutable magnetic moment in Cartesian space.\n"
                 "\n"
                 "Represents a magnetic moment with its position, direction, and material properties.\n"
                 "\n"
                 "@param coordinates Position in Cartesian space (m).\n"
                 "@param direction Unit vector representing the spin direction.\n"
                 "@param material Material properties associated with the moment.\n"
                 "\n"
                 "Provides methods to access and modify the moment's properties.";

    // ---------------------------------------------------------------------
    //  Spherical submodule bindings
    // ---------------------------------------------------------------------

    py::module_ spherical = types_module.def_submodule("spherical");

    spherical.doc() =
        "@brief  Spherical coordinate system types.\n"
        "\n"
        "This submodule contains marker types for the Spherical coordinate system, "
        "  including `Coordinates`, `Direction`, and `Moment` classes.\n"
        "\n"
        "These types are intended for future use in representing points, directions, and magnetic "
        "  moments in a spherical coordinate system.";
}

#endif // ! __SPINDYNAPY_TYPES_HPP__
