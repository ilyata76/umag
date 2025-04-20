#ifndef __INTERACTIONS_HPP__
#define __INTERACTIONS_HPP__

/**
 * Интерфейсы взаимодействий между элементами системы.
 * Таковыми могут быть потенциальные поля, силы, etc.
 */

#include "constants.hpp"
#include "geometries.hpp"
#include "registries.hpp"
#include "types.hpp"

#include <format>
#include <memory>
#include <pybind11/pybind11.h>
#include <stdexcept>
#include <string>

namespace py = pybind11;

namespace spindynapy {

/**
 * Базовый интерфейс взаимодействий (сил, полей, etc.)
 */
template <CoordSystemConcept CoordSystem> class IInteraction {
  protected:
    IInteraction() = default;

  public:
    virtual ~IInteraction() = default;

    virtual void saveStepBuffer(std::vector<EffectiveField> buffer_from) {
        throw std::logic_error("Method saveStepBuffer Not implemented");
    }

    virtual EffectiveField calculateFieldContribution(
        size_t moment_index, IGeometry<CoordSystem> &geometry, MaterialRegistry &material_registry
    ) const {
        throw std::logic_error("Method calculateFieldContribution Not implemented");
    };

    virtual std::string __str__() const { return nullptr; };
    virtual std::string __repr__() const { return nullptr; };
};

/**
 * Обменное взаимодействие
 */

using CartesianAbstractInteraction = IInteraction<CartesianCoordSystem>;

class CartesianExchangeInteraction : public CartesianAbstractInteraction {
  protected:
    double _cutoff_radius;

  public:
    CartesianExchangeInteraction(double cutoff_radius) : _cutoff_radius(cutoff_radius) {};

    virtual EffectiveField calculateFieldContribution(
        size_t moment_index, IGeometry<CartesianCoordSystem> &geometry, MaterialRegistry &material_registry
    ) const override {
        // H = (J over abs(magnetic_moment)) * SUM on S_neighbors
        // S = mu / |mu| normalized (!)
        EffectiveField exchange_field = EffectiveField::Zero();
        auto &current_moment = geometry[moment_index];
        auto &current_material = current_moment.getMaterial();

        // величина магнитного момента
        auto atomic_magnetic_moments_norm =
            current_material.atomic_magnetic_saturation_magnetization * constants::BOHR_MAGNETON;
        const auto GENERALIZED_PREFIX = current_material.exchange_constant_J / atomic_magnetic_moments_norm;

        auto neighbor_indices = geometry.getNeighbors(moment_index, this->_cutoff_radius);
        for (size_t neighbor_index : neighbor_indices) {
            exchange_field += (geometry[neighbor_index].getMaterial() == current_material
                                   ? GENERALIZED_PREFIX
                                   : GENERALIZED_PREFIX /* тут можно будет интерфейс использовать */) *
                              geometry[neighbor_index].getDirection().asVector();
        }
        return exchange_field;
    }

    virtual std::string __str__() const override {
        return std::format("CartesianExchangeInteraction(r={})", _cutoff_radius);
    };
    virtual std::string __repr__() const override {
        return std::format("CartesianExchangeInteraction(cutoff_radius={})", _cutoff_radius);
    };
};

class CartesianExternalInteraction : public CartesianAbstractInteraction {
  protected:
    Eigen::Vector3d external_field;

  public:
    CartesianExternalInteraction(const Eigen::Vector3d &external_field) : external_field(external_field) {};
    CartesianExternalInteraction(double sx, double sy, double sz) : external_field(sx, sy, sz) {};

    virtual EffectiveField calculateFieldContribution(
        size_t moment_index, IGeometry<CartesianCoordSystem> &geometry, MaterialRegistry &material_registry
    ) const override {
        return this->external_field;
    }

    virtual std::string __str__() const override {
        return std::format(
            "CartesianExternalInteraction(field=({}, {}, {}))",
            this->external_field.x(),
            this->external_field.y(),
            this->external_field.z()
        );
    };

    virtual std::string __repr__() const override {
        return std::format(
            "CartesianExternalInteraction(sx={}, sy={}, sz={})",
            this->external_field.x(),
            this->external_field.y(),
            this->external_field.z()
        );
    };
};

class CartesianAnisotropyInteraction : public CartesianAbstractInteraction {
  public:
    CartesianAnisotropyInteraction() {};

    virtual EffectiveField calculateFieldContribution(
        size_t moment_index, IGeometry<CartesianCoordSystem> &geometry, MaterialRegistry &material_registry
    ) const override {

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
};

template <CoordSystemConcept CoordSystem> using InteractionRegistry = Registry<IInteraction<CoordSystem>>;
using CartesianInteractionRegistry = InteractionRegistry<CartesianCoordSystem>;

}; // namespace spindynapy

inline void pyBindInteractions(py::module_ &module) {
    using namespace spindynapy;

    // clang-format off

    module.doc() = "Интерфейсы взаимодействий между элементами системы.\n"
                   "Таковыми могут быть потенциальные поля, силы, etc. \n";

    py::class_<CartesianAbstractInteraction, std::shared_ptr<CartesianAbstractInteraction>>(module, "CartesianAbstractInteraction")
        .def("__str__", &CartesianAbstractInteraction::__str__)
        .def("__repr__", &CartesianAbstractInteraction::__repr__)
        .def(
            "calculate_field_contribution",
            &CartesianAbstractInteraction::calculateFieldContribution,
            py::arg("moment_index"), py::arg("geometry"), py::arg("material_registry")
        )
        .doc() = "Базовый интерфейс взаимодействий (сил, полей, etc.)";

    py::class_<CartesianExchangeInteraction, CartesianAbstractInteraction, std::shared_ptr<CartesianExchangeInteraction>>(module, "CartesianExchangeInteraction")
        .def(py::init<double>(), py::arg("cutoff_radius"))
        .doc() = "Обменное взаимодействие";

    py::class_<CartesianExternalInteraction, CartesianAbstractInteraction, std::shared_ptr<CartesianExternalInteraction>>(module, "CartesianExternalInteraction")
        .def(py::init<const Eigen::Vector3d&>(), py::arg("external_field"))
        .def(py::init<double, double, double>(), py::arg("sx"), py::arg("sy"), py::arg("sz"))
        .doc() = "Взаимодействие с внешним полем";

    py::class_<CartesianAnisotropyInteraction, CartesianAbstractInteraction, std::shared_ptr<CartesianAnisotropyInteraction>>(module, "CartesianAnisotropyInteraction")
        .def(py::init())
        .doc() = "оси магнитной анизотропии";

    py::class_<CartesianInteractionRegistry, std::shared_ptr<CartesianInteractionRegistry>>(module, "CartesianInteractionRegistry")
        .def(py::init<>())
        .def(py::init<RegistryContainer<CartesianAbstractInteraction>>(), py::arg("container"))
        BUILD_REGISTRY_TEMPLATE_METHODS(CartesianInteractionRegistry)
        .doc() = "Интерфейс для регистра-контейнера для разных взаимодействий";

    // clang-format on
}

#endif // ! __INTERACTIONS_HPP__
