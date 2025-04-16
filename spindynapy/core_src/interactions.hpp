#ifndef __INTERACTIONS_HPP__
#define __INTERACTIONS_HPP__

/**
 * Интерфейсы взаимодействий между элементами системы.
 * Таковыми могут быть потенциальные поля, силы, etc.
 */

#include "geometries.hpp"
#include "registries.hpp"
#include "types.hpp"

#include <format>
#include <iostream>
#include <memory>
#include <pybind11/pybind11.h>
#include <string>

namespace py = pybind11;

namespace spindynapy {

// Определим тип для эффективного поля
using EffectiveField = Eigen::Vector3d;

/**
 * Базовый интерфейс взаимодействий (сил, полей, etc.)
 */
template <CoordSystemConcept CoordSystem> class IInteraction {
  protected:
    IInteraction() = default;

  public:
    virtual ~IInteraction() = default;

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
        EffectiveField exchange_field = EffectiveField::Zero();
        auto &current_moment = geometry[moment_index];
        auto &current_material = current_moment.getMaterial();
        auto neighbor_indices = geometry.getNeighbors(moment_index, this->_cutoff_radius);
        for (size_t neighbor_index : neighbor_indices) {
            geometry[neighbor_index].getDirection();
        }
        std::cout << "CartesianExchangeInteraction" << std::endl;
        return {1, 2, 3};
    }

    virtual std::string __str__() const override {
        return std::format("CartesianExchangeInteraction(r={})", _cutoff_radius);
    };
    virtual std::string __repr__() const override {
        return std::format("CartesianExchangeInteraction(cutoff_radius={})", _cutoff_radius);
    };
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

    py::class_<CartesianInteractionRegistry, std::shared_ptr<CartesianInteractionRegistry>>(module, "CartesianInteractionRegistry")
        .def(py::init<>())
        .def(py::init<RegistryContainer<CartesianAbstractInteraction>>(), py::arg("container"))
        BUILD_REGISTRY_TEMPLATE_METHODS(CartesianInteractionRegistry)
        .doc() = "Интерфейс для регистра-контейнера для разных взаимодействий";

    // clang-format on
}

#endif // ! __INTERACTIONS_HPP__
