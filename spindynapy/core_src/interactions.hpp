#ifndef __INTERACTIONS_HPP__
#define __INTERACTIONS_HPP__

/**
 * Интерфейсы взаимодействий между элементами системы.
 * Таковыми могут быть потенциальные поля, силы, etc.
 */

#include "geometries.hpp"
#include "registries.hpp"
#include "types.hpp"

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
        size_t moment_index, const IGeometry<CoordSystem> &geometry, const MaterialRegistry &material_registry
    ) const {
        throw std::logic_error("Method calculateFieldContribution Not implemented");
    };

    virtual std::string __str__() const { return nullptr; };
    virtual std::string __repr__() const { return nullptr; };
};

/**
 * Обменное взаимодействие
 */

using CartesianInteraction = IInteraction<CartesianCoordSystem>;

class CartesianExchangeInteraction : public CartesianInteraction {
  public:
    CartesianExchangeInteraction() = default;

    virtual EffectiveField calculateFieldContribution(
        size_t moment_index,
        const IGeometry<CartesianCoordSystem> &geometry,
        const MaterialRegistry &material_registry
    ) const override {
        throw std::logic_error("Method calculateFieldContribution Not implemented");
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

    py::class_<CartesianInteraction, std::shared_ptr<CartesianInteraction>>(module, "CartesianInteraction")
        .def("__str__", &CartesianInteraction::__str__)
        .def("__repr__", &CartesianInteraction::__repr__)
        .def(
            "calculate_field_contribution",
            &CartesianInteraction::calculateFieldContribution,
            py::arg("moment_index"), py::arg("geometry"), py::arg("material_registry")
        )
        .doc() = "Базовый интерфейс взаимодействий (сил, полей, etc.)";

    py::class_<CartesianExchangeInteraction, CartesianInteraction, std::shared_ptr<CartesianExchangeInteraction>>(module, "CartesianExchangeInteraction")
        .def(py::init<>())
        .doc() = "Обменное взаимодействие";

    py::class_<CartesianInteractionRegistry, std::shared_ptr<CartesianInteractionRegistry>>(module, "CartesianInteractionRegistry")
        .def(py::init<>())
        .def(py::init<RegistryContainer<CartesianInteraction>>(), py::arg("container"))
        BUILD_REGISTRY_TEMPLATE_METHODS(CartesianInteractionRegistry)
        .doc() = "Интерфейс для регистра-контейнера для разных взаимодействий";

    // clang-format on
}

#endif // ! __INTERACTIONS_HPP__
