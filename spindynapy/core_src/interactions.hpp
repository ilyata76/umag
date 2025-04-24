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

#include <cmath>
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

    virtual void
    prepareData(size_t moment_index, IGeometry<CoordSystem> &geometry, MaterialRegistry &material_registry) {
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

template <CoordSystemConcept CoordSystem> using InteractionRegistry = Registry<IInteraction<CoordSystem>>;

namespace cartesian {

using AbstractInteraction = IInteraction<NamespaceCoordSystem>;

class ExchangeInteraction : public AbstractInteraction {
  protected:
    double _cutoff_radius;

  public:
    ExchangeInteraction(double cutoff_radius) : _cutoff_radius(cutoff_radius) {};

    virtual void prepareData(
        size_t moment_index, IGeometry<NamespaceCoordSystem> &geometry, MaterialRegistry &material_registry
    ) override {
        // обновить кэш по соседям
        geometry.getNeighbors(moment_index, this->_cutoff_radius);
    }

    virtual EffectiveField calculateFieldContribution(
        size_t moment_index, IGeometry<NamespaceCoordSystem> &geometry, MaterialRegistry &material_registry
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
        return std::format("ExchangeInteraction(r={})", _cutoff_radius);
    };
    virtual std::string __repr__() const override {
        return std::format("ExchangeInteraction(cutoff_radius={})", _cutoff_radius);
    };
};

class ExternalInteraction : public AbstractInteraction {
  protected:
    Eigen::Vector3d external_field;

  public:
    ExternalInteraction(const Eigen::Vector3d &external_field) : external_field(external_field) {};
    ExternalInteraction(double sx, double sy, double sz) : external_field(sx, sy, sz) {};

    virtual void prepareData(
        size_t moment_index, IGeometry<NamespaceCoordSystem> &geometry, MaterialRegistry &material_registry
    ) override {
        return;
    }

    virtual EffectiveField calculateFieldContribution(
        size_t moment_index, IGeometry<NamespaceCoordSystem> &geometry, MaterialRegistry &material_registry
    ) const override {
        return this->external_field;
    }

    virtual std::string __str__() const override {
        return std::format(
            "ExternalInteraction(field=({}, {}, {}))",
            this->external_field.x(),
            this->external_field.y(),
            this->external_field.z()
        );
    };

    virtual std::string __repr__() const override {
        return std::format(
            "ExternalInteraction(sx={}, sy={}, sz={})",
            this->external_field.x(),
            this->external_field.y(),
            this->external_field.z()
        );
    };
};

class AnisotropyInteraction : public AbstractInteraction {
  public:
    AnisotropyInteraction() {};

    virtual void prepareData(
        size_t moment_index, IGeometry<NamespaceCoordSystem> &geometry, MaterialRegistry &material_registry
    ) override {
        return;
    }

    virtual EffectiveField calculateFieldContribution(
        size_t moment_index, IGeometry<NamespaceCoordSystem> &geometry, MaterialRegistry &material_registry
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

class DemagnetizationInteraction : public AbstractInteraction {
  protected:
    std::string _strategy;
    double _cutoff_radius;

  public:
    DemagnetizationInteraction(double cutoff_radius, std::string strategy = "cutoff")
        : _cutoff_radius(cutoff_radius), _strategy(strategy) {
        if (strategy != "cutoff") {
            throw std::invalid_argument("Invalid strategy string");
        }
    };

    virtual void prepareData(
        size_t moment_index, IGeometry<NamespaceCoordSystem> &geometry, MaterialRegistry &material_registry
    ) override {
        // if cutoff
        // обновить кэш по соседям
        geometry.getNeighbors(moment_index, this->_cutoff_radius);
    }

    virtual EffectiveField calculateFieldContribution(
        size_t moment_index, IGeometry<NamespaceCoordSystem> &geometry, MaterialRegistry &material_registry
    ) const override {
        // if cutoff
        EffectiveField demagnetization_field = EffectiveField::Zero();
        auto &current_moment = geometry[moment_index];
        auto &current_material = current_moment.getMaterial();
        auto neighbor_indices = geometry.getNeighbors(moment_index, _cutoff_radius);
        // auto neighbor_moments = geometry.getFromIndexes(neighbor_indices)

        auto atomic_magnetic_moments_norm =
            current_material.atomic_magnetic_saturation_magnetization * constants::BOHR_MAGNETON;

        for (size_t neighbor_index : neighbor_indices) {
            auto &neighbor_moment = geometry[neighbor_index];
            auto distance_vector =
                current_moment.getCoordinates().asVector() - neighbor_moment.getCoordinates().asVector();
            auto neighbor_atomistic_moment =
                neighbor_moment.getMaterial().atomic_magnetic_saturation_magnetization *
                constants::BOHR_MAGNETON;
            auto distance_norm = distance_vector.norm();

            demagnetization_field -=
                neighbor_atomistic_moment *
                (neighbor_moment.getDirection().asVector() / pow(distance_norm, 3) -
                 (3 * neighbor_moment.getDirection().asVector().dot(distance_vector) * distance_vector) /
                     pow(distance_norm, 5));
        }
        demagnetization_field =
            (constants::VACUUM_MAGNETIC_PERMEABILITY / 4 / constants::NUMBER_PI) * demagnetization_field;

        return demagnetization_field;
    }
};

using InteractionRegistry = InteractionRegistry<NamespaceCoordSystem>;

}; // namespace cartesian

}; // namespace spindynapy

inline void pyBindInteractions(py::module_ &module) {
    using namespace spindynapy;

    // -------- | INTERACTIONS | --------
    py::module_ interaction_module = module.def_submodule("interactions");

    interaction_module.doc() = "Интерфейсы взаимодействий между элементами системы.\n"
                               "Таковыми могут быть потенциальные поля, силы, etc. \n";

    // -------- | CARTESIAN INTERACTIONS | --------
    py::module_ cartesian = interaction_module.def_submodule("cartesian");

    using cartesian::AbstractInteraction;
    using cartesian::AnisotropyInteraction;
    using cartesian::DemagnetizationInteraction;
    using cartesian::ExchangeInteraction;
    using cartesian::ExternalInteraction;
    using cartesian::InteractionRegistry;

    py::class_<AbstractInteraction, std::shared_ptr<AbstractInteraction>>(cartesian, "AbstractInteraction")
        .def("__str__", &AbstractInteraction::__str__)
        .def("__repr__", &AbstractInteraction::__repr__)
        .def(
            "calculate_field_contribution",
            &AbstractInteraction::calculateFieldContribution,
            py::arg("moment_index"),
            py::arg("geometry"),
            py::arg("material_registry")
        )
        .doc() = "Базовый интерфейс взаимодействий (сил, полей, etc.)";

    py::class_<ExchangeInteraction, AbstractInteraction, std::shared_ptr<ExchangeInteraction>>(
        cartesian, "ExchangeInteraction"
    )
        .def(py::init<double>(), py::arg("cutoff_radius"))
        .doc() = "Обменное взаимодействие";

    py::class_<ExternalInteraction, AbstractInteraction, std::shared_ptr<ExternalInteraction>>(
        cartesian, "ExternalInteraction"
    )
        .def(py::init<const Eigen::Vector3d &>(), py::arg("external_field"))
        .def(py::init<double, double, double>(), py::arg("sx"), py::arg("sy"), py::arg("sz"))
        .doc() = "Взаимодействие с внешним полем";

    py::class_<AnisotropyInteraction, AbstractInteraction, std::shared_ptr<AnisotropyInteraction>>(
        cartesian, "AnisotropyInteraction"
    )
        .def(py::init())
        .doc() = "оси магнитной анизотропии";

    py::class_<DemagnetizationInteraction, AbstractInteraction, std::shared_ptr<DemagnetizationInteraction>>(
        cartesian, "DemagnetizationInteraction"
    )
        .def(py::init<double, std::string>(), py::arg("cutoff_radius"), py::arg("strategy") = "cutoff")
        .doc() = "demag";

    py::class_<InteractionRegistry, std::shared_ptr<InteractionRegistry>>(cartesian, "InteractionRegistry")
        .def(py::init<>())
        .def(py::init<RegistryContainer<AbstractInteraction>>(), py::arg("container"))
            BUILD_REGISTRY_TEMPLATE_METHODS(InteractionRegistry)
        .doc() = "Интерфейс для регистра-контейнера для разных взаимодействий";
}

#endif // ! __INTERACTIONS_HPP__
