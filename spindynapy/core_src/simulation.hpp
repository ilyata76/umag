#ifndef __SIMULATION_HPP__
#define __SIMULATION_HPP__

/**
 * Интерфейсы и классы, отвечающие за главный функционал - симуляцию,
 * управление симуляцией.
 */

#include "geometries.hpp"
#include "registries.hpp"
#include "solvers.hpp"
#include "types.hpp"

#include <memory>
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <stdexcept>
#include <string>

namespace py = pybind11;

namespace spindynapy {

/**
 * Базовый интерфейс управлятора симуляцией. Вход в программу.
 *
 * Через управлятор, содержащий в себе всю информацию о проводимом
 * эксперименте (богатое состояние), можно развивать систему согласно заданным настройкам,
 * а также извне, вызывая специальные для этого методы (API симулятора)
 */
template <CoordSystemConcept CoordSystem> class Simulation {
  protected:
    // спиновое состояние (изменяемое, подменяемое извне) в разных системах координат
    std::shared_ptr<IGeometry<CoordSystem>> _geometry;
    // полиморфный решатель в разных системах координат И с разными методами решения (стратегия расчётов)
    std::shared_ptr<ISolver<CoordSystem>> _solver;

    // регистр, который хранит мета-информацию о материалах, на которые ссылаются классы моментов
    std::shared_ptr<MaterialRegistry> _material_registry;
    // регистр, хранящий используемые рассчитываемые величины (взаимодействия между моментами)
    std::shared_ptr<InteractionRegistry<CoordSystem>> _interaction_registry;

    // буфер эффективных полей на каждый элемент геометрии (поиндексная связка)
    std::vector<EffectiveField> _effective_fields;

    // текущее время
    double _current_time = .0;

    // временной шаг
    double _dt;

    virtual void calculateEffectiveFields() {
        size_t moments_size = this->_geometry->size();

        // корректировка размера буфера (если изменилась геометрия!)
        if (this->_effective_fields.size() != moments_size) {
            this->_effective_fields.resize(moments_size, EffectiveField::Zero());
        }

        // обсчёт полей от зарегистрированных взаимодействий
        for (size_t i; i < moments_size; ++i) {
            this->_effective_fields[i] = EffectiveField::Zero();
            for (auto &[_, interaction] : *_interaction_registry) {
                this->_effective_fields[i] +=
                    interaction->calculateFieldContribution(i, *this->_geometry, *this->_material_registry);
            }
        }
    }

  public:
    /**
     * ATTENTION!
     *
     * Изменение геометрии или регистров может привести к изменению процесса симуляции.
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
        //
        if (!_geometry) throw std::invalid_argument("Геометрия не может быть None");
        if (!_solver) throw std::invalid_argument("Решатель не может быть None");
        if (!_interaction_registry || _material_registry->isEmpty())
            throw std::invalid_argument("Регистр взаимодействий не может быть None");
        if (!_material_registry || _material_registry->isEmpty())
            throw std::invalid_argument("Регистр материалов не может быть None");
        if (_dt <= 0) throw std::invalid_argument("Time step dt must be positive.");

        // Инициализируем буфер эффективных полей нужного размера
        this->_effective_fields.resize(geometry->size(), EffectiveField::Zero());
    };

    virtual std::string __str__() const { return _geometry->__str__(); };
    virtual std::string __repr__() const { return _geometry->__repr__(); };

    virtual void simulateOneStep() {
        this->calculateEffectiveFields();
        this->_solver->updateMoments(*this->_geometry, this->_effective_fields, this->_dt);
        this->_current_time += this->_dt;
    };

    virtual void simulateManySteps(uint steps) {
        for (uint i = 0; i < steps; ++i) {
            this->_solver->updateMoments(*this->_geometry, this->_effective_fields, this->_dt);
        }
    };
};

using CartesianSimulation = Simulation<CartesianCoordSystem>;

}; // namespace spindynapy

inline void pyBindSimulation(py::module_ &module) {
    using namespace spindynapy;

    // clang-format off

    module.doc() = "Интерфейсы и классы, отвечающие за главный функционал - симуляцию.";

    py::class_<CartesianSimulation>(module, "CartesianSimulation")
        .def("__str__", &CartesianSimulation::__str__)
        .def("__repr__", &CartesianSimulation::__repr__)
        .def("simulate_one_step", &CartesianSimulation::simulateOneStep)
        .def("simulate_many_steps", &CartesianSimulation::simulateManySteps, py::arg("steps"))
        .def(
            py::init<
                std::shared_ptr<CartesianAbstractGeometry>,
                std::shared_ptr<CartesianAbstractSolver>,
                std::shared_ptr<MaterialRegistry>,
                std::shared_ptr<CartesianInteractionRegistry>,
                double>(),
            py::arg("geometry"),
            py::arg("solver"),
            py::arg("material_registry"),
            py::arg("interaction_registry"),
            py::arg("dt") = 1e-13
        )
        .doc() = "TODO";

    // clang-format on
}

#endif // ! __SIMULATION_HPP__
