#ifndef __SIMULATION_HPP__
#define __SIMULATION_HPP__

/**
 * Интерфейсы и классы, отвечающие за главный функционал - симуляцию,
 * управление симуляцией.
 */

#include "geometries.hpp"
#include "interactions.hpp"
#include "registries.hpp"
#include "solvers.hpp"
#include "types.hpp"

#include <algorithm>
#include <format>
#include <memory>
#include <pybind11/detail/common.h>
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <stdexcept>
#include <string>

namespace py = pybind11;

namespace spindynapy {

using EffectiveFieldVector = std::vector<EffectiveField>;

template <CoordSystemConcept CoordSystem> struct SimulationStepData {
    double time;
    std::unique_ptr<IGeometry<CoordSystem>> geometry;
    // буфер эффективных полей на каждый элемент геометрии (поиндексная связка)
    EffectiveFieldVector total_fields;
    // STEP BUFFER <INTERACTION REGISTRY NUMBER, FIELD MASSIVE> связки через ИНДЕКСЫ
    // т.е. inter_fields[0] - список конкретного поля по индексам.
    // inter_fields[0][100] - конкретное поле на 100 атоме
    std::unordered_map<regnum, EffectiveFieldVector> interaction_fields;

    // Конструктор
    SimulationStepData(
        double time,
        std::unique_ptr<IGeometry<CoordSystem>> geometry,
        const EffectiveFieldVector &total_fields,
        const std::unordered_map<regnum, EffectiveFieldVector> &interaction_fields
    )
        : time(time),
          geometry(std::move(geometry)),
          total_fields(total_fields),
          interaction_fields(interaction_fields) {}

    std::string __str__() const { return this->asString(false); };

    std::string asString(bool use_descriptions = true) const {
        std::string result = std::format(
            "# Simulation Step Data for time {:.5e}. \n"
            "# <i> <...Geometry[i]> <total_field value> <...interaction_fields values> \n",
            time,
            geometry->size()
        );
        for (size_t i = 0; i < geometry->size(); ++i) {
            result +=
                use_descriptions
                    ? std::format(
                          "[{}]\tmoment[x, y, z, sx, sy, sz, material]: {},\ttotal_field_norm: {:.5e} x: "
                          "{:.5e} y: {:.5e} z: {:.5e}",
                          i,
                          (*geometry)[i].__str__(),
                          total_fields[i].norm(),
                          total_fields[i].x(),
                          total_fields[i].y(),
                          total_fields[i].z()
                      )
                    : std::format("{}\t{}\t{:.5e}", i, (*geometry)[i].__str__(), total_fields[i].norm());
            for (const auto &[number, field_vector] : interaction_fields) {
                result += use_descriptions
                              ? std::format(
                                    ", interaction[{}]_norm: {:.5e} x: {:.5e} y: {:.5e} z: {:.5e}",
                                    number,
                                    field_vector[i].norm(),
                                    field_vector[i].x(),
                                    field_vector[i].y(),
                                    field_vector[i].z()
                                )
                              : std::format(", {}: {:.5e}", number, field_vector[i].norm());
            }
            result += "\n";
        }
        return result;
    }
};

namespace cartesian {

using SimulationStepData = SimulationStepData<NamespaceCoordSystem>;

};

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
    EffectiveFieldVector _effective_fields;
    // STEP BUFFER <INTERACTION REGISTRY NUMBER, FIELD MASSIVE>
    std::unordered_map<regnum, EffectiveFieldVector> _interaction_effective_fields;

    std::vector<SimulationStepData<CoordSystem>> steps;

    // текущее время
    double _current_time = .0;

    // временной шаг
    double _dt;

    void correctBuffers(size_t moments_size) {
        // корректировка размера буфера (если изменилась геометрия!)
        if (this->_effective_fields.size() != moments_size) {
            this->_effective_fields.resize(moments_size, EffectiveField::Zero());

            for (auto &[interaction_regnum, interaction] : *_interaction_registry) {
                this->_interaction_effective_fields[interaction_regnum].resize(
                    moments_size, EffectiveField::Zero()
                );
            }
        }
    }

    void clearBuffers(size_t moments_size) {
        this->_effective_fields = EffectiveFieldVector(moments_size, EffectiveField::Zero());
        for (auto &[interaction_regnum, interaction] : *_interaction_registry) {
            this->_interaction_effective_fields[interaction_regnum] =
                EffectiveFieldVector(moments_size, EffectiveField::Zero());
        }
    };

    void calculateEffectiveFields() {
        size_t moments_size = this->_geometry->size();
        // если изменилась геометрия
        this->correctBuffers(moments_size);
        // подготовка (чистка буфферов)
        this->clearBuffers(moments_size);

        // обсчёт полей от зарегистрированных взаимодействий
        for (auto &[interaction_regnum, interaction] : *_interaction_registry) {
            EffectiveField calc_field;
            for (size_t i = 0; i < moments_size; ++i) {
                calc_field =
                    interaction->calculateFieldContribution(i, *this->_geometry, *this->_material_registry);
                this->_effective_fields[i] += calc_field; // это суммарное поле
                this->_interaction_effective_fields[interaction_regnum][i] =
                    calc_field; // это - единственное поле
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

    std::string __str__() const { return _geometry->__str__(); };
    std::string __repr__() const { return _geometry->__repr__(); };

    void simulateOneStep(bool save_step = false) {
        this->_current_time += this->_dt;
        this->calculateEffectiveFields();
        this->_solver->updateMoments(*this->_geometry, this->_effective_fields, this->_dt);
        this->steps.push_back(SimulationStepData<CoordSystem>(
            this->_current_time,
            this->_geometry->clone(),
            this->_effective_fields,
            this->_interaction_effective_fields
        ));
    };

    void simulateManySteps(uint steps, uint save_every_step = 1) {
        for (uint i = 0; i < steps; ++i) {
            simulateOneStep(i % save_every_step == 0);
        }
    };

    std::vector<SimulationStepData<CoordSystem>> &getSteps() { return this->steps; }
    void clearSteps() { return this->steps.clear(); }
};

namespace cartesian {

using Simulation = Simulation<NamespaceCoordSystem>;

};

}; // namespace spindynapy

inline void pyBindSimulation(py::module_ &module) {
    using namespace spindynapy;

    // -------- | SIMULATION | --------
    py::module_ simulation_module = module.def_submodule("simulation");

    simulation_module.doc() = "Интерфейсы и классы, отвечающие за главный функционал - симуляцию.";

    // -------- | CARTESIAN SIMULATION | --------
    py::module_ cartesian = simulation_module.def_submodule("cartesian");

    using cartesian::AbstractGeometry;
    using cartesian::AbstractSolver;
    using cartesian::InteractionRegistry;
    using cartesian::Simulation;
    using cartesian::SimulationStepData;

    py::class_<SimulationStepData>(cartesian, "SimulationStepData")
        .def_readonly("time", &SimulationStepData::time)
        .def_property_readonly(
            "geometry",
            [](const SimulationStepData &self) { return self.geometry.get(); },
            py::return_value_policy::reference_internal
        )
        .def("__str__", &SimulationStepData::__str__)
        .def("as_string", &SimulationStepData::asString, py::arg("use_descriptions") = true)
        .def_readonly("total_fields", &SimulationStepData::total_fields)
        .def_readonly("interaction_fields", &SimulationStepData::interaction_fields)
        .doc() = "Data structure holding simulation step information";

    py::class_<Simulation>(cartesian, "Simulation")
        .def("__str__", &Simulation::__str__)
        .def("__repr__", &Simulation::__repr__)
        .def("simulate_one_step", &Simulation::simulateOneStep, py::arg("save_step") = false)
        .def(
            "simulate_many_steps",
            &Simulation::simulateManySteps,
            py::arg("steps"),
            py::arg("save_every_step") = 1
        )
        .def("get_steps", &Simulation::getSteps, py::return_value_policy::reference)
        .def("clear_steps", &Simulation::clearSteps)
        .def(
            py::init<
                std::shared_ptr<AbstractGeometry>,
                std::shared_ptr<AbstractSolver>,
                std::shared_ptr<MaterialRegistry>,
                std::shared_ptr<InteractionRegistry>,
                double>(),
            py::arg("geometry"),
            py::arg("solver"),
            py::arg("material_registry"),
            py::arg("interaction_registry"),
            py::arg("dt") = 1e-13
        )
        .doc() = "TODO";
}

#endif // ! __SIMULATION_HPP__
