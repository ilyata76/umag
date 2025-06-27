#ifndef __SIMULATION_HPP__
#define __SIMULATION_HPP__

/**
 * Интерфейсы и классы, отвечающие за главный функционал - симуляцию,
 * управление симуляцией.
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

template <CoordSystemConcept CoordSystem> struct SimulationStepData : public SolverData<CoordSystem> {
  public:
    double time;
    uint step;
    std::unique_ptr<IGeometry<CoordSystem>> geometry;

    SimulationStepData() {}; // TODO не забыть про nullptr и нулевые массивы...

    // Конструктор
    SimulationStepData(
        double time,
        uint step,
        std::unique_ptr<IGeometry<CoordSystem>> geometry,
        SolverData<CoordSystem> solver_data
    )
        : SolverData<CoordSystem>(solver_data), time(time), step(step), geometry(std::move(geometry)) {};

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

    PYTHON_API double getMeanMagnetizationNorm() const { return this->getMeanMagnetization().norm(); }

    PYTHON_API double getEnergy() const {
        return std::accumulate(this->energies.begin(), this->energies.end(), 0.0);
    }

    PYTHON_API std::unordered_map<regnum, double> getEnergyByInteraction() const {
        std::unordered_map<regnum, double> energy_by_interaction;
        for (const auto &[reg, energies] : this->interaction_energies) {
            energy_by_interaction[reg] = std::accumulate(energies.begin(), energies.end(), 0.0);
        }
        return energy_by_interaction;
    }
};

namespace cartesian {

using AbstractSimulationStepData = SimulationStepData<NamespaceCoordSystem>;

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
    // полиморфный решатель в разных системах координат И с разными методами решения (стратегия
    // расчётов)
    std::shared_ptr<ISolver<CoordSystem>> _solver;

    // регистр, который хранит мета-информацию о материалах, на которые ссылаются классы моментов
    std::shared_ptr<MaterialRegistry> _material_registry;
    // регистр, хранящий используемые рассчитываемые величины (взаимодействия между моментами)
    std::shared_ptr<InteractionRegistry<CoordSystem>> _interaction_registry;

    // текущий шаг (буфер)
    SolverData<CoordSystem> _step_solver_data;

    // сохранёненые шаги (снимки)
    std::vector<SimulationStepData<CoordSystem>> steps;

    // счётчик шагов
    uint _step = 0;

    // текущее время
    double _current_time = 0.0;

    // временной шаг
    double _dt;

    // была ли предпоготовка для обсчёта вкладов взаимодействий
    bool _system_prepared = false;

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

    void saveStep() {
        this->steps.push_back(SimulationStepData<CoordSystem>(
            this->_current_time, this->_step, this->_geometry->clone(false, true), this->_step_solver_data
        ));
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

    std::string __str__() const { return _geometry->__str__(); };

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

    void simulateManySteps(uint sim_steps, uint save_every_step = 1, uint update_macrocells_every_step = 1) {
        for (uint i = 0; i < sim_steps; ++i) {
            simulateOneStep(i % save_every_step == 0, i % update_macrocells_every_step == 0);
        }
    };

    std::vector<SimulationStepData<CoordSystem>> &getSteps() { return this->steps; }
    void clearSteps() { return this->steps.clear(); }
    // TODO pop_back last step чтобы не хранить все...
};

namespace cartesian {

using AbstractSimulation = Simulation<NamespaceCoordSystem>;

};

}; // namespace PYTHON_API spindynapy

inline void pyBindSimulation(py::module_ &module) {
    using namespace spindynapy;

    // -------- | SIMULATION | --------
    py::module_ simulation_module = module.def_submodule("simulation");

    simulation_module.doc() = "Интерфейсы и классы, отвечающие за главный функционал - симуляцию.";

    // -------- | CARTESIAN SIMULATION | --------
    py::module_ cartesian = simulation_module.def_submodule("cartesian");

    using cartesian::AbstractGeometry;
    using cartesian::AbstractInteractionRegistry;
    using cartesian::AbstractSimulation;
    using cartesian::AbstractSimulationStepData;
    using cartesian::AbstractSolver;
    using cartesian::AbstractSolverData;

    py::class_<AbstractSimulationStepData, AbstractSolverData>(cartesian, "SimulationStepData")
        .def_readonly("time", &AbstractSimulationStepData::time)
        .def_readonly("step", &AbstractSimulationStepData::step)
        .def_property_readonly(
            "geometry",
            [](const AbstractSimulationStepData &self) { return self.geometry.get(); },
            py::return_value_policy::reference_internal
        )
        .def(
            "get_mean_magnetization",
            &AbstractSimulationStepData::getMeanMagnetization,
            py::doc("Get the mean magnetization of the system at this step.")
        )
        .def(
            "get_mean_magnetization_norm",
            &AbstractSimulationStepData::getMeanMagnetizationNorm,
            py::doc("Get the norm of the mean magnetization of the system at this step.")
        )
        .def(
            "get_energy",
            &AbstractSimulationStepData::getEnergy,
            py::doc("Get the total energy of the system at this step.")
        )
        .def(
            "get_energy_by_interaction",
            &AbstractSimulationStepData::getEnergyByInteraction,
            py::doc("Get the energy contributions by interaction type at this step.")
        )
        .doc() = "Data structure holding simulation step information";

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
        )
        .def("__str__", &AbstractSimulation::__str__)
        .def(
            "simulate_one_step",
            &AbstractSimulation::simulateOneStep,
            py::call_guard<py::gil_scoped_release>(),
            py::arg("save_step") = false,
            py::arg("update_macrocells") = true
        )
        .def(
            "simulate_many_steps",
            &AbstractSimulation::simulateManySteps,
            py::call_guard<py::gil_scoped_release>(),
            py::arg("steps"),
            py::arg("save_every_step") = 1,
            py::arg("update_macrocells_every_step") = 1
        )
        .def("get_steps", &AbstractSimulation::getSteps, py::return_value_policy::reference)
        .def("clear_steps", &AbstractSimulation::clearSteps)
        .doc() = "TODO";
}

#endif // ! __SIMULATION_HPP__
