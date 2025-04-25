#ifndef __SOLVERS_HPP__
#define __SOLVERS_HPP__

/**
 * Решатели (интеграторы) берут на себя задачу эволюционирования
 * или любого другого прогрессирования системы в зависимости от
 * предоставленных внешних параметров.
 * Абстракция решателя - именно метод. Фиксируется внутри системы координат (коих немного).
 */

// полиморфный решатель: шаблонный метод + IF PREPROCESSOR ELSE.
// решатель может реализовать сколь угодно много подвидов своих алгоритмов
// через функции, вызов которых разрешается при раскрытии <темплейта CoordSystem>

#include "geometries.hpp"
#include "types.hpp"

#include <Eigen/Dense>
#include <cassert>
#include <cmath>
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <stdexcept>

namespace py = pybind11;

namespace spindynapy {

enum class SolverStrategy { EULER = 0, HEUN = 1 };

/**
 * Базовый интерфейс решателя-интегратора
 */

template <CoordSystemConcept CoordSystem> class ISolver {
  protected:
    ISolver() = default;

  public:
    virtual ~ISolver() = default;

    virtual void
    updateMoments(IGeometry<CoordSystem> &geometry, std::vector<EffectiveField> effective_fields, double dt) {
        throw std::logic_error("The solving Not implemented");
    }
};

namespace cartesian {

using AbstractSolver = ISolver<NamespaceCoordSystem>;

class LLGSolver : public AbstractSolver {
  protected:
    SolverStrategy _strategy;

    Eigen::Vector3d calculateLLGMomentChange(
        Material &material, Eigen::Vector3d moment_vector, Eigen::Vector3d effective_field
    ) {
        auto prefix_term = material.gyromagnetic_ratio / (1.0 + pow(material.damping_constant, 2));

        Eigen::Vector3d moment_change =
            -(prefix_term * moment_vector.cross(effective_field) +
              prefix_term * material.damping_constant *
                  moment_vector.cross(moment_vector.cross(effective_field)));

        return moment_change;
    }

    virtual void updateMomentsViaEuler(
        IGeometry<NamespaceCoordSystem> &geometry, std::vector<EffectiveField> effective_fields, double dt
    ) {
        auto moments_size = geometry.size();
        for (size_t i = 0; i < moments_size; ++i) {
            auto &moment = geometry[i];
            auto &base_moment_vector = moment.getDirection().asVector(); // нормализованный

            auto moment_change =
                this->calculateLLGMomentChange(moment.getMaterial(), base_moment_vector, effective_fields[i]);

            // МЕТОД ЭЙЛЕРА
            Eigen::Vector3d new_moment = base_moment_vector + moment_change * dt;
            moment.setDirection(new_moment); // normalize вызывается внутри
        }
    }

    virtual void updateMomentsViaHeun(
        IGeometry<NamespaceCoordSystem> &geometry, std::vector<EffectiveField> effective_fields, double dt
    ) {
        auto moments_size = geometry.size();
        for (size_t i = 0; i < moments_size; ++i) {
            auto &moment = geometry[i];
            auto &base_moment_vector = moment.getDirection().asVector();

            // 1. Предиктор
            auto moment_change_1 =
                this->calculateLLGMomentChange(moment.getMaterial(), base_moment_vector, effective_fields[i]);
            Eigen::Vector3d moment_predictor = base_moment_vector + moment_change_1 * dt;
            moment_predictor.normalize(); // надо руками

            // 2. Корректор
            auto moment_change_2 =
                this->calculateLLGMomentChange(moment.getMaterial(), moment_predictor, effective_fields[i]);

            // ТЕПЕРЬ БЕРЁМ СРЕДНЕЕ
            Eigen::Vector3d new_moment = base_moment_vector + (moment_change_1 + moment_change_2) * dt / 2;
            // меняем момент
            moment.setDirection(new_moment); // normalize вызывается внутри
        }
    }

  public:
    LLGSolver(SolverStrategy strategy = SolverStrategy::EULER) : _strategy(strategy) {};

    virtual void updateMoments(
        IGeometry<NamespaceCoordSystem> &geometry, std::vector<EffectiveField> effective_fields, double dt
    ) override {
        // geometry и effective_fields связаны через индексную привязку
        // effective_fields[i] AT geometry[i]
        assert(geometry.size() == effective_fields.size());

        if (this->_strategy == SolverStrategy::EULER) {
            return this->updateMomentsViaEuler(geometry, effective_fields, dt);
        } else if (this->_strategy == SolverStrategy::HEUN) {
            return this->updateMomentsViaHeun(geometry, effective_fields, dt);
        };
    };
};

}; // namespace cartesian

}; // namespace spindynapy

inline void pyBindSolvers(py::module_ &module) {
    using namespace spindynapy;

    // -------- | SOLVERS | --------
    py::module_ solvers_module = module.def_submodule("solvers");

    solvers_module.doc() = "Модуль решателей\n"
                           "Решатели (интеграторы) берут на себя задачу эволюционирования\n"
                           "или любого другого прогрессирования системы в зависимости от\n"
                           "предоставленных внешних параметров.";

    // Биндим enum
    py::enum_<SolverStrategy>(solvers_module, "SolverStrategy")
        .value("EULER", SolverStrategy::EULER)
        .value("HEUN", SolverStrategy::HEUN);

    // -------- | CARTESIAN SOLVERS | --------
    py::module_ cartesian = solvers_module.def_submodule("cartesian");

    using cartesian::AbstractSolver;
    using cartesian::LLGSolver;

    py::class_<AbstractSolver, std::shared_ptr<AbstractSolver>>(cartesian, "AbstractSolver")
        .def(
            "update_moments",
            &AbstractSolver::updateMoments,
            py::arg("geometry"),
            py::arg("effective_fields"),
            py::arg("dt")
        )
        .doc() = "по сути - базовый абстрактный солвер, в декартовых координатах";

    py::class_<LLGSolver, AbstractSolver, std::shared_ptr<LLGSolver>>(cartesian, "LLGSolver")
        .def(py::init<SolverStrategy>(), py::arg("strategy") = SolverStrategy::EULER)
        .doc() = "солвер Ландау-Лифшица-Гильберта в декартовых координатах (классический алгоритм)";
}

#endif // ! __SOLVERS_HPP__
