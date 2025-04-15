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

#include <iostream>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <stdexcept>
#include <string>
#include <type_traits>

namespace py = pybind11;

namespace spindynapy {

/**
 * Базовый интерфейс решателя-интегратора
 */

class ISolver {
  protected:
    virtual std::string getSolverName() { throw std::logic_error("Method getSolverName Not implemented"); };
    virtual void updateMomentsCartesian(IGeometry &geometry) {
        throw std::logic_error("Method updateMomentsCartesian Not implemented");
    };
    virtual void updateMomentsSpherical(IGeometry &geometry) {
        throw std::logic_error("Method updateMomentsSpherical Not implemented");
    };

  public:
    ISolver() = default;
    virtual ~ISolver() {};

    virtual std::string __str__() const { return nullptr; };
    virtual std::string __repr__() const { return nullptr; };

    virtual void updateMoments(IGeometry &geometry) {
        throw std::logic_error("Method step Not implemented");
    };
};

/**
 * SOLVER может реализовываться в разных координатных системах.
 * При построении из шаблона, через if constexpr, будет разрешаться,
 * какой метод будет вызывать солвер (для разных координатных систем = разные алгоритмы и разные методы
 * типов).
 *
 * Зато на уровне симулятора неважно, какой шаблон подставляется, если есть метод step().
 */
template <typename CoordSystem> class AbstractSolver : public ISolver {
  protected:
    virtual std::string getSolverName() override {
        throw std::logic_error("Method getSolverName Not implemented");
    };
    virtual void updateMomentsCartesian(IGeometry &geometry) override {
        throw std::logic_error(
            "The solving via " + getSolverName() + " is not implemented in the Cartesian coordinate system"
        );
    };
    virtual void updateMomentsSpherical(IGeometry &geometry) override {
        throw std::logic_error(
            "The solving via " + getSolverName() + " is not implemented in the Spherical coordinate system"
        );
    };

  public:
    AbstractSolver() = default;
    virtual ~AbstractSolver() = default;

    virtual void updateMoments(IGeometry &geometry) override { // метод разрешает ещё при препроцессинге какой
                                                               // метод вызвать для текущей реализации шаблона
        if constexpr (std::is_same_v<CoordSystem, CartesianCoordSystem>) {
            return this->updateMomentsCartesian(geometry);
        } else if constexpr (std::is_same_v<CoordSystem, SphericalCoordSystem>) {
            return this->updateMomentsSpherical(geometry);
        } else {
            throw std::logic_error("The solving is not implemented with this coords system");
        }
    };
};

// по сути - базовый абстрактный солвер, в декартовых координатах
using CartesianAbstractSolver = AbstractSolver<CartesianCoordSystem>;
// не реализован, создан для демонстрации абстракции солвера
using SphericalAbstractSolver = AbstractSolver<SphericalCoordSystem>;

template <typename CoordSystem> class LLGSolver : public AbstractSolver<CoordSystem> {
  protected: // todo: стратегия численного исчисления (эйлер, хойн и т.д.)
    virtual std::string getSolverName() override { return "LLG Solver"; };
    virtual void updateMomentsCartesian(IGeometry &geometry) override {
        std::cout << "AMOGUS\0";
        return;
    }
    virtual void updateMomentsSpherical(IGeometry &geometry) override {
        std::cout << "SPHERICAL AMOGUS\0";
        return;
    }

  public:
    LLGSolver() {};
};

// солвер Ландау-Лифшица-Гильберта в декартовых координатах (классический алгоритм)
using CartesianLLGSolver = LLGSolver<CartesianCoordSystem>;
// не реализован, создан для демонстрации абстракции солвера
using SphericalLLGSolver = LLGSolver<SphericalCoordSystem>;

}; // namespace spindynapy

inline void pyBindSolvers(py::module_ &module) {
    using namespace spindynapy;

    // clang-format off

    module.doc() = "Модуль решателей\n"
                   "Решатели (интеграторы) берут на себя задачу эволюционирования\n"
                   "или любого другого прогрессирования системы в зависимости от\n"
                   "предоставленных внешних параметров.";

    py::class_<ISolver, std::shared_ptr<ISolver>>(module, "ISolver")
        .def("__str__", &ISolver::__str__)
        .def("__repr__", &ISolver::__repr__)
        .def("update_moments", &ISolver::updateMoments, py::arg("geometry"))
        .doc() = "Базовый интерфейс решателя-интегратора";

    py::class_<CartesianAbstractSolver, ISolver, std::shared_ptr<CartesianAbstractSolver>>(module, "CartesianAbstractSolver")
        .doc() = "по сути - базовый абстрактный солвер, в декартовых координатах";

    py::class_<SphericalAbstractSolver, ISolver, std::shared_ptr<SphericalAbstractSolver>>(module, "SphericalAbstractSolver")
        .doc() = "не реализован, создан для демонстрации абстракции солвера";

    py::class_<CartesianLLGSolver, CartesianAbstractSolver, std::shared_ptr<CartesianLLGSolver>>(module, "CartesianLLGSolver")
        .def(py::init<>())
        .doc() = "солвер Ландау-Лифшица-Гильберта в декартовых координатах (классический алгоритм)";

    py::class_<SphericalLLGSolver, SphericalAbstractSolver, std::shared_ptr<SphericalLLGSolver>>(module, "SphericalLLGSolver")
        .def(py::init<>())
        .doc() = "не реализован, создан для демонстрации абстракции солвера";

    // clang-format on
}

#endif // ! __SOLVERS_HPP__
