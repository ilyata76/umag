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
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <stdexcept>
#include <string>

namespace py = pybind11;

namespace spindynapy {

/**
 * Базовый интерфейс решателя-интегратора
 */

typedef std::vector<Eigen::Vector3d> eff_f;

template <CoordSystemConcept CoordSystem> class ISolver {
  protected:
    ISolver() = default;
  public:
    virtual ~ISolver() = default;
    virtual void updateMoments(IGeometry<CoordSystem> &geometry, eff_f effective_fields, double dt) {
        throw std::logic_error("The solving Not implemented");
    }
};

class CartesianSolver : public ISolver<CartesianCoordSystem> {
  protected:
    CartesianSolver() = default;
  public:
    virtual void
    updateMoments(IGeometry<CartesianCoordSystem> &geometry, eff_f effective_fields, double dt) override {
        throw std::logic_error("The solving Not implemented");
    };
};

class CartesianLLGSolver : public CartesianSolver {
  public:
    virtual void
    updateMoments(IGeometry<CartesianCoordSystem> &geometry, eff_f effective_fields, double dt) override {
        std::cout << "AMOGUS\n";
    };
};

}; // namespace spindynapy

inline void pyBindSolvers(py::module_ &module) {
    using namespace spindynapy;

    // clang-format off

    module.doc() = "Модуль решателей\n"
                   "Решатели (интеграторы) берут на себя задачу эволюционирования\n"
                   "или любого другого прогрессирования системы в зависимости от\n"
                   "предоставленных внешних параметров.";

    py::class_<CartesianSolver, std::shared_ptr<CartesianSolver>>(module, "CartesianAbstractSolver")
        .def("update_moments", &CartesianSolver::updateMoments, py::arg("geometry"), py::arg("effective_fields"), py::arg("dt"))
        .doc() = "по сути - базовый абстрактный солвер, в декартовых координатах";

    py::class_<CartesianLLGSolver, CartesianSolver, std::shared_ptr<CartesianLLGSolver>>(module, "CartesianLLGSolver")
        .def(py::init<>())
        .doc() = "солвер Ландау-Лифшица-Гильберта в декартовых координатах (классический алгоритм)";

    // clang-format on
}

#endif // ! __SOLVERS_HPP__
