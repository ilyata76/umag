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
#include "interactions.hpp"
#include "types.hpp"

#include <iostream>
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <stdexcept>

namespace py = pybind11;

namespace spindynapy {

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

using CartesianAbstractSolver = ISolver<CartesianCoordSystem>;

class CartesianLLGSolver : public CartesianAbstractSolver {
  public:
    CartesianLLGSolver() = default;

    virtual void updateMoments(
        IGeometry<CartesianCoordSystem> &geometry, std::vector<EffectiveField> effective_fields, double dt
    ) override {
        size_t num_moments = geometry.size();
        std::cout << "\n\nUPDATE MOMENTS\n\n";

        // проверка на изменение геометрии
        if (effective_fields.size() != num_moments) {
            throw std::runtime_error("Mismatch between geometry size and effective fields size in LLG solver."
            );
        }
        // для каждого момента обсчёт уравнения
        for (size_t i = 0; i < num_moments; ++i) {
            //
        }
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

    py::class_<CartesianAbstractSolver, std::shared_ptr<CartesianAbstractSolver>>(module, "CartesianAbstractSolver")
        .def("update_moments", &CartesianAbstractSolver::updateMoments, py::arg("geometry"), py::arg("effective_fields"), py::arg("dt"))
        .doc() = "по сути - базовый абстрактный солвер, в декартовых координатах";

    py::class_<CartesianLLGSolver, CartesianAbstractSolver, std::shared_ptr<CartesianLLGSolver>>(module, "CartesianLLGSolver")
        .def(py::init<>())
        .doc() = "солвер Ландау-Лифшица-Гильберта в декартовых координатах (классический алгоритм)";

    // clang-format on
}

#endif // ! __SOLVERS_HPP__
