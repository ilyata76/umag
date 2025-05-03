#ifndef __SOLVERS_HPP__
#define __SOLVERS_HPP__

/**
 * Решатель - это абстракция, которая берёт на себя задачу
 *  эволюционирования или любого другого прогрессирования системы в зависимости от
 *  предоставленных внешних параметров.
 *
 * Стратегия решателя (интегратора) - метод интегрирования: EULER, HEUN
 */

#include "constants.hpp"
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

namespace PYTHON_API spindynapy {

/**
 * Стратегия решателя (интегратора) - метод интегрирования.
 *
 * EULER - метод Эйлера (самый простой)
 * HEUN - метод Хьюна (метод предиктор-корректор, улучшенный метод Эйлера)
 */
enum class PYTHON_API SolverStrategy { EULER = 0, HEUN = 1 };

/**
 * Базовый интерфейс решателя (интегратора) в выбранной системе координат.
 *
 * Решатель - это абстракция, которая берёт на себя задачу
 *  эволюционирования или любого другого прогрессирования системы в зависимости от
 *  предоставленных внешних параметров.
 */
template <CoordSystemConcept CoordSystem> class PYTHON_API ISolver {
  protected:
    // конструктор только для наследников
    ISolver() = default;

  public:
    // деструктор
    virtual ~ISolver() = default;

    // обновить моменты в геометрии исходя из приложенных к ним эффективных полей
    PYTHON_API virtual void updateMoments(
        IGeometry<CoordSystem> &geometry, std::vector<EffectiveField> effective_fields, double dt
    ) = 0;
};

namespace cartesian {

/**
 * Базовый интерфейс решателя (интегратора) в выбранной (ДЕКАРТОВОЙ) системе координат.
 *
 * Решатель - это абстракция, которая берёт на себя задачу
 *  эволюционирования или любого другого прогрессирования системы в зависимости от
 *  предоставленных внешних параметров.
 */
using AbstractSolver = PYTHON_API ISolver<NamespaceCoordSystem>;

/**
 * Решатель ЛЛГ (Ландау-Лифшица-Гильберта) - это решатель, который использует
 *  уравнение ЛЛГ спиновой динамики для обновления моментов в геометрии.
 *
 * Решатель использует метод Эйлера или метод Хойна для интегрирования уравнения ЛЛГ с шагом dt.
 */
class LLGSolver : public AbstractSolver {
  protected:
    // выбранная стратегия решателя (интегратора) - метод интегрирования
    SolverStrategy _strategy;

    // вычисление изменения момента по уравнению ЛЛГ (правая часть уравнения) (TODO: команда + стратегия)
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

    // решение по методу Эйлера: сделать шаг с конечными разностями
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

    // решение по методу Хойна: сделать шаг с конечными разностями, но дважды (+корректор с подсчитанными
    // полями)
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
    // базовый конструктор
    LLGSolver(SolverStrategy strategy = SolverStrategy::EULER) : _strategy(strategy) {};

    // обновить моменты в геометрии исходя из приложенных к ним эффективных полей
    virtual void updateMoments(
        IGeometry<NamespaceCoordSystem> &geometry, std::vector<EffectiveField> effective_fields, double dt
    ) override {
        // geometry и effective_fields связаны через индексную привязку effective_fields[i] ON geometry[i]
        assert(geometry.size() == effective_fields.size());

        if (this->_strategy == SolverStrategy::EULER) {
            return this->updateMomentsViaEuler(geometry, effective_fields, dt);
        } else if (this->_strategy == SolverStrategy::HEUN) {
            return this->updateMomentsViaHeun(geometry, effective_fields, dt);
        };

        throw std::invalid_argument("Invalid solver strategy enum value / doesnt support");
    };
};

}; // namespace cartesian

}; // namespace PYTHON_API spindynapy

#define ABSTRACTSOLVER_TEMPLATE_BINDINGS(cls)                                                                \
    .def(                                                                                                    \
        "update_moments",                                                                                    \
        &cls::updateMoments,                                                                                 \
        py::arg("geometry"),                                                                                 \
        py::arg("effective_fields"),                                                                         \
        py::arg("dt"),                                                                                       \
        py::doc("обновить моменты в геометрии исходя из приложенных к ним эффективных полей")                \
    )

inline void pyBindSolvers(py::module_ &module) {
    using namespace spindynapy;

    // -------- | SOLVERS | --------
    py::module_ solvers_module = module.def_submodule("solvers");

    solvers_module.doc() = "Решатель - это абстракция, которая берёт на себя задачу\n"
                           " эволюционирования или любого другого прогрессирования системы в зависимости от\n"
                           " предоставленных внешних параметров.\n"
                           "\n"
                           "Стратегия решателя (интегратора) - метод интегрирования: EULER, HEUN, ...";

    // Биндим enum
    py::enum_<SolverStrategy>(solvers_module, "SolverStrategy")
        .value("EULER", SolverStrategy::EULER, "Метод Эйлера (самый простой)")
        .value(
            "HEUN",
            SolverStrategy::HEUN,
            "Метод Хойна (предиктор-корректор (2 расчёта), улучшенный метод Эйлера)"
        );

    // -------- | CARTESIAN SOLVERS | --------
    py::module_ cartesian = solvers_module.def_submodule("cartesian");

    cartesian.doc() = "Решатель - это абстракция, которая берёт на себя задачу\n"
                      " эволюционирования или любого другого прогрессирования системы в зависимости от\n"
                      " предоставленных внешних параметров.\n"
                      "\n"
                      "Стратегия решателя (интегратора) - метод интегрирования: EULER, HEUN, ...";

    using cartesian::AbstractSolver;
    using cartesian::LLGSolver;

    py::class_<AbstractSolver, std::shared_ptr<AbstractSolver>>(cartesian, "AbstractSolver")
        ABSTRACTSOLVER_TEMPLATE_BINDINGS(AbstractSolver)
            .doc() = "Базовый интерфейс решателя (интегратора) в выбранной (ДЕКАРТОВОЙ) системе координат.\n"
                     "  Решатель - это абстракция, которая берёт на себя задачу\n"
                     "  эволюционирования или любого другого прогрессирования системы в зависимости от\n"
                     "  предоставленных внешних параметров.";

    py::class_<LLGSolver, AbstractSolver, std::shared_ptr<LLGSolver>>(cartesian, "LLGSolver")
        .def(py::init<SolverStrategy>(), py::arg("strategy") = SolverStrategy::EULER)
        .doc() =
        "Решатель ЛЛГ (Ландау-Лифшица-Гильберта) - это решатель, который использует\n"
        "  уравнение ЛЛГ спиновой динамики для обновления моментов в геометрии.\n"
        "\n"
        "  Решатель использует метод Эйлера или метод Хойна для интегрирования уравнения ЛЛГ с шагом dt.";
}

#endif // ! __SOLVERS_HPP__
