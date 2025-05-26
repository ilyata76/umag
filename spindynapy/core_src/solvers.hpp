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
#include "interactions.hpp"
#include "registries.hpp"
#include "types.hpp"

#include <Eigen/Dense>
#include <cassert>
#include <cmath>
#include <memory>
#ifdef _OPENMP // если OpenMP доступен (флаг компиляции, заголовки)
#include <omp.h>
#endif
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
 * MIDPOINT - метод Эйлера (частный случай Хойна без расчёта полей)
 * HEUN - метод Хойна (метод предиктор-корректор)
 *
 */
enum class PYTHON_API SolverStrategy { EULER = 0, HEUN = 1, MIDPOINT = 2 };

/**
 * Структура для организации данных, которые возвращает и использует решатель (интегратор).
 *
 * Содержит эффективные поля и энергии, в т.ч. по каждому из взаимодействий.
 * Может использоваться как буфер для хранения промежуточных данных.
 */
template <CoordSystemConcept CoordSystem> struct PYTHON_API SolverData {
  public:
    // массив эффективных полей на каждый спин [index]
    PYTHON_API EffectiveFieldVector effective_fields;
    // массив энергий по каждому из эффективных на каждый спин [index]
    PYTHON_API std::vector<double> energies;
    // карта взаимодействий и массива его эффективных полей на каждый спин [index]
    PYTHON_API std::unordered_map<regnum, EffectiveFieldVector> interaction_effective_fields;
    // карта взаимодействий и массива его потенциалов-энергий на каждый спин [index]
    PYTHON_API std::unordered_map<regnum, std::vector<double>> interaction_energies;

    // Конструктор по умолчанию
    SolverData() {}; // TODO не забыть про nullptr и нулевые массивы...

    // Конструктор
    SolverData(
        const EffectiveFieldVector &effective_fields,
        std::vector<double> energies,
        const std::unordered_map<regnum, EffectiveFieldVector> &interaction_effective_fields,
        const std::unordered_map<regnum, std::vector<double>> &interaction_energies
    )
        : effective_fields(effective_fields),
          energies(energies),
          interaction_effective_fields(interaction_effective_fields),
          interaction_energies(interaction_energies) {}

    // Очистить и занулить состояние на основе регистра взаимодействий
    PYTHON_API void clear(size_t moments_size, InteractionRegistry<CoordSystem> &interaction_registry) {
        this->effective_fields = EffectiveFieldVector(moments_size, EffectiveField::Zero());
        this->energies = std::vector<double>(moments_size, 0.0);
        this->interaction_effective_fields.clear();
        this->interaction_energies.clear();
        for (auto &[interaction_regnum, interaction] : interaction_registry) {
            this->interaction_effective_fields[interaction_regnum] =
                EffectiveFieldVector(moments_size, EffectiveField::Zero());
            this->interaction_energies[interaction_regnum] = std::vector<double>(moments_size, 0.0);
        }
    };

    // В случае, если изменилась геометрия, то нужно обновить размеры массивов
    PYTHON_API void correct(size_t moments_size, InteractionRegistry<CoordSystem> &interaction_registry) {
        if (moments_size == this->effective_fields.size())
            return;
        this->effective_fields.clear();
        this->effective_fields.resize(moments_size);
        this->energies.clear();
        this->energies.resize(moments_size);
        this->interaction_effective_fields.clear();
        this->interaction_energies.clear();
        for (auto &[interaction_regnum, interaction] : interaction_registry) {
            this->interaction_effective_fields[interaction_regnum].resize(moments_size);
            this->interaction_energies[interaction_regnum].resize(moments_size);
        }
    }
};

/**
 * Базовый интерфейс решателя (интегратора) в выбранной системе координат.
 *
 * Решатель - это абстракция, которая берёт на себя задачу
 *  эволюционирования или любого другого прогрессирования системы в зависимости от
 *  предоставленных внешних параметров.
 *
 * Типичный решатель использует стратегию интегрирования (метод Эйлера, Хьюна и т.д.).
 * А также рассчитывает самостоятельно эффективные поля и энергии взаимодействий
 *     на основе предоставленной ему стратегии расчёта (IFieldUpdater).
 */
template <CoordSystemConcept CoordSystem> class PYTHON_API ISolver {
  protected:
    // конструктор только для наследников
    ISolver() = default;

  public:
    // деструктор
    virtual ~ISolver() = default;

    // обновить моменты в геометрии исходя из приложенных к ним эффективных полей
    PYTHON_API virtual SolverData<CoordSystem> updateMoments(
        IGeometry<CoordSystem> &geometry,
        InteractionRegistry<CoordSystem> &interaction_registry,
        MaterialRegistry &material_registry,
        double dt
    ) = 0;
};

/**
 * Базовый интерфейс стратегии обновления полей в геометрии в выбранной системе координат.
 *
 * Стратегия обновления полей - это абстракция, которая берёт на себя задачу
 *      обновления и рассчёта эффективных полей и энергий в геометрии на каждый магнитный момент
 */
template <CoordSystemConcept CoordSystem> class PYTHON_API IFieldUpdater {
  protected:
    // конструктор только для наследников
    IFieldUpdater() = default;

  public:
    // деструктор
    virtual ~IFieldUpdater() = default;

    // Рассчитать эффективные поля и энергии согласно заданной геометрии
    PYTHON_API virtual SolverData<CoordSystem> calculateFields(
        IGeometry<CoordSystem> &geometry,
        InteractionRegistry<CoordSystem> &interaction_registry,
        MaterialRegistry &material_registry
    ) = 0;
};

namespace cartesian {

/**
 * Базовый интерфейс решателя (интегратора) в выбранной (ДЕКАРТОВОЙ) системе координат.
 *
 * Решатель - это абстракция, которая берёт на себя задачу
 *  эволюционирования или любого другого прогрессирования системы в зависимости от
 *  предоставленных внешних параметров.
 * Типичный решатель использует стратегию интегрирования (метод Эйлера, Хойна и т.д.).
 * А также рассчитывает самостоятельно эффективные поля и энергии взаимодействий
 *     на основе предоставленной ему стратегии расчёта (IFieldUpdater).
 */
using AbstractSolver = PYTHON_API ISolver<NamespaceCoordSystem>;

/**
 * Базовый интерфейс стратегии обновления полей в геометрии в выбранной (ДЕКАРТОВОЙ) системе координат.
 *
 * Стратегия обновления полей - это абстракция, которая берёт на себя задачу
 *      обновления и рассчёта эффективных полей и энергий в геометрии на каждый магнитный момент
 */
using AbstractFieldUpdater = PYTHON_API IFieldUpdater<NamespaceCoordSystem>;

/**
 * Структура для организации данных, которые возвращает и использует решатель (интегратор).
 *      В декартовой системе координат.
 *
 * Содержит эффективные поля и энергии, в т.ч. по каждому из взаимодействий.
 * Может использоваться как буфер для хранения промежуточных данных.
 */
using AbstractSolverData = PYTHON_API SolverData<NamespaceCoordSystem>;

/**
 * Стратегия обновления поле в декаторой геометрии.
 *
 * Использует OpenMP для параллельного расчёта эффективных полей и энергий взаимодействий.
 */
class PYTHON_API OMPFieldUpdater : public AbstractFieldUpdater {
  protected:
    AbstractSolverData _step_data;

    // Рассчитать эффективные поля и энергии согласно заданной геометрии для конкретного взаимодействия
    void calculateFieldContribution(
        AbstractGeometry &geometry,
        regnum interaction_regnum,
        AbstractInteraction *interaction,
        MaterialRegistry &material_registry
    ) {
        // --- предсоздать буфферы ---

        size_t moments_size = geometry.size(); // Размер системы
        EffectiveFieldVector contribution_vector(moments_size
        ); // Вектор для полевого вклада текущего interaction
        std::vector<double> energy(moments_size); // Вектор для энергетич. вклада текущего interaction

        // --- посчитать энергии и эфф. поля в CPU-параллели ---

        // clang-format off
        #pragma omp parallel for schedule(dynamic)
        // clang-format on
        for (size_t i = 0; i < moments_size; ++i) {
            contribution_vector[i] = interaction->calculateFieldContribution(i, geometry, material_registry);
            energy[i] = interaction->calculateEnergy(geometry[i], contribution_vector[i]);
        }

        // --- проверить на правильность результатов ---

        for (size_t i = 0; i < moments_size; ++i) {
            if (contribution_vector[i].hasNaN()) {
                throw std::runtime_error("Invalid contribution_vector data (NaN detected)");
            }
        }

        // --- записать в текущее шаговое состояние (гарантировать последовательность записи) ---

        // clang-format off
        #pragma omp critical
        // clang-format on
        for (size_t i = 0; i < moments_size; ++i) {
            this->_step_data.effective_fields[i] += contribution_vector[i];
            this->_step_data.interaction_effective_fields[interaction_regnum][i] = contribution_vector[i];
            // и энергии
            this->_step_data.energies[i] += energy[i];
            this->_step_data.interaction_energies[interaction_regnum][i] = energy[i];
        }
    }

  public:
    // Конструктор
    PYTHON_API OMPFieldUpdater() {};

    // Рассчитать эффективные поля и энергии согласно заданной геометрии
    PYTHON_API virtual AbstractSolverData calculateFields(
        AbstractGeometry &geometry,
        AbstractInteractionRegistry &interaction_registry,
        MaterialRegistry &material_registry
    ) override {
        size_t moments_size = geometry.size();                      // Размер системы
        this->_step_data.clear(moments_size, interaction_registry); // Занулить буфер
        this->_step_data.correct(moments_size, interaction_registry); // Не потерять входящие изменения

        // обсчёт полей и энергий от зарегистрированных взаимодействий с записью во внутреннее состояние
        for (auto &[interaction_regnum, interaction] : interaction_registry) {
            this->calculateFieldContribution(
                geometry, interaction_regnum, interaction.get(), material_registry
            );
        }

        return this->_step_data;
    }
};

/**
 * Решатель ЛЛГ (Ландау-Лифшица-Гильберта) - это решатель, который использует
 *  уравнение ЛЛГ спиновой динамики для обновления моментов в геометрии.
 *
 * Решатель использует метод Эйлера (+ midpoint) или метод Хойна для интегрирования уравнения ЛЛГ с шагом dt.
 * А также рассчитывает самостоятельно эффективные поля и энергии взаимодействий
 *     на основе предоставленной ему стратегии расчёта (IFieldUpdater).
 */
class PYTHON_API LLGSolver : public AbstractSolver {
  protected:
    // выбранная стратегия решателя (интегратора) - метод интегрирования
    SolverStrategy _strategy;
    // стратегия вычисления эффективного поля на каждый спин
    std::unique_ptr<AbstractFieldUpdater> _field_updater;

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

    // вычисление изменения моментов по уравнению ЛЛГ (правые части системы уравнений)
    virtual std::vector<Eigen::Vector3d> calculateLLGMomentsChange(
        IGeometry<NamespaceCoordSystem> &geometry, AbstractSolverData &data
    ) {
        auto moments_size = geometry.size();
        std::vector<Eigen::Vector3d> result;
        result.resize(moments_size, Eigen::Vector3d::Zero());
        for (size_t i = 0; i < moments_size; ++i) {
            auto &moment = geometry[i];
            auto &base_moment_vector = moment.getDirection().asVector();

            result[i] = this->calculateLLGMomentChange(
                moment.getMaterial(), base_moment_vector, data.effective_fields[i]
            );
        }
        return result;
    }

    // обновить моменты в геометрии согласно правым частям уравнения ЛЛГ
    virtual void updateMomentsInGeometry(
        IGeometry<NamespaceCoordSystem> &geometry, std::vector<Eigen::Vector3d> &moment_changes, double dt
    ) {
        auto moments_size = geometry.size();
        if (moment_changes.size() != moments_size) {
            throw std::invalid_argument("Invalid moment changes size");
        }
        for (size_t i = 0; i < moments_size; ++i) {
            auto &moment = geometry[i];
            Eigen::Vector3d new_moment = moment.getDirection().asVector() + moment_changes[i] * dt;
            moment.setDirection(new_moment);
        }
    }

    // решение по методу Эйлера: сделать шаг с конечными разностями
    virtual AbstractSolverData updateMomentsViaEuler(
        IGeometry<NamespaceCoordSystem> &geometry,
        InteractionRegistry<NamespaceCoordSystem> &interaction_registry,
        MaterialRegistry &material_registry,
        double dt
    ) {
        auto solver_data =
            this->_field_updater->calculateFields(geometry, interaction_registry, material_registry);
        auto moments_changes = this->calculateLLGMomentsChange(geometry, solver_data);
        this->updateMomentsInGeometry(geometry, moments_changes, dt);
        return solver_data;
    }

    // решение по методу Хойна: предиктор-корректор
    virtual AbstractSolverData updateMomentsViaHeun(
        IGeometry<NamespaceCoordSystem> &geometry,
        InteractionRegistry<NamespaceCoordSystem> &interaction_registry,
        MaterialRegistry &material_registry,
        double dt
    ) {
        // рассчитать поля и изменения для предиктора
        auto predictor_data =
            this->_field_updater->calculateFields(geometry, interaction_registry, material_registry);
        auto predictor_geometry = geometry.clone();
        auto predictor_moments_changes = this->calculateLLGMomentsChange(geometry, predictor_data);

        // рассчитать новое направление спина (predictor_geometry изменилась)
        this->updateMomentsInGeometry(*predictor_geometry, predictor_moments_changes, dt);
        // рассчитать новые поля для корректора на основе изменённой геометрии предиктором
        auto corrector_data = this->_field_updater->calculateFields(
            *predictor_geometry, interaction_registry, material_registry
        );
        auto corrector_moments_changes = this->calculateLLGMomentsChange(*predictor_geometry, corrector_data);

        // усреднить изменения dS
        std::vector<Eigen::Vector3d> avg_moments_changes(geometry.size());
        for (size_t i = 0; i < geometry.size(); ++i) {
            avg_moments_changes[i] = (predictor_moments_changes[i] + corrector_moments_changes[i]) * 0.5;
        }

        // посчитать изменение на обобщённом dS и начальной геометрии
        this->updateMomentsInGeometry(geometry, avg_moments_changes, dt);

        // вернуть данные, которые привели к текущей ориентации (предиктор) TODO: усреднить пред+корр
        return predictor_data;
    }

  public:
    // базовый конструктор
    PYTHON_API LLGSolver(SolverStrategy strategy = SolverStrategy::EULER)
        : _strategy(strategy), _field_updater(std::make_unique<OMPFieldUpdater>()) {};

    // обновить моменты в геометрии, сделать "шаг"
    PYTHON_API virtual AbstractSolverData updateMoments(
        IGeometry<NamespaceCoordSystem> &geometry,
        InteractionRegistry<NamespaceCoordSystem> &interaction_registry,
        MaterialRegistry &material_registry,
        double dt
    ) override {
        //
        if (this->_strategy == SolverStrategy::EULER) {
            return this->updateMomentsViaEuler(geometry, interaction_registry, material_registry, dt);
        } else if (this->_strategy == SolverStrategy::HEUN) {
            return this->updateMomentsViaHeun(geometry, interaction_registry, material_registry, dt);
        };

        throw std::invalid_argument("Invalid solver strategy enum value / doesnt support");
    };
};

}; // namespace cartesian

}; // namespace PYTHON_API spindynapy

#define SOLVER_TEMPLATE_BINDINGS(cls)                                                                        \
    .def(                                                                                                    \
        "update_moments",                                                                                    \
        &cls::updateMoments,                                                                                 \
        py::arg("geometry"),                                                                                 \
        py::arg("interaction_registry"),                                                                     \
        py::arg("material_registry"),                                                                        \
        py::arg("dt"),                                                                                       \
        py::doc("обновить моменты в геометрии, сделать \"шаг\"")                                             \
    )

#define SOLVERDATA_TEMPLATE_BINDINGS(cls)                                                                    \
    .def_readwrite("effective_fields", &cls::effective_fields)                                               \
        .def_readwrite("energies", &cls::energies)                                                           \
        .def_readwrite("interaction_effective_fields", &cls::interaction_effective_fields)                   \
        .def_readwrite("interaction_energies", &cls::interaction_energies)                                   \
        .def(                                                                                                \
            "clear",                                                                                         \
            &cls::clear,                                                                                     \
            py::arg("moments_size"),                                                                         \
            py::arg("interaction_registry"),                                                                 \
            py::doc("очистить и занулить состояние на основе регистра взаимодействий")                       \
        )                                                                                                    \
        .def(                                                                                                \
            "correct",                                                                                       \
            &cls::correct,                                                                                   \
            py::arg("moments_size"),                                                                         \
            py::arg("interaction_registry"),                                                                 \
            py::doc("в случае, если изменилась геометрия, то нужно обновить размеры массивов")               \
        )

#define FIELDUPDATER_TEMPLATE_BINDINGS(cls)                                                                  \
    .def(                                                                                                    \
        "calculate_fields",                                                                                  \
        &cls::calculateFields,                                                                               \
        py::arg("geometry"),                                                                                 \
        py::arg("interaction_registry"),                                                                     \
        py::arg("material_registry"),                                                                        \
        py::doc("рассчитать эффективные поля и энергии согласно заданной геометрии")                         \
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
        .value("HEUN", SolverStrategy::HEUN, "Метод Хойна (предиктор-корректор (2 расчёта))");

    // -------- | CARTESIAN SOLVERS | --------
    py::module_ cartesian = solvers_module.def_submodule("cartesian");

    cartesian.doc() = "Решатель - это абстракция, которая берёт на себя задачу\n"
                      " эволюционирования или любого другого прогрессирования системы в зависимости от\n"
                      " предоставленных внешних параметров.\n"
                      "\n"
                      "Стратегия решателя (интегратора) - метод интегрирования: EULER, HEUN, ...";
    {
        using cartesian::AbstractFieldUpdater;
        using cartesian::AbstractSolver;
        using cartesian::AbstractSolverData;
        using cartesian::LLGSolver;
        using cartesian::OMPFieldUpdater;

        py::class_<AbstractSolverData>(cartesian, "SolverData")
            SOLVERDATA_TEMPLATE_BINDINGS(AbstractSolverData)
                .doc() =
            "Структура для организации данных, которые возвращает и использует решатель (интегратор).\n"
            "     В декартовой системе координат.\n"
            "\n"
            "Содержит эффективные поля и энергии, в т.ч. по каждому из взаимодействий.\n"
            "Может использоваться как буфер для хранения промежуточных данных.\n";

        py::class_<AbstractFieldUpdater, std::shared_ptr<AbstractFieldUpdater>>(
            cartesian, "AbstractFieldUpdater"
        ) FIELDUPDATER_TEMPLATE_BINDINGS(AbstractFieldUpdater)
            .doc() =
            "Базовый интерфейс стратегии обновления полей в геометрии в выбранной (ДЕКАРТОВОЙ) системе "
            "координат.\n"
            "\n"
            "Стратегия обновления полей - это абстракция, которая берёт на себя задачу\n"
            "  обновления и рассчёта эффективных полей и энергий в геометрии на каждый магнитный момент.";

        py::class_<OMPFieldUpdater, AbstractFieldUpdater, std::shared_ptr<OMPFieldUpdater>>(
            cartesian, "OMPFieldUpdater"
        )
            .def(py::init<>())
            .doc() =
            "Стратегия обновления поле в декаторой геометрии.\n"
            "\n"
            "  Использует OpenMP для параллельного расчёта эффективных полей и энергий взаимодействий.";

        py::class_<AbstractSolver, std::shared_ptr<AbstractSolver>>(cartesian, "AbstractSolver")
            SOLVER_TEMPLATE_BINDINGS(AbstractSolver)
                .doc() =
            "Базовый интерфейс решателя (интегратора) в выбранной (ДЕКАРТОВОЙ) системе координат.\n"
            "Решатель - это абстракция, которая берёт на себя задачу\n"
            "    эволюционирования или любого другого прогрессирования системы в зависимости от\n"
            "    предоставленных внешних параметров.\n"
            "\n"
            "Типичный решатель использует стратегию интегрирования (метод Эйлера, Хойна и т.д.).\n"
            "А также рассчитывает самостоятельно эффективные поля и энергии взаимодействий\n"
            "    на основе предоставленной ему стратегии расчёта (IFieldUpdater).";

        py::class_<LLGSolver, AbstractSolver, std::shared_ptr<LLGSolver>>(cartesian, "LLGSolver")
            .def(py::init<SolverStrategy>(), py::arg("strategy") = SolverStrategy::EULER)
            .doc() = "Решатель ЛЛГ (Ландау-Лифшица-Гильберта) - это решатель, который использует\n"
                     "  уравнение ЛЛГ спиновой динамики для обновления моментов в геометрии.\n"
                     "\n"
                     "Решатель использует метод Эйлера (+ midpoint) или метод Хойна для интегрирования "
                     "   уравнения ЛЛГ с шагом dt.\n"
                     "А также рассчитывает самостоятельно эффективные поля и энергии взаимодействий\n"
                     "   на основе предоставленной ему стратегии расчёта (IFieldUpdater).";
    }
}

#endif // ! __SOLVERS_HPP__
