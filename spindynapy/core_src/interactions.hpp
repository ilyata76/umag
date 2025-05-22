#ifndef __INTERACTIONS_HPP__
#define __INTERACTIONS_HPP__

/**
 * Интерфейсы и реализации взаимодействий в выбранной системе координат в симуляции.
 *   Абстрагируют взаимодействие между моментами в геометрии.
 *   Предоставляют интерфейс для расчёта эффективного поля на моменте и энергии взаимодействия.
 *
 * Взаимодействия могут быть с внешнием полем, между моментами и т.д.
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

namespace PYTHON_API spindynapy {

/**
 * Базовый интерфейс взаимодействий в выбранной системе координат.
 *
 * Абстрагирует взаимодействие между моментами в геометрии.
 *   Предоставляет интерфейс для расчёта эффективного поля на моменте и энергии взаимодействия.
 */
template <CoordSystemConcept CoordSystem> class PYTHON_API IInteraction {
  protected:
    // конструктор только для наследников
    IInteraction() = default;

  public:
    // деструктор
    virtual ~IInteraction() = default;

    // посчитать вклад в эффективное поле на моменте от взаимодействия в геометрии
    PYTHON_API virtual EffectiveField calculateFieldContribution(
        size_t moment_index, IGeometry<CoordSystem> &geometry, MaterialRegistry &material_registry
    ) const = 0;

    // посчитать энергию взаимодействия между моментом и приложенным эффективным полем
    PYTHON_API virtual double calculateEnergy(CoordSystem::Moment &moment, EffectiveField &field) const = 0;

    // получить название взаимодействия
    PYTHON_API virtual std::string getName() const = 0;

    // строковое представление для принтинга (TODO заменить на ostream/stringstream)
    PYTHON_API virtual std::string __str__() const = 0;
};

/*
 * Регистр взаимодействий (см. registires.hpp).
 *   Позволяет хранить и управлять различными взаимодействиями в системе.
 *
 * Основной цикл симуляции может использовать этот регистр для получения всех взаимодействий,
 *    которые необходимо учитывать при расчёте эффективного поля на моменте.
 */
template <CoordSystemConcept CoordSystem>
using InteractionRegistry = PYTHON_API Registry<IInteraction<CoordSystem>>;

// ^
// | base template interfaces
// |
// ================= NEW_BLOCK ===================
// |
// | cartesian realization of interfaces (namespace cartesian)
// v

namespace PYTHON_API cartesian {

/**
 * Базовый интерфейс взаимодействий в выбранной (декартовой) системе координат.
 *
 * Абстрагирует взаимодействие между моментами в декартовой геометрии.
 *   Предоставляет интерфейс для расчёта эффективного поля на моменте и энергии взаимодействия.
 */
using AbstractInteraction = PYTHON_API IInteraction<NamespaceCoordSystem>;

/*
 * Обменное взаимодействие между моментами.
 *
 * Выбранная стратегия - радиус обрезки: считаются все соседи в радиусе обрезки.
 *   Модель Гейзенберга.
 *
 * TODO: формулу сюда вставить
 */
class PYTHON_API ExchangeInteraction : public AbstractInteraction {
  protected:
    // раидус обрезки соседей (внутри которого будем считать)
    double _cutoff_radius;

  public:
    // базовый конструктор
    PYTHON_API ExchangeInteraction(double cutoff_radius) : _cutoff_radius(cutoff_radius) {};

    // посчитать энергию взаимодействия между моментом и приложенным эффективным полем
    PYTHON_API double calculateEnergy(NamespaceCoordSystem::Moment &moment, EffectiveField &field)
        const override {
        const auto &current_material = moment.getMaterial();
        const auto &direction = moment.getDirection().asVector();
        const auto current_moment_norm =
            current_material.atomic_magnetic_saturation_magnetization * constants::BOHR_MAGNETON;
        return -current_moment_norm / 2 *
               direction.dot(field); // TODO сделать через флаг pairwise на вызывающей стороне
    }

    // посчитать энергию взаимодействия между моментом и приложенным эффективным полем
    PYTHON_API virtual EffectiveField
    calculateFieldContribution(size_t moment_index, IGeometry<NamespaceCoordSystem> &geometry, MaterialRegistry &)
        const override {
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
                                   : GENERALIZED_PREFIX /* тут можно будет интерфейс использовать */
                              ) *
                              geometry[neighbor_index].getDirection().asVector();
        }
        return exchange_field;
    }

    // получить название взаимодействия
    PYTHON_API virtual std::string getName() const override { return "EXCHANGE"; }

    // строковое представление для принтинга
    PYTHON_API virtual std::string __str__() const override {
        return std::format("ExchangeInteraction(r={})", _cutoff_radius);
    };
};

/*
 * Взаимодействие с внешним полем.
 *
 * TODO: формулу сюда вставить
 */
class PYTHON_API ExternalInteraction : public AbstractInteraction {
  protected:
    Eigen::Vector3d external_field;

  public:
    // конструктор по умолчанию
    PYTHON_API ExternalInteraction(const Eigen::Vector3d &external_field) : external_field(external_field) {};
    // конструктор по умолчанию
    PYTHON_API ExternalInteraction(double sx, double sy, double sz) : external_field(sx, sy, sz) {};

    // посчитать энергию взаимодействия между моментом и приложенным эффективным полем
    PYTHON_API double calculateEnergy(NamespaceCoordSystem::Moment &moment, EffectiveField &field)
        const override {
        const auto &current_material = moment.getMaterial();
        const auto &direction = moment.getDirection().asVector();
        const auto current_moment_norm =
            current_material.atomic_magnetic_saturation_magnetization * constants::BOHR_MAGNETON;
        return -current_moment_norm * direction.dot(field);
    }

    // посчитать энергию взаимодействия между моментом и приложенным эффективным полем
    PYTHON_API virtual EffectiveField
    calculateFieldContribution(size_t, IGeometry<NamespaceCoordSystem> &, MaterialRegistry &) const override {
        return this->external_field;
    }

    // получить название взаимодействия
    PYTHON_API virtual std::string getName() const override { return "EXTERNAL"; }

    // строковое представление для принтинга
    PYTHON_API virtual std::string __str__() const override {
        return std::format(
            "ExternalInteraction(field=({}, {}, {}))",
            this->external_field.x(),
            this->external_field.y(),
            this->external_field.z()
        );
    };
};

/*
 * Воздействие на моменты в силу магнитной анизотропии.
 *
 * Выбранная стратегия: в зависимости от типа анизотропии (односторонняя, кубическая и т.д.)
 *   рассчитывается вклад в эффективное поле на моменте.
 *
 * TODO: формулу сюда вставить
 */
class PYTHON_API AnisotropyInteraction : public AbstractInteraction {
  public:
    // конструктор по умолчанию
    PYTHON_API AnisotropyInteraction() {};

    // посчитать энергию взаимодействия между моментом и приложенным эффективным полем
    PYTHON_API double calculateEnergy(NamespaceCoordSystem::Moment &moment, EffectiveField &field)
        const override {
        const auto &current_material = moment.getMaterial();
        const auto &direction = moment.getDirection().asVector();
        const auto current_moment_norm =
            current_material.atomic_magnetic_saturation_magnetization * constants::BOHR_MAGNETON;

        // из-за взятия производной от квадрата => / 2
        return -current_moment_norm / 2 * direction.dot(field);
    }

    // посчитать энергию взаимодействия между моментом и приложенным эффективным полем
    PYTHON_API virtual EffectiveField
    calculateFieldContribution(size_t moment_index, IGeometry<NamespaceCoordSystem> &geometry, MaterialRegistry &)
        const override {
        //
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

    // получить название взаимодействия
    PYTHON_API virtual std::string getName() const override { return "ANISOTROPY"; }

    // строковое представление для принтинга
    PYTHON_API virtual std::string __str__() const override { return "AnisotropyInteraction()"; };
};

/*
 * Магнитостатическое (диполь-дипольное) взаимодействие между моментами.
 *
 * В зависимости от выбранной стратегии (cutoff или macrocells)
 *   рассчитывается вклад в эффективное поле на моменте.
 *
 * 1. Макроячейки - это группы моментов, которые считаются как единое целое с усреднением по величинам.
 * 2. cutoff - это радиус, в пределах которого считаются все соседи или соседние макроячейки.
 *
 * TODO: формулу сюда вставить
 */
class PYTHON_API DipoleDipoleInteraction : public AbstractInteraction {
  protected:
    std::string _strategy;
    double _cutoff_radius;

    // посчитать энергию взаимодействия между моментом и другими, полученных по стратегии
    EffectiveField calculate(
        Moment &current_moment, Material &, MomentsContainer<NamespaceCoordSystem> calculation_moments
    ) const {
        EffectiveField dipole_field = EffectiveField::Zero();

        for (auto moment : calculation_moments) {
            auto distance_vector =
                moment->getCoordinates().asVector() - current_moment.getCoordinates().asVector();

            auto neighbor_atomistic_moment = moment->amplitude *
                                             moment->getMaterial().atomic_magnetic_saturation_magnetization *
                                             constants::BOHR_MAGNETON;
            auto distance_norm = distance_vector.norm();

            if (distance_norm < 1e-30) {
                // предотвратить деление на ноль (превращение в NaN)
                continue;
            }

            dipole_field += neighbor_atomistic_moment *
                            ((3 * moment->getDirection().asVector().dot(distance_vector) * distance_vector) /
                                 pow(distance_norm, 5) -
                             moment->getDirection().asVector() / pow(distance_norm, 3));
        }

        dipole_field = (constants::VACUUM_MAGNETIC_PERMEABILITY / 4 / constants::NUMBER_PI) * dipole_field;

        return dipole_field;
    }

  public:
    // базовый конструктор
    PYTHON_API DipoleDipoleInteraction(double cutoff_radius, std::string strategy = "cutoff")
        : _strategy(strategy), _cutoff_radius(cutoff_radius) {
        if (strategy != "cutoff" && strategy != "macrocells") {
            throw std::invalid_argument("Invalid strategy string");
        }
    };

    // посчитать энергию взаимодействия между моментом и приложенным эффективным полем
    PYTHON_API double calculateEnergy(NamespaceCoordSystem::Moment &moment, EffectiveField &field)
        const override {
        const auto &current_material = moment.getMaterial();
        const auto &direction = moment.getDirection().asVector();
        const auto current_moment_norm =
            current_material.atomic_magnetic_saturation_magnetization * constants::BOHR_MAGNETON;
        return -0.5 * current_moment_norm * moment.amplitude *
               direction.dot(field); // TODO сделать через флаг pairwise на вызывающей стороне
    }

    // посчитать энергию взаимодействия между моментом и приложенным эффективным полем
    PYTHON_API virtual EffectiveField
    calculateFieldContribution(size_t moment_index, IGeometry<NamespaceCoordSystem> &geometry, MaterialRegistry &)
        const override {
        auto &current_moment = geometry[moment_index];
        auto &current_material = current_moment.getMaterial();

        if (this->_strategy == "cutoff") {
            auto neighbor_indices = geometry.getNeighbors(moment_index, this->_cutoff_radius);
            auto calculation_moments = geometry.getFromIndexes(neighbor_indices);
            return this->calculate(current_moment, current_material, calculation_moments);
        } else if (this->_strategy == "macrocells") {
            // TODO как учесть внутри? получать список внутри ячейки находящихся
            auto calculation_moments = geometry.getMomentsFromMacrocells(moment_index, this->_cutoff_radius);
            return this->calculate(current_moment, current_material, calculation_moments);
        }
        // Если не удалось найти подходящую стратегию
        throw std::invalid_argument("Invalid strategy string");
    }

    // получить название взаимодействия
    PYTHON_API virtual std::string getName() const override { return "DIPOLE-DIPOLE"; }

    // строковое представление для принтинга
    PYTHON_API virtual std::string __str__() const override {
        return std::format(
            "DipoleDipoleInteraction(cutoff_radius={}, strategy={})", this->_cutoff_radius, this->_strategy
        );
    };
};

/**
 * Считается, что взаимодействие между моментами в макроячейках
 *   обосновывает т.н. энергию демагнетизации.
 *
 * Текущее взаимодействие - это диполь-дипольное взаимодействие между моментами в макроячейках,
 *    с учётом демагнетизирующего характера.
 *
 * В зависимости от выбранной стратегии (cutoff или macrocells) рассчитывается вклад в эффективное поле на
 * моменте.
 *
 * 1. Макроячейки - это группы моментов, которые считаются как единое целое с усреднением по величинам.
 *          Для макроячеек предусмотрен т.н. self-term, который считается через кубический фактор
 *                  размагничивания. (см. VAMPIRE)
 *
 * 2. cutoff - это радиус, в пределах которого считаются все соседи или соседние макроячейки.
 *
 * TODO поправить описание... cutoff подходит для всех . а здесь просто ШТРАФ т.е. -(-E)
 * и селф терм. а выше другой учёт у макроячеек
 */
class DemagnetizationInteraction : public DipoleDipoleInteraction {
  public:
    // базовый конструктор
    PYTHON_API DemagnetizationInteraction(double cutoff_radius, std::string strategy = "macrocells")
        : DipoleDipoleInteraction(cutoff_radius, strategy) {};

    // посчитать энергию взаимодействия между моментом и приложенным эффективным полем
    PYTHON_API virtual EffectiveField
    calculateFieldContribution(size_t moment_index, IGeometry<NamespaceCoordSystem> &geometry, MaterialRegistry &)
        const override {
        auto &current_moment = geometry[moment_index];
        auto &current_material = current_moment.getMaterial();

        if (this->_strategy == "cutoff") {
            auto neighbor_indices = geometry.getNeighbors(moment_index, this->_cutoff_radius);
            auto calculation_moments = geometry.getFromIndexes(neighbor_indices);
            return this->calculate(current_moment, current_material, calculation_moments);
        } else if (this->_strategy == "macrocells") {

            auto calculation_moments = geometry.getMomentsFromMacrocells(moment_index, this->_cutoff_radius);

            auto macrocell_moment =
                geometry.getMacrocells()[geometry.getMacrocellIndexBySpin(moment_index)].avg_moment;

            auto moment_term = macrocell_moment->amplitude * macrocell_moment->getDirection().asVector() *
                               macrocell_moment->getMaterial().atomic_magnetic_saturation_magnetization *
                               constants::BOHR_MAGNETON;

            // TODO с VAMPIRE не сходится, почему?
            auto self_term = (constants::VACUUM_MAGNETIC_PERMEABILITY / 3.0) *
                             (moment_term / pow(geometry.getMacrocellSize(), 3));

            return this->calculate(current_moment, current_material, calculation_moments) + self_term;
        }
        // Если не удалось найти подходящую стратегию
        throw std::invalid_argument("Invalid strategy string");
    }

    // получить название взаимодействия
    PYTHON_API virtual std::string getName() const override { return "DEMAGNETIZATION"; }

    // строковое представление для принтинга
    PYTHON_API virtual std::string __str__() const override {
        return std::format(
            "DemagnetizationInteraction(cutoff_radius={}, strategy={})", this->_cutoff_radius, this->_strategy
        );
    };
};

/*
 * Регистр взаимодействий (см. registires.hpp).
 *   Позволяет хранить и управлять различными взаимодействиями в системе.
 *
 * Основной цикл симуляции может использовать этот регистр для получения всех взаимодействий,
 *    которые необходимо учитывать при расчёте эффективного поля на моменте.
 */
using AbstractInteractionRegistry = PYTHON_API InteractionRegistry<NamespaceCoordSystem>;

}; // namespace PYTHON_API cartesian

}; // namespace PYTHON_API spindynapy

// ^
// | cartesian realization of interfaces (namespace cartesian)
// |
// ================= NEW_BLOCK ===================
// |
// | MACROSES and BINDINGS for PYTHON
// v

#define INTERACTION_TEMPLATE_BINDINGS(cls)                                                                   \
    .def("__str__", &cls::__str__, py::doc("строковое представление для принтинга"))                         \
        .def(                                                                                                \
            "calculate_field_contribution",                                                                  \
            &cls::calculateFieldContribution,                                                                \
            py::arg("moment_index"),                                                                         \
            py::arg("geometry"),                                                                             \
            py::arg("material_registry"),                                                                    \
            py::doc("посчитать вклад в эффективное поле на моменте от взаимодействия в геометрии")           \
        )                                                                                                    \
        .def(                                                                                                \
            "calculate_energy",                                                                              \
            &cls::calculateEnergy,                                                                           \
            py::arg("moment"),                                                                               \
            py::arg("field"),                                                                                \
            py::doc("посчитать энергию взаимодействия между моментом и приложенным эффективным полем")       \
        )                                                                                                    \
        .def("get_name", &cls::getName, py::doc("получить название взаимодействия"))

// функция для связывания взаимодействий с Python
inline void pyBindInteractions(py::module_ &module) {
    using namespace spindynapy;

    // -------- | INTERACTIONS | --------
    py::module_ interaction_module = module.def_submodule("interactions");

    interaction_module.doc() =
        "Интерфейсы и реализации взаимодействий в выбранной системе координат в симуляции.\n"
        "  Абстрагируют взаимодействие между моментами в геометрии.\n"
        "  Предоставляют интерфейс для расчёта эффективного поля на моменте и энергии взаимодействия.\n"
        "\n"
        "Взаимодействия могут быть с внешнием полем, между моментами и т.д.";

    // -------- | CARTESIAN INTERACTIONS | --------
    py::module_ cartesian = interaction_module.def_submodule("cartesian");

    cartesian.doc() =
        "Интерфейсы и реализации взаимодействий в выбранной системе (ДЕКАРТОВОЙ) координат в симуляции.\n"
        "  Абстрагируют взаимодействие между моментами в геометрии.\n"
        "  Предоставляют интерфейс для расчёта эффективного поля на моменте и энергии взаимодействия.\n"
        "\n"
        "Взаимодействия могут быть с внешнием полем, между моментами и т.д.";

    { // TODO везде такую "защиту" расставить
        using cartesian::AbstractInteraction;
        using cartesian::AbstractInteractionRegistry;
        using cartesian::AnisotropyInteraction;
        using cartesian::DemagnetizationInteraction;
        using cartesian::DipoleDipoleInteraction;
        using cartesian::ExchangeInteraction;
        using cartesian::ExternalInteraction;

        py::class_<AbstractInteraction, std::shared_ptr<AbstractInteraction>>(
            cartesian, "AbstractInteraction"
        ) INTERACTION_TEMPLATE_BINDINGS(AbstractInteraction)
            .doc() =
            "Базовый интерфейс взаимодействий в выбранной (декартовой) системе координат.\n"
            "\n"
            "Абстрагирует взаимодействие между моментами в геометрии.\n"
            "  Предоставляет интерфейс для расчёта эффективного поля на моменте и энергии взаимодействия.";

        py::class_<ExchangeInteraction, AbstractInteraction, std::shared_ptr<ExchangeInteraction>>(
            cartesian, "ExchangeInteraction"
        )
            .def(py::init<double>(), py::arg("cutoff_radius"))
            .doc() = "Обменное взаимодействие между моментами.\n"
                     "\n"
                     "Выбранная стратегия - радиус обрезки: считаются все соседи в радиусе обрезки.\n"
                     "  Модель Гейзенберга.";

        py::class_<ExternalInteraction, AbstractInteraction, std::shared_ptr<ExternalInteraction>>(
            cartesian, "ExternalInteraction"
        )
            .def(py::init<const Eigen::Vector3d &>(), py::arg("external_field"))
            .def(py::init<double, double, double>(), py::arg("sx"), py::arg("sy"), py::arg("sz"))
            .doc() = "Взаимодействие с внешним полем.\n"
                     "\n"
                     "Выбранная стратегия: просто поле в Тл";

        py::class_<AnisotropyInteraction, AbstractInteraction, std::shared_ptr<AnisotropyInteraction>>(
            cartesian, "AnisotropyInteraction"
        )
            .def(py::init())
            .doc() =
            "Воздействие на моменты в силу магнитной анизотропии.\n"
            "\n"
            "Выбранная стратегия: в зависимости от типа анизотропии (односторонняя, кубическая и т.д.)\n"
            "  рассчитывается вклад в эффективное поле на моменте.";

        py::class_<DipoleDipoleInteraction, AbstractInteraction, std::shared_ptr<DipoleDipoleInteraction>>(
            cartesian, "DipoleDipoleInteraction"
        )
            .def(py::init<double, std::string>(), py::arg("cutoff_radius"), py::arg("strategy") = "cutoff")
            .doc() =
            "Магнитостатическое взаимодействие между моментами.\n"
            "\n"
            "В зависимости от выбранной стратегии (cutoff или macrocells)\n"
            "  рассчитывается вклад в эффективное поле на моменте.\n"
            "\n"
            "1. Макроячейки - это группы моментов, которые считаются как единое целое с усреднением по "
            "величинам.\n"
            "2. cutoff - это радиус, в пределах которого считаются все соседи или соседние макроячейки.";

        py::class_<
            DemagnetizationInteraction,
            DipoleDipoleInteraction,
            std::shared_ptr<DemagnetizationInteraction>>(cartesian, "DemagnetizationInteraction")
            .def(py::init<double, std::string>(), py::arg("cutoff_radius"), py::arg("strategy") = "cutoff")
            .doc() =
            ("Считается, что взаимодействие между моментами в макроячейках\n"
             "  обосновывает т.н. энергию демагнетизации.\n"
             "\n"
             "Текущее взаимодействие - это диполь-дипольное взаимодействие между моментами в макроячейках,\n"
             "   с учётом демагнетизирующего характера.\n"
             "\n"
             "В зависимости от выбранной стратегии (cutoff или macrocells) рассчитывается вклад в "
             "эффективное "
             "поле на моменте.\n"
             "\n"
             "1. Макроячейки - это группы моментов, которые считаются как единое целое с усреднением по "
             "величинам.\n"
             "         Для макроячеек предусмотрен т.н. self-term, который считается через кубический фактор "
             "размагничивания.\n"
             "          (см. VAMPIRE)\n"
             "\n"
             "2. cutoff - это радиус, в пределах которого считаются все соседи или соседние макроячейки.");

        py::class_<AbstractInteractionRegistry, std::shared_ptr<AbstractInteractionRegistry>>(
            cartesian, "InteractionRegistry"
        )
            .def(py::init<>())
            .def(py::init<RegistryContainer<AbstractInteraction>>(), py::arg("container"))
                REGISTRY_TEMPLATE_BINDINGS(AbstractInteractionRegistry)
            .doc() =
            "Регистр взаимодействий (см. registires.hpp).\n"
            "  Позволяет хранить и управлять различными взаимодействиями в системе.\n"
            "\n"
            "Основной цикл симуляции может использовать этот регистр для получения всех взаимодействий,\n"
            "   которые необходимо учитывать при расчёте эффективного поля на моменте.\n"
            "БЕЗ ЗАПИСИ => ПОТОКОБЕЗОПАСНЫЙ";
    }
}

#endif // ! __INTERACTIONS_HPP__
