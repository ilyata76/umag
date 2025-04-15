#ifndef __TYPES_HPP__
#define __TYPES_HPP__

/**
 * Заголовки с базовыми интерфейсами для базовых типов/классов:
 *    - стандартные миксины
 *    - стандартные единицы абстракции (материал, спин, координаты, etc.)
 *
 * Запрещают конструктор по умолчанию.
 * Все процессы должны зависеть от этих интерфейсов.
 */

#include <format>
#include <memory>
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <stdexcept>
#include <string>

namespace py = pybind11;

#define BIND_STR_REPR(cls) .def("__str__", &cls::__str__).def("__repr__", &cls::__repr__)

namespace spindynapy {

struct CartesianCoordSystem {
    const char *name = "CartesianCoordSystem";
};
struct SphericalCoordSystem {
    const char *name = "SphericalCoordSystem";
};

typedef short regnum;

/**
 * (Интерфейс) Базовая единица абстракции - координаты.
 * Описывает расположение в пространстве объекта.
 */
class ICoordinates {
  protected:
    ICoordinates() = default;

  public:
    virtual ~ICoordinates() {};

    virtual std::string __str__() const { return nullptr; };
    virtual std::string __repr__() const { return nullptr; };
};

/**
 * Декартовые координаты
 */
class CartesianCoordinates : public ICoordinates {
  protected:
    Eigen::Vector3d _coords;

  public:
    CartesianCoordinates(double x, double y, double z) noexcept : _coords(x, y, z) {};
    CartesianCoordinates(const Eigen::Vector3d &coords) noexcept : _coords(coords) {};

    virtual std::string __str__() const override {
        return std::format("{:.3f}\t{:.3f}\t{:.3f}", _coords(0), _coords(1), _coords(2));
    };
    virtual std::string __repr__() const override {
        return std::format(
            "CartesianCoordinates(x={:.3f}, y={:.3f}, z={:.3f})", _coords(0), _coords(1), _coords(2)
        );
    };

    /**
     * КОМПРОМИСС (нарушение полиморфности):
     *      в солверах с темплейтами используется через if constexpr
     *   Иначе - сложно поддерживать множество ветвей наследования.
     *   + невиртуальный - быстрее.
     */
    void setCoordinates(double x, double y, double z) {
        // TODO: рассмотреть возможность setCoordinates(params), а дальше идентифицировать по номерам
        _coords[0] = x;
        _coords[1] = y;
        _coords[2] = z;
        return;
    }
    /**
     * КОМПРОМИСС (нарушение полиморфности):
     *      в солверах с темплейтами используется через if constexpr
     *   Иначе - сложно поддерживать множество ветвей наследования.
     *   + невиртуальный - быстрее.
     */
    void setCoordinates(Eigen::Vector3d coords) {
        _coords[0] = coords[0];
        _coords[1] = coords[1];
        _coords[2] = coords[2];
        return;
    }
};

/**
 * (Интерфейс) Базовая единица абстракции - направление, вектор.
 * Описывает направление сил, движения, etc.
 */
class IDirection {
  protected:
    IDirection() = default;

  public:
    virtual ~IDirection() {};

    virtual std::string __str__() const { return nullptr; };
    virtual std::string __repr__() const { return nullptr; };
};

/**
 * Вектор в декартовых координатах
 */
class CartesianDirection : public IDirection {
  protected:
    Eigen::Vector3d _vector;

  public:
    CartesianDirection(double x, double y, double z) noexcept : _vector(x, y, z) {};
    CartesianDirection(const Eigen::Vector3d &vector) noexcept : _vector(vector) {};

    virtual std::string __str__() const override {
        return std::format("{:.3f}\t{:.3f}\t{:.3f}", _vector(0), _vector(1), _vector(2));
    };
    virtual std::string __repr__() const override {
        return std::format(
            "CartesianDirection(x={:.3f}, y={:.3f}, z={:.3f})", _vector(0), _vector(1), _vector(2)
        );
    };

    /**
     * КОМПРОМИСС (нарушение полиморфности):
     *      в солверах с темплейтами используется через if constexpr
     *   Иначе - сложно поддерживать множество ветвей наследования.
     *   + невиртуальный - быстрее.
     */
    void setDirection(double x, double y, double z) {
        _vector[0] = x;
        _vector[1] = y;
        _vector[2] = z;
        return;
    }
    /**
     * КОМПРОМИСС (нарушение полиморфности):
     *      в солверах с темплейтами используется через if constexpr
     *   Иначе - сложно поддерживать множество ветвей наследования.
     *   + невиртуальный - быстрее.
     */
    void setDirection(Eigen::Vector3d vector) {
        _vector[0] = vector[0];
        _vector[1] = vector[1];
        _vector[2] = vector[2];
        return;
    }
};

/**
 * (Интерфейс) Структурная базовая единица абстракции - момент.
 * Магнитный, квантово-механический, etc.
 * ХРАНИТ В СЕБЕ КООРДИНАТЫ И ВЕКТОР.
 */
class IMoment {
  protected:
    IMoment() = default;

  public:
    virtual ~IMoment() {};

    /**
     * Получить вектор момента
     */
    virtual IDirection &getDirection() { throw std::logic_error("Not implemented"); };

    /**
     * Получить расположение момента
     */
    virtual ICoordinates &getCoordinates() { throw std::logic_error("Not implemented"); };

    virtual std::string __str__() const { return nullptr; };
    virtual std::string __repr__() const { return nullptr; };
};

/**
 * Момент в декартовых координатах.
 * ХРАНИТ В СЕБЕ КООРДИНАТЫ И ВЕКТОР.
 */
class CartesianMoment : virtual public IMoment {
  protected:
    std::unique_ptr<CartesianCoordinates> _coordinates;
    std::unique_ptr<CartesianDirection> _direction;
    regnum _material_number;

  public:
    CartesianMoment(
        const CartesianCoordinates &coordinates, const CartesianDirection &direction, regnum material_number
    ) noexcept
        : _direction(std::make_unique<CartesianDirection>(direction)),
          _coordinates(std::make_unique<CartesianCoordinates>(coordinates)),
          _material_number(material_number) {};

    CartesianMoment(
        CartesianCoordinates &&coordinates, CartesianDirection &&direction, regnum material_number
    ) noexcept
        : _direction(std::make_unique<CartesianDirection>(std::move(direction))),
          _coordinates(std::make_unique<CartesianCoordinates>(std::move(coordinates))),
          _material_number(material_number) {};

    /**
     * Получить вектор момента
     */
    virtual CartesianDirection &getDirection() override { return *_direction; };

    /**
     * Получить расположение момента
     */
    virtual CartesianCoordinates &getCoordinates() override { return *_coordinates; };

    virtual std::string __str__() const override {
        return _direction->__str__() + "\t" + _coordinates->__str__() + "\t" +
               std::to_string(_material_number);
    };
    virtual std::string __repr__() const override {
        return std::format(
            "CartesianMoment(coordinates={}, direction={}, material_number={})",
            _coordinates->__repr__(),
            _direction->__repr__(),
            _material_number
        );
    };
};

/**
 * (Интерфейс) Структурная единица абстракции
 * Спиновый момент. ХРАНИТ В СЕБЕ КООРДИНАТЫ И ВЕКТОР.
 */
class ISpin : virtual public IMoment {
  protected:
    ISpin() = default;

  public:
    virtual ~ISpin() = default;
};

/**
 * Спиновый момент в декартовых координатах.
 * ХРАНИТ В СЕБЕ КООРДИНАТЫ И ВЕКТОР.
 */
class CartesianSpin : public CartesianMoment, public ISpin { // CartesianSpin + ISpin diamond
  public:
    CartesianSpin(
        const CartesianCoordinates &coordinates, const CartesianDirection &direction, regnum material_number
    ) noexcept
        : CartesianMoment(coordinates, direction, material_number) {};

    CartesianSpin(
        CartesianCoordinates &&coordinates, CartesianDirection &&direction, regnum material_number
    ) noexcept
        : CartesianMoment(coordinates, direction, material_number) {};
};

/**
 * (Интерфейс) Базовая единица абстракции - свойства материала.
 * Хранится в регистре, в остальных - в виде ссылки, добываемой из регистра.
 */
class IMaterial {
  protected:
    IMaterial() = default;

  public:
    virtual ~IMaterial() = default;

    virtual std::string __str__() const { return nullptr; };
    virtual std::string __repr__() const { return nullptr; };
};

/**
 * Свойства магнитного материала (ферромагнетиков, etc.)
 */
class MagneticMaterial : public IMaterial {
  protected:
    double exchange_constant_J;

  public:
    MagneticMaterial(double _exchange_constant_J) : exchange_constant_J(_exchange_constant_J) {};

    virtual std::string __str__() const override {
        return std::format("(exchange: {:.3f})", this->exchange_constant_J);
    };
    virtual std::string __repr__() const override {
        return std::format("MagneticMaterial(exchange_constant_J={:.3f}", this->exchange_constant_J);
    };
};

class IRegion {
  protected:
    IRegion() = default;

  public:
    virtual ~IRegion() = default;
};

class AccuracyRegion : public IRegion {
  protected:
    double _accuracy;

  public:
    AccuracyRegion(double accuracy) : _accuracy(accuracy) {};
};

}; // namespace spindynapy

inline void pyBindTypes(py::module_ &module) {
    using namespace spindynapy;

    // clang-format off

    module.doc() = "Модуль предоставляет базовые абстракные классы и интерфейсы \n"
                   "  - базовые единицы абстракции \n"
                   "  - миксины \n\n"
                   "и их реализации";

    py::class_<CartesianCoordSystem, std::shared_ptr<CartesianCoordSystem>>(module, "CartesianCoordSystem")
        .def(py::init<>())
        .def_readonly("name", &CartesianCoordSystem::name);

    py::class_<SphericalCoordSystem, std::shared_ptr<SphericalCoordSystem>>(module, "SphericalCoordSystem")
        .def(py::init<>())
        .def_readonly("name", &SphericalCoordSystem::name);

    py::class_<IMaterial, std::shared_ptr<IMaterial>>(module, "IMaterial")
        BIND_STR_REPR(IMaterial)
        .doc() = "(Интерфейс) Базовая единица абстракции - свойства материала.\n"
                 "Хранится в регистре, в остальных - в виде ссылки, добываемой из регистра.\n";

    py::class_<ICoordinates>(module, "ICoordinates")
        BIND_STR_REPR(ICoordinates)
        .doc() = "(Интерфейс) Базовая единица абстракции - координаты.\n"
                 "Описывает расположение в пространстве объекта.";

    py::class_<IDirection>(module, "IDirection")
        BIND_STR_REPR(IDirection)
        .doc() = "(Интерфейс) Базовая единица абстракции - направление, вектор.\n"
                 "Описывает направление сил, движения, etc.";

    py::class_<IMoment, std::shared_ptr<IMoment>>(module, "IMoment") BIND_STR_REPR(IMoment)
        .def(
            "get_direction",
            &IMoment::getDirection,
            py::return_value_policy::reference,
            py::doc("Получить вектор момента")
        )
        .def(
            "get_coordinates",
            &IMoment::getCoordinates,
            py::return_value_policy::reference,
            py::doc("Получить расположение момента")
        )
        .doc() = "(Интерфейс) Структурная базовая единица абстракции - момент.\n"
                 "Магнитный, квантово-механический, etc."
                 "ХРАНИТ В СЕБЕ КООРДИНАТЫ И ВЕКТОР.";

    py::class_<ISpin, IMoment, std::shared_ptr<ISpin>>(module, "ISpin")
        .doc() = "(Интерфейс) Структурная единица абстракции. Спиновый момент.\n"
                 "ХРАНИТ В СЕБЕ КООРДИНАТЫ И ВЕКТОР.";

    py::class_<MagneticMaterial, IMaterial, std::shared_ptr<MagneticMaterial>>(module, "MagneticMaterial")
        .def(py::init<double>(), py::arg("_exchange_constant_J"))
        .doc() = "Свойства магнитного материала (ферромагнетиков, etc.)";

    py::class_<CartesianCoordinates, ICoordinates>(module, "CartesianCoordinates")
        .def(
            py::init<double, double, double>(),
            py::arg("x"), py::arg("y"), py::arg("z")
        )
        .def(py::init<const Eigen::Vector3d &>(), py::arg("coords"))
        .doc() = "Декартовые координаты";

    py::class_<CartesianDirection, IDirection>(module, "CartesianDirection")
        .def(
            py::init<double, double, double>(),
            py::arg("x"), py::arg("y"), py::arg("z")
        )
        .def(py::init<const Eigen::Vector3d &>(), py::arg("vector"))
        .doc() = "TODO";

    py::class_<CartesianMoment, IMoment, std::shared_ptr<CartesianMoment>>(module, "CartesianMoment")
        .def(
            py::init<const CartesianCoordinates &, const CartesianDirection &, regnum>(),
            py::arg("coordinates"), py::arg("direction"), py::arg("material_number")
        )
        .def("get_direction", &CartesianMoment::getDirection, py::return_value_policy::reference)
        .def("get_coordinates", &CartesianMoment::getCoordinates, py::return_value_policy::reference)
        .doc() = "TODO";

    py::class_<CartesianSpin, CartesianMoment, ISpin, std::shared_ptr<CartesianSpin>>(module, "CartesianSpin")
        .def(
            py::init<const CartesianCoordinates &, const CartesianDirection &, regnum>(),
            py::arg("coordinates"), py::arg("direction"), py::arg("material_number")
        )
        .doc() = "Спиновый момент в декартовых координатах. ХРАНИТ В СЕБЕ КООРДИНАТЫ И ВЕКТОР.";

    py::class_<IRegion>(module, "IRegion")
        .doc() = "TODO";

    py::class_<AccuracyRegion, IRegion>(module, "AccuracyRegion")
        .def(py::init<double>(), py::arg("accuracy"))
        .doc() = "TODO";

    // clang-format on
}
#endif // ! __TYPES_HPP__

/*

Реализуем пока что следующие связи:

- Координаты + Направления через наследование от базового интерфейсного абстрактного класса,
    причём нарушая полную полиморфность в качестве компромисса
    (т.е. у координат и направлений по-разному интерпретируются методы get/set)
  - через шаблоны не получится, т.к. РАЗНЫЕ КОНСТРУКТОРЫ

- Геометрия реализует:
  - зависимость от системы координат (делаем через обычное наследование геометрии)
  - простая VS сложная (сложная геометрия НЕ зависит от типа вложенных в неё других геометрий),
    т.е. она просто комплексная - сшивалка
  - через шаблоны не получится, т.к. РАЗНЫЕ КОНСТРУКТОРЫ

- Решатели, солверы
  - реализуют собственное наследование по типу солвера ("сделай шаг" - а какой уже солвер решает)
  - зависят от системы координат (делаем через шаблонные методы и if constexpr elsif else throw exception)
  - в питон будут всё равно ретранслироваться по типу CartesianLLGSolver...

*/
