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
#include <string>

namespace py = pybind11;

#define BIND_STR_REPR(cls) .def("__str__", &cls::__str__).def("__repr__", &cls::__repr__)

namespace spindynapy {

typedef short regnum;

/**
 * (Интерфейс) Базовая единица абстракции - свойства материала.
 * Хранится в регистре, в остальных - в виде ссылки, добываемой из регистра.

 * компромисс: полноценный класс с набором параметров.
 */
class Material {
  protected:
    double _exchange_constant_J;
    regnum _material_number;

  public:
    Material(regnum material_number, double exchange_constant_J)
        : _material_number(material_number), _exchange_constant_J(exchange_constant_J) {};

    regnum getNumber() const { return this->_material_number; };
    std::string __str__() const { return std::to_string(this->_material_number); };
    std::string __repr__() const {
        return std::format(
            "Material(material_number={}, exchange_constant_J={})",
            this->_material_number,
            this->_exchange_constant_J
        );
    }
};

class Region {
  protected:
    double _accuracy;

  public:
    Region(double accuracy) : _accuracy(accuracy) {};
};

// Собираем концепт шаблона:
// классы будут как бы независимые. Склейка по координатным системам только через шаблоны

template <typename CoordSystem>
concept CoordSystemConcept = requires {
    typename CoordSystem::Moment;
    typename CoordSystem::Coordinates;
    typename CoordSystem::Direction;
    { CoordSystem {}.name } -> std::convertible_to<const char *>;
};

struct CoordSystem {};

/**
 * Декартовые координаты
 */
class CartesianCoordinates {
  protected:
    Eigen::Vector3d _coords;

  public:
    CartesianCoordinates(double x, double y, double z) noexcept : _coords(x, y, z) {};
    CartesianCoordinates(const Eigen::Vector3d &coords) noexcept : _coords(coords) {};

    std::string __str__() const {
        return std::format("{:.3f}\t{:.3f}\t{:.3f}", _coords(0), _coords(1), _coords(2));
    };
    std::string __repr__() const {
        return std::format(
            "CartesianCoordinates(x={:.3f}, y={:.3f}, z={:.3f})", _coords(0), _coords(1), _coords(2)
        );
    };

    double getDistanceFrom(const CartesianCoordinates &other) const {
        return (this->_coords - other._coords).norm();
    };

    double getDistanceSquareFrom(const CartesianCoordinates &other) const {
        return (this->_coords - other._coords).squaredNorm();
    };

    void setCoordinates(double x, double y, double z) {
        // TODO: рассмотреть возможность setCoordinates(params), а дальше идентифицировать по номерам
        _coords[0] = x;
        _coords[1] = y;
        _coords[2] = z;
        return;
    }

    void setCoordinates(Eigen::Vector3d coords) {
        _coords[0] = coords[0];
        _coords[1] = coords[1];
        _coords[2] = coords[2];
        return;
    }
};

/**
 * Вектор в декартовых координатах
 */
class CartesianDirection {
  protected:
    Eigen::Vector3d _vector;

  public:
    CartesianDirection(double x, double y, double z) noexcept : _vector(x, y, z) { _vector.normalize(); };
    CartesianDirection(const Eigen::Vector3d &vector) noexcept : _vector(vector) { _vector.normalize(); };

    std::string __str__() const {
        return std::format("{:.3f}\t{:.3f}\t{:.3f}", _vector(0), _vector(1), _vector(2));
    };
    std::string __repr__() const {
        return std::format(
            "CartesianDirection(x={:.3f}, y={:.3f}, z={:.3f})", _vector(0), _vector(1), _vector(2)
        );
    };

    void setDirection(double x, double y, double z) {
        _vector[0] = x;
        _vector[1] = y;
        _vector[2] = z;
        _vector.normalize();
        return;
    }

    void setDirection(Eigen::Vector3d vector) {
        _vector = vector;
        _vector.normalize();
        return;
    }

    void normalize() { return _vector.normalize(); }
};

class CartesianMoment {
  protected:
    std::unique_ptr<CartesianCoordinates> _coordinates;
    std::unique_ptr<CartesianDirection> _direction;
    std::shared_ptr<Material> _material;

  public:
    CartesianMoment(
        const CartesianCoordinates &coordinates,
        const CartesianDirection &direction,
        std::shared_ptr<Material> material
    ) noexcept
        : _coordinates(std::make_unique<CartesianCoordinates>(coordinates)),
          _direction(std::make_unique<CartesianDirection>(direction)),
          _material(material) {};

    CartesianMoment(
        CartesianCoordinates &&coordinates, CartesianDirection &&direction, std::shared_ptr<Material> material
    ) noexcept
        : _coordinates(std::make_unique<CartesianCoordinates>(std::move(coordinates))),
          _direction(std::make_unique<CartesianDirection>(std::move(direction))),
          _material(std::move(material)) {};

    /**
     * Получить вектор момента
     */
    CartesianDirection &getDirection() { return *_direction; };

    /**
     * Получить расположение момента
     */
    CartesianCoordinates &getCoordinates() { return *_coordinates; };

    /**
     * Получить ссылку на метериал
     */
    Material &getMaterial() { return *_material; };

    std::string __str__() const {
        return _coordinates->__str__() + "\t" + _direction->__str__() + "\t" + _material->__str__();
    };
    std::string __repr__() const {
        return std::format(
            "CartesianMoment(coordinates={}, direction={}, material={})",
            _coordinates->__repr__(),
            _direction->__repr__(),
            _material->__repr__()
        );
    };
};

struct CartesianCoordSystem : public CoordSystem {
    using Moment = CartesianMoment;
    using Coordinates = CartesianCoordinates;
    using Direction = CartesianDirection;
    const char *name = "CartesianCoordSystem";
};

class SphericalMoment {};
class SphericalCoordinates {};
class SphericalDirection {};

struct SphericalCoordSystem : public CoordSystem {
    using Moment = SphericalMoment;
    using Coordinates = SphericalCoordinates;
    using Direction = SphericalDirection;
    const char *name = "SphericalCoordSystem";
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

    py::class_<Material, std::shared_ptr<Material>>(module, "Material")
        BIND_STR_REPR(Material)
        .def(py::init<regnum, double>(), py::arg("material_number"), py::arg("exchange_constant_J"))
        .def("get_number", &Material::getNumber)
        .doc() = "(Интерфейс) Базовая единица абстракции - свойства материала.\n"
                 "Хранится в регистре, в остальных - в виде ссылки, добываемой из регистра.\n";

    py::class_<CartesianCoordinates>(module, "CartesianCoordinates")
        .def(
            py::init<double, double, double>(),
            py::arg("x"), py::arg("y"), py::arg("z")
        )
        .def(py::init<const Eigen::Vector3d &>(), py::arg("coords"))
        .def("get_distance_from", &CartesianCoordinates::getDistanceFrom, py::arg("other"))
        .def("get_distance_square_from", &CartesianCoordinates::getDistanceSquareFrom, py::arg("other"))
        BIND_STR_REPR(CartesianCoordinates)
        .doc() = "Декартовые координаты";

    py::class_<CartesianDirection>(module, "CartesianDirection")
        .def(
            py::init<double, double, double>(),
            py::arg("x"), py::arg("y"), py::arg("z")
        )
        .def(py::init<const Eigen::Vector3d &>(), py::arg("vector"))
        BIND_STR_REPR(CartesianDirection)
        .doc() = "Направление в декартовых координатах";

    py::class_<CartesianMoment, std::shared_ptr<CartesianMoment>>(module, "CartesianMoment")
        .def(
            py::init([](const CartesianCoordinates &coords, const CartesianDirection &dir, Material& mat) {
                return new CartesianMoment(coords, dir, std::make_shared<Material>(mat));
            }),
            py::arg("coordinates"), py::arg("direction"), py::arg("material")
        )
        BIND_STR_REPR(CartesianMoment)
        .def("get_material", &CartesianMoment::getMaterial, py::return_value_policy::reference)
        .def("get_coordinates", &CartesianMoment::getCoordinates, py::return_value_policy::reference)
        .def("get_direction", &CartesianMoment::getDirection, py::return_value_policy::reference)
        .doc() = "TODO";

    py::class_<Region>(module, "Region")
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
