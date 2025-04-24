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

#include "constants.hpp"

#include <format>
#include <memory>
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>
#include <stdexcept>
#include <string>

namespace py = pybind11;

#define BIND_STR_REPR(cls) .def("__str__", &cls::__str__).def("__repr__", &cls::__repr__)
// TODO переделать через ostream и lambda: ostringstream.str()

namespace spindynapy {

typedef short regnum;

// Определим тип для эффективного поля
using EffectiveField = Eigen::Vector3d;

// Базовый класс для анизотропии (НЕ понятно: как делить анизотропию по CoordSystem?)
class Anisotropy {
  protected:
    Anisotropy() = default;

  public:
    virtual ~Anisotropy() = default;
    virtual std::string __str__() const { throw std::logic_error("Method __str__ Not implemented"); };
};

class UniaxialAnisotropy : public Anisotropy {
  public:
    Eigen::Vector3d axis; // Ось анизотропии (нормализованная)
    double constant;      // Константа анизотропии K (Дж/м³)

    UniaxialAnisotropy(const Eigen::Vector3d &axis, double constant)
        : axis(axis.normalized()), constant(constant) {}

    std::string __str__() const override {
        return std::format(
            "UniaxialAnisotropy(axis=({}, {}, {}), constant={:.2e})", axis.x(), axis.y(), axis.z(), constant
        );
    }
};

/**
 * (Интерфейс) Базовая единица абстракции - свойства материала.
 * Хранится в регистре, в остальных - в виде ссылки, добываемой из регистра.

 * компромисс: полноценный класс с набором параметров.
 */
class Material {
  public:
    regnum material_number;     // порядковый номер
    double exchange_constant_J; // обменный интеграл
    double atomic_magnetic_saturation_magnetization; // намагниченность насыщения в магнетонах бора
    std::shared_ptr<Anisotropy> anisotropy; // присущая материалу магнитная анизотропия
    double gyromagnetic_ratio; // гиромагнитное отношение
    double damping_constant;   // демпфирующая константа Гильберта

    /** Гиромагнитное отношение зависит от эффективного фактора Ланде */
    Material(
        regnum material_number,
        double exchange_constant_J,
        double atomic_magnetic_saturation_magnetization,
        double gyromagnetic_ratio = constants::FREE_SPIN_GYROMAGNETIC_RATIO, // ~=~ Co Ni Fe
        double damping_constant = 0.1,
        std::shared_ptr<Anisotropy> anisotropy = nullptr
    )
        : material_number(material_number),
          exchange_constant_J(exchange_constant_J),
          atomic_magnetic_saturation_magnetization(atomic_magnetic_saturation_magnetization),
          damping_constant(damping_constant),
          gyromagnetic_ratio(gyromagnetic_ratio),
          anisotropy(anisotropy) {};

    regnum getNumber() const { return this->material_number; };
    std::string __str__() const { return std::to_string(this->material_number); };
    std::string __repr__() const {
        return std::format(
            "Material(material_number={}, exchange_constant_J={}, "
            "atomic_magnetic_saturation_magnetization={})",
            this->material_number,
            this->exchange_constant_J,
            this->atomic_magnetic_saturation_magnetization
        );
    }

    bool operator==(const Material &other) const { return material_number == other.material_number; }
};

class MaterialInterface {};

class Region {
  protected:
    double _accuracy;

  public:
    Region(double accuracy) : _accuracy(accuracy) {};
};

namespace cartesian {

/**
 * Декартовые координаты
 */
class Coordinates {
  protected:
    Eigen::Vector3d _coords;

  public:
    Coordinates(double x, double y, double z) noexcept : _coords(x, y, z) {};
    Coordinates(const Eigen::Vector3d &coords) noexcept : _coords(coords) {};

    std::string __str__() const {
        return std::format("{:.3f}\t{:.3f}\t{:.3f}", _coords(0) * 1e10, _coords(1) * 1e10, _coords(2) * 1e10);
    };
    std::string __repr__() const {
        return std::format("Coordinates(x={:.3f}, y={:.3f}, z={:.3f})", _coords(0), _coords(1), _coords(2));
    };

    double getDistanceFrom(const Coordinates &other) const { return (this->_coords - other._coords).norm(); };

    double getDistanceSquareFrom(const Coordinates &other) const {
        return (this->_coords - other._coords).squaredNorm();
    };

    void setCoordinates(double x, double y, double z) {
        // TODO: рассмотреть возможность setCoordinates(params), а дальше идентифицировать по номерам
        _coords[0] = x;
        _coords[1] = y;
        _coords[2] = z;
        return;
    }

    Eigen::Vector3d &asVector() { return this->_coords; }

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
class Direction {
  protected:
    Eigen::Vector3d _vector;

  public:
    Direction(double x, double y, double z) noexcept : _vector(x, y, z) { _vector.normalize(); };
    Direction(const Eigen::Vector3d &vector) noexcept : _vector(vector) { _vector.normalize(); };

    std::string __str__() const {
        return std::format("{:.3f}\t{:.3f}\t{:.3f}", _vector(0), _vector(1), _vector(2));
    };
    std::string __repr__() const {
        return std::format("Direction(x={:.3f}, y={:.3f}, z={:.3f})", _vector(0), _vector(1), _vector(2));
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

    Eigen::Vector3d &asVector() { return this->_vector; }

    void normalize() { return _vector.normalize(); }
};

class Moment {
  protected:
    std::unique_ptr<Coordinates> _coordinates;
    std::unique_ptr<Direction> _direction;
    std::shared_ptr<Material> _material;

  public:
    Moment(
        const Coordinates &coordinates, const Direction &direction, std::shared_ptr<Material> material
    ) noexcept
        : _coordinates(std::make_unique<Coordinates>(coordinates)),
          _direction(std::make_unique<Direction>(direction)),
          _material(material) {};

    Moment(Coordinates &&coordinates, Direction &&direction, std::shared_ptr<Material> material) noexcept
        : _coordinates(std::make_unique<Coordinates>(std::move(coordinates))),
          _direction(std::make_unique<Direction>(std::move(direction))),
          _material(std::move(material)) {};

    /**
     * Получить вектор момента
     */
    Direction &getDirection() { return *_direction; };

    void setDirection(double sx, double sy, double sz) { this->_direction->setDirection(sx, sy, sz); }

    void setDirection(Eigen::Vector3d new_direction_vector) {
        this->_direction->setDirection(new_direction_vector);
    }

    /**
     * Получить расположение момента
     */
    Coordinates &getCoordinates() { return *_coordinates; };

    /**
     * Получить ссылку на метериал
     */
    Material &getMaterial() { return *_material; };

    std::shared_ptr<Moment> clone() const {
        return std::make_shared<Moment>(*_coordinates, *_direction, _material);
    }

    std::string __str__() const {
        return _material->__str__() + "\t0\t" + _coordinates->__str__() + "\t" + _direction->__str__();
    };
    std::string __repr__() const {
        return std::format(
            "Moment(coordinates={}, direction={}, material={})",
            _coordinates->__repr__(),
            _direction->__repr__(),
            _material->__repr__()
        );
    };
};

} // namespace cartesian

namespace spherical {

class Moment {};
class Coordinates {};
class Direction {};

} // namespace spherical

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

struct SphericalCoordSystem : public CoordSystem {
    using Moment = spherical::Moment;
    using Coordinates = spherical::Coordinates;
    using Direction = spherical::Direction;
    const char *name = "SphericalCoordSystem";
};

struct CartesianCoordSystem : public CoordSystem {
    using Moment = cartesian::Moment;
    using Coordinates = cartesian::Coordinates;
    using Direction = cartesian::Direction;
    const char *name = "CoordSystem";
};

namespace cartesian {

using NamespaceCoordSystem = CartesianCoordSystem;

};

namespace spherical {

using NamespaceCoordSystem = SphericalCoordSystem;

};

}; // namespace spindynapy

inline void pyBindTypes(py::module_ &module) {
    using namespace spindynapy;

    // -------- | TYPES | --------
    py::module_ types_module = module.def_submodule("types");
    types_module.doc() = "Модуль предоставляет базовые абстракные классы и интерфейсы \n"
                         "  - базовые единицы абстракции \n"
                         "  - миксины \n\n"
                         "и их реализации";

    py::class_<Anisotropy, std::shared_ptr<Anisotropy>>(types_module, "Anisotropy")
        .def("__str__", &Anisotropy::__str__)
        .doc() = "TODO";

    py::class_<UniaxialAnisotropy, Anisotropy, std::shared_ptr<UniaxialAnisotropy>>(
        types_module, "UniaxialAnisotropy"
    )
        .def(py::init<const Eigen::Vector3d &, double>(), py::arg("axis"), py::arg("constant"))
        .def_readwrite("axis", &UniaxialAnisotropy::axis)
        .def_readwrite("constant", &UniaxialAnisotropy::constant)
        .doc() = "Uniaxial anisotropy data";

    py::class_<Material, std::shared_ptr<Material>>(types_module, "Material") BIND_STR_REPR(Material)
        .def(
            py::init<regnum, double, double, double, double, std::shared_ptr<Anisotropy>>(),
            py::arg("material_number"),
            py::arg("exchange_constant_J"),
            py::arg("atomic_magnetic_saturation_magnetization"),
            py::arg("gyromagnetic_ratio") = constants::FREE_SPIN_GYROMAGNETIC_RATIO,
            py::arg("damping_constant") = 0.1,
            py::arg("anisotropy").none(true) = pybind11::none()
        )
        .def("get_number", &Material::getNumber)
        .doc() = "(Интерфейс) Базовая единица абстракции - свойства материала.\n"
                 "Хранится в регистре, в остальных - в виде ссылки, добываемой из регистра.\n";

    py::class_<Region>(module, "Region").def(py::init<double>(), py::arg("accuracy")).doc() = "TODO";

    py::class_<CartesianCoordSystem, std::shared_ptr<CartesianCoordSystem>>(module, "CoordSystem")
        .def(py::init<>())
        .def_readonly("name", &CartesianCoordSystem::name);

    py::class_<SphericalCoordSystem, std::shared_ptr<SphericalCoordSystem>>(module, "SphericalCoordSystem")
        .def(py::init<>())
        .def_readonly("name", &SphericalCoordSystem::name);

    // -------- | CARTESIAN TYPES | --------
    py::module_ cartesian = types_module.def_submodule("cartesian");

    using cartesian::Coordinates;
    using cartesian::Direction;
    using cartesian::Moment;

    py::class_<Coordinates>(cartesian, "Coordinates")
        .def(py::init<double, double, double>(), py::arg("x"), py::arg("y"), py::arg("z"))
        .def(py::init<const Eigen::Vector3d &>(), py::arg("coords"))
        .def("get_distance_from", &Coordinates::getDistanceFrom, py::arg("other"))
        .def("get_distance_square_from", &Coordinates::getDistanceSquareFrom, py::arg("other"))
            BIND_STR_REPR(Coordinates)
        .doc() = "Декартовые координаты";

    py::class_<Direction>(cartesian, "Direction")
        .def(py::init<double, double, double>(), py::arg("x"), py::arg("y"), py::arg("z"))
        .def(py::init<const Eigen::Vector3d &>(), py::arg("vector")) BIND_STR_REPR(Direction)
        .doc() = "Направление в декартовых координатах";

    py::class_<Moment, std::shared_ptr<Moment>>(cartesian, "Moment")
        .def(
            py::init([](const Coordinates &coords, const Direction &dir, Material &mat) {
                return new Moment(coords, dir, std::make_shared<Material>(mat));
            }),
            py::arg("coordinates"),
            py::arg("direction"),
            py::arg("material")
        ) BIND_STR_REPR(Moment)
        .def("get_material", &Moment::getMaterial, py::return_value_policy::reference)
        .def("get_coordinates", &Moment::getCoordinates, py::return_value_policy::reference)
        .def("get_direction", &Moment::getDirection, py::return_value_policy::reference)
        .doc() = "TODO";

    // -------- | SPHERICAL TYPES | --------
    py::module_ spherical = types_module.def_submodule("spherical");
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
  - в питон будут всё равно ретранслироваться по типу LLGSolver...

*/
