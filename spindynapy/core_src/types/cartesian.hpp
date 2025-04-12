#ifndef __TYPES_CARTESIAN_HPP__
#define __TYPES_CARTESIAN_HPP__

/**
 * Реализации интерфейсов базовых единиц абстракций из base.hpp,
 * зависящих от системы координат, в декартовой системе координат
 */

#include "base.hpp"

#include <pybind11/eigen.h>

namespace spindynapy {

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
        return std::format("({:.3f}, {:.3f}, {:.3f})", _coords(0), _coords(1), _coords(2));
    };
    virtual std::string __repr__() const override {
        return std::format("CartesianCoordinates(x={:.3f}, y={:.3f}, z={:.3f})", _coords(0), _coords(1), _coords(2));
    };
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
        return std::format("({:.3f}, {:.3f}, {:.3f})", _vector(0), _vector(1), _vector(2));
    };
    virtual std::string __repr__() const override {
        return std::format("CartesianDirection(sx={:.3f}, sy={:.3f}, sz={:.3f})", _vector(0), _vector(1), _vector(2));
    };
};

/**
 * Момент в декартовых координатах.
 * ХРАНИТ В СЕБЕ КООРДИНАТЫ И ВЕКТОР.
 */
class CartesianMoment : virtual public IMoment {
  protected:
    std::unique_ptr<CartesianDirection> _direction;
    std::unique_ptr<CartesianCoordinates> _coordinates;

  public:
    CartesianMoment(const CartesianCoordinates &coordinates, const CartesianDirection &direction) noexcept
        : _direction(std::make_unique<CartesianDirection>(direction)),
          _coordinates(std::make_unique<CartesianCoordinates>(coordinates)) {};

    CartesianMoment(CartesianCoordinates &&coordinates, CartesianDirection &&direction) noexcept
        : _direction(std::make_unique<CartesianDirection>(std::move(direction))),
          _coordinates(std::make_unique<CartesianCoordinates>(std::move(coordinates))) {};

    /**
     * Получить вектор момента
     */
    virtual CartesianDirection &getDirection() override { return *_direction; };

    /**
     * Получить расположение момента
     */
    virtual CartesianCoordinates &getCoordinates() override { return *_coordinates; };

    virtual std::string __str__() const override { return _direction->__str__() + " v:" + _coordinates->__str__(); };
    virtual std::string __repr__() const override { return _direction->__str__() + " v:" + _coordinates->__str__(); };
};

/**
 * Спиновый момент в декартовых координатах.
 * ХРАНИТ В СЕБЕ КООРДИНАТЫ И ВЕКТОР.
 */
class CartesianSpin : public CartesianMoment, public ISpin { // CartesianSpin + ISpin diamond
  public:
    CartesianSpin(const CartesianCoordinates &_coordinates, const CartesianDirection &_direction) noexcept
        : CartesianMoment(_coordinates, _direction) {};

    CartesianSpin(CartesianCoordinates &&_coordinates, CartesianDirection &&_direction) noexcept
        : CartesianMoment(_coordinates, _direction) {};
};

}; // namespace spindynapy

namespace spindynapy::doc {

constexpr char module_types_cartesian[] = "Модуль предоставляет реализации базовых абстракных единиц \n"
                                          ", описанные в контексте декартовой системы координат";

constexpr char CartesianCoordinates[] = "Декартовые координаты";
constexpr char CartesianDirection[] = "Вектор в декартовых координатах";
constexpr char CartesianMoment[] = "Момент в декартовых координатах. ХРАНИТ В СЕБЕ КООРДИНАТЫ И ВЕКТОР.";
constexpr char CartesianSpin[] = "Спиновый момент в декартовых координатах. ХРАНИТ В СЕБЕ КООРДИНАТЫ И ВЕКТОР.";

}; // namespace spindynapy::doc

#endif // ! __TYPES_CARTESIAN_HPP__
