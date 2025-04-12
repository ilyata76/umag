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
struct CartesianCoordinates : public ICoordinates {
    Eigen::Vector3d _coords;

    CartesianCoordinates(double x, double y, double z) : _coords(x, y, z) {};
    CartesianCoordinates(const Eigen::Vector3d &coords) : _coords(coords) {};

    virtual std::string __str__() const override;
    virtual std::string __repr__() const override;
};

/**
 * Вектор в декартовых координатах
 */
struct CartesianDirection : public IDirection {
    Eigen::Vector3d _vector;

    CartesianDirection(double x, double y, double z) : _vector(x, y, z) {};
    CartesianDirection(const Eigen::Vector3d &vector) : _vector(vector) {};

    virtual std::string __str__() const override;
    virtual std::string __repr__() const override;
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
    CartesianMoment(const CartesianCoordinates &coordinates, const CartesianDirection &direction)
        : _direction(std::make_unique<CartesianDirection>(direction)),
          _coordinates(std::make_unique<CartesianCoordinates>(coordinates)) {};

    CartesianMoment(CartesianCoordinates &&coordinates, CartesianDirection &&direction)
        : _direction(std::make_unique<CartesianDirection>(std::move(direction))),
          _coordinates(std::make_unique<CartesianCoordinates>(std::move(coordinates))) {};

    virtual std::string __str__() const override;
    virtual std::string __repr__() const override;

    /**
     * Получить вектор момента
     */
    virtual CartesianDirection &getDirection() override;

    /**
     * Получить расположение момента
     */
    virtual CartesianCoordinates &getCoordinates() override;
};

/**
 * Спиновый момент в декартовых координатах.
 * ХРАНИТ В СЕБЕ КООРДИНАТЫ И ВЕКТОР.
 */
class CartesianSpin : public CartesianMoment, public ISpin { // CartesianSpin + ISpin diamond
  public:
    CartesianSpin(const CartesianCoordinates &_coordinates, const CartesianDirection &_direction)
        : CartesianMoment(_coordinates, _direction) {};

    CartesianSpin(CartesianCoordinates &&_coordinates, CartesianDirection &&_direction)
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
