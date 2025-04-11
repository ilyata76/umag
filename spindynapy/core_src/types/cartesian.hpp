#ifndef __TYPES_CARTESIAN_HPP__
#define __TYPES_CARTESIAN_HPP__

/**
 * Реализации интерфейсов базовых единиц абстракций из base.hpp,
 * зависящих от системы координат, в декартовой системе координат
 */

#include "base.hpp"

namespace spindynapy {

/**
 * Декартовые координаты
 */
struct CartesianCoordinates : public ICoordinates {
    double x, y, z;

    CartesianCoordinates(double _x, double _y, double _z) : x(_x), y(_y), z(_z) {};

    virtual std::string __str__() const override;
    virtual std::string __repr__() const override;
};

/**
 * Вектор в декартовых координатах
 */
struct CartesianDirection : public IDirection {
    double sx, sy, sz;

    CartesianDirection(double _sx, double _sy, double _sz) : sx(_sx), sy(_sy), sz(_sz) {};

    virtual std::string __str__() const override;
    virtual std::string __repr__() const override;
};

/**
 * Момент в декартовых координатах.
 * ХРАНИТ В СЕБЕ КООРДИНАТЫ И ВЕКТОР.
 */
class CartesianMoment : virtual public IMoment {
  protected:
    std::unique_ptr<CartesianDirection> direction;
    std::unique_ptr<CartesianCoordinates> coordinates;

  public:
    CartesianMoment(CartesianCoordinates &_coordinates, CartesianDirection &_direction)
        : direction(std::make_unique<CartesianDirection>(_direction)),
          coordinates(std::make_unique<CartesianCoordinates>(_coordinates)) {};

    CartesianMoment(CartesianCoordinates &&_coordinates, CartesianDirection &&_direction)
        : direction(std::make_unique<CartesianDirection>(std::move(_direction))),
          coordinates(std::make_unique<CartesianCoordinates>(std::move(_coordinates))) {};

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
    CartesianSpin(CartesianCoordinates &_coordinates, CartesianDirection &_direction)
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
