#ifndef __SIMULATION_HPP__
#define __SIMULATION_HPP__

/**
 * Интерфейсы и классы, отвечающие за главный функционал - симуляцию,
 * управление симуляцией.
 */

#include "geometries/base.hpp"
#include "types/base.hpp"

#include <stdexcept>

namespace spindynapy {

/**
 * Базовый интерфейс управлятора симуляцией. Вход в программу.
 *
 * Через управлятор, содержащий в себе всю информацию о проводимом
 * эксперименте (богатое состояние), можно развивать систему согласно заданным настройкам,
 * а также извне, вызывая специальные для этого методы (API симулятора)
 */
class ISimulation {
  public:
    ISimulation() = default;
    virtual ~ISimulation() {};

    virtual std::string __str__() const { return nullptr; };
    virtual std::string __repr__() const { return nullptr; };
};

class Simulation : public ISimulation {
    std::shared_ptr<IGeometry> _geometry;

  public:
    Simulation(std::shared_ptr<IGeometry> geometry) : _geometry(geometry) {
        if (!geometry) {
            throw std::invalid_argument("Геометрия не может быть None");
        }
    };

    virtual std::string __str__() const override { return _geometry->__str__(); };
    virtual std::string __repr__() const override { return _geometry->__repr__(); };
};

}; // namespace spindynapy

namespace spindynapy::doc {

constexpr char module_simulation[] = "Интерфейсы и классы, отвечающие за главный функционал - симуляцию.";

constexpr char ISimulation[] =
    "Базовый интерфейс управлятора симуляцией. Вход в программу.\n\n"
    "Через управлятор, содержащий в себе всю информацию о проводимом\n"
    "эксперименте (богатое состояние), можно развивать систему согласно заданным настройкам,\n"
    "а также извне, вызывая специальные для этого методы (API симулятора)\n";

constexpr char Simulation[] = "TODO";

}; // namespace spindynapy::doc

#endif // ! __SIMULATION_HPP__
