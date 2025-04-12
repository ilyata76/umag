#ifndef __TYPES_MATERIALS_HPP__
#define __TYPES_MATERIALS_HPP__

/**
 * Реализация интерфейсов базовых единиц абстракций из base.hpp,
 * связанных с представлением различных материалов и их свойств.
 */

#include "base.hpp"

namespace spindynapy {

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

}; // namespace spindynapy

namespace spindynapy::doc {

constexpr char module_types_material[] = "Модуль предоставляет реализации базовых абстракных единиц \n"
                                         ", связанных с материалами и их свойствами.";

constexpr char MagneticMaterial[] = "Свойства магнитного материала (ферромагнетиков, etc.)";

}; // namespace spindynapy::doc

#endif // ! __TYPES_MATERIALS_HPP__
