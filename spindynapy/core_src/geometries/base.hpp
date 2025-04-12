#ifndef __GEOMETRIES_BASE_HPP__
#define __GEOMETRIES_BASE_HPP__

/**
 * Заголовки с базовыми интерфейсами для геометрий образцов.
 * Геометрии задают устройство образца, сохраняя в себе
 * состояние расположения моментов и их свойств.
 *
 * Геометрия также может поставлять ограниченное количество
 * обсчитываемых параметров. Может состоять из разных участков разной точности.
 */

#include "../types/base.hpp"

namespace spindynapy {

/**
 * Базовый интерфейс геометрии
 */
class IGeometry {
  protected:
    IGeometry() = default;

  public:
    virtual ~IGeometry() = default;

    virtual std::string __str__() const { return nullptr; };
    virtual std::string __repr__() const { return nullptr; };
};

}; // namespace spindynapy

namespace spindynapy::doc {

constexpr char module_geometries[] = "Модуль геометрий\n"
                                     "Геометрии задают устройство образца, сохраняя в себе\n"
                                     "состояние расположения моментов и их свойств.\n"
                                     "Геометрия также может поставлять ограниченное количество\n"
                                     "обсчитываемых параметров. Может состоять из разных участков разной точности.\n";

constexpr const char *module_geometries_base = module_geometries;

constexpr char IGeometry[] = "Базовый интерфейс геометрии";

}; // namespace spindynapy::doc

#endif // ! __GEOMETRIES_BASE_HPP__
