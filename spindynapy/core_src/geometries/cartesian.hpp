#ifndef __GEOMETRIES_CARTESIAN_HPP__
#define __GEOMETRIES_CARTESIAN_HPP__

/**
 * Геометрии, задаваемые в декартовой системе координат
 */

#include "base.hpp"

namespace spindynapy {

/**
 * Простая геометрия, задаваемая в декартовой системе координат
 */
class CartesianGeometry : public IGeometry {

  public:
    int _x;
    CartesianGeometry(int x) : _x(x) {};

    virtual std::string __str__() const override { return std::to_string(_x); };
    virtual std::string __repr__() const override { return std::to_string(_x); };
};

}; // namespace spindynapy

namespace spindynapy::doc {

constexpr char module_geometries_cartesian[] = "Геометрии, задаваемые в декартовой системе координат";

constexpr char CartesianGeometry[] = "Простая геометрия, задаваемая в декартовой системе координат";

}; // namespace spindynapy::doc

#endif // ! __GEOMETRIES_CARTESIAN_HPP__
