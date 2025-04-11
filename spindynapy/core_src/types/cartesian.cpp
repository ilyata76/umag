#include "cartesian.hpp"

/**
 * Реализации методов классов и функций из cartesian.hpp
 * (единицы абстракции в декартовой системе координат)
 */

namespace spindynapy {

// Координаты

std::string CartesianCoordinates::__str__() const {
    return std::format("({:.3f}, {:.3f}, {:.3f})", this->x, this->y, this->z);
};

std::string CartesianCoordinates::__repr__() const {
    return std::format("CartesianCoordinates(x={:.3f}, y={:.3f}, z={:.3f})", this->x, this->y, this->z);
};

// Вектор, направление

std::string CartesianDirection::__str__() const {
    return std::format("({:.3f}, {:.3f}, {:.3f})", this->sx, this->sy, this->sz);
}

std::string CartesianDirection::__repr__() const {
    return std::format("CartesianDirection(sx={:.3f}, sy={:.3f}, sz={:.3f})", this->sx, this->sy, this->sz);
}

// Спин

std::string CartesianMoment::__str__() const {
    return this->direction->__str__() + " v:" + this->coordinates->__str__();
}

std::string CartesianMoment::__repr__() const {
    return this->direction->__str__() + " v:" + this->coordinates->__str__();
}

CartesianDirection &CartesianMoment::getDirection() {
    return *this->direction;
};

CartesianCoordinates &CartesianMoment::getCoordinates() {
    return *this->coordinates;
};

} // namespace spindynapy
