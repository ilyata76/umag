#include "cartesian.hpp"

/**
 * Реализации методов классов и функций из cartesian.hpp
 * (единицы абстракции в декартовой системе координат)
 */

namespace spindynapy {

// Координаты

std::string CartesianCoordinates::__str__() const {
    return std::format("({:.3f}, {:.3f}, {:.3f})", this->_coords(0), this->_coords(1), this->_coords(2));
};

std::string CartesianCoordinates::__repr__() const {
    return std::format(
        "CartesianCoordinates(x={:.3f}, y={:.3f}, z={:.3f})", this->_coords(0), this->_coords(1), this->_coords(2)
    );
};

// Вектор, направление

std::string CartesianDirection::__str__() const {
    return std::format("({:.3f}, {:.3f}, {:.3f})", this->_vector(0), this->_vector(1), this->_vector(2));
}

std::string CartesianDirection::__repr__() const {
    return std::format(
        "CartesianDirection(sx={:.3f}, sy={:.3f}, sz={:.3f})", this->_vector(0), this->_vector(1), this->_vector(2)
    );
}

// Спин

std::string CartesianMoment::__str__() const {
    return this->_direction->__str__() + " v:" + this->_coordinates->__str__();
}

std::string CartesianMoment::__repr__() const {
    return this->_direction->__str__() + " v:" + this->_coordinates->__str__();
}

CartesianDirection &CartesianMoment::getDirection() {
    return *this->_direction;
};

CartesianCoordinates &CartesianMoment::getCoordinates() {
    return *this->_coordinates;
};

} // namespace spindynapy
