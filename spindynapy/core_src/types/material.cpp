#include "material.hpp"

namespace spindynapy {

// Magnetics

std::string MagneticMaterial::__str__() const {
    return std::format("(exchange: {:.3f})", this->exchange_constant_J);
};

std::string MagneticMaterial::__repr__() const {
    return std::format("MagneticMaterial(exchange_constant_J={:.3f}", this->exchange_constant_J);
};

} // namespace spindynapy
