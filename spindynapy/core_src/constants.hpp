#ifndef __CONSTANTS_HPP__
#define __CONSTANTS_HPP__

/**
 * Константы, используемые в приложении, и их "короткие" версии - алиасы (::sci)
 * Фундаментальные и не только. Просто константы.
 */

#include <pybind11/pybind11.h>

namespace py = pybind11;

#define PYTHON_API [[]]

namespace PYTHON_API spindynapy {
namespace PYTHON_API constants {

PYTHON_API constexpr auto NUMBER_PI = std::numbers::pi_v<double>;
PYTHON_API constexpr auto VACUUM_MAGNETIC_PERMEABILITY = 4.0 * NUMBER_PI * 1.0e-7; // Н/А^2 или Гн/м
PYTHON_API constexpr auto BOHR_MAGNETON = 9.2740100783e-24; // Дж/Тл (или А*м^2)
PYTHON_API constexpr auto PLANCK_CONSTANT = 6.62607015e-34; // Дж/Гц
PYTHON_API constexpr auto REDUCED_PLANCK_CONSTANT = PLANCK_CONSTANT / 2 / NUMBER_PI; // Дж * сек
PYTHON_API constexpr auto FREE_SPIN_LANDE_FACTOR = 2.00231930436256;
PYTHON_API constexpr auto FREE_SPIN_GYROMAGNETIC_RATIO =
    FREE_SPIN_LANDE_FACTOR * BOHR_MAGNETON / REDUCED_PLANCK_CONSTANT; // rad/(s*T) ~ 1.76e11
PYTHON_API constexpr auto BOLTZMANN_CONSTANT = 1.380649e-23;          // Дж/К

namespace PYTHON_API sci {

PYTHON_API constexpr auto mu0 = VACUUM_MAGNETIC_PERMEABILITY;
PYTHON_API constexpr auto pi = NUMBER_PI;
PYTHON_API constexpr auto mu_B = BOHR_MAGNETON;
PYTHON_API constexpr auto h = PLANCK_CONSTANT;
PYTHON_API constexpr auto h_rad = REDUCED_PLANCK_CONSTANT;
PYTHON_API constexpr auto fe_lande = FREE_SPIN_LANDE_FACTOR;
PYTHON_API constexpr auto fe_gamma = FREE_SPIN_GYROMAGNETIC_RATIO;
PYTHON_API constexpr auto k_b = BOLTZMANN_CONSTANT;

}; // namespace PYTHON_API sci

}; // namespace PYTHON_API constants
}; // namespace PYTHON_API spindynapy

// Привязка констант к Python
inline void pyBindConstants(py::module_ &module) {
    using namespace spindynapy::constants;

    py::module_ constants_module = module.def_submodule("constants");

    constants_module.doc() =
        "Модуль, отвечающий за предоставление фундаментальных констант \n"
        "и любых других, используемых в приложении. \n"
        "Здесь - константы через полные названия. \n"
        "В .sci предоставлены \"научные\" короткие версии наименований констант (mu0, pi, etc...)";

    constants_module.attr("VACUUM_MAGNETIC_PERMEABILITY") = VACUUM_MAGNETIC_PERMEABILITY;
    constants_module.attr("NUMBER_PI") = NUMBER_PI;
    constants_module.attr("BOHR_MAGNETON") = BOHR_MAGNETON;
    constants_module.attr("PLANCK_CONSTANT") = PLANCK_CONSTANT;
    constants_module.attr("REDUCED_PLANCK_CONSTANT") = REDUCED_PLANCK_CONSTANT;
    constants_module.attr("FREE_SPIN_LANDE_FACTOR") = FREE_SPIN_LANDE_FACTOR;
    constants_module.attr("FREE_SPIN_GYROMAGNETIC_RATIO") = FREE_SPIN_GYROMAGNETIC_RATIO;
    constants_module.attr("BOLTZMANN_CONSTANT") = BOLTZMANN_CONSTANT;

    constants_module.doc() =
        "Модуль, отвечающий за предоставление фундаментальных констант \n"
        "и любых других, используемых в приложении. \n"
        "Здесь предоставлены \"научные\" короткие версии наименований констант (mu0, pi, etc...)";

    constants_module.attr("mu0") = sci::mu0;
    constants_module.attr("pi") = sci::pi;
    constants_module.attr("mu_B") = sci::mu_B;
    constants_module.attr("h") = sci::h;
    constants_module.attr("h_rad") = sci::h_rad;
    constants_module.attr("fe_lande") = sci::fe_lande;
    constants_module.attr("fe_gamma") = sci::fe_gamma;
    constants_module.attr("k_b") = sci::k_b;
}

#endif // ! __CONSTANTS_HPP__
