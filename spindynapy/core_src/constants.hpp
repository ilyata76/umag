#ifndef __CONSTANTS_HPP__
#define __CONSTANTS_HPP__

/**
 * Константы, используемые в приложении, и их "короткие" версии - алиасы (::sci)
 * Фундаментальные и не только. Просто константы.
 */

#include <pybind11/pybind11.h>

namespace py = pybind11;

namespace spindynapy::constants {

constexpr auto VACUUM_MAGNETIC_PERMEABILITY = 1.256637061e-6; // Н/А^2 или Гн/м
constexpr auto NUMBER_PI = 3.14159;
constexpr auto BOHR_MAGNETON = 9.2740100783e-24; // Дж/Тл (или А*м^2)
constexpr auto PLANK_CONSTANT = 6.62607015e-34; // Дж/Гц
constexpr auto REDUCED_PLANK_CONSTANT = PLANK_CONSTANT / 2 / NUMBER_PI; // Дж * сек

namespace sci {

constexpr auto mu0 = VACUUM_MAGNETIC_PERMEABILITY;
constexpr auto pi = NUMBER_PI;
constexpr auto mu_B = BOHR_MAGNETON;
constexpr auto h = PLANK_CONSTANT;
constexpr auto h_rad = REDUCED_PLANK_CONSTANT;

}; // namespace sci

}; // namespace spindynapy::constants

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
    constants_module.attr("PLANK_CONSTANT") = PLANK_CONSTANT;
    constants_module.attr("REDUCED_PLANK_CONSTANT") = REDUCED_PLANK_CONSTANT;

    constants_module.doc() =
        "Модуль, отвечающий за предоставление фундаментальных констант \n"
        "и любых других, используемых в приложении. \n"
        "Здесь предоставлены \"научные\" короткие версии наименований констант (mu0, pi, etc...)";

    constants_module.attr("mu0") = sci::mu0;
    constants_module.attr("pi") = sci::pi;
    constants_module.attr("mu_B") = sci::mu_B;
    constants_module.attr("h") = sci::h;
    constants_module.attr("h_rad") = sci::h_rad;
}

#endif // ! __CONSTANTS_HPP__
