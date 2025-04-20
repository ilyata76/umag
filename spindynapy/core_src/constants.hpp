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

namespace sci {

constexpr auto mu0 = VACUUM_MAGNETIC_PERMEABILITY;
constexpr auto pi = NUMBER_PI;
constexpr auto mu_B = BOHR_MAGNETON;

}; // namespace sci

}; // namespace spindynapy::constants

inline void pyBindConstants(py::module_ &module) {
    using namespace spindynapy::constants;

    // clang-format off

    module.doc() = "Модуль, отвечающий за предоставление фундаментальных констант \n"
                   "и любых других, используемых в приложении. \n"
                   "Здесь - константы через полные названия. \n"
                   "В .sci предоставлены \"научные\" короткие версии наименований констант (mu0, pi, etc...)";

    module.attr("VACUUM_MAGNETIC_PERMEABILITY") = VACUUM_MAGNETIC_PERMEABILITY;
    module.attr("NUMBER_PI") = NUMBER_PI;
    module.attr("BOHR_MAGNETON") = BOHR_MAGNETON;

    module.doc() = "Модуль, отвечающий за предоставление фундаментальных констант \n"
                   "и любых других, используемых в приложении. \n"
                   "Здесь предоставлены \"научные\" короткие версии наименований констант (mu0, pi, etc...)";

    module.attr("mu0") = sci::mu0;
    module.attr("pi") = sci::pi;
    module.attr("mu_B") = sci::mu_B;

    // clang-format on
}

#endif // ! __CONSTANTS_HPP__
