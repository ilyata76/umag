#ifndef __CONSTANTS_HPP__
#define __CONSTANTS_HPP__

/**
 * Константы, используемые в приложении, и их "короткие" версии - алиасы (::sci)
 * Фундаментальные и не только. Просто константы.
 */

namespace spindynapy::constants {

constexpr auto VACUUM_MAGNETIC_PERMEABILITY = 1e-6;
constexpr auto NUMBER_PI = 3.14159;

namespace sci {

constexpr auto mu0 = VACUUM_MAGNETIC_PERMEABILITY;
constexpr auto pi = NUMBER_PI;

}; // namespace sci

}; // namespace spindynapy::constants

namespace spindynapy::doc {

constexpr char module_constants[] =
    ("Модуль, отвечающий за предоставление фундаментальных констант \n"
     "и любых других, используемых в приложении. \n"
     "Здесь - константы через полные названия. \n"
     "В .sci предоставлены \"научные\" короткие версии наименований констант (mu0, pi, etc...)");

constexpr char module_constants_sci[] =
    ("Модуль, отвечающий за предоставление фундаментальных констант \n"
     "и любых других, используемых в приложении. \n"
     "Здесь предоставлены \"научные\" короткие версии наименований констант (mu0, pi, etc...)");

}; // namespace spindynapy::doc

#endif // ! __CONSTANTS_HPP__
