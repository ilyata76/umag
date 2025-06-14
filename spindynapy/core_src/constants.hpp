#ifndef __SPINDYNAPY_CONSTANTS_HPP__
#define __SPINDYNAPY_CONSTANTS_HPP__

/**
 * @file   constants.hpp
 * @brief  Compile‑time definitions of fundamental physical constants.
 *
 * This header provides a centralised, constexpr, type‑safe repository of the
 *   fundamental constants required by the SpinDynaPy project.
 * All values are expressed in SI units.
 *
 * Two nested namespaces are provided:
 *   - `spindynapy::constants` — descriptive identifiers.
 *   - `spindynapy::constants::sci` — concise scientific aliases (mu0, pi, …).
 *
 * The helper function ::pyBindConstants exports both namespaces to Python,
 * creating the module layout:
 *     spindynapy.constants         # long and short scientific names
 *
 * where every constant is mirrored as a `float` (Python double).
 *
 * @copyright 2025 SpinDynaPy
 */

#include <numbers>
#include <pybind11/pybind11.h>

namespace py = pybind11;

#define PYTHON_API [[]]

namespace PYTHON_API spindynapy {

// ===========================================================================
//  Constants and aliases
// ===========================================================================

namespace PYTHON_API constants {

/** π — ratio of circumference to diameter (dimensionless). */
PYTHON_API constexpr double NUMBER_PI = std::numbers::pi_v<double>;

/** μ₀ — magnetic permeability of vacuum [N·A⁻²] (a.k.a. H·m⁻¹). */
PYTHON_API constexpr double VACUUM_MAGNETIC_PERMEABILITY = 4.0 * NUMBER_PI * 1.0e-7;

/** μ_B — Bohr magneton [J·T⁻¹] (A·m²). */
PYTHON_API constexpr double BOHR_MAGNETON = 9.2740100783e-24; // Дж/Тл (или А*м^2)

/** h — Planck constant [J·s]. */
PYTHON_API constexpr double PLANCK_CONSTANT = 6.626'070'15e-34;

/** ħ — reduced Planck constant [J·s]. */
PYTHON_API constexpr double REDUCED_PLANCK_CONSTANT = PLANCK_CONSTANT / 2 / NUMBER_PI;

/** g_e — free‑electron Landé g‑factor (dimensionless). */
PYTHON_API constexpr double FREE_SPIN_LANDE_FACTOR = 2.002'319'304'362'56;

/** γ_e — free‑electron gyromagnetic ratio [rad·s⁻¹·T⁻¹]. */
PYTHON_API constexpr double FREE_SPIN_GYROMAGNETIC_RATIO =
    (FREE_SPIN_LANDE_FACTOR * BOHR_MAGNETON) / REDUCED_PLANCK_CONSTANT;

/** k_B — Boltzmann constant [J·K⁻¹]. */
PYTHON_API constexpr double BOLTZMANN_CONSTANT = 1.380'649e-23;

namespace PYTHON_API sci {
PYTHON_API constexpr double mu0 = VACUUM_MAGNETIC_PERMEABILITY;
PYTHON_API constexpr double pi = NUMBER_PI;
PYTHON_API constexpr double mu_B = BOHR_MAGNETON;
PYTHON_API constexpr double h = PLANCK_CONSTANT;
PYTHON_API constexpr double h_rad = REDUCED_PLANCK_CONSTANT;
PYTHON_API constexpr double fe_lande = FREE_SPIN_LANDE_FACTOR;
PYTHON_API constexpr double fe_gamma = FREE_SPIN_GYROMAGNETIC_RATIO;
PYTHON_API constexpr double k_b = BOLTZMANN_CONSTANT;
}; // namespace PYTHON_API sci

}; // namespace PYTHON_API constants
}; // namespace PYTHON_API spindynapy

// ===========================================================================
//  Python bindings
// ===========================================================================

/**
 * @brief Bind all constants into a Python sub‑module.
 *
 * @param module Parent PyBind11 module (usually the core extension module).
 * @returns void – extends the parent module.
 */
inline void pyBindConstants(py::module_ &module) {
    using namespace spindynapy::constants;

    // ---------------------------------------------------------------------
    //  Create submodule
    // ---------------------------------------------------------------------
    py::module_ constants_module = module.def_submodule("constants");

    constants_module.doc() =
        "Модуль, отвечающий за предоставление фундаментальных констант \n"
        "и любых других, используемых в приложении. \n"
        "Здесь - константы через полные названия. \n"
        "В .sci предоставлены \"научные\" короткие версии наименований констант (mu0, pi, etc...)";

    // ---------------------------------------------------------------------
    //  Long descriptive names
    // ---------------------------------------------------------------------
    constants_module.attr("VACUUM_MAGNETIC_PERMEABILITY") = VACUUM_MAGNETIC_PERMEABILITY;
    constants_module.attr("NUMBER_PI") = NUMBER_PI;
    constants_module.attr("BOHR_MAGNETON") = BOHR_MAGNETON;
    constants_module.attr("PLANCK_CONSTANT") = PLANCK_CONSTANT;
    constants_module.attr("REDUCED_PLANCK_CONSTANT") = REDUCED_PLANCK_CONSTANT;
    constants_module.attr("FREE_SPIN_LANDE_FACTOR") = FREE_SPIN_LANDE_FACTOR;
    constants_module.attr("FREE_SPIN_GYROMAGNETIC_RATIO") = FREE_SPIN_GYROMAGNETIC_RATIO;
    constants_module.attr("BOLTZMANN_CONSTANT") = BOLTZMANN_CONSTANT;

    // ---------------------------------------------------------------------
    //  Scientific aliases
    // ---------------------------------------------------------------------
    constants_module.attr("mu0") = sci::mu0;
    constants_module.attr("pi") = sci::pi;
    constants_module.attr("mu_B") = sci::mu_B;
    constants_module.attr("h") = sci::h;
    constants_module.attr("h_rad") = sci::h_rad;
    constants_module.attr("fe_lande") = sci::fe_lande;
    constants_module.attr("fe_gamma") = sci::fe_gamma;
    constants_module.attr("k_b") = sci::k_b;
}

#endif // ! __SPINDYNAPY_CONSTANTS_HPP__
