#ifndef __PRINTER_HPP__
#define __PRINTER_HPP__

/**
 * @file   printer.hpp
 * @brief  Textual output utilities for atomistic spin dynamics simulations.
 *
 * This header defines interfaces and implementations for emitting human- and machine-readable
 *   representations of simulation state, compatible with debugging, logging, and external
 *   visualisation tools.
 *
 * The printing framework operates on simulation steps, geometries, and moment data,
 *   and provides formatting hooks to emit rich output in some formats.
 *
 * Exposed entities:
 *   - `IPrinter<CoordSystem>` — abstract printer interface
 *   - `cartesian::SimulationPrinter` — concrete printer for Cartesian simulations
 *
 * Output formats:
 *   - SHOT — detailed, tabular dump of moment and macrocell data (with energy, fields, etc.)
 *   - vvis — line-per-spin format for visualisation tools (position, direction, material, etc.)
 *   - timeSeries  — tabular time series with per-step energy, magnetisation summary, etc.
 *
 * All formats are fully customisable through format strings and precision parameters.
 *
 * Format default definitions:
 *   - `vvis_format_default` — default template for `.vvis` lines
 *   - `time_series_format_default` — default template for `.tsv` lines
 *
 * Python binding macros:
 *   - `ABSTRACTPRINTER_TEMPLATE_BINDINGS(cls)`
 *
 * Python submodule binder:
 *   - `pyBindPrinter()` — exports the `core.printer.cartesian` submodule
 *
 * Design notes:
 *   - Format strings use `{key}` placeholders (compatible with `fmt::format`).
 *   - Position output is in ångströms (Å), all other quantities are in SI units.
 *   - Precision is tunable per-component (e.g., coordinates, field, energy).
 *
 * @warning This printer system is not thread-safe. Concurrent formatting must be externally synchronised.
 *
 * @copyright 2025 SpinDynaPy
 */

#include "constants.hpp"
#include "fmt/base.h"
#include "geometries.hpp"
#include "interactions.hpp"
#include "registries.hpp"
#include "simulation.hpp"
#include "types.hpp"

#include <fmt/args.h>
#include <fmt/format.h>
#include <iomanip>
#include <ios>
#include <sstream>
#include <string>

namespace PYTHON_API spindynapy {

// ============================================================================
//  Formatting constants
// ============================================================================

/**
 * @brief  Default line format for `.vvis` visualization output.
 *
 * This format string controls how each magnetic moment is printed to the `.vvis`
 *   text representation, which is typically used for visualisation or debugging.
 *
 * Each line corresponds to a single moment and includes:
 *   - material id
 *   - position (x, y, z), (Å)
 *   - direction vector (x, y, z), unitless
 *
 * @see SimulationPrinter::vvis
 */
constexpr auto vvis_format_default =
    "{material}\t{coord_x:.3f}\t{coord_y:.3f}\t{coord_z:.3f}\t{dir_x:.3f}\t{dir_y:.3f}\t{dir_z:.3f}";

/**
 * @brief  Default format string for time-series tabular output (`.shot`-like).
 *
 * This format string is used to control the output of a time-series report,
 *   which is typically used for logging or analysis of simulation results.
 *
 * This format string controls the output line per simulation step in a
 *   time-series report. It contains:
 *   - step number
 *   - simulation time (s)
 *   - total energy (J)
 *   - net magnetization norm (unitless)
 *   - net magnetization components (x, y, z) – unitless
 *
 * @see SimulationPrinter::timeSeries
 */
constexpr auto time_series_format_default = "{step}\t{sim_time:.3e}\t{full_energy:+.3e}\t{magnetization:.3f}"
                                            "\t{magnetization_x:+.3f}\t{magnetization_y:"
                                            "+.3f}\t{magnetization_z:+.3f}";

// ============================================================================
//  Printer
// ============================================================================

/**
 * @class  IPrinter
 * @brief  Interface for textual/binary output of simulation state in arbitrary coordinate systems.
 *
 * @tparam CoordSystem  Any type satisfying the coordinate system concept (e.g., Cartesian).
 *
 * This interface defines methods for exporting the state of a spin dynamics simulation
 *   (individual steps, etc., system snapshots) into human\machine-readable text formats.
 *
 * All derived classes are responsible for formatting data such as moment directions,
 *   positions, energies, and fields in a user-configurable style. Typical uses include
 *   debugging, scientific logging, and visualisation preparation.
 *
 * Concrete implementations may customise precision, units, and formatting templates.
 *
 * Functional overview:
 * - `shot(...)`        → Produces a rich full-state dump for a single step (human readable, `.shot`).
 * - `vvis(...)`        → Exports per-moment data for visualisation (e.g. `.vvis` format).
 * - `timeSeries(...)`  → Prints tabular summary per step (e.g. energy and net magnetisation).
 *
 * @note Units are assumed to follow the SI convention unless overridden. Position output
 *       may use angstroms (Å) for human readability, but this is handled explicitly in formatters.
 */
template <typename CoordSystem> class PYTHON_API IPrinter {
  protected:
    /**
     * @brief Protected constructor (non-instantiable base).
     *
     * Only derived printers may be constructed.
     */
    IPrinter() = default;

  public:
    /**
     * @brief Virtual destructor (default).
     */
    virtual ~IPrinter() = default;

    /**
     * @brief Emit a full diagnostic snapshot of a simulation step.
     *
     * This method is intended to produce a complete textual representation of a system state
     *   at a particular step, including energy decomposition, spin configurations,
     *   moment fields, and optionally macrocell statistics, etc.
     *
     * Format is implementation-defined but human-readable (typically multi-line).
     *
     * @param step_data  Simulation step to be printed.
     * @returns Formatted string containing full state snapshot.
     */
    PYTHON_API virtual std::string shot(SimulationStepData<CoordSystem> &step_data) const = 0;

    /**
     * @brief Export spin configuration in `.vvis` (vis_mag-compatible but configurable) format.
     *
     * This output is suitable for downstream visualisation tools and includes:
     *   - position in ångströms (x, y, z)
     *   - spin direction (sx, sy, sz)
     *   - material number
     *
     * Output format can be customized via the format string. Header lines is optional.
     *
     * @param step_data     Simulation step to export.
     * @param format        Format string of line construction.
     * @param print_header  If true, include metadata and column specifier.
     *
     * @returns Text string in specified format (line per moment).
     */
    PYTHON_API virtual std::string vvis(
        SimulationStepData<CoordSystem> &step_data,
        std::string format = vvis_format_default,
        bool print_header = true
    ) const = 0;

    /**
     * @brief Export a time-series trace over simulation steps.
     *
     * The output is a tabular format containing key summary quantities per step,
     *   such as time, total energy, and net magnetisation components.
     *
     * Interaction-specific energy fields may also be embedded if format string includes them.
     *
     * @param simulation     Reference to simulation whose history will be printed.
     * @param format         Line format string.
     * @param print_header   If true, prepend header lines with metadata.
     *
     * @returns Formate specified time-series string (multi-line, table-like).
     */
    PYTHON_API virtual std::string timeSeries(
        Simulation<CoordSystem> &simulation,
        std::string format = time_series_format_default,
        bool print_header = true
    ) const = 0;
};

// ============================================================================
//  artesian output implementation
// ============================================================================

namespace PYTHON_API cartesian {

/**
 * @class  AbstractPrinter
 * @brief  Interface for textual/binary output of simulation state in arbitrary (cartesian) coordinate
 *   systems.
 *
 * @tparam CoordSystem  Any type satisfying the coordinate system concept (e.g., Cartesian).
 *
 * This interface defines methods for exporting the state of a spin dynamics simulation
 *   (individual steps, etc., system snapshots) into human\machine-readable text formats.
 *
 * All derived classes are responsible for formatting data such as moment directions,
 *   positions, energies, and fields in a user-configurable style. Typical uses include
 *   debugging, scientific logging, and visualisation preparation.
 *
 * Concrete implementations may customise precision, units, and formatting templates.
 *
 * Functional overview:
 * - `shot(...)`        → Produces a rich full-state dump for a single step (human readable, `.shot`).
 * - `vvis(...)`        → Exports per-moment data for visualisation (e.g. `.vvis` format).
 * - `timeSeries(...)`  → Prints tabular summary per step (e.g. energy and net magnetisation).
 *
 * @note Units are assumed to follow the SI convention unless overridden. Position output
 *       may use angstroms (Å) for human readability, but this is handled explicitly in formatters.
 */
using AbstractPrinter = PYTHON_API IPrinter<NamespaceCoordSystem>;

/**
 * @class  SimulationPrinter
 * @brief  Concrete printer for simulation output in the Cartesian coordinate system.
 *
 * This printer produces formatted text representations of simulation states in various
 *   human\machine-readable formats. It supports single-step diagnostics (`shot()`), moment-wise dumps
 *   (`vvis()`), and stepwise summaries (`timeSeries()`).
 *
 * Intended uses include:
 *   - debugging simulation dynamics,
 *   - visualisation preprocessing (e.g., .vvis),
 *   - saving intermediate or final results in textual form.
 *
 * The output formats are customisable in both structure and numeric precision. Units are
 *   expressed in SI by default, with select exceptions (positions printed in Å for readability).
 *
 * Internally, this class delegates formatting into modular helpers.
 *
 * @note This class is not thread-safe and assumes external synchronisation
 *       if used in multi-threaded logging.
 */
class PYTHON_API SimulationPrinter : public AbstractPrinter {
  private:
    /** @brief Material registry (reference only). Used to print material details. */
    [[maybe_unused]] MaterialRegistry &_material_registry;
    /** @brief Interaction registry (reference only). Used to print interaction details. */
    [[maybe_unused]] InteractionRegistry<NamespaceCoordSystem> &_interaction_registry;

  protected:
    int _time_precision;          ///< Decimal digits for time values.
    int _energy_precision;        ///< Decimal digits for energies (J).
    int _magnetization_precision; ///< Decimal digits for magnetisation components.
    int _coordinates_precision;   ///< Decimal digits for coordinates (Å).
    int _direction_precision;     ///< Decimal digits for spin directions (unitless).
    int _field_precision;         ///< Decimal digits for effective fields (T).
    int _precision;               ///< Reserved (global override, unused).

    /**
     * @brief Generate a single `.vvis` line for a magnetic moment.
     *
     * @param moment Moment to format.
     * @param format Format string (must include {coord_x}, {dir_x}, etc.).
     *
     * @returns Formatted line.
     */
    std::string vvisMomentLine(Moment &moment, std::string format = vvis_format_default) const {
        auto coords = moment.getCoordinates().asVector() * 1e10; // Å
        auto &dir = moment.getDirection().asVector();
        return fmt::format(
            fmt::runtime(format),
            fmt::arg("material", moment.getMaterial().getNumber()),
            fmt::arg("coord_x", coords.x()),
            fmt::arg("coord_y", coords.y()),
            fmt::arg("coord_z", coords.z()),
            fmt::arg("dir_x", dir.x()),
            fmt::arg("dir_y", dir.y()),
            fmt::arg("dir_z", dir.z())
        );
    }

    /**
     * @brief Format a single time-series line (per step).
     *
     * Embeds global values like energy and magnetisation; may include
     *   interaction-resolved energies depending on the format string.
     *
     * @param step_data Step to export.
     * @param format    Format string (default: time_series_format_default).
     *
     * @returns Formatted line.
     */
    std::string timeSeriesLine(
        const SimulationStepData<NamespaceCoordSystem> &step_data,
        std::string format = time_series_format_default
    ) const {
        const auto mean_magnetization = step_data.getMeanMagnetization();
        fmt::dynamic_format_arg_store<fmt::format_context> store;
        store.push_back(fmt::arg("step", step_data.step));
        store.push_back(fmt::arg("sim_time", step_data.time));
        store.push_back(fmt::arg("full_energy", step_data.getEnergy()));
        store.push_back(fmt::arg("magnetization", mean_magnetization.norm()));
        store.push_back(fmt::arg("magnetization_x", mean_magnetization.x()));
        store.push_back(fmt::arg("magnetization_y", mean_magnetization.y()));
        store.push_back(fmt::arg("magnetization_z", mean_magnetization.z()));
        for (auto &[interaction_number, interaction_energy] : step_data.getEnergyByInteraction()) {
            store.push_back(
                fmt::arg(("energy_" + std::to_string(interaction_number)).c_str(), interaction_energy)
            );
        }
        return fmt::vformat(format, store);
    }

    /**
     * @brief Emit metadata and energy breakdown header for `.shot` format.
     *
     * Includes step number, time, atom/macrocell count, full system energy,
     *   and per-interaction energy contributions.
     *
     * @param step_data Step being printed.
     *
     * @returns SHOT Header string.
     */
    std::string shot_header(AbstractSimulationStepData &step_data) const {
        std::stringstream ss;
        ss << "STEP: " << step_data.step << "\n"
           << "TIME (s): " << std::scientific << std::setprecision(this->_time_precision) << step_data.time
           << "\n"
           << "ATOMS_COUNT: " << step_data.geometry->size() << "\n"
           << "MACROCELLS_COUNT: " << step_data.geometry->getMacrocells().size() << "\n\n";

        // --- Энергии по системе ---
        ss << "Energies by whole system (J):\n";
        // общая
        ss << "ENERGY[-1|FULL]: " << std::showpos << std::scientific
           << std::setprecision(this->_energy_precision) << step_data.getEnergy() << "\n";
        // и по каждому потенциалу
        for (auto &[reg, energy] : step_data.getEnergyByInteraction()) {
            auto interaction = this->_interaction_registry.getElement(reg);
            ss << std::noshowpos << "ENERGY[" << reg << "|" << interaction->getName() << "]: " << std::showpos
               << std::scientific << std::setprecision(this->_energy_precision) << energy << "\n"
               << std::noshowpos;
        }
        ss << "\n";
        Magnetization mean_magnetization = step_data.getMeanMagnetization();
        ss << "Mean magnetization vector (normalized):\n"
           << "M = (" << std::showpos << std::fixed << std::setprecision(this->_magnetization_precision)
           << mean_magnetization.x() << ", " << mean_magnetization.y() << ", " << mean_magnetization.z()
           << ") |M| = " << std::noshowpos << step_data.getMeanMagnetizationNorm() << "\n";
        return ss.str();
    }

    /**
     * @brief Emit per-moment full dump in tabular format.
     *
     * Each moment includes:
     *  - id, material, position, direction
     *  - net field and per-interaction fields
     *  - energy and per-interaction energies
     *
     * Units: position [Å], field [T], energy [J]
     *
     * @param step_data Step to export.
     *
     * @returns Tabular string (multi-line).
     */
    std::string shot_atoms(AbstractSimulationStepData &step_data) const {
        std::stringstream ss;
        auto moments_size = step_data.geometry->size();

        int w_id = std::max(2, static_cast<int>(std::to_string(moments_size - 1).length()));
        int w_mat = 4;
        int w_coord = this->_coordinates_precision + 5;
        int w_spin = this->_direction_precision + 5;
        int w_field = this->_field_precision + 10;
        int w_energy = this->_energy_precision + 10;

        // clang-format off

        ss << std::setw(w_id) << "id" << " |" 
           << std::setw(w_mat) << "mat" << " |"
           << std::setw(w_coord) << "x[A]" << std::setw(w_coord) << "y[A]" << std::setw(w_coord) << "z[A]" << " |"
           << std::setw(w_spin) << "sx" << std::setw(w_spin) << "sy" << std::setw(w_spin) << "sz" << " |"
           << std::setw(w_field) << "|B|[T]" << std::setw(w_field) << "Bx[T]" << std::setw(w_field) << "By[T]" << std::setw(w_field) << "Hz[T]" << " |";

        // поля по каждому interaction
        for (auto & [reg, _] : step_data.interaction_effective_fields) {
            ss << std::setw(w_field) << ("|B[" + std::to_string(reg) + "]|[T]")
               << std::setw(w_field) << ("Bx[" + std::to_string(reg) + "][T]")
               << std::setw(w_field) << ("By[" + std::to_string(reg) + "][T]")
               << std::setw(w_field) << ("Bz[" + std::to_string(reg) + "][T]") << " |";
        }

        // энергии по каждому interaction
        for (auto & [reg, _] : step_data.interaction_energies) {
            ss << std::setw(w_energy) << ("E[" + std::to_string(reg) + "][J]") << " |";
        }

        ss << "\n";

        // сами спины и атомы
        for (size_t i = 0; i < moments_size; ++i) {
            auto &m = (*step_data.geometry)[i];
            auto &c3 = m.getCoordinates();
            auto &s3 = m.getDirection();
            auto &field = step_data.effective_fields[i];

            ss << std::setw(w_id) << i << " |"
               << std::setw(w_mat) << m.getMaterial().getNumber() << " |"
               << std::fixed << std::setprecision(this->_coordinates_precision)
               << std::setw(w_coord) << c3[0] * 1e10 << std::setw(w_coord) << c3[1] * 1e10 << std::setw(w_coord) << c3[2] * 1e10 << " |" 
               << std::setprecision(this->_direction_precision)
               << std::setw(w_spin) <<  s3[0] << std::setw(w_spin) << s3[1] << std::setw(w_spin) << s3[2] << " |"
               << std::scientific << std::setprecision(this->_field_precision)
               << std::setw(w_field) << field.norm()
               << std::setw(w_field) << field[0] << std::setw(w_field) << field[1] << std::setw(w_field) << field[2] << " |";

            // поля по interaction
            for (auto &[reg, field] : step_data.interaction_effective_fields) {
                ss << std::setprecision(this->_field_precision) 
                   << std::setw(w_field) << field[i].norm() 
                   << std::setw(w_field) << field[i].x() << std::setw(w_field) << field[i].y() << std::setw(w_field) << field[i].z() << " |";
            }

            // энергии по interaction
            for (auto &[reg, energy] : step_data.interaction_energies) {
                ss << std::setprecision(this->_energy_precision) << std::setw(w_energy) << energy[i] << " |";
            }

            ss << "\n";
        }

        // clang-format on

        return ss.str();
    }

    /**
     * @brief Emit macrocell summary section.
     *
     * Each macrocell includes:
     *  - index, spin count, representative material
     *  - average position (Å) and direction
     *
     * Only present if geometry supports macrocells.
     *
     * @param step_data Step to export.
     *
     * @returns Formatted macrocell section.
     */
    std::string shot_macrocells(AbstractSimulationStepData &step_data) const {
        const auto &macrocells = step_data.geometry->getMacrocells();
        std::stringstream ss;
        size_t macrocells_size = macrocells.size();
        auto moments_size = step_data.geometry->size();

        int w_id = std::max(2, static_cast<int>(std::to_string(moments_size - 1).length()));
        int w_mat = 4;
        int w_coord = this->_coordinates_precision + 5;
        int w_spin = this->_direction_precision + 5;

        // clang-format off

        if (macrocells_size > 0) {
            // 
            int w_id_mc = std::max(2, static_cast<int>(std::to_string(macrocells_size - 1).length()));

            ss << "\nMACROCELLS:\n"
               << std::setw(w_id) << "id"  << " |"
               << std::setw(w_coord) << "spins [N]" << " |"
               << std::setw(w_mat) << "mat" << " |"
               << std::setw(w_coord) << "x[A]" << std::setw(w_coord) << "y[A]" << std::setw(w_coord) << "z[A]" << " |"
               << std::setw(w_spin) << "sx"   << std::setw(w_spin)  << "sy"   << std::setw(w_spin)  << "sz" << " |";

            ss << "\n";

            // сами макроячейки
            for (size_t j = 0; j < macrocells_size; ++j) {
                auto &cell = macrocells[j];
                auto &m = *cell.avg_moment;
                auto const &coord = m.getCoordinates().asVector();
                auto const  dir = m.getDirection().asVector();

                ss << std::setw(w_id_mc) << j << " |"
                   << std::setw(w_coord) << cell.moment_indices.size() << " |"
                   << std::setw(w_mat) << m.getMaterial().getNumber() << " |"
                   << std::fixed << std::setprecision(this->_coordinates_precision)
                   << std::setw(w_coord) << coord.x() * 1e10 << std::setw(w_coord) << coord.y() * 1e10 << std::setw(w_coord) << coord.z() * 1e10 << " |"
                   << std::setprecision(this->_direction_precision)
                   << std::setw(w_spin) << dir.x() << std::setw(w_spin) << dir.y() << std::setw(w_spin) << dir.z() << " |";

                ss << "\n";
            }

            ss << "\n";
        }

        // clang-format on

        return ss.str();
    }

    /**
     * @brief Emit material list section.
     *
     * Uses the material registry to list registered materials by id and its properties.
     *
     * @returns Material section string.
     */
    std::string shot_materials() const {
        std::stringstream ss;

        // clang-format off

        ss << "Materials:\n";
        for (auto &[reg, material] : this->_material_registry) {
            ss << "Material[" << reg << "]" << "\n"; 
        }
        ss << "\n";

        // clang-format on

        return ss.str();
    }

  public:
    /**
     * @brief Construct a printer for Cartesian simulation output.
     *
     * @param material_registry       Material REGISTRY reference.
     * @param interaction_registry    Interaction REGISTRY reference.
     * @param time_precision          Decimal precision for time values.
     * @param energy_precision        Decimal precision for energies.
     * @param magnetization_precision Decimal precision for magnetisation.
     * @param coordinates_precision   Decimal precision for coordinates.
     * @param direction_precision     Decimal precision for spin vectors.
     * @param field_precision         Decimal precision for effective fields.
     * @param precision               Global precision override (unused).
     */
    SimulationPrinter(
        MaterialRegistry &material_registry,
        InteractionRegistry<NamespaceCoordSystem> &interaction_registry,
        int time_precision = 5,
        int energy_precision = 5,
        int magnetization_precision = 3,
        int coordinates_precision = 5,
        int direction_precision = 5,
        int field_precision = 5,
        int precision = 5
    )
        : _material_registry(material_registry),
          _interaction_registry(interaction_registry),
          _time_precision(time_precision),
          _energy_precision(energy_precision),
          _magnetization_precision(magnetization_precision),
          _coordinates_precision(coordinates_precision),
          _direction_precision(direction_precision),
          _field_precision(field_precision),
          _precision(precision) {};

    /**
     * @brief Produce a full `.shot`-like diagnostic dump for a simulation step.
     *
     * The output is a multiline string containing all relevant simulation state information.
     *
     * Combines:
     *   - header that includes metadata and energy breakdown,
     *   - materials that lists all materials and their properties,
     *   - atoms that provides a detailed moment dump with positions, directions, fields and energies,
     *   - macrocells that summarises macrocell statistics (if applicable).
     *
     * @param step_data Step data to export.
     *
     * @returns Multiline full snapshot string.
     */
    PYTHON_API virtual std::string shot(AbstractSimulationStepData &step_data) const override {
        std::stringstream ss;
        // clang-format off
        ss << this->shot_header(step_data) << "\n"
           << this->shot_materials() << "\n"
           << this->shot_atoms(step_data) << "\n"
           << this->shot_macrocells(step_data) << "\n";
        // clang-format on
        return ss.str();
    }

    /**
     * @brief Emit text vvis-compatible output for the current geometry state using format string.
     *
     * Includes one line per magnetic moment in the geometry.
     *
     * @param step_data     Step to export.
     * @param format        Format string for each moment line (default: vvis_format_default).
     * @param print_header  If true, includes geometry count and header lines.
     *
     * @returns "Raw"-compatible string.
     */
    PYTHON_API virtual std::string vvis(
        AbstractSimulationStepData &step_data,
        std::string format = vvis_format_default,
        bool print_header = true
    ) const override {
        std::stringstream out;
        auto &geometry = *step_data.geometry;
        if (print_header) {
            out << "# count\n" << geometry.size() << "\n# " << format << "\n";
        }
        for (auto &moment : geometry) {
            out << this->vvisMomentLine(*moment, format) << '\n';
        }
        return out.str();
    }

    /**
     * @brief Emit a time-series report over all steps in a simulation.
     *
     * Includes:
     *   - header with interaction mapping
     *   - per-step formatted lines
     *
     * @param simulation     Simulation whose steps are exported.
     * @param format         Format string for each line (see `time_series_format_default`).
     * @param print_header   If true, emits header before data.
     *
     * @returns Multi-line string containing time-series data.
     */
    PYTHON_API virtual std::string timeSeries(
        AbstractSimulation &simulation,
        std::string format = time_series_format_default,
        bool print_header = true
    ) const override {
        std::stringstream ss;
        if (print_header) {
            ss << "# INTERACTIONS: ";
            for (auto &[regnum, interaction] : this->_interaction_registry) {
                ss << interaction->getName() << "=" << regnum << ", ";
            }
            ss << "\n# " << format << "\n";
        }
        for (auto &step_data : simulation.getSteps()) {
            ss << this->timeSeriesLine(step_data, format) << "\n";
        }
        return ss.str();
    }
};

}; // namespace PYTHON_API cartesian

}; // namespace PYTHON_API spindynapy

// ============================================================================
//  Python bindings
// ============================================================================

/**
 * @def   ABSTRACTPRINTER_TEMPLATE_BINDINGS
 * @brief PyBind11 boilerplate for exposing Printer methods.
 */
#define ABSTRACTPRINTER_TEMPLATE_BINDINGS(cls)                                                               \
    .def(                                                                                                    \
        "shot",                                                                                              \
        &cls::shot,                                                                                          \
        py::arg("step_data"),                                                                                \
        py::doc("@brief Produce a full `.shot`-like diagnostic dump for a simulation step.\n"                \
                "\n"                                                                                         \
                "The output is a multiline string containing all relevant simulation state information.\n"   \
                "\n"                                                                                         \
                "Combines:\n"                                                                                \
                "  - header that includes metadata and energy breakdown,\n"                                  \
                "  - materials that lists all materials and their properties,\n"                             \
                "  - atoms that provides a detailed moment dump with positions, directions, fields and "     \
                "energies,\n"                                                                                \
                "  - macrocells that summarises macrocell statistics (if applicable).\n"                     \
                "\n"                                                                                         \
                "@param step_data Step data to export.\n"                                                    \
                "\n"                                                                                         \
                "@returns Multiline full snapshot string.")                                                  \
    )                                                                                                        \
        .def(                                                                                                \
            "vvis",                                                                                          \
            &cls::vvis,                                                                                      \
            py::arg("step_data"),                                                                            \
            py::arg("format") = vvis_format_default,                                                         \
            py::arg("print_header") = true,                                                                  \
            py::doc(                                                                                         \
                "@brief Emit text vvis-compatible output for the current geometry state using format "       \
                "string.\n"                                                                                  \
                "\n"                                                                                         \
                "Includes one line per magnetic moment in the geometry.\n"                                   \
                "\n"                                                                                         \
                "@param step_data     Step to export.\n"                                                     \
                "@param format        Format string for each moment line (default: vvis_format_default).\n"  \
                "@param print_header  If true, includes geometry count and header lines.\n"                  \
                "\n"                                                                                         \
                "@returns \"Raw\"-compatible string."                                                        \
            )                                                                                                \
        )                                                                                                    \
        .def(                                                                                                \
            "time_series",                                                                                   \
            &cls::timeSeries,                                                                                \
            py::arg("simulation"),                                                                           \
            py::arg("format") = time_series_format_default,                                                  \
            py::arg("print_header") = true,                                                                  \
            py::doc(                                                                                         \
                "@brief Emit a time-series report over all steps in a simulation.\n"                         \
                "\n"                                                                                         \
                "Includes:\n"                                                                                \
                "  - header with interaction mapping\n"                                                      \
                "  - per-step formatted lines\n"                                                             \
                "\n"                                                                                         \
                "@param simulation     Simulation whose steps are exported.\n"                               \
                "@param format         Format string for each line (see `time_series_format_default`).\n"    \
                "@param print_header   If true, emits header before data.\n"                                 \
                "\n"                                                                                         \
                "@returns Multi-line string containing time-series data."                                    \
            )                                                                                                \
        )

/**
 * @brief Bind the logging utilities to a Python sub‑module.
 *
 * @param module Parent PyBind11 module (usually the core extension module).
 * @returns void – extends the parent module.
 */
inline void pyBindPrinter(py::module_ &module) {
    using namespace spindynapy;

    // ---------------------------------------------------------------------
    //  Create submodule
    // ---------------------------------------------------------------------

    py::module_ printer_module = module.def_submodule("printer");

    printer_module.doc() =
        "@brief  Textual output utilities for atomistic spin dynamics simulations.\n"
        "\n"
        "This header defines interfaces and implementations for emitting human- and machine-readable\n"
        "  representations of simulation state, compatible with debugging, logging, and external\n"
        "  visualisation tools.\n"
        "\n"
        "The printing framework operates on simulation steps, geometries, and moment data,\n"
        "  and provides formatting hooks to emit rich output in some formats.\n"
        "\n"
        "Exposed entities:\n"
        "  - `IPrinter<CoordSystem>` — abstract printer interface\n"
        "  - `cartesian::SimulationPrinter` — concrete printer for Cartesian simulations\n"
        "\n"
        "Output formats:\n"
        "  - SHOT — detailed, tabular dump of moment and macrocell data (with energy, fields, etc.)\n"
        "  - vvis — line-per-spin format for visualisation tools (position, direction, material, etc.)\n"
        "  - timeSeries  — tabular time series with per-step energy, magnetisation summary, etc.\n"
        "All formats are fully customisable through format strings and precision parameters.\n"
        "\n"
        "Python submodule binder:\n"
        "  - `pyBindPrinter()` — exports the `core.printer.cartesian` submodule\n"
        "\n"
        "Design notes:\n"
        "  - Format strings use `{key}` placeholders (compatible with `fmt::format`).\n"
        "  - Position output is in ångströms (Å), all other quantities are in SI units.\n"
        "  - Precision is tunable per-component (e.g., coordinates, field, energy).\n"
        "\n"
        "@warning This printer system is not thread-safe. Concurrent formatting must be externally "
        "synchronised.";

    // ---------------------------------------------------------------------
    //  Cartesian submodule bindings
    // ---------------------------------------------------------------------

    py::module_ cartesian = printer_module.def_submodule("cartesian");
    cartesian.doc() = printer_module.doc();

    {
        using cartesian::AbstractInteractionRegistry;
        using cartesian::AbstractPrinter;
        using cartesian::SimulationPrinter;

        py::class_<AbstractPrinter>(cartesian, "AbstractPrinter")
            ABSTRACTPRINTER_TEMPLATE_BINDINGS(AbstractPrinter)
                .doc() =
            "@class  AbstractPrinter\n"
            "@brief  Interface for textual/binary output of simulation state in arbitrary (cartesian) "
            "coordinate\n"
            "  systems.\n"
            "\n"
            "@tparam CoordSystem  Any type satisfying the coordinate system concept (e.g., Cartesian).\n"
            "\n"
            "This interface defines methods for exporting the state of a spin dynamics simulation\n"
            "  (individual steps, etc., system snapshots) into human\\machine-readable text formats.\n"
            "\n"
            "All derived classes are responsible for formatting data such as moment directions,\n"
            "  positions, energies, and fields in a user-configurable style. Typical uses include\n"
            "  debugging, scientific logging, and visualisation preparation.\n"
            "\n"
            "Concrete implementations may customise precision, units, and formatting templates.\n"
            "\n"
            "Functional overview:\n"
            "- `shot(...)`        → Produces a rich full-state dump for a single step (human readable, "
            "`.shot`).\n"
            "- `vvis(...)`        → Exports per-moment data for visualisation (e.g. `.vvis` format).\n"
            "- `timeSeries(...)`  → Prints tabular summary per step (e.g. energy and net magnetisation).\n"
            "\n"
            "@note Units are assumed to follow the SI convention unless overridden. Position output\n"
            "      may use angstroms (Å) for human readability, but this is handled explicitly in "
            "formatters.";

        py::class_<SimulationPrinter, AbstractPrinter>(cartesian, "SimulationPrinter")
            .def(
                py::init<
                    MaterialRegistry &,
                    AbstractInteractionRegistry &,
                    int,
                    int,
                    int,
                    int,
                    int,
                    int,
                    int>(),
                py::arg("material_registry"),
                py::arg("interaction_registry"),
                py::arg("time_precision") = 5,
                py::arg("energy_precision") = 5,
                py::arg("magnetization_precision") = 3,
                py::arg("coordinates_precision") = 5,
                py::arg("direction_precision") = 5,
                py::arg("field_precision") = 5,
                py::arg("precision") = 5
            )
            .doc() =
            "@class  SimulationPrinter\n"
            "@brief  Concrete printer for simulation output in the Cartesian coordinate system.\n"
            "\n"
            "This printer produces formatted text representations of simulation states in various\n"
            "  human\\machine-readable formats. It supports single-step diagnostics (`shot()`), moment-wise "
            "dumps\n"
            "  (`vvis()`), and stepwise summaries (`timeSeries()`).\n"
            "\n"
            "Intended uses include:\n"
            "  - debugging simulation dynamics,\n"
            "  - visualisation preprocessing (e.g., .vvis),\n"
            "  - saving intermediate or final results in textual form.\n"
            "\n"
            "The output formats are customisable in both structure and numeric precision. Units are\n"
            "  expressed in SI by default, with select exceptions (positions printed in Å for readability).\n"
            "\n"
            "Internally, this class delegates formatting into modular helpers.\n"
            "\n"
            "@note This class is not thread-safe and assumes external synchronisation\n"
            "      if used in multi-threaded logging.";
    }
}

#endif