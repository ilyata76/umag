#ifndef __PRINTER_HPP__
#define __PRINTER_HPP__

/**
 * Принтер выводит в виде текста
 *   состояние на основе получаемых данных.
 * Может использоваться для отладки и визуализации.
 */

#include "constants.hpp"
#include "geometries.hpp"
#include "interactions.hpp"
#include "registries.hpp"
#include "simulation.hpp"
#include "types.hpp"

#include <iomanip>
#include <ios>
#include <numeric>
#include <sstream>
#include <string>

namespace PYTHON_API spindynapy {

/**
 * Базовый интерфейс принтера для получения состояния шага/симуляции/etc.
 *      в виде текста в выбранной системе координат.
 */
template <typename CoordSystem> class PYTHON_API IPrinter {
  protected:
    // конструктор только для наследников
    IPrinter() = default;

  public:
    // деструктор
    virtual ~IPrinter() = default;

    // почти полный дамп состояния на шаге step_data
    PYTHON_API virtual std::string shot(
        SimulationStepData<CoordSystem> &step_data,
        MaterialRegistry &material_registry,
        InteractionRegistry<CoordSystem> &interaction_registry
    ) const = 0;

    // вывод в формате vis_mag .vvis
    PYTHON_API virtual std::string vvis(SimulationStepData<CoordSystem> &step_data) const = 0;
};

namespace PYTHON_API cartesian {

/**
 * Базовый интерфейс принтера для получения состояния шага/симуляции/etc.
 *      в виде текста в выбранной (ДЕКАРТОВОЙ) системе координат.
 */
using AbstractPrinter = PYTHON_API IPrinter<NamespaceCoordSystem>;

/**
 * Принтер для симуляции для получения состояния шага/симуляции/etc.
 *      в виде текста в выбранной (ДЕКАРТОВОЙ) системе координат.
 *
 * Может настраиваться на точность вывода данных.
 */
class PYTHON_API SimulationPrinter : public AbstractPrinter {
  protected:
    int _time_precision = 5;
    int _energy_precision = 5;
    int _magnetization_precision = 3;
    int _coordinates_precision = 5;
    int _direction_precision = 5;
    int _field_precision = 5;
    int _precision = 5;

  public:
    // конструктор по умолчанию
    SimulationPrinter(
        int time_precision = 5,
        int energy_precision = 5,
        int magnetization_precision = 3,
        int coordinates_precision = 5,
        int direction_precision = 5,
        int field_precision = 5,
        int precision = 5
    )
        : _time_precision(time_precision),
          _energy_precision(energy_precision),
          _magnetization_precision(magnetization_precision),
          _coordinates_precision(coordinates_precision),
          _direction_precision(direction_precision),
          _field_precision(field_precision),
          _precision(precision) {};

    // почти полный дамп состояния на шаге step_data
    PYTHON_API virtual std::string shot(
        AbstractSimulationStepData &step_data,
        MaterialRegistry &,
        AbstractInteractionRegistry &interaction_registry
    ) const override {
        std::stringstream ss;
        // --- Общие данные ---
        ss << "STEP: " << step_data.step << "\n"
           << "TIME (s): " << std::scientific << std::setprecision(this->_time_precision) << step_data.time
           << "\n"
           << "ATOMS_COUNT: " << step_data.geometry->size() << "\n"
           << "MACROCELLS_COUNT: " << step_data.geometry->getMacrocells().size() << "\n\n";

        // --- Энергии по системе ---
        ss << "Energies by whole system (J):\n";
        // общая
        ss << "ENERGY[-1|FULL]: " << std::showpos << std::scientific
           << std::setprecision(this->_energy_precision)
           << std::accumulate(step_data.energies.begin(), step_data.energies.end(), 0.0) << "\n";
        // и по каждому потенциалу
        for (auto &[reg, energies] : step_data.interaction_energies) {
            auto interaction = interaction_registry.getElement(reg);
            double sumE = std::accumulate(energies.begin(), energies.end(), 0.0);
            ss << std::noshowpos << "ENERGY[" << reg << "|" << interaction->getName() << "]: " << std::showpos
               << std::scientific << std::setprecision(this->_energy_precision) << sumE << "\n"
               << std::noshowpos;
        }
        ss << "\n";

        // --- Намагниченность образца ---
        size_t moments_size = step_data.geometry->size();
        Eigen::Vector3d mean_magnetization = Eigen::Vector3d::Zero();
        for (size_t i = 0; i < moments_size; ++i) {
            mean_magnetization += step_data.geometry->operator[](i).getDirection().asVector();
        }
        mean_magnetization /= double(moments_size);
        ss << "Mean magnetization vector (normalized):\n"
           << "M = (" << std::showpos << std::fixed << std::setprecision(this->_magnetization_precision)
           << mean_magnetization.x() << ", " << mean_magnetization.y() << ", " << mean_magnetization.z()
           << ") |M| = " << std::noshowpos << mean_magnetization.norm() << "\n\n";

        // TODO и материалы здесь выводить

        // АТОМЫ И МОМЕНТЫ

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

        // --- МАКРОЯЧЕЙКИ ---

        const auto &macrocells = step_data.geometry->getMacrocells();
        size_t macrocells_size = macrocells.size();

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
        }

        // clang-format on

        return ss.str();
    }

    // вывод в формате vis_mag .vvis
    PYTHON_API virtual std::string vvis(AbstractSimulationStepData &step_data) const override {
        std::stringstream ss;
        size_t moments_size = step_data.geometry->size();
        auto &geometry = *step_data.geometry;

        ss << "# count\n"
           << moments_size << "\n"
           << "#M\tL\tX\tY\tZ\tSX\tSY\tSZ\n";

        // clang-format off
        for (size_t i = 0; i < moments_size; ++i) {
            auto &moment = geometry[i];
            auto coords = moment.getCoordinates().asVector();
            auto dir = moment.getDirection().asVector();
            ss << moment.getMaterial().getNumber() 
               << "\t0\t"
               << std::fixed << std::setprecision(this->_coordinates_precision) 
               << (coords.x() * 1e10) << "\t" << (coords.y() * 1e10) << "\t" << (coords.z() * 1e10) << "\t"
               << std::setprecision(this->_direction_precision)
               << dir.x() << "\t" << dir.y() << "\t" << dir.z();
            ss << "\n";
        }
        // clang-format on
        return ss.str();
    }
};

}; // namespace PYTHON_API cartesian

}; // namespace PYTHON_API spindynapy

#define ABSTRACTPRINTER_TEMPLATE_BINDINGS(cls)                                                               \
    .def(                                                                                                    \
        "shot",                                                                                              \
        &cls::shot,                                                                                          \
        py::arg("step_data"),                                                                                \
        py::arg("material_registry"),                                                                        \
        py::arg("interaction_registry"),                                                                     \
        py::doc("вывести SHOT-дамп о состоянии системы для шага step_data")                                  \
    )                                                                                                        \
        .def(                                                                                                \
            "vvis",                                                                                          \
            &cls::vvis,                                                                                      \
            py::arg("step_data"),                                                                            \
            py::doc("вывести данные о конфигурации системы в формате .vvis")                                 \
        )

// функция для связывания взаимодействий с Python
inline void pyBindPrinter(py::module_ &module) {
    using namespace spindynapy;

    // -------- | PRINTER | --------
    py::module_ printer_module = module.def_submodule("printer");

    printer_module.doc() = "Принтер выводит в виде текста\n"
                           "  состояние на основе получаемых данных.\n"
                           "Может использоваться для отладки и визуализации.";

    // -------- | CARTESIAN PRINTER | --------
    py::module_ cartesian = printer_module.def_submodule("cartesian");

    cartesian.doc() = "Принтер выводит в виде текста\n"
                      "  состояние на основе получаемых данных.\n"
                      "Может использоваться для отладки и визуализации.";

    {
        using cartesian::AbstractPrinter;
        using cartesian::SimulationPrinter;

        py::class_<AbstractPrinter>(cartesian, "AbstractPrinter")
            ABSTRACTPRINTER_TEMPLATE_BINDINGS(AbstractPrinter)
                .doc() = "Базовый интерфейс принтера для получения состояния шага/симуляции/etc.\n"
                         "     в виде текста в выбранной (ДЕКАРТОВОЙ) системе координат.";

        py::class_<SimulationPrinter, AbstractPrinter>(cartesian, "SimulationPrinter")
            .def(
                py::init<double, double, double, double, double, double, double>(),
                py::arg("time_precision") = 5,
                py::arg("energy_precision") = 5,
                py::arg("magnetization_precision") = 3,
                py::arg("coordinates_precision") = 5,
                py::arg("direction_precision") = 5,
                py::arg("field_precision") = 5,
                py::arg("precision") = 5
            )
            .doc() = "Принтер для симуляции для получения состояния шага/симуляции/etc.\n"
                     "     в виде текста в выбранной (ДЕКАРТОВОЙ) системе координат.\n"
                     "\n"
                     "Может настраиваться на точность вывода данных.";
    }
}

#endif