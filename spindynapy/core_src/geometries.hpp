#ifndef __GEOMETRIES_BASE_HPP__
#define __GEOMETRIES_BASE_HPP__

/**
 * Заголовки с базовыми интерфейсами для геометрий образцов.
 * Геометрии задают устройство образца, сохраняя в себе
 * состояние расположения моментов и их свойств.
 *
 * Геометрия также может поставлять ограниченное количество
 * обсчитываемых параметров. Может состоять из разных участков разной точности.
 */

#include "registries.hpp"
#include "types.hpp"

#include <memory>
#include <object.h>
#include <pybind11/detail/common.h>
#include <pybind11/pybind11.h>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace py = pybind11;

using index_radius_key = std::pair<size_t, double>;

template <typename T> void hash_combine(std::size_t &seed, T const &key) {
    std::hash<T> hasher;
    seed ^= hasher(key) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

template <typename T1, typename T2> struct std::hash<std::pair<T1, T2>> {
    std::size_t operator()(std::pair<T1, T2> const &p) const {
        std::size_t seed1(0);
        ::hash_combine(seed1, p.first);
        ::hash_combine(seed1, p.second);

        std::size_t seed2(0);
        ::hash_combine(seed2, p.second);
        ::hash_combine(seed2, p.first);

        return std::min(seed1, seed2);
    }
};

namespace spindynapy {

template <CoordSystemConcept CoordSystem>
using MomentsContainer = std::vector<std::shared_ptr<typename CoordSystem::Moment>>;

using NeighborsContainers = std::unordered_map<index_radius_key, std::vector<size_t>>;

template <CoordSystemConcept CoordSystem> struct MacroCell {
    std::vector<size_t> moment_indices;
    std::shared_ptr<typename CoordSystem::Moment> avg_moment;

    MacroCell() : moment_indices(0), avg_moment(nullptr) {};

    MacroCell(const std::vector<size_t> &indices, std::shared_ptr<typename CoordSystem::Moment> moment)
        : moment_indices(indices), avg_moment(moment) {};
};

/**
 * Базовый интерфейс геометрии
 */
template <CoordSystemConcept CoordSystem> class IGeometry {
  protected:
    IGeometry() = default;
    virtual bool hasInNeighborsCanche(size_t index, double cutoff_radius) const {
        throw std::logic_error("Method hasInNeighborsCanche Not implemented");
    };
    virtual std::vector<size_t> getFromNeighborsCache(size_t index, double cutoff_radius) const {
        throw std::logic_error("Method getFromNeighborsCache Not implemented");
    };
    virtual void updateNeighborsCache(size_t index, double cutoff_radius, std::vector<size_t> neighbors) {
        throw std::logic_error("Method updateNeighborsCache Not implemented");
    };
    virtual void clearNeighborsCache(size_t index, double cutoff_radius) {
        throw std::logic_error("Method clearNeighborsCache Not implemented");
    };
    virtual void clearAllNeighborsCache() {
        throw std::logic_error("Method clearAllNeighborsCache Not implemented");
    };

  public:
    virtual ~IGeometry() = default;

    virtual std::string __str__() const { return nullptr; };

    virtual void prepareData() { throw std::logic_error("Method hasInNeighborsCanche Not implemented"); };

    virtual const std::vector<MacroCell<CoordSystem>> &getMacrocells() {
        throw std::logic_error("Method getMacrocells Not implemented");
    };
    virtual size_t getMacrocellIndexBySpin(size_t spin_index) {
        throw std::logic_error("Method getMacrocells Not implemented");
    };
    virtual void updateMacrocells() { throw std::logic_error("Method updateMacrocells Not implemented"); };

    virtual void createMacrocellsIfNotCreated(bool recreate = false) {
        throw std::logic_error("Method createMacrocells Not implemented");
    }
    virtual void createMacrocells() { throw std::logic_error("Method createMacrocells Not implemented"); };
    virtual void clearMacrocells() { throw std::logic_error("Method clearMacrocells Not implemented"); };
    virtual MomentsContainer<CoordSystem> getMomentsFromMacrocells(size_t index, double cutoff_radius) {
        throw std::logic_error("Method getMomentsFromMacrocells Not implemented");
    };

    virtual bool hasInNeighborsMacrocellsCache(size_t index, double cutoff_radius) const = 0;
    virtual std::vector<size_t> getFromNeighborsMacrocellsCache(size_t index, double cutoff_radius) const = 0;
    virtual void
    updateNeighborsMacrocellsCache(size_t index, double cutoff_radius, std::vector<size_t> neighbors) = 0;
    virtual void clearNeighborsMacrocellsCache(size_t index, double cutoff_radius) = 0;
    virtual void clearAllNeighborsMacrocellsCache() = 0;

    virtual CoordSystem::Moment &operator[](size_t index) {
        throw std::logic_error("Operator [index] Not implemented");
    };
    virtual MomentsContainer<CoordSystem> getFromIndexes(std::vector<size_t> indexes) {
        throw std::logic_error("Operator getFromIndexes Not implemented");
    };
    virtual std::vector<size_t> getNeighbors(size_t index, double cutoff_radius) {
        throw std::logic_error("Method getNeighbors Not implemented");
    };
    virtual size_t size() const { throw std::logic_error("Method size Not implemented"); };

    using iterator = MomentsContainer<CoordSystem>::iterator;
    virtual iterator begin() { throw std::logic_error("Method iterator.begin Not implemented"); }
    virtual iterator end() { throw std::logic_error("Method iterator.end Not implemented"); }

    virtual std::unique_ptr<IGeometry<CoordSystem>> clone() const {
        throw std::logic_error("Method clone Not implemented");
    };

    virtual Eigen::MatrixXd asNumpy() const { throw std::logic_error("Method asNumPY Not implemented"); }
};

namespace cartesian {

using AbstractGeometry = IGeometry<NamespaceCoordSystem>;
using NamespaceMacroCell = MacroCell<NamespaceCoordSystem>;

/**
 * Простая геометрия, задаваемая в декартовой системе координат
 */
class Geometry : public AbstractGeometry {
  public:
    MomentsContainer<NamespaceCoordSystem> _moments;
    NeighborsContainers _neighbors_cache;
    NeighborsContainers _neighbors_macrocells_cache;

    std::vector<size_t> _spin2cell;              // Карта спин -> ячейка
    std::vector<NamespaceMacroCell> _macrocells; // Макроячейки
    double _macrocell_size = 1e-9;               // Размер макроячейки

    virtual void createMacrocells() override {
        if (_moments.empty())
            return;
        size_t num_spins = _moments.size();

        // Определение границ геометрии ("что" будем разбивать на ячейки)
        auto min_coords = _moments[0]->getCoordinates().asVector();
        auto max_coords = min_coords;

        for (const auto &moment : _moments) {
            Eigen::Vector3d coords = moment->getCoordinates().asVector();
            min_coords = min_coords.cwiseMin(coords); // покомпонентно минимум
            max_coords = max_coords.cwiseMax(coords); // покомпонентно максимум
        }

        // отступ
        min_coords -= Eigen::Vector3d::Constant(this->_macrocell_size / 2);
        max_coords += Eigen::Vector3d::Constant(this->_macrocell_size / 2);

        // размер системы в декартовых измерениях
        Eigen::Vector3d system_size = max_coords - min_coords;

        // Определение размеров сетки (по количеству макроячеек)
        size_t nx = std::max(1UL, static_cast<size_t>(std::ceil(system_size.x() / this->_macrocell_size)));
        size_t ny = std::max(1UL, static_cast<size_t>(std::ceil(system_size.y() / this->_macrocell_size)));
        size_t nz = std::max(1UL, static_cast<size_t>(std::ceil(system_size.z() / this->_macrocell_size)));
        size_t num_cells = nx * ny * nz;

        size_t ix = 0, iy = 0, iz = 0;
        size_t cell_idx = 0;

        // this->_macrocells.assign(num_cells, NamespaceMacroCell());
        this->_spin2cell.assign(num_spins, -1);

        std::map<size_t, std::vector<size_t>>
            spins_by_cell_indicies; // Временная карта для группировки индексов

        // Распределение спинов по ячейкам (лийнерезация)
        for (size_t spin_idx = 0; spin_idx < num_spins; ++spin_idx) {
            auto &pos = _moments[spin_idx]->getCoordinates().asVector();
            ix = std::min(
                nx - 1, static_cast<size_t>(std::floor((pos.x() - min_coords.x()) / this->_macrocell_size))
            );
            iy = std::min(
                ny - 1, static_cast<size_t>(std::floor((pos.y() - min_coords.y()) / this->_macrocell_size))
            );
            iz = std::min(
                nz - 1, static_cast<size_t>(std::floor((pos.z() - min_coords.z()) / this->_macrocell_size))
            );
            cell_idx = ix + iy * nx + iz * nx * ny;
            this->_spin2cell[spin_idx] = cell_idx;
            spins_by_cell_indicies[cell_idx].push_back(spin_idx);
        }

        this->_macrocells = {};               // очищаем макроячейки
        this->_macrocells.reserve(num_cells); // резервируем место под макроячейки

        // для каждой макроячейки k
        for (size_t k = 0; k < num_cells; ++k) {
            // проверяем, что есть хотя бы один спин
            if (spins_by_cell_indicies.contains(k) && !spins_by_cell_indicies[k].empty()) {
                // и записываем (передаём полностью) предсозданную ячейку в финальный список
                MacroCell<NamespaceCoordSystem> new_cell;
                new_cell.moment_indices = spins_by_cell_indicies[k];
                this->_macrocells.push_back(std::move(new_cell));

                // в конечном итоге, записываем за каждым из спинов из карты его ячейку
                auto cell_index = this->_macrocells.size() - 1;
                for (size_t spin_idx : spins_by_cell_indicies[k]) {
                    this->_spin2cell[spin_idx] = cell_index;
                }
            }
        };
    };

    //

    virtual bool hasInNeighborsCanche(size_t index, double cutoff_radius) const override {
        return _neighbors_cache.contains(std::make_pair(index, cutoff_radius));
    };
    virtual std::vector<size_t> getFromNeighborsCache(size_t index, double cutoff_radius) const override {
        return _neighbors_cache.at(std::make_pair(index, cutoff_radius));
    };
    virtual void
    updateNeighborsCache(size_t index, double cutoff_radius, std::vector<size_t> neighbors) override {
        _neighbors_cache[std::make_pair(index, cutoff_radius)] = neighbors;
    };
    virtual void clearNeighborsCache(size_t index, double cutoff_radius) override {
        _neighbors_cache.erase(std::make_pair(index, cutoff_radius));
    };
    virtual void clearAllNeighborsCache() override { _neighbors_cache.clear(); };

    //

    virtual bool hasInNeighborsMacrocellsCache(size_t index, double cutoff_radius) const override {
        return _neighbors_macrocells_cache.contains(std::make_pair(index, cutoff_radius));
    };
    virtual std::vector<size_t>
    getFromNeighborsMacrocellsCache(size_t index, double cutoff_radius) const override {
        return _neighbors_macrocells_cache.at(std::make_pair(index, cutoff_radius));
    };
    virtual void updateNeighborsMacrocellsCache(
        size_t index, double cutoff_radius, std::vector<size_t> neighbors
    ) override {
        _neighbors_macrocells_cache[std::make_pair(index, cutoff_radius)] = neighbors;
    };
    virtual void clearNeighborsMacrocellsCache(size_t index, double cutoff_radius) override {
        _neighbors_macrocells_cache.erase(std::make_pair(index, cutoff_radius));
    };
    virtual void clearAllNeighborsMacrocellsCache() override { _neighbors_macrocells_cache.clear(); };

  public:
    Geometry(const std::vector<std::shared_ptr<Moment>> &moments, double macrocell_size = 1e-9)
        : _moments(moments), _macrocell_size(macrocell_size) {};
    Geometry(std::vector<std::shared_ptr<Moment>> &&moments, double macrocell_size = 1e-9)
        : _moments(std::move(moments)), _macrocell_size(macrocell_size) {};
    Geometry(
        const Eigen::MatrixXd &moments, MaterialRegistry &material_registry, double macrocell_size = 1e-9
    )
        : _macrocell_size(macrocell_size) {
        //
        if (moments.cols() < 7) {
            throw std::invalid_argument("Expected 7 columns: [x, y, z, sx, sy, sz, material]");
        };
        _moments.reserve(moments.rows());
        for (Eigen::Index i = 0; i < moments.rows(); ++i) {
            _moments.emplace_back(std::make_shared<Moment>(
                Coordinates(moments(i, 0), moments(i, 1), moments(i, 2)),
                Direction(moments(i, 3), moments(i, 4), moments(i, 5)),
                material_registry.getElementShared(static_cast<regnum>(moments(i, 6)))
            ));
        };
    };

    virtual MomentsContainer<NamespaceCoordSystem> getFromIndexes(std::vector<size_t> indexes) override {
        MomentsContainer<NamespaceCoordSystem> result(indexes.size(), nullptr);
        for (size_t i = 0; i < indexes.size(); ++i) {
            result[i] = this->_moments[indexes[i]];
        }
        return result;
    };

    virtual void createMacrocellsIfNotCreated(bool recreate = false) override {
        if (this->_macrocells.empty() || recreate) {
            this->createMacrocells(); // создать макроячейки, поделив спины
        }
    };

    virtual void prepareData() override {
        this->createMacrocells(); // создать макроячейки, поделив спины
        this->updateMacrocells(); //
    };

    virtual const std::vector<NamespaceMacroCell> &getMacrocells() override {
        if (this->_macrocells.empty()) {
            throw std::runtime_error("Macrocells are not created. Call prepareData() first.");
        }
        return this->_macrocells;
    };

    virtual size_t getMacrocellIndexBySpin(size_t spin_index) override {
        if (spin_index >= this->_spin2cell.size()) {
            throw std::out_of_range("Spin index out of range (для него нет макроячейки...)");
        }
        return this->_spin2cell[spin_index];
    };

    virtual MomentsContainer<CartesianCoordSystem>
    getMomentsFromMacrocells(size_t index, double cutoff_radius) override {
        // проверяем кэш
        if (this->hasInNeighborsMacrocellsCache(index, cutoff_radius)) {
            auto macro_indexes = this->getFromNeighborsMacrocellsCache(index, cutoff_radius);
            MomentsContainer<CartesianCoordSystem> moments_from_macrocells(macro_indexes.size(), nullptr);
            for (size_t macro_index = 0; macro_index < macro_indexes.size(); ++macro_index) {
                moments_from_macrocells[macro_index] =
                    this->_macrocells[macro_indexes[macro_index]].avg_moment;
            }

            return moments_from_macrocells;
        }

        // иначе считаем
        std::vector<size_t> neighbors;
        const auto &target_coords = _moments[index]->getCoordinates();

        for (size_t i = 0; i < _macrocells.size(); ++i) {
            const double distance =
                target_coords.getDistanceFrom(_macrocells[i].avg_moment->getCoordinates());
            if (distance <= cutoff_radius) {
                neighbors.push_back(i);
            }
        }

        // сохраняем в кэш
        this->updateNeighborsMacrocellsCache(index, cutoff_radius, neighbors);

        MomentsContainer<CartesianCoordSystem> moments_from_macrocells(neighbors.size(), nullptr);
        for (size_t macro_index = 0; macro_index < neighbors.size(); ++macro_index) {
            moments_from_macrocells[macro_index] = this->_macrocells[neighbors[macro_index]].avg_moment;
        }

        return moments_from_macrocells;
    };

    virtual void updateMacrocells() override {
        for (auto &cell : _macrocells) { // Итерируем по существующим ячейкам
            size_t spin_count = cell.moment_indices.size();
            if (spin_count == 0)
                continue; // Пропускаем пустые ячейки (их быть так-то не должно)

            Eigen::Vector3d avg_coordinates_vector(0, 0, 0);
            Eigen::Vector3d avg_direction_vector(0, 0, 0);
            std::map<std::shared_ptr<Material>, size_t> material_counts;

            // обход по всем

            std::shared_ptr<Moment> moment = nullptr;
            std::shared_ptr<Material> material = nullptr;
            for (size_t spin_idx : cell.moment_indices) {
                moment = this->_moments[spin_idx];
                avg_coordinates_vector += moment->getCoordinates().asVector();
                avg_direction_vector += moment->getDirection().asVector();
                if ((material = moment->getMaterialSharedPtr())) {
                    material_counts[material]++;
                }
            }

            // усреднение
            avg_coordinates_vector /= static_cast<double>(spin_count);
            avg_direction_vector /= static_cast<double>(spin_count);

            // Поиск наиболее часто встречающегося материала
            std::shared_ptr<Material> predominant_material = nullptr;
            size_t max_count = 0;
            for (const auto &[material, count] : material_counts) {
                if (count > max_count) {
                    max_count = count;
                    predominant_material = material;
                }
            }

            cell.avg_moment =
                std::make_shared<Moment>(avg_coordinates_vector, avg_direction_vector, predominant_material);
        }
    };

    virtual void clearMacrocells() override {
        this->_macrocells.clear();
        this->_spin2cell.clear();
    };

    virtual Moment &operator[](size_t index) override { return *this->_moments[index]; };

    virtual std::vector<size_t> getNeighbors(size_t index, double cutoff_radius) override {
        if (cutoff_radius <= 0)
            return {};

        // проверяем кэш
        if (this->hasInNeighborsCanche(index, cutoff_radius)) {
            return this->getFromNeighborsCache(index, cutoff_radius);
        }

        // иначе считаем
        std::vector<size_t> neighbors;
        const auto &target_coords = _moments[index]->getCoordinates();

        for (size_t i = 0; i < _moments.size(); ++i) {
            if (i == index) {
                continue;
            }
            const double distance_sq = target_coords.getDistanceFrom(_moments[i]->getCoordinates());
            if (distance_sq <= cutoff_radius) {
                neighbors.push_back(i);
            }
        }

        // сохраняем в кэш
        this->updateNeighborsCache(index, cutoff_radius, neighbors);

        return neighbors;
    };

    virtual std::string __str__() const override {
        std::stringstream ss;
        for (const auto &elem : _moments) {
            ss << "\n" << elem->__str__();
        }
        return ss.str();
    };
    virtual size_t size() const override { return _moments.size(); }

    using iterator = std::vector<std::shared_ptr<Moment>>::iterator;
    virtual iterator begin() override { return _moments.begin(); }
    virtual iterator end() override { return _moments.end(); }

    virtual std::unique_ptr<IGeometry<NamespaceCoordSystem>> clone() const override {
        std::vector<std::shared_ptr<Moment>> cloned_moments;
        cloned_moments.reserve(_moments.size());
        for (const auto &moment : _moments) {
            cloned_moments.push_back(moment->clone());
        }
        return std::make_unique<Geometry>(cloned_moments, this->_macrocell_size);
    }

    virtual Eigen::MatrixXd asNumpy() const override {
        Eigen::MatrixXd numpy_matrix(_moments.size(), 7);
        for (size_t i = 0; i < _moments.size(); ++i) {
            numpy_matrix(i, 0) = _moments[i]->getCoordinates().asVector().x();
            numpy_matrix(i, 1) = _moments[i]->getCoordinates().asVector().y();
            numpy_matrix(i, 2) = _moments[i]->getCoordinates().asVector().z();
            numpy_matrix(i, 3) = _moments[i]->getDirection().asVector().x();
            numpy_matrix(i, 4) = _moments[i]->getDirection().asVector().y();
            numpy_matrix(i, 5) = _moments[i]->getDirection().asVector().z();
            numpy_matrix(i, 6) = _moments[i]->getMaterial().getNumber();
        }
        return numpy_matrix;
    };
};

}; // namespace cartesian

}; // namespace spindynapy

inline void pyBindGeometries(py::module_ &module) {
    using namespace spindynapy;

    // -------- | GEOMETRIES | --------
    py::module_ geometries_module = module.def_submodule("geometries");

    geometries_module.doc() =
        "Модуль геометрий\n"
        "Геометрии задают устройство образца, сохраняя в себе\n"
        "состояние расположения моментов и их свойств.\n"
        "Геометрия также может поставлять ограниченное количество\n"
        "обсчитываемых параметров. Может состоять из разных участков разной точности.\n";

    // -------- | CARTESIAN GEOMETRIES | --------
    py::module_ cartesian = geometries_module.def_submodule("cartesian");

    using cartesian::AbstractGeometry;
    using cartesian::Geometry;
    using cartesian::Moment;
    using cartesian::NamespaceMacroCell;

    py::class_<NamespaceMacroCell, std::shared_ptr<NamespaceMacroCell>>(cartesian, "MacroCell")
        .def_readonly("moment_indices", &NamespaceMacroCell::moment_indices)
        .def_property_readonly(
            "avg_moment",
            [](const NamespaceMacroCell &self) { return self.avg_moment.get(); },
            py::return_value_policy::reference_internal
        )
        .doc() = "Макроячейка";

    py::class_<AbstractGeometry, std::shared_ptr<AbstractGeometry>>(cartesian, "AbstractGeometry")
        .def("__str__", &AbstractGeometry::__str__)
        .def("__len__", &AbstractGeometry::size)
        .def("get_neighbors", &AbstractGeometry::getNeighbors, py::arg("index"), py::arg("cutoff_radius"))
        .def("as_numpy", &AbstractGeometry::asNumpy)
        .def("update_macrocells", &AbstractGeometry::updateMacrocells)
        .def("get_macrocell_index_by_spin", &AbstractGeometry::getMacrocellIndexBySpin, py::arg("spin_index"))
        .def("get_macrocells", &AbstractGeometry::getMacrocells)
        .def(
            "create_macrocells_if_not_created",
            &AbstractGeometry::createMacrocellsIfNotCreated,
            py::arg("recreate") = false
        )
        .def("clear_macrocells", &AbstractGeometry::clearMacrocells)
        .def(
            "__getitem__", &AbstractGeometry::operator[], py::arg("index"), py::return_value_policy::reference
        )
        .doc() = "TODO";

    py::class_<Geometry, AbstractGeometry, std::shared_ptr<Geometry>>(cartesian, "Geometry")
        .def(
            py::init<const std::vector<std::shared_ptr<Moment>> &, double>(),
            py::arg("moments"),
            py::arg("macrocell_size") = 1e-9
        )
        .def(
            py::init<const Eigen::MatrixXd &, MaterialRegistry &, double>(),
            py::arg("moments"),
            py::arg("material_registry"),
            py::arg("macrocell_size") = 1e-9
        )
        .doc() = "Простая геометрия, задаваемая в декартовой системе координат";
}

#endif // ! __GEOMETRIES_BASE_HPP__
