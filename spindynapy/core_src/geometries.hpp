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

    virtual CoordSystem::Moment &operator[](size_t index) {
        throw std::logic_error("Operator [index] Not implemented");
    };
    virtual CoordSystem::Moment getFromIndexes(std::vector<size_t> indexes) {
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
};

namespace cartesian {

using AbstractGeometry = IGeometry<NamespaceCoordSystem>;

/**
 * Простая геометрия, задаваемая в декартовой системе координат
 */
class Geometry : public AbstractGeometry {
  protected:
    MomentsContainer<NamespaceCoordSystem> _moments;
    NeighborsContainers _neighbors_cache;

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

  public:
    Geometry(const std::vector<std::shared_ptr<Moment>> &moments) : _moments(moments) {};
    Geometry(std::vector<std::shared_ptr<Moment>> &&moments) : _moments(std::move(moments)) {}
    Geometry(const Eigen::MatrixXd &moments, MaterialRegistry &material_registry) {
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
        }
    };

    virtual Moment &operator[](size_t index) override { return *this->_moments[index]; };
    virtual Moment getFromIndexes(std::vector<size_t> indexes) override {
        throw std::logic_error("Operator getFromIndecies Not implemented");
    };

    virtual std::vector<size_t> getNeighbors(size_t index, double cutoff_radius) override {
        if (cutoff_radius <= 0) return {};

        // проверяем кэш
        auto cache_key = std::make_pair(index, cutoff_radius);
        if (this->hasInNeighborsCanche(index, cutoff_radius)) {
            return this->getFromNeighborsCache(index, cutoff_radius);
        }

        // иначе считаем
        std::vector<size_t> neighbors;
        const auto &target_coords = _moments[index]->getCoordinates();

        for (size_t i = 0; i < _moments.size(); ++i) {
            if (i == index) continue;
            const double distance_sq = target_coords.getDistanceFrom(_moments[i]->getCoordinates());
            if (distance_sq <= cutoff_radius) neighbors.push_back(i);
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
        return std::make_unique<Geometry>(cloned_moments);
    }
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

    py::class_<AbstractGeometry, std::shared_ptr<AbstractGeometry>>(cartesian, "AbstractGeometry")
        .def("__str__", &AbstractGeometry::__str__)
        .def("__len__", &AbstractGeometry::size)
        .def("get_neighbors", &AbstractGeometry::getNeighbors, py::arg("index"), py::arg("cutoff_radius"))
        .def(
            "__getitem__", &AbstractGeometry::operator[], py::arg("index"), py::return_value_policy::reference
        )
        .doc() = "TODO";

    py::class_<Geometry, AbstractGeometry, std::shared_ptr<Geometry>>(cartesian, "Geometry")
        .def(py::init<const std::vector<std::shared_ptr<Moment>> &>(), py::arg("moments"))
        .def(
            py::init<const Eigen::MatrixXd &, MaterialRegistry &>(),
            py::arg("moments"),
            py::arg("material_registry")
        )
        .doc() = "Простая геометрия, задаваемая в декартовой системе координат";
}

#endif // ! __GEOMETRIES_BASE_HPP__
