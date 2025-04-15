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
#include <pybind11/detail/common.h>
#include <pybind11/pybind11.h>
#include <string>
#include <utility>
#include <vector>

namespace py = pybind11;

namespace spindynapy {

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
    virtual std::string __repr__() const { return nullptr; };

    virtual CoordSystem::Moment &operator[](size_t index) {
        throw std::logic_error("Operator [index] Not implemented");
    };
    virtual std::vector<size_t> getNeighbors(size_t index, double cutoff_radius) {
        throw std::logic_error("Method getNeighbors Not implemented");
    };
    virtual size_t __len__() const { throw std::logic_error("Method __len__ Not implemented"); };
};

/**
 * Простая геометрия, задаваемая в декартовой системе координат
 */
class CartesianGeometry : public IGeometry<CartesianCoordSystem> {
  protected:
    std::vector<std::shared_ptr<CartesianMoment>> _moments;
    std::map<std::pair<size_t, double>, std::vector<size_t>> _neighbors_cache;

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
    CartesianGeometry(const std::vector<std::shared_ptr<CartesianMoment>> &moments) : _moments(moments) {};
    CartesianGeometry(std::vector<std::shared_ptr<CartesianMoment>> &&moments)
        : _moments(std::move(moments)) {}
    CartesianGeometry(const Eigen::MatrixXd &moments, MaterialRegistry &material_registry) {
        if (moments.cols() < 7) {
            throw std::invalid_argument("Expected 7 columns: [x, y, z, sx, sy, sz, material]");
        };
        _moments.reserve(moments.rows());
        for (Eigen::Index i = 0; i < moments.rows(); ++i) {
            _moments.emplace_back(std::make_shared<CartesianMoment>(
                CartesianCoordinates(moments(i, 0), moments(i, 1), moments(i, 2)),
                CartesianDirection(moments(i, 3), moments(i, 4), moments(i, 5)),
                material_registry.getElementShared(static_cast<regnum>(moments(i, 6)))
            ));
        }
    };

    virtual CartesianMoment &operator[](size_t index) override { return *this->_moments[index]; };

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
            const double distance_sq = target_coords.getDistanceFrom(_moments[i]->getCoordinates());
            if (distance_sq <= cutoff_radius) neighbors.push_back(i);
        }

        // сохраняем в кэш
        this->updateNeighborsCache(index, cutoff_radius, neighbors);

        return neighbors;
    };

    virtual std::string __str__() const override {
        std::string result = "CartesianGeometry : size = " + std::to_string(_moments.size()) + "\n";
        for (const auto &elem : _moments) {
            result += elem->__str__() + "\n";
        }
        return result;
    };
    virtual std::string __repr__() const override {
        std::string result = "CartesianGeometry(moments=[";
        for (const auto &elem : _moments) {
            result += elem->__repr__() + ", ";
        }
        result.pop_back();
        result.pop_back();
        return result + "])";
    };
    virtual size_t __len__() const override { return _moments.size(); }
};

}; // namespace spindynapy

inline void pyBindGeometries(py::module_ &module) {
    using namespace spindynapy;

    // clang-format off

    module.doc() = "Модуль геометрий\n"
                   "Геометрии задают устройство образца, сохраняя в себе\n"
                   "состояние расположения моментов и их свойств.\n"
                   "Геометрия также может поставлять ограниченное количество\n"
                   "обсчитываемых параметров. Может состоять из разных участков разной точности.\n";

    py::class_<CartesianGeometry, std::shared_ptr<CartesianGeometry>>(module, "CartesianGeometry")
        .def("__str__", &CartesianGeometry::__str__)
        .def("__repr__", &CartesianGeometry::__repr__)
        .def("__len__", &CartesianGeometry::__len__)
        .def("get_neighbors", &CartesianGeometry::getNeighbors, py::arg("index"), py::arg("cutoff_radius"))
        .def("__getitem__", &CartesianGeometry::operator[], py::arg("index"), py::return_value_policy::reference)
        .def(
            py::init<const std::vector<std::shared_ptr<CartesianMoment>>&>(),
            py::arg("moments")
        )
        .def(
            py::init<const Eigen::MatrixXd &, MaterialRegistry&>(),
            py::arg("moments"), py::arg("material_registry")
        )
        .doc() = "Простая геометрия, задаваемая в декартовой системе координат";

    // clang-format on
}

#endif // ! __GEOMETRIES_BASE_HPP__
