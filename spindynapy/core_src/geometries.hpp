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

#include "types.hpp"

#include <memory>
#include <pybind11/pybind11.h>
#include <string>

namespace py = pybind11;

namespace spindynapy {

/**
 * Базовый интерфейс геометрии
 */
class IGeometry {
  protected:
    IGeometry() = default;

  public:
    virtual ~IGeometry() = default;

    virtual std::string __str__() const { return nullptr; };
    virtual std::string __repr__() const { return nullptr; };
    virtual size_t __len__() const { throw std::logic_error("Method __len__ Not implemented"); };
};

/**
 * Простая геометрия, задаваемая в декартовой системе координат
 */
class CartesianGeometry : public IGeometry {
  protected:
    std::vector<std::shared_ptr<CartesianMoment>> _moments;

  public:
    CartesianGeometry(const std::vector<std::shared_ptr<CartesianMoment>> &moments) : _moments(moments) {};
    CartesianGeometry(std::vector<std::shared_ptr<CartesianMoment>> &&moments)
        : _moments(std::move(moments)) {}
    CartesianGeometry(const Eigen::MatrixXd &moments) {
        if (moments.cols() < 7) {
            throw std::invalid_argument("Expected 7 columns: [x, y, z, sx, sy, sz, material]");
        };
        _moments.reserve(moments.rows());
        for (Eigen::Index i = 0; i < moments.rows(); ++i) {
            _moments.emplace_back(std::make_shared<CartesianMoment>(
                CartesianCoordinates(moments(i, 0), moments(i, 1), moments(i, 2)),
                CartesianDirection(moments(i, 3), moments(i, 4), moments(i, 5)),
                regnum(moments(i, 6))
            ));
        }
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

    py::class_<IGeometry, std::shared_ptr<IGeometry>>(module, "IGeometry")
        .def("__str__", &IGeometry::__str__)
        .def("__repr__", &IGeometry::__repr__)
        .def("__len__", &IGeometry::__len__)
        .doc() = "Базовый интерфейс геометрии";

    py::class_<CartesianGeometry, IGeometry, std::shared_ptr<CartesianGeometry>>(module, "CartesianGeometry")
        .def(
            py::init<const std::vector<std::shared_ptr<CartesianMoment>>&>(),
            py::arg("moments")
        )
        .def(
            py::init<const Eigen::MatrixXd &>(),
            py::arg("moments")
        )
        .doc() = "Простая геометрия, задаваемая в декартовой системе координат";

    // clang-format on
}

#endif // ! __GEOMETRIES_BASE_HPP__
