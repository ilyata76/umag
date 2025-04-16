#ifndef __REGISTRIES_HPP__
#define __REGISTRIES_HPP__

/**
 * Регистры - контейнеры, хранящие в себе экземпляры переиспользуемых классов,
 * которые задают базовые настройки симуляции (например, свойства конкретного материала).
 * У регистров может быть расширенный функционал.
 * Реализуются для паттерна "легковес". Живут всю программу.
 */

#include "types.hpp"

#include <map>
#include <memory>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <stdexcept>
#include <string>

namespace py = pybind11;

namespace spindynapy {

template <typename Element> using RegistryContainer = std::map<regnum, std::shared_ptr<Element>>;

/**
 * Регистр-контейнер
 */
template <typename Element> class Registry {
  protected:
    RegistryContainer<Element> _container;

  public:
    Registry() : _container() {};
    Registry(RegistryContainer<Element> container) : _container(container) {}
    ~Registry() = default;

    std::string __str__() const { throw std::logic_error("Method __Str__ Not implemented"); };
    std::string __repr__() const { throw std::logic_error("Method __Repr__ Not implemented"); };

    Element &getElement(regnum number) {
        if (this->_container.contains(number))
            return *this->_container.at(number);
        else
            throw std::invalid_argument(std::format(
                "Такой элемент (<{}> {}) не был зарегистрирован в регистре", typeid(this).name(), number
            ));
    };

    std::shared_ptr<Element> getElementShared(regnum number) {
        if (this->_container.contains(number))
            return this->_container.at(number);
        else
            throw std::invalid_argument(std::format(
                "Такой элемент (<{}> {}) не был зарегистрирован в регистре", typeid(this).name(), number
            ));
    }

    bool isEmpty() { return _container.empty(); };

    using iterator = RegistryContainer<Element>::iterator;
    iterator begin() { return iterator(_container.begin()); }
    iterator end() { return iterator(_container.end()); }
    iterator cbegin() const { return iterator(_container.cbegin()); }
    iterator cend() const { return iterator(_container.cend()); }
};

using MaterialRegistry = Registry<Material>;
using RegionRegistry = Registry<Region>;

}; // namespace spindynapy

#define BUILD_REGISTRY_TEMPLATE_METHODS(cls)                                                                 \
    .def("__str__", &cls::__str__)                                                                           \
        .def("__repr__", &cls::__repr__)                                                                     \
        .def("get_element", &cls::getElement, py::arg("number"), py::return_value_policy::reference)         \
        .def("is_empty", &cls::isEmpty, "true, если в регистре нет ни одного элемента")

inline void pyBindRegistries(py::module_ &module) {
    using namespace spindynapy;

    // clang-format off

    module.doc() = " Регистры - контейнеры, хранящие в себе экземпляры переиспользуемых классов,\n"
                   " которые задают базовые настройки симуляции (например, свойства конкретного материала).\n"
                   " У регистров может быть расширенный функционал.\n"
                   " Реализуются для паттерна \"легковес\". Живут всю программу.";

    py::class_<MaterialRegistry, std::shared_ptr<MaterialRegistry>>(module, "MaterialRegistry")
        .def(py::init<>())
        .def(py::init<RegistryContainer<Material>>(), py::arg("container"))
        BUILD_REGISTRY_TEMPLATE_METHODS(MaterialRegistry)
        .doc() = "Интерфейс для регистра-контейнера для материалов";

    py::class_<RegionRegistry, std::shared_ptr<RegionRegistry>>(module, "RegionRegistry")
        .def(py::init<>())
        .def(py::init<RegistryContainer<Region>>(), py::arg("container"))
        BUILD_REGISTRY_TEMPLATE_METHODS(RegionRegistry)
        .doc() = "Интерфейс для регистра-контейнера для разных областей материалов";

    // clang-format on
}

#endif // ! __REGISTRIES_HPP__
