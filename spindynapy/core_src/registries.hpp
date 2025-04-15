#ifndef __REGISTRIES_HPP__
#define __REGISTRIES_HPP__

/**
 * Регистры - контейнеры, хранящие в себе экземпляры переиспользуемых классов,
 * которые задают базовые настройки симуляции (например, свойства конкретного материала).
 * У регистров может быть расширенный функционал.
 * Реализуются для паттерна "легковес". Живут всю программу.
 */

#include "interactions.hpp"
#include "types.hpp"

#include <map>
#include <memory>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <string>

namespace py = pybind11;

namespace spindynapy {

template <typename Element> using RegistryContainer = std::map<int, std::shared_ptr<Element>>;

/**
 * Регистр-контейнер
 */
template <typename Element> class Registry {
  protected:
    RegistryContainer<Element> _container;

  public:
    Registry() : _container() {};
    Registry(RegistryContainer<Element> container) : _container(container) {}
    virtual ~Registry() = default;

    virtual std::string __str__() const { throw std::logic_error("Method __Str__ Not implemented"); };
    virtual std::string __repr__() const { throw std::logic_error("Method __Repr__ Not implemented");};

    // virtual void registerElement(regnum, ...) { throw std::logic_error("Method registerMaterial Not implemented"); };
    // virtual void changeElement(regnum, ...) { throw std::logic_error("Method changeMaterial Not implemented"); };
    // virtual Element &getElement(regnum) { throw std::logic_error("Method getMaterial Not implemented"); };
    // virtual void removeElement(regnum) { throw std::logic_error("Method removeMaterial Not implemented"); };

    virtual bool isEmpty() { return _container.empty(); };
};

using MaterialRegistry = Registry<IMaterial>;
using RegionRegistry = Registry<IRegion>;
using InteractionRegistry = Registry<IInteraction>;

}; // namespace spindynapy

inline void pyBindRegistries(py::module_ &module) {
    using namespace spindynapy;

    // clang-format off

    module.doc() = " Регистры - контейнеры, хранящие в себе экземпляры переиспользуемых классов,\n"
                   " которые задают базовые настройки симуляции (например, свойства конкретного материала).\n"
                   " У регистров может быть расширенный функционал.\n"
                   " Реализуются для паттерна \"легковес\". Живут всю программу.";

    py::class_<MaterialRegistry, std::shared_ptr<MaterialRegistry>>(module, "MaterialRegistry")
        .def(py::init<RegistryContainer<IMaterial>>(), py::arg("container"))
        .def("__str__", &MaterialRegistry::__str__)
        .def("__repr__", &MaterialRegistry::__repr__)
        .def("is_empty", &MaterialRegistry::isEmpty, "true, если в регистре нет ни одного элемента")
        .doc() = "Интерфейс для регистра-контейнера для материалов";

    py::class_<RegionRegistry, std::shared_ptr<RegionRegistry>>(module, "RegionRegistry")
        .def(py::init<RegistryContainer<IRegion>>(), py::arg("container"))
        .def("__str__", &RegionRegistry::__str__)
        .def("__repr__", &RegionRegistry::__repr__)
        .def("is_empty", &RegionRegistry::isEmpty, "true, если в регистре нет ни одного элемента")
        .doc() = "Интерфейс для регистра-контейнера для разных областей материалов";

    py::class_<InteractionRegistry, std::shared_ptr<InteractionRegistry>>(module, "InteractionRegistry")
        .def(py::init<RegistryContainer<IInteraction>>(), py::arg("container"))
        .def("__str__", &InteractionRegistry::__str__)
        .def("__repr__", &InteractionRegistry::__repr__)
        .def("is_empty", &InteractionRegistry::isEmpty, "true, если в регистре нет ни одного элемента")
        .doc() = "Интерфейс для регистра-контейнера для разных взаимодействий";

    // clang-format on
}

#endif // ! __REGISTRIES_HPP__
