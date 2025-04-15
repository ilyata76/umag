#ifndef __INTERACTIONS_HPP__
#define __INTERACTIONS_HPP__

/**
 * Интерфейсы взаимодействий между элементами системы.
 * Таковыми могут быть потенциальные поля, силы, etc.
 */

#include <pybind11/pybind11.h>
#include <string>

namespace py = pybind11;

namespace spindynapy {

/**
 * Базовый интерфейс взаимодействий (сил, полей, etc.)
 */
class IInteraction {
  protected:
    IInteraction() = default;

  public:
    virtual ~IInteraction() {};

    virtual std::string __str__() const { return nullptr; };
    virtual std::string __repr__() const { return nullptr; };
};

/**
 * Обменное взаимодействие
 */
class ExchangeInteraction : public IInteraction {
  public:
    ExchangeInteraction() {};
};

}; // namespace spindynapy

inline void pyBindInteractions(py::module_ &module) {
    using namespace spindynapy;

    // clang-format off

    module.doc() = "Интерфейсы взаимодействий между элементами системы.\n"
                   "Таковыми могут быть потенциальные поля, силы, etc. \n";

    py::class_<IInteraction, std::shared_ptr<IInteraction>>(module, "IInteraction")
        .def("__str__", &IInteraction::__str__)
        .def("__repr__", &IInteraction::__repr__)
        .doc() = "Базовый интерфейс взаимодействий (сил, полей, etc.)";

    py::class_<ExchangeInteraction, IInteraction, std::shared_ptr<ExchangeInteraction>>(module, "ExchangeInteraction")
        .def(py::init<>())
        .doc() = "Обменное взаимодействие";

    // clang-format on
}

#endif // ! __INTERACTIONS_HPP__
