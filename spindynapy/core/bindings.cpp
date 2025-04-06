/*
 * Точка входа в библиотеку, сборка с PyBind11
 */

#include <pybind11/pybind11.h>

#include "hello.hpp"

PYBIND11_MODULE(spindynapy_core, module) {
    pybind11::class_<HelloClass>(module, "HelloClass")
        .def(pybind11::init<const std::string&>())
        .def("greet", &HelloClass::greet);

    module.def("hello_function", &hello_function, pybind11::arg("name"), "AMOGUS");
};
