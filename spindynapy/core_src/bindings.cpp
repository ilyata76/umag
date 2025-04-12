/*
 * Точка входа в библиотеку, сборка с PyBind11 всех модулей, подмодулей,
 * функций и классов.
 */

#include "constants.hpp"
#include "geometries/base.hpp"
#include "registries/base.hpp"
#include "types/base.hpp"
#include "types/cartesian.hpp"
#include "types/material.hpp"

#include <pybind11/pybind11.h>

namespace py = pybind11;
namespace sd = spindynapy;

PYBIND11_MODULE(core, module) {

    // types
    py::module_ types_module = module.def_submodule("types");
    types_module.doc() = sd::doc::module_types;

    // types/base
    py::module_ types_base_module = types_module.def_submodule("base");
    types_base_module.doc() = sd::doc::module_types_base;

    auto StrPresentationMixin =
        py::class_<sd::StrPresentationMixin>(types_base_module, "StrPresentationMixin")
            .def("__str__", &sd::StrPresentationMixin::__str__, py::doc(sd::doc::StrPresentationMixin___str__))
            .def("__repr__", &sd::StrPresentationMixin::__repr__, py::doc(sd::doc::StrPresentationMixin___repr__));
    StrPresentationMixin.doc() = sd::doc::StrPresentationMixin;

    auto IMaterial = py::class_<sd::IMaterial, sd::StrPresentationMixin>(types_base_module, "IMaterial");
    IMaterial.doc() = sd::doc::IMaterial;

    auto ICoordinates = py::class_<sd::ICoordinates, sd::StrPresentationMixin>(types_base_module, "ICoordinates");
    ICoordinates.doc() = sd::doc::ICoordinates;

    auto IDirection = py::class_<sd::IDirection, sd::StrPresentationMixin>(types_base_module, "IDirection");
    IDirection.doc() = sd::doc::IDirection;

    auto IMoment = py::class_<sd::IMoment, sd::StrPresentationMixin>(types_base_module, "IMoment")
                       .def("getDirection", &sd::IMoment::getDirection, py::return_value_policy::reference,
                            py::doc(sd::doc::IMoment_getDirection))
                       .def("getCoordinates", &sd::IMoment::getCoordinates, py::return_value_policy::reference,
                            py::doc(sd::doc::IMoment_getCoordinates));
    IMoment.doc() = sd::doc::IMoment;

    auto ISpin = py::class_<sd::ISpin, sd::IMoment>(types_base_module, "ISpin");
    ISpin.doc() = sd::doc::ISpin;

    // types/material
    py::module_ types_material_module = types_module.def_submodule("material");
    types_material_module.doc() = sd::doc::module_types_material;

    auto MagneticMaterial = py::class_<sd::MagneticMaterial, sd::IMaterial>(types_material_module, "MagneticMaterial")
                                .def(py::init<double>(), py::arg("_exchange_constant_J"));
    MagneticMaterial.doc() = sd::doc::MagneticMaterial;

    // types/cartesian
    py::module_ types_cartesian_module = types_module.def_submodule("cartesian");
    types_cartesian_module.doc() = sd::doc::module_types_cartesian;

    auto CartesianCoordinates =
        py::class_<sd::CartesianCoordinates, sd::ICoordinates>(types_cartesian_module, "CartesianCoordinates")
            .def(py::init<double, double, double>(), py::arg("_x"), py::arg("_y"), py::arg("_z"))
            .def_readwrite("x", &sd::CartesianCoordinates::x)
            .def_readwrite("y", &sd::CartesianCoordinates::y)
            .def_readwrite("z", &sd::CartesianCoordinates::z);
    CartesianCoordinates.doc() = sd::doc::CartesianCoordinates;

    auto CartesianDirection =
        py::class_<sd::CartesianDirection, sd::IDirection>(types_cartesian_module, "CartesianDirection")
            .def(py::init<double, double, double>(), py::arg("_sx"), py::arg("_sy"), py::arg("_sz"))
            .def_readwrite("sx", &sd::CartesianDirection::sx)
            .def_readwrite("sy", &sd::CartesianDirection::sy)
            .def_readwrite("sz", &sd::CartesianDirection::sz);
    CartesianDirection.doc() = sd::doc::CartesianDirection;

    auto CartesianMoment =
        py::class_<sd::CartesianMoment, sd::IMoment>(types_cartesian_module, "CartesianMoment")
            .def(py::init<sd::CartesianCoordinates &, sd::CartesianDirection &>(), py::arg("_coordinates"),
                 py::arg("_direction"))
            .def("getDirection", &sd::CartesianMoment::getDirection, py::return_value_policy::reference)
            .def("getCoordinates", &sd::CartesianMoment::getCoordinates, py::return_value_policy::reference);
    CartesianMoment.doc() = sd::doc::CartesianMoment;

    auto CartesianSpin =
        py::class_<sd::CartesianSpin, sd::CartesianMoment, sd::ISpin>(types_cartesian_module, "CartesianSpin")
            .def(py::init<sd::CartesianCoordinates &, sd::CartesianDirection &>(), py::arg("_coordinates"),
                 py::arg("_direction"));
    CartesianSpin.doc() = sd::doc::CartesianSpin;

    // geometries
    py::module_ geometries_module = module.def_submodule("geometries");
    geometries_module.doc() = sd::doc::module_geometries;

    // geometries/base
    py::module_ geometries_base_module = geometries_module.def_submodule("base");
    geometries_base_module.doc() = sd::doc::module_geometries_base;

    auto IGeometry = py::class_<sd::IGeometry, sd::StrPresentationMixin>(geometries_base_module, "IGeometry");
    IGeometry.doc() = sd::doc::IGeometry;

    // registries
    py::module_ registries_module = module.def_submodule("registries");
    registries_module.doc() = sd::doc::module_registries;

    // registries/base
    py::module_ registries_base_module = registries_module.def_submodule("base");
    registries_base_module.doc() = sd::doc::module_registries_base;

    auto IRegistry = py::class_<sd::IRegistry, sd::StrPresentationMixin>(registries_base_module, "IRegistry");
    IRegistry.doc() = sd::doc::IRegistry;

    auto IMaterialRegistry =
        py::class_<sd::IMaterialRegistry, sd::IRegistry>(registries_base_module, "IMaterialRegistry");
    IMaterialRegistry.doc() = sd::doc::IMaterialRegistry;

    auto IRegionRegistry = py::class_<sd::IRegionRegistry, sd::IRegistry>(registries_base_module, "IRegionRegistry");
    IRegionRegistry.doc() = sd::doc::IRegionRegistry;

    // constants
    py::module_ constants_module = module.def_submodule("constants");
    constants_module.doc() = sd::doc::module_constants;

    constants_module.doc() = sd::doc::module_constants;
    constants_module.attr("VACUUM_MAGNETIC_PERMEABILITY") = sd::constants::VACUUM_MAGNETIC_PERMEABILITY;
    constants_module.attr("NUMBER_PI") = sd::constants::NUMBER_PI;

    // constants/sci
    py::module_ constants_sci_module = constants_module.def_submodule("sci");
    constants_sci_module.doc() = sd::doc::module_constants_sci;

    constants_sci_module.doc() = sd::doc::module_constants_sci;
    constants_sci_module.attr("mu0") = sd::constants::sci::mu0;
    constants_sci_module.attr("pi") = sd::constants::sci::pi;

}; // ! PYBIND11_MODULE
