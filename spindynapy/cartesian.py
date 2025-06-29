"""
Cartesian coordinate system shortcut imports.

This module re-exports all core simulation components implemented
  for the Cartesian coordinate system. It serves as a unified
  entry point for accessing geometry, solvers, interactions,
  registries, simulation drivers, and printer utilities.

All imported components are pulled from the compiled C++ backend via pybind11 bindings.

Usage:
    from spindynapy.cartesian import Geometry, Anisotropy, LLGSolver, ...

Note:
    This module exists purely for ergonomic convenience and static typing.
"""

from .core.types import Material, Anisotropy, UniaxialAnisotropy  # type: ignore # noqa
from .core.types.cartesian import *  # type: ignore # noqa

from .core.geometries.cartesian import *  # type: ignore # noqa

from .core.interactions.cartesian import *  # type: ignore # noqa

from .core.solvers import SolverStrategy  # type: ignore # noqa
from .core.solvers.cartesian import *  # type: ignore # noqa

from .core.simulation.cartesian import *  # type: ignore # noqa

from .core.registries import *  # type: ignore # noqa

from .core.printer.cartesian import *  # type: ignore # noqa

from .core import constants  # type: ignore # noqa
