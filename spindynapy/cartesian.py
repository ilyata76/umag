"""
Шорткат для импорта всех модулей в spindynapy.cartesian для декартовой системы координат.
И для удобства, конечно. Из CORE
"""

from .core.types import Material, Anisotropy, UniaxialAnisotropy  # type: ignore # noqa
from .core.types.cartesian import *  # type: ignore # noqa

from .core.geometries.cartesian import *  # type: ignore # noqa

from .core.interactions.cartesian import *  # type: ignore # noqa

from .core.solvers import SolverStrategy  # type: ignore # noqa
from .core.solvers.cartesian import *  # type: ignore # noqa

from .core.simulation.cartesian import *  # type: ignore # noqa

from .core.registries import *  # type: ignore # noqa

from .core import constants  # type: ignore # noqa
