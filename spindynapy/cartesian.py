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

from .core.printer.cartesian import *  # type: ignore # noqa

from .core import constants  # type: ignore # noqa


def save_data(data: str, path: str) -> None:
    """Сохранить данные в файл по пути path"""
    with open(path, "w") as f:
        f.write(data)


def get_shot_data(
    step_data: SimulationStepData,  # noqa
    material_registry: MaterialRegistry | None = None,  # noqa
    interaction_registry: InteractionRegistry | None = None,  # noqa
    printer: AbstractPrinter | None = None,  # noqa
) -> str:
    """Получить данные в формате .shot для анализа (текстовый квази-дамп-файл)"""
    if printer is None:
        printer = SimulationPrinter()  # noqa
    return printer.shot(step_data, material_registry, interaction_registry)


def get_vvis_data(
    step_data: SimulationStepData,  # noqa
    printer: AbstractPrinter | None = None,  # noqa
) -> str:
    """Получить данные в формате .vvis для vis_mag(visualis)"""
    if printer is None:
        printer = SimulationPrinter(coordinates_precision=3, direction_precision=3)  # noqa
    return printer.vvis(step_data)
