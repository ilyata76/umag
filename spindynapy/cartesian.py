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


def save_data(data: str, path: str) -> None:
    """TODO"""
    with open(path, "w") as f:
        f.write(data)


def get_shot_data(
    step_data: SimulationStepData,  # noqa
    material_registry: MaterialRegistry | None = None,  # noqa
    interaction_registry: InteractionRegistry | None = None,  # noqa
) -> str:
    """TODO

    Args:
        step_data (SimulationStepData): _description_

    Returns:
        str: _description_
    """
    n_atoms = len(step_data.geometry)
    id_width = max(len(str(n_atoms - 1)), len("id"))

    # === GENERALIZED DATA ===

    result: str = f"STEP: {step_data.step}\n" f"TIME (s): {step_data.time:.5e}\n" f"ATOM_COUNT: {n_atoms}\n"

    result += "\nEnergies by whole system (J):\n"
    for regnum, energy_vec in step_data.interaction_energies.items():
        result += (
            f"ENERGY[{regnum}|{interaction_registry.get_element(regnum).get_name() if interaction_registry else ''}]: "
            f"{sum(energy_vec):+.6e}\n"
        )

    # === MAGNETIZATION VECTOR ===

    total_mx, total_my, total_mz = 0.0, 0.0, 0.0
    for i in range(n_atoms):
        s = step_data.geometry[i].get_direction()
        total_mx += s.x
        total_my += s.y
        total_mz += s.z

    total_mx /= n_atoms
    total_my /= n_atoms
    total_mz /= n_atoms

    mag_norm = (total_mx**2 + total_my**2 + total_mz**2) ** 0.5

    result += "\nMean magnetization vector (normalized):\n"
    result += f"M = ({total_mx:+.5f}, {total_my:+.5f}, {total_mz:+.5f}) |M| = {mag_norm:.5f}\n"

    # === ATOM CONCRETE DATA ===

    result += "\nConcrete atoms:\n"

    # шапка для строк
    result += (
        f" {'id':>{id_width}}"
        f" |{'mat':>4}"
        f" |{'x[A]':>7}{'y[A]':>7}{'z[A]':>7}"
        f" |{'sx':>10}{'sy':>10}{'sz':>10}"
        f" |{'|H|':>15}{'Hx':>15}{'Hy':>15}{'Hz':>15}"
    )

    for regnum in step_data.interaction_fields:
        result += (
            f" |{f'|H[{regnum}]|':>15}" f"{f'H[{regnum}]_x':>15}{f'H[{regnum}]_y':>15}{f'H[{regnum}]_z':>15}"
        )

    for regnum in step_data.interaction_energies:
        result += (
            f" |{f'|E[{regnum}]|':>15}"
        )

    result += "\n"

    for i in range(n_atoms):
        moment = step_data.geometry[i]
        mat_number = moment.get_material().get_number()
        c = moment.get_coordinates()
        s = moment.get_direction()
        H = step_data.total_fields[i]
        H_norm = sum(H[j] ** 2 for j in range(3)) ** 0.5

        line = (
            f" {i:>{id_width}d}"
            f" |{mat_number:>4d}"
            f" |{c.x * 1e10:7.4f}{c.y * 1e10:7.4f}{c.z * 1e10:7.4f}"
            f" |{s.x:+10.5f}{s.y:+10.5f}{s.z:+10.5f}"
            f" |{H_norm:+15.5e}{H[0]:+15.5e}{H[1]:+15.5e}{H[2]:+15.5e}"
        )

        for regnum, field_vec in step_data.interaction_fields.items():
            F = field_vec[i]
            F_norm = sum(F[j] ** 2 for j in range(3)) ** 0.5
            line += f" |{F_norm:+15.5e}{F[0]:+15.5e}{F[1]:+15.5e}{F[2]:+15.5e}"

        for regnum, energy in step_data.interaction_energies.items():
            line += f" |{energy[i]:+15.5e}"

        result += line + "\n"

    return result + "\n"


def get_vvis_data(
    step_data: SimulationStepData,  # noqa
    material_registry: MaterialRegistry | None = None,  # noqa
    interaction_registry: InteractionRegistry | None = None,  # noqa
) -> str:
    """
    Генерирует строку в формате .vvis из SimulationStepData.

    Формат:
    # count
    N
    #M	L	X	Y	Z	SX	SY	SZ
    ...
    """
    n_atoms = len(step_data.geometry)
    result = "# count\n"
    result += f"{n_atoms}\n"
    result += "#M\tL\tX\tY\tZ\tSX\tSY\tSZ\n"

    for i in range(n_atoms):
        moment = step_data.geometry[i]
        m = moment.get_material().get_number()
        c = moment.get_coordinates()
        s = moment.get_direction()

        line = (
            f"{m}\t0\t"
            f"{c.x * 1e10:.3f}\t{c.y * 1e10:.3f}\t{c.z * 1e10:.3f}\t"
            f"{s.x:.3f}\t{s.y:.3f}\t{s.z:.3f}"
        )
        result += line + "\n"

    return result
