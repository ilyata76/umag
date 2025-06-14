"""
Модуль с материалами для спиновой динамики
"""

from enum import Enum
from typing import NamedTuple
import numpy as np

from .unit import LatticeConstant, nano
from .core import constants  # type: ignore # noqa
from .core.types import Material, UniaxialAnisotropy  # type: ignore # noqa


class MaterialEnum(Enum):
    """Перечисление для материалов (индекс для регистра)"""

    COBALT = 1
    HYPOTHETIC = 100


class CellSize(NamedTuple):
    """Размеры ячейки"""

    unit: float
    atom: float


def _hcp_unit_cell_size(lattice: LatticeConstant) -> float:
    """гексогональная плотная упаковка"""
    return 3 * np.sqrt(3) / 2 * lattice.a**2 * lattice.c


lattice_lib: dict[MaterialEnum, LatticeConstant] = {
    MaterialEnum.COBALT: LatticeConstant(nano(0.2507), nano(0.408)),
}

cell_size_lib: dict[MaterialEnum, CellSize] = {
    MaterialEnum.COBALT: CellSize(
        unit=_hcp_unit_cell_size(lattice_lib[MaterialEnum.COBALT]),
        atom=_hcp_unit_cell_size(lattice_lib[MaterialEnum.COBALT]) / 6,
    ),
}

mat_lib: dict[MaterialEnum, Material] = {
    MaterialEnum.COBALT: Material(
        material_number=MaterialEnum.COBALT.value,
        exchange_constant_J=5e-21,  # Дж
        atomic_magnetic_saturation_magnetization=1.72,  # в му_B
        damping_constant=0.12,
        unit_cell_size=cell_size_lib[MaterialEnum.COBALT].unit,
        atom_cell_size=cell_size_lib[MaterialEnum.COBALT].atom,
        gyromagnetic_ratio=constants.FREE_SPIN_GYROMAGNETIC_RATIO,
        anisotropy=UniaxialAnisotropy(np.array([0.0, 0.0, 1.0]), 6.69e-24),  # Дж/atom
    )
}
"""
Словарь материалов, где ключ - материал из перечисления, значение - объект Material.
Может использоваться в том числе и для регистрации материалов в MaterialRegistry.
"""
