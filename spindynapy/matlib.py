"""
Модуль с материалами для спиновой динамики
"""

from enum import Enum
import numpy as np

from .core import constants  # type: ignore # noqa
from .core.types import Material, UniaxialAnisotropy  # type: ignore # noqa


class MaterialEnum(Enum):
    """Перечисление для материалов (индекс для регистра)"""

    COBALT = 1


mat_lib: dict[str, Material] = {
    "Co": Material(
        material_number=MaterialEnum.COBALT.value,
        exchange_constant_J=6.064e-21,  # Дж
        atomic_magnetic_saturation_magnetization=1.72,  # в му_B
        damping_constant=0.2,
        gyromagnetic_ratio=constants.FREE_SPIN_GYROMAGNETIC_RATIO,
        anisotropy=UniaxialAnisotropy(np.array([0.0, 0.0, 1.0]), 6.69e-24),  # Дж
    )
}
"""
Словарь материалов, где ключ - название материала, значение - объект Material.
Может использоваться в том числе и для регистрации материалов в MaterialRegistry.
"""
