"""
Модуль для работы с единицами измерения и векторами.
"""


class Unit:
    """Класс для определения единиц измерения (порядок)."""

    def __init__(self, factor):
        self.factor = factor

    def __call__(self, value):
        return value * self.factor


nano = Unit(1e-9)
"""Нано, 10^-9"""

pico = Unit(1e-12)
"""Пико, 10^-9"""

femto = Unit(1e-15)
"""Фемто, 10^-9"""


class XYZ:
    """Класс для определения векторов/координат/etc. в 3D пространстве."""

    def __init__(self, x: float, y: float, z: float) -> None:
        self.x = x
        self.y = y
        self.z = z


class ThetaPhi:
    """Класс для определения векторов в пространстве через углы."""

    def __init__(self, theta: float, phi: float) -> None:
        self.theta = theta
        self.phi = phi
