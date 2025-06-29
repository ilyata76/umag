"""
Units and vector structures for simulation parameters.

This module defines basic physical units (nano, pico, femto)
and helper data structures for describing positions, directions,
and crystallographic constants.
"""


class Unit:
    """
    Scalable unit prefix wrapper.

    Provides a callable object representing a physical unit prefix.
    Multiplying a numeric value by a `Unit` instance scales it accordingly.

    Example:
        nano = Unit(1e-9)
        length = nano(5.0)  # 5 nanometers → 5e-9

    Multiplier for the unit (e.g., 1e-9 for nano).
    """

    def __init__(self, factor):
        self.factor = factor

    def __call__(self, value):
        return value * self.factor


nano = Unit(1e-9)
"""10⁻⁹ multiplier – nano-scale unit."""

pico = Unit(1e-12)
"""10⁻¹² multiplier – pico-scale unit."""

femto = Unit(1e-15)
"""10⁻¹⁵ multiplier – femto-scale unit."""


class XYZ:
    """
    3D Cartesian vector structure.

    Represents a vector or coordinate in three-dimensional space.

    Used for spatial positions, moment directions, or general geometry.

    :param x: X component.
    :param y: Y component.
    :param z: Z component.
    """

    def __init__(self, x: float, y: float, z: float) -> None:
        self.x = x
        self.y = y
        self.z = z


class LatticeConstant:
    """
    Crystallographic lattice constants for hexagonal systems.

    Stores the 'a' and 'c' parameters of a hexagonal unit cell,
      used in volume calculations and geometry generation.

    :param a: In-plane lattice spacing [m].
    :param c: Vertical lattice spacing [m].
    """

    def __init__(self, a: float, c: float) -> None:
        self.a = a
        self.c = c


class ThetaPhi:
    """
    Angular coordinate pair (spherical system).

    Describes a direction in 3D space using spherical angles:
      - `theta`: polar angle from z-axis [0, π],
      - `phi`: azimuthal angle in xy-plane [0, 2π].

    :param theta: Polar angle θ [rad].
    :param phi: Azimuthal angle φ [rad].
    """

    def __init__(self, theta: float, phi: float) -> None:
        self.theta = theta
        self.phi = phi


__all__ = [
    "Unit",
    "nano",
    "pico",
    "femto",
    "XYZ",
    "LatticeConstant",
    "ThetaPhi",
]
