"""
Material definitions for spin dynamics simulations.

This module defines predefined materials with crystallographic, magnetic,
  and simulation-specific properties used in the spin dynamics core.

It includes:
- Lattice constants and unit cell volumes,
- Mappings to the C++ backend `Material` and `Anisotropy` classes,
- Convenient enums for registry keys and physical identifiers.

Usage:
    from spindynapy.materials import mat_lib, MaterialEnum
    mat_lib[MaterialEnum.COBALT].exchange_constant_J = 5.5e-21
    material_registry = MaterialRegistry({MaterialEnum.COBALT.value: mat_lib[MaterialEnum.COBALT]})
"""

from enum import Enum
from typing import NamedTuple
import numpy as np

from .unit import LatticeConstant, nano
from .core import constants  # type: ignore # noqa
from .core.types import Material, UniaxialAnisotropy  # type: ignore # noqa


class MaterialEnum(Enum):
    """Enumerated identifiers for built-in materials (used as registry keys)."""

    COBALT = 1
    HYPOTHETIC = 100


class CellSize(NamedTuple):
    """Effective volume estimates of unit and atomic cells (in m³)."""

    unit: float
    atom: float


def _hcp_unit_cell_size(lattice: LatticeConstant) -> float:
    """
    Compute the crystallographic unit cell volume for HCP structure.

    Formula:
        V_HCP = (3√3 / 2) · a² · c

    :param lattice: Lattice constants: a, c [SI units]

    :returns: Unit cell volume in cubic meters
    """
    return 3 * np.sqrt(3) / 2 * lattice.a**2 * lattice.c


lattice_lib: dict[MaterialEnum, LatticeConstant] = {
    MaterialEnum.COBALT: LatticeConstant(nano(0.2507), nano(0.408)),
}
"""
Predefined lattice constants for built-in materials.

Each entry maps a `MaterialEnum` identifier to its crystallographic constants
  in the hexagonal close-packed (HCP) system.

@note All distances are in meters (SI). These values are used to compute unit cell volumes.
"""

cell_size_lib: dict[MaterialEnum, CellSize] = {
    MaterialEnum.COBALT: CellSize(
        unit=_hcp_unit_cell_size(lattice_lib[MaterialEnum.COBALT]),
        atom=_hcp_unit_cell_size(lattice_lib[MaterialEnum.COBALT]) / 6,
    ),
}
"""
Precomputed physical volumes for each material.

Maps a `MaterialEnum` to its CellSize, which includes:
- full crystallographic unit cell volume,
- atomic cell volume (unit cell divided by atoms per cell).

@note Used to compute volume-dependent interaction scaling (e.g. dipole coupling, demagnetization).
@note All in SI (m³)
"""

mat_lib: dict[MaterialEnum, Material] = {
    MaterialEnum.COBALT: Material(
        material_number=MaterialEnum.COBALT.value,  # [int]
        exchange_constant_J=5e-21,  # [J]
        atomic_magnetic_saturation_magnetization=1.72,  # [μB]
        damping_constant=0.5,  # [0-1 float]
        unit_cell_size=cell_size_lib[MaterialEnum.COBALT].unit,
        atom_cell_size=cell_size_lib[MaterialEnum.COBALT].atom,
        gyromagnetic_ratio=constants.FREE_SPIN_GYROMAGNETIC_RATIO,  # [rad/s/T]
        anisotropy=UniaxialAnisotropy(np.array([0.0, 0.0, 1.0]), 6.69e-24),  # [J/atom]
    )
}
"""
Material database with simulation-ready physical parameters.

Maps a `MaterialEnum` identifier to a fully constructed `Material` object,
  containing all required attributes for computations.

@note Can be registered into a `MaterialRegistry` object:
>>> material_registry = MaterialRegistry({MaterialEnum.COBALT.value: mat_lib[MaterialEnum.COBALT]})
"""

__all__ = ["MaterialEnum", "lattice_lib", "cell_size_lib", "mat_lib"]
