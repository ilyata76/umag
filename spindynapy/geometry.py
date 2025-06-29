"""
Geometry generation and serialization utilities based on NumPy.

This module provides a geometry manager class for:
  - saving and loading NumPy-based geometry arrays,
  - generating structured 3D atomic lattices (SC HCP, etc.),

Each geometry is represented as a NumPy array of shape (N, 7):
    [x, y, z, Sx, Sy, Sz, material_id]
Where:
    - (x, y, z)       — coordinates in SI units [m],
    - (Sx, Sy, Sz)    — normalized spin vector components,
    - material_id     — numeric identifier linked to a Material.

Supported lattice types:
    - SC (simple cubic) lattice with rectangular XY base,
    - HCP (hexagonal close-packed) with ABAB stacking.

Example:
    manager = NumpyGeometryManager()
    geom = manager.generate_hcp_monomaterial_parallelepiped(...)
    manager.save_geometry("geom_hcp", geom)

Usage example:
    numpy_geometry = NumpyGeometryManager.generate_hcp_monomaterial_parallelepiped(
        lattice_constant=lattice_lib[MaterialEnum.COBALT],
        size=XYZ(nano(5.0), nano(5.0), nano(5.0)),
        material_number=mat_lib[MaterialEnum.COBALT].get_number(),
        initial_direction=None,
        base_shift=None,
    )
    material_registry = MaterialRegistry({MaterialEnum.COBALT.value: mat_lib[MaterialEnum.COBALT]})
    geometry = Geometry(moments=numpy_geometry, material_registry=material_registry)

Note:
    All dimensions are in SI units. Lattice spacing is interpreted in meters.
    Initial spin direction can be random or fixed.
    Seeded RNG ensures reproducible geometries.
"""

import numpy as np

from .unit import XYZ, LatticeConstant


class NumpyGeometryManager:
    """
    Utilities for managing geometry creation and storage using NumPy arrays.

    Provides methods for:
      - Saving and loading spin geometries to/from `.npy` files.
      - Procedurally generating atomic lattices (SC, HCP, etc.).
      - Initializing spins and assigning material identifiers.

    All generated geometries follow the format: [x, y, z, Sx, Sy, Sz, material_id].
    """

    @staticmethod
    def save_geometry(path: str, geometry: np.typing.ArrayLike) -> None:
        """
        Save a NumPy array geometry to disk as a .npy file.

        :param path:      Output path (without extension).
        :param geometry:  NumPy-compatible array to save (shape: [N, 7]).

        :return: void – writes file to "{path}.npy".
        """
        np.save(f"{path}.npy", geometry)

    @staticmethod
    def load_geometry(name: str) -> np.typing.ArrayLike | None:
        """
        Load a previously saved NumPy geometry file.

        :param name: Path to file (without extension).

        :returns: Array-like geometry if found and valid, else None.
        """
        try:
            loaded_numpy = np.load(f"{name}.npy")
        except Exception:
            loaded_numpy = None
        return loaded_numpy

    @staticmethod
    def generate_sc_monomaterial_parallelepiped(
        lattice_constant: LatticeConstant,
        size: XYZ,
        material_number: int,
        initial_direction: XYZ | None = None,
        base_shift: XYZ | None = None,
        seed: int | None = None,
    ) -> np.typing.ArrayLike:
        """
        Generate a simple cubic (SC) monomaterial parallelepiped geometry.

        Each layer forms a rectangular grid in XY, stacked vertically in Z.

        @note
            The number of lattice points is determined by integer division: floor(size / spacing).
            The SC unit cell uses (a, a, c) spacing.

        :param lattice_constant:   Lattice spacing object (a, c) in meters.
        :param size:               Total system size in meters (XYZ).
        :param material_number:    Integer identifier for material.
        :param initial_direction:  Optional fixed spin direction (XYZ). If None, spins are randomized.
        :param base_shift:         Optional XYZ shift for the entire geometries spins.
        :param seed:               Optional RNG seed for reproducibility (random directions).

        :returns: Array of shape (N, 7) representing [x, y, z, Sx, Sy, Sz, material_id].
        """
        # 0) Local RNG —— deterministic if *seed* is provided.
        rng = np.random.default_rng(seed)

        # 1) Unpack lattice spacings.  In an SC cell a == b, but we allow a≠c.
        a, c = lattice_constant.a, lattice_constant.c

        # 2) Integer node counts per axis, at least one node each.
        nx = max(int(size.x / a), 1)
        ny = max(int(size.y / a), 1)
        nz = max(int(size.z / c), 1)

        # 3) Build 3‑D index grids without explicit Python loops.
        #    *ix*, *iy*, *iz* are integer ndarrays of identical shape
        #    (nx, ny, nz) that encode the Cartesian product of indices.
        #
        #   Get coordinates:
        #        ix[i, j, k] == <x>    -- along X (first dimension)
        #        iy[i, j, k] == <y>    -- along Y (second dimension)
        #        iz[i, j, k] == <z>    -- along Z (third dimension)
        #
        #    ────────────────────────────────────────────────────────────
        #    2‑D ILLUSTRATION
        #    -----------------------------------------------------------
        #    Suppose we have 1‑D vectors (blue = X index, red = Y index):
        #        x = [1 2 3 4]          # four X lattice points
        #        y = [7 6 5]            # three Y lattice points (top→bottom)
        #
        #    Call:  XX, YY = np.meshgrid(x, y, indexing="ij")
        #
        #    Combined logical grid  (pairs "x,y"):
        #
        #        +-----+-----+-----+-----+
        #  y=7 → | 1,7 | 2,7 | 3,7 | 4,7 |
        #        +-----+-----+-----+-----+
        #  y=6 → | 1,6 | 2,6 | 3,6 | 4,6 |
        #        +-----+-----+-----+-----+
        #  y=5 → | 1,5 | 2,5 | 3,5 | 4,5 |
        #        +-----+-----+-----+-----+
        #             x=1   x=2   x=3   x=4
        #
        #    Separate arrays returned by *meshgrid*:
        #        XX (=ix) = [[1 1 1]
        #                     [2 2 2]
        #                     [3 3 3]
        #                     [4 4 4]]
        #
        #        YY (=iy) = [[7 6 5]
        #                     [7 6 5]
        #                     [7 6 5]
        #                     [7 6 5]]
        #
        #    In 3‑D the pattern extends so that *iz* repeats planes
        #    0 … nz‑1 across every X,Y combination.
        #    ────────────────────────────────────────────────────────────
        ix, iy, iz = np.meshgrid(np.arange(nx), np.arange(ny), np.arange(nz), indexing="ij")

        # 4) Convert indices to physical coordinates R = (i·a, j·a, k·c).
        #    After stacking and reshaping we get a flat array POS of shape
        #    (N, 3) where N = nx·ny·nz.
        poses = np.stack((ix * a, iy * a, iz * c), axis=-1).reshape(-1, 3)  # combine into last axis + flat

        # 5) Optional global shift every spin.
        if base_shift is not None:
            poses += np.array([base_shift.x, base_shift.y, base_shift.z])

        # 6) Spin initialisation.
        if initial_direction is None:
            # Random components in U(‑1,1);
            dir = rng.uniform(-1.0, 1.0, size=(poses.shape[0], 3))
        else:
            # Copy the same vector for every node.
            dir = np.repeat(
                [[initial_direction.x, initial_direction.y, initial_direction.z]],
                poses.shape[0],
                axis=0,
            )

        # 7) Material ID column.
        mat = np.full((poses.shape[0], 1), material_number, dtype=np.int32)

        # 8) Concatenate into final array: [x y z Sx Sy Sz material].
        coords = np.hstack((poses, dir, mat))  # [N,7]

        return coords

    @staticmethod
    def generate_hcp_monomaterial_parallelepiped(
        lattice_constant: LatticeConstant,
        size: XYZ,
        material_number: int,
        initial_direction: XYZ | None = None,
        base_shift: XYZ | None = None,
        seed: int | None = None,
    ) -> np.typing.ArrayLike:
        """
        Generate a hexagonal close-packed (HCP) monomaterial parallelepiped geometry.

        Each layer forms a triangular lattice in XY with ABAB stacking in Z.

        @note
            The lattice consists of AB layers, where each B layer is offset by half a lattice constant
            in both X and Y directions to form a close-packed structure.

        :param lattice_constant:   Lattice spacing object (a, c) in meters.
        :param size:               Total system size in meters (XYZ).
        :param material_number:    Integer identifier for material.
        :param initial_direction:  Optional fixed spin direction (XYZ). If None, spins are randomized.
        :param base_shift:         Optional XYZ shift for the entire geometries spins.
        :param seed:               Optional RNG seed for reproducibility (random directions)

        :returns: Array of shape (N, 7) representing [x, y, z, Sx, Sy, Sz, material_id].
        """
        # 0) Local RNG —— deterministic if *seed* is provided.
        rng = np.random.default_rng(seed)

        # 1) Unpack lattice parameters.
        a = lattice_constant.a  # in‑plane lattice spacing [m]
        c = lattice_constant.c  # vertical period [m]
        dy = a * np.sqrt(3) / 2  # row spacing within plane [m]
        layer_h = c / 2  # height per half‑layer (A→B)

        # 2) Integer node counts along each axis.
        nx = max(int(size.x / a), 1)
        ny = max(int(size.y / dy), 1)
        nz = max(int(size.z / layer_h), 1)  # half‑layers (k‑index)

        # 3) Build 3‑D grids of indices (vectorised).
        ix, iy, iz = np.meshgrid(np.arange(nx), np.arange(ny), np.arange(nz), indexing="ij")

        # 4) Flatten index arrays for element‑wise arithmetic.
        ix = ix.ravel()
        iy = iy.ravel()
        iz = iz.ravel()

        # 5) Compute AB‑layer boolean mask and per‑node shifts.
        is_B = iz % 2 == 1  # True for B‑type half‑layers

        # 5a) X‑shift due to B layer.
        x_shift_layer = is_B * (a / 2)  # ½a for B, 0 for A
        # 5b) Y‑shift due to B layer.
        y_shift_layer = is_B * (dy / 3)  # dy/3 for B, 0 for A
        # 5c) X‑shift for every second in‑plane row (triangular lattice).
        x_shift_row = (iy % 2 == 1) * (a / 2)

        # 6) Convert indices + shifts to real coordinates.
        x = ix * a + x_shift_layer + x_shift_row  # ndarray (N,)
        y = iy * dy + y_shift_layer  # ndarray (N,)
        z = iz * layer_h  # ndarray (N,)
        pos = np.column_stack((x, y, z))  # shape (N, 3)

        # 7) Optional global base shift.
        if base_shift is not None:
            pos += np.array([base_shift.x, base_shift.y, base_shift.z])

        # 8) Spin initialisation
        if initial_direction is None:
            dir = rng.uniform(-1.0, 1.0, size=(pos.shape[0], 3))
        else:
            dir = np.repeat(
                [[initial_direction.x, initial_direction.y, initial_direction.z]],
                repeats=pos.shape[0],
                axis=0,
            )

        # 9) Material ID column.
        mat = np.full((pos.shape[0], 1), material_number, dtype=np.float64)

        # 10) Concatenate into final array: [x y z Sx Sy Sz material].
        coords = np.hstack((pos, dir, mat))

        return coords


__all__ = ["NumpyGeometryManager"]
