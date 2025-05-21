import numpy as np

from .unit import XYZ, LatticeConstant


class NumpyGeometryManager:
    """
    Менеджер геометрии, использующий numpy для хранения и загрузки геометрии.
    А также для генерации геометрии (разные методы).
    """

    @staticmethod
    def save_geometry(name: str, geometry: np.typing.ArrayLike) -> None:
        """Сохранить геометрию в файл по пути name"""
        np.save(f"{name}.npy", geometry)

    @staticmethod
    def load_geometry(name: str) -> np.typing.ArrayLike | None:
        """Загрузить геометрию из файла по пути name"""
        try:
            loaded_numpy = np.load(f"{name}.npy")
        except Exception:
            loaded_numpy = None
        return loaded_numpy

    @staticmethod
    def generate_sc_monomaterial_parallelepiped(
        lattice_constant: XYZ,
        size: XYZ,
        material_number: int,
        initial_direction: XYZ | None = None,
        base_shift: XYZ | None = None,
    ) -> np.typing.ArrayLike:
        """
        Сгенерировать параллелепипед из одного материала с заданными параметрами.
        В каждом слое: прямоугольная подрешётка по XY

         - SC решётка
         - параллельные слои
        """
        nx, ny, nz = (
            int(size.x / lattice_constant.x) or 1,
            int(size.y / lattice_constant.y) or 1,
            int(size.z / lattice_constant.z) or 1,
        )
        total_nodes = nx * ny * nz
        coords, idx = np.zeros((total_nodes, 7)), 0  # [x, y, z, sx, sy, sz, material]
        for k in range(nz):
            for i in range(nx):
                for j in range(ny):
                    coords[idx, 0] = i * lattice_constant.x
                    coords[idx, 1] = j * lattice_constant.y
                    coords[idx, 2] = k * lattice_constant.z  # todo HCP, cubic etc.. сложнее
                    coords[idx, 3] = initial_direction.x if initial_direction else np.random.uniform(-1, 1)
                    coords[idx, 4] = initial_direction.y if initial_direction else np.random.uniform(-1, 1)
                    coords[idx, 5] = initial_direction.z if initial_direction else np.random.uniform(-1, 1)
                    coords[idx, 6] = material_number
                    idx += 1
        return coords

    @staticmethod
    def generate_hcp_monomaterial_parallelepiped(
        lattice_constant: LatticeConstant,
        size: XYZ,
        material_number: int,
        initial_direction: XYZ | None = None,
        base_shift: XYZ | None = None,
    ) -> np.typing.ArrayLike:
        """
        Сгенерировать параллелепипед из одного материала с заданными параметрами.
        В каждом слое: треугольная подрешётка по XY.

         - HCP решётка
         - ABAB слоистость (со смещением)
        """
        a = lattice_constant.a
        dy = a * np.sqrt(3) / 2
        layer_height = lattice_constant.c / 2  # половина кристаллографического c

        nx = int(size.x / a) or 1
        ny = int(size.y / dy) or 1
        nz = int(size.z / layer_height) or 1   # пересчитываем число слоёв

        coords_list = []

        for k in range(nz):
            z = k * layer_height
            shift_x = a / 2 if k % 2 == 1 else 0
            shift_y = dy / 3 if k % 2 == 1 else 0

            for j in range(ny):
                y = j * dy + shift_y
                for i in range(nx):
                    x = i * a + (a / 2 if j % 2 else 0) + shift_x

                    sx = initial_direction.x if initial_direction else np.random.uniform(-1, 1)
                    sy = initial_direction.y if initial_direction else np.random.uniform(-1, 1)
                    sz = initial_direction.z if initial_direction else np.random.uniform(-1, 1)

                    coords_list.append([x, y, z, sx, sy, sz, material_number])

        return np.array(coords_list)
