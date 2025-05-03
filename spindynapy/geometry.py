import numpy as np

from .unit import XYZ


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
    def generate_parallelepiped_monomaterial_geometry(
        lattice_constant: XYZ,
        size: XYZ,
        material_number: int,
        initial_direction: XYZ | None = None,
    ) -> np.typing.ArrayLike:
        """Сгенерировать параллелепипед из одного материала с заданными параметрами"""
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
