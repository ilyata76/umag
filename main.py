import numpy as np
import time
import os
from enum import Enum

from spindynapy.core import constants  # type: ignore
from spindynapy.core.geometries.cartesian import Geometry  # type: ignore
from spindynapy.core.simulation.cartesian import Simulation  # type: ignore
from spindynapy.core.solvers.cartesian import AbstractSolver, LLGSolver  # type: ignore
from spindynapy.core.solvers import SolverStrategy  # type: ignore
from spindynapy.core.registries import MaterialRegistry  # type: ignore
from spindynapy.core.interactions.cartesian import (  # type: ignore # noqa
    InteractionRegistry,
    ExchangeInteraction,
    ExternalInteraction,  # noqa
    AnisotropyInteraction,
    DemagnetizationInteraction,
)
from spindynapy.core.types import Material, UniaxialAnisotropy  # type: ignore


class Unit:
    """TODO"""

    def __init__(self, factor):
        self.factor = factor

    def __call__(self, value):
        return value * self.factor


nano = Unit(1e-9)
pico = Unit(1e-12)
femto = Unit(1e-15)


class MaterialEnum(Enum):
    COBALT = 1


mat_lib: dict[str, Material] = {
    "Co": Material(
        material_number=MaterialEnum.COBALT.value,
        exchange_constant_J=6.064e-21,  # Дж
        atomic_magnetic_saturation_magnetization=1.72,  # в му_B
        damping_constant=0.5,
        gyromagnetic_ratio=constants.FREE_SPIN_GYROMAGNETIC_RATIO,
        anisotropy=UniaxialAnisotropy(np.array([-1.0, 0.0, 0.0]), 6.69e-24),  # Дж
    )
}


class XYZ:
    """TODO"""

    def __init__(self, x: float, y: float, z: float) -> None:
        self.x = x
        self.y = y
        self.z = z


class ThetaPhi:
    """TODO"""

    def __init__(self, theta: float, phi: float) -> None:
        self.theta = theta
        self.phi = phi


class NumpyGeometryManager:
    """TODO"""

    @staticmethod
    def save_geometry(name: str, geometry: np.typing.ArrayLike) -> None:
        """TODO"""
        np.save(f"{name}.npy", geometry)

    @staticmethod
    def load_geometry(name: str) -> np.typing.ArrayLike | None:
        """TODO"""
        try:
            loaded_numpy = np.load(f"{name}.npy")
        except Exception:
            loaded_numpy = None
        return loaded_numpy

    @staticmethod
    def generate_parallelepiped_monomaterial_geometry(
        lattice_constant: XYZ,
        size: XYZ,
        material: Material,
        initial_direction: XYZ | None = None,
    ) -> np.typing.ArrayLike:
        """TODO"""
        nx, ny, nz = (
            int(size.x / lattice_constant.x) or 1,
            int(size.y / lattice_constant.y) or 1,
            int(size.z / lattice_constant.z) or 1,
        )
        total_nodes = nx * ny * nz
        coords, idx = np.zeros((total_nodes, 7)), 0  # [x, y, z, sx, sy, sz, material]
        matnumber = material.get_number()
        for k in range(nz):
            for i in range(nx):
                for j in range(ny):
                    coords[idx, 0] = i * lattice_constant.x
                    coords[idx, 1] = j * lattice_constant.y
                    coords[idx, 2] = k * lattice_constant.z  # todo HCP, cubic etc.. сложнее
                    coords[idx, 3] = initial_direction.x if initial_direction else np.random.uniform(-1, 1)
                    coords[idx, 4] = initial_direction.y if initial_direction else np.random.uniform(-1, 1)
                    coords[idx, 5] = initial_direction.z if initial_direction else np.random.uniform(-1, 1)
                    coords[idx, 6] = matnumber
                    idx += 1
        return coords


class InteractionEnum(Enum):
    EXCHANGE = 0
    DEMAGNETIZATION = 1
    ANISOTROPY = 2
    EXTERNAL = 3


def simulate(
    simulation: Simulation, steps: int = 100, save_every_step: int = 2, update_macrocells_every_step: int = 2
) -> None:
    print("[PYTHON] Simulation started at", time.strftime("%H:%M:%S"))
    try:
        start_time = time.time()
        # шаги симуляции
        for i in range(steps):
            simulation.simulate_one_step(i % save_every_step == 0, i % update_macrocells_every_step == 0)
    finally:
        print(f"[PYTHON] Simulation time: {time.time() - start_time:.2f} seconds")


def simulate_and_save_every_step(
    simulation: Simulation,
    output_dir: str,
    steps: int = 100,
    save_every_step: int = 2,
    update_macrocells_every_step: int = 2,
) -> None:
    print("[PYTHON] Simulation started at", time.strftime("%H:%M:%S"))
    try:
        start_time = time.time()
        # шаги симуляции
        for i in range(steps):
            simulation.simulate_one_step(i % save_every_step == 0, i % update_macrocells_every_step == 0)
            if i % save_every_step == 0:
                with open(f"{output_dir}/sconfiguration-{i:08d}.vvis", mode="w") as file:
                    file.writelines(simulation.get_steps()[-1].vvisString())
    finally:
        print(f"[PYTHON] Simulation time: {time.time() - start_time:.2f} seconds")


def process(
    material_registry: MaterialRegistry,
    interaction_registry: InteractionRegistry,
    codename: str,
    solver: AbstractSolver,
    time_step: float,
    generate_geometry_args: dict | None = None,
    macrocell_size: float = nano(1),
    steps: int = 1000,
    save_every_step: int = 100,
    update_macrocells_every_step: int = 100,
    load_geometry_path: str | None = None,  # можно загрузить геометрию после другой симуляции
    save_geometry_path: str | None = None,  # сохранить геометрию в файл
) -> None:
    """TODO это - фасадная функция так сказать"""

    # --- В КАКУЮ ПАПКУ ПИСАТЬ ---

    output_dir = f"./temp/{codename}"
    os.makedirs(output_dir, exist_ok=True)

    # --- ГЕОМЕТРИЯ ---

    print(f"[PYTHON] Reading geometry from {load_geometry_path or f'{output_dir}/INITIAL'}")
    numpy_geometry = NumpyGeometryManager.load_geometry(load_geometry_path or f"{output_dir}/INITIAL")
    if numpy_geometry is None or not numpy_geometry.any():  # type: ignore
        print("[PYTHON] Geometry not found, generating new one")
        if not generate_geometry_args:
            raise ValueError("[PYTHON] No geometry to load and no geometry generation parameters provided.")
        numpy_geometry = NumpyGeometryManager.generate_parallelepiped_monomaterial_geometry(
            **generate_geometry_args
        )
        NumpyGeometryManager.save_geometry(save_geometry_path or f"{output_dir}/INITIAL", numpy_geometry)

    geometry = Geometry(numpy_geometry, material_registry, macrocell_size=macrocell_size)  # type: ignore

    # --- СИМУЛЯЦИЯ ---

    simulation = Simulation(
        geometry=geometry,
        solver=solver,
        material_registry=material_registry,
        interaction_registry=interaction_registry,
        dt=time_step
    )

    try:
        simulate_and_save_every_step(
            simulation,
            output_dir,
            steps=steps,
            save_every_step=save_every_step,
            update_macrocells_every_step=update_macrocells_every_step,
        )
    finally:
        # --- СОХРАНЕНИЕ РЕЗУЛЬТАТОВ ---
        # print(f"[PYTHON] Saving results for {simulation.get_steps().__len__()} steps")
        # for index, elem in enumerate(simulation.get_steps()):
        #     with open(f"{output_dir}/sconfiguration-{index:08d}.vvis", mode="w") as file:
        #         file.writelines(elem.vvisString())

        print(f"[PYTHON] Saving geometry in {output_dir}/RESULT as numpy array")
        NumpyGeometryManager.save_geometry(
            f"{output_dir}/RESULT",
            simulation.get_steps()[-1].geometry.as_numpy(),
        )

    # --- КОНЕЦ ---


if __name__ == "__main__":

    # --- ВЗАИМОДЕЙСТВИЯ ---

    interactions = {
        InteractionEnum.EXCHANGE.value: ExchangeInteraction(cutoff_radius=nano(0.51)),
        InteractionEnum.DEMAGNETIZATION.value: DemagnetizationInteraction(
            cutoff_radius=nano(20), strategy="macrocells"
        ),
        InteractionEnum.ANISOTROPY.value: AnisotropyInteraction(),
        # InteractionEnum.EXTERNAL.value: ExternalInteraction(0.5, 0.05, 0.0),
    }

    # --- МАТЕРИАЛЫ ---

    material_registry = MaterialRegistry({MaterialEnum.COBALT.value: mat_lib["Co"]})

    # --- РЕШАТЕЛЬ ---

    llg_solver = LLGSolver(strategy=SolverStrategy.HEUN)

    # --- ЗАПУСК ---

    print("[PYTHON] Starting simulation with MACROCELLS demag")
    # process(
    #     material_registry=material_registry,
    #     interaction_registry=InteractionRegistry(interactions),
    #     codename="_120x120x0_Co_Relaxation_macrocells_CONTINUES_(4)",
    #     solver=llg_solver,
    #     time_step=femto(1),
    #     steps=1_000_000,
    #     save_every_step=5_000,
    #     update_macrocells_every_step=100,
    #     macrocell_size=nano(1),
    #     generate_geometry_args={
    #         "lattice_constant": XYZ(nano(0.2507), nano(0.2507), nano(0.2507)),
    #         "size": XYZ(nano(120), nano(120), 0),
    #         "material": mat_lib["Co"],
    #         "initial_direction": None,
    #     },
    #     load_geometry_path="./temp/_120x120x0_Co_Relaxation_macrocells_(3)/RESULT"
    # )

    process(
        material_registry=material_registry,
        interaction_registry=InteractionRegistry(interactions),
        codename="100x100x0_Co_Relaxation_macrocells_30nmcutoff___",
        solver=llg_solver,
        time_step=femto(1),
        steps=1_000_000,
        save_every_step=500,
        update_macrocells_every_step=10,
        macrocell_size=nano(2),
        generate_geometry_args={
            "lattice_constant": XYZ(nano(0.2507), nano(0.2507), nano(0.2507)),
            "size": XYZ(nano(100), nano(100), 0),
            "material": mat_lib["Co"],
            "initial_direction": None,
        },
        # load_geometry_path="./temp/____TEST4____40x40x0_Co_Relaxation_macrocells_____/INITIAL",
    )
