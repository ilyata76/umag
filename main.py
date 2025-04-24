import numpy as np
import time
from enum import Enum

from spindynapy.core.geometries.cartesian import Geometry  # type: ignore
from spindynapy.core.simulation.cartesian import Simulation  # type: ignore
from spindynapy.core.solvers.cartesian import LLGSolver  # type: ignore
from spindynapy.core.solvers import SolverStrategy
from spindynapy.core.registries import MaterialRegistry  # type: ignore
from spindynapy.core.interactions.cartesian import (  # type: ignore
    InteractionRegistry,
    ExchangeInteraction,
    ExternalInteraction,
    AnisotropyInteraction,
    DemagnetizationInteraction
)
from spindynapy.core.types import Material, UniaxialAnisotropy  # type: ignore

# Константы
nm = 1e-9
lattice_constant = 0.2507 * nm  # Размер ячейки (2.507 Å)


class MaterialEnum(Enum):
    COBALT = 1


co_mat = Material(
    material_number=MaterialEnum.COBALT.value,
    exchange_constant_J=6.064e-21,  # Дж
    atomic_magnetic_saturation_magnetization=1.72,  # в му_B
    damping_constant=0.5,
    anisotropy=UniaxialAnisotropy(np.array([1.0, 0.0, 0.0]), 6.69e-24),  # Дж
)

material_registry = MaterialRegistry({MaterialEnum.COBALT.value: co_mat})

nx = int(10 * nm / lattice_constant)
ny = int(10 * nm / lattice_constant)
total_nodes = nx * ny

try:
    coords = np.load("coords.npy")
except Exception:
    coords = None

if coords is None or not coords.any():
    # Генерация координат и случайных направлений спинов
    coords = np.zeros((total_nodes, 7))  # [x, y, z, sx, sy, sz, material]
    idx = 0
    for i in range(nx):
        for j in range(ny):
            coords[idx, 0] = i * lattice_constant  # x
            coords[idx, 1] = j * lattice_constant  # y
            coords[idx, 2] = 0.0  # z (однослойная пластина)
            # Случайное направление спина
            theta = np.random.uniform(0, np.pi)
            phi = np.random.uniform(0, 2 * np.pi)
            coords[idx, 3] = np.sin(theta) * np.cos(phi)  # sx
            coords[idx, 4] = np.sin(theta) * np.sin(phi)  # sy
            coords[idx, 5] = np.cos(theta)  # sz
            coords[idx, 6] = MaterialEnum.COBALT.value  # material ID
            idx += 1

np.save("coords.npy", coords)


geometry = Geometry(coords, material_registry)


class InteractionEnum(Enum):
    EXCHANGE = 0
    DEMAGNETIZATION = 1
    ANISOTROPY = 2
    EXTERNAL = 3


exchange_interaction = ExchangeInteraction(cutoff_radius=0.51 * nm)  # 2 шага решетки
demagnetization_interaction = DemagnetizationInteraction(cutoff_radius=0.51 * nm, strategy="cutoff")
anisotropy_interaction = AnisotropyInteraction()
external_interaction = ExternalInteraction(1000.0, 0.0, 0.0)

interaction_registry = InteractionRegistry(
    {
        InteractionEnum.EXCHANGE.value: exchange_interaction,
        InteractionEnum.DEMAGNETIZATION.value: demagnetization_interaction,
        InteractionEnum.ANISOTROPY.value: anisotropy_interaction,
        # InteractionEnum.EXTERNAL.value: external_interaction,
    }
)
llg_solver = LLGSolver(strategy=SolverStrategy.HEUN)
dt = 1e-16

simulation = Simulation(geometry, llg_solver, material_registry, interaction_registry, dt, use_openmp=True)


def main():
    try:
        # Запуск симуляции
        start_time = time.time()
        for i in range(50000):
            simulation.simulate_one_step(i % 1000 == 0)
            # elem = simulation.get_steps()[-1]
            # if i % 200 == 0:
            #     with open(f"./temp/sconfiguration-{i:08d}.vvis", mode="w") as file:
            #         file.writelines(elem.vvisString())
        # simulation.simulate_many_steps(100000, save_every_step=1000)  # 10 пс, снимки каждые 0.1 пс
        print(f"Simulation time: {time.time() - start_time:.2f} seconds")
        np.save("coords.npy", simulation.get_steps()[-1].geometry.as_numpy())
    finally:
        for index, elem in enumerate(simulation.get_steps()):
            with open(f"./temp/sconfiguration-{index:08d}.vvis", mode="w") as file:
                file.writelines(elem.vvisString())
            # with open(f"./temp/SHOT-{index}", mode="w") as file:
            #     file.writelines(elem.as_string(True)) 

# import signal, sys, os


# def signal_handler(sig, frame):
#     """Обработчик сигналов SIGINT и SIGTSTP."""
#     signal_name = "SIGTSTP" if sig == signal.SIGTSTP else "SIGINT"
#     print(f"\nCaught {signal_name}, saving simulation steps...", file=sys.stderr)
#     if sig == signal.SIGTSTP:
#         # Возвращаем стандартный обработчик SIGTSTP для приостановки процесса
#         signal.signal(signal.SIGTSTP, signal.SIG_DFL)
#         os.kill(os.getpid(), signal.SIGTSTP)
#     else:
#         sys.exit(1)

# # Установка обработчиков сигналов
# signal.signal(signal.SIGINT, signal_handler)   # Для Ctrl+C
# signal.signal(signal.SIGTSTP, signal_handler)  # Для Ctrl+Z


if __name__ == "__main__":
    main()
