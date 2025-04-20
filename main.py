import numpy as np
import time
from enum import Enum

from spindynapy.core.geometries.cartesian import Geometry  # type: ignore
from spindynapy.core.simulation.cartesian import Simulation  # type: ignore
from spindynapy.core.solvers.cartesian import LLGSolver  # type: ignore
from spindynapy.core.registries import MaterialRegistry  # type: ignore
from spindynapy.core.interactions.cartesian import (  # type: ignore
    InteractionRegistry,
    ExchangeInteraction,
    ExternalInteraction,
    AnisotropyInteraction,
)
from spindynapy.core.types import Material, UniaxialAnisotropy  # type: ignore

linear_chain = np.array(
    [
        [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1],
        [1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1],
        [2.0, 0.0, 0.0, -1.0, 0.0, 0.0, 1],
        [3.0, 0.0, 0.0, -1.0, 0.0, 0.0, 1],
        [4.0, 0.0, 0.0, -1.0, 0.0, 0.0, 1],
        [5.0, 0.0, 0.0, -1.0, 0.0, 0.0, 1],
        # *[[i, 0.0, 0.0, 1.0, 0.0, 0.0, 1] for i in range(5, 1000)],
    ]
)

print("O!")

cubic_lattice = np.array(
    [
        [0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 2],
        [0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 2],
        [0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 2],
        [0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 2],
        [1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 2],
        [1.0, 0.0, 1.0, 0.0, 0.0, 1.0, 2],
        [1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 2],
        [1.0, 1.0, 1.0, 0.0, 0.0, 1.0, 2],
    ]
)


class MaterialEnum(Enum):
    COBALT = 1


class InteractionEnum(Enum):
    EXCHANGE = 0
    EXTERNAL = 1
    ANISOTROPY = 2


co_mat = Material(
    MaterialEnum.COBALT.value,  # должен соответствовать с ID из MaterialRegistry!
    6.064e-21,
    1.72,
    anisotropy=UniaxialAnisotropy(np.array([1, 0, 0]), 1e-24),
)

exchange_interaction = ExchangeInteraction(cutoff_radius=5)
external_interaction = ExternalInteraction(100, 100, 100)
anisotropy_interaction = AnisotropyInteraction()

material_registry = MaterialRegistry({MaterialEnum.COBALT.value: co_mat})
llg_solver = LLGSolver()
interaction_registry = InteractionRegistry(
    {
        InteractionEnum.EXCHANGE.value: exchange_interaction,
        InteractionEnum.EXTERNAL.value: external_interaction,
        InteractionEnum.ANISOTROPY.value: anisotropy_interaction,
    }
)
geometry = Geometry(linear_chain, material_registry)

simulation = Simulation(geometry, llg_solver, material_registry, interaction_registry)

old = time.time()
simulation.simulate_one_step()

print(time.time() - old)
old = time.time()
simulation.simulate_one_step()

print(time.time() - old)
old = time.time()
simulation.simulate_one_step()

print(time.time() - old)
old = time.time()
simulation.simulate_many_steps(1000)

print(time.time() - old)

print(simulation.get_steps()[-1].as_string(True))

# with open("./a.txt", mode="w") as file:
#     file.write()

# print(simulation.get_steps()[2].total_fields)
