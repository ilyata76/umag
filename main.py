from spindynapy.core.types import CartesianDirection, CartesianMoment, CartesianCoordinates  # type: ignore

import numpy as np
import sys
from enum import Enum

from spindynapy.core.geometries import CartesianSimpleGeometry
from spindynapy.core.simulation import CartesianSimulation
from spindynapy.core.solvers import CartesianLLGSolver
from spindynapy.core.registries import MaterialRegistry
from spindynapy.core.interactions import CartesianInteractionRegistry, CartesianExchangeInteraction, CartesianExternalInteraction, CartesianAnisotropyInteraction
from spindynapy.core.types import Material, UniaxialAnisotropy

linear_chain = np.array([
    [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1],
    [1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1],
    [2.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1],
    [3.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1],
    [4.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1],
    *[
        [i, 0.0, 0.0, 1.0, 0.0, 0.0, 1] for i in range(5, 1000)
    ]
])

print("O!")

cubic_lattice = np.array([
    [0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 2],
    [0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 2],
    [0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 2],
    [0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 2],
    [1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 2],
    [1.0, 0.0, 1.0, 0.0, 0.0, 1.0, 2],
    [1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 2],
    [1.0, 1.0, 1.0, 0.0, 0.0, 1.0, 2],
])


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
    anisotropy=UniaxialAnisotropy(np.array([1, 0, 0]), 5e-24)
)

exchange_interaction = CartesianExchangeInteraction(cutoff_radius=5)
external_interaction = CartesianExternalInteraction(100, 100, 100)
anisotropy_interaction = CartesianAnisotropyInteraction()

material_registry = MaterialRegistry({MaterialEnum.COBALT.value: co_mat})
llg_solver = CartesianLLGSolver()
interaction_registry = CartesianInteractionRegistry(
    {InteractionEnum.EXCHANGE.value: exchange_interaction,
     InteractionEnum.EXTERNAL.value: external_interaction,
     InteractionEnum.ANISOTROPY.value: anisotropy_interaction
    }
)
geometry = CartesianSimpleGeometry(linear_chain, material_registry)

simulation = CartesianSimulation(
    geometry,
    llg_solver,
    material_registry,
    interaction_registry
)

import time

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

print(simulation.get_steps()[-1].as_string(False))

# with open("./a.txt", mode="w") as file:
#     file.write()

##print(simulation.get_steps()[2].total_fields)
