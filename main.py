from spindynapy.core.types import CartesianDirection, CartesianMoment, CartesianCoordinates  # type: ignore

import numpy as np
import sys

from spindynapy.core.geometries import CartesianGeometry
from spindynapy.core.simulation import Simulation
from spindynapy.core.solvers import CartesianLLGSolver, SphericalLLGSolver
from spindynapy.core.registries import MaterialRegistry, InteractionRegistry
from spindynapy.core.interactions import ExchangeInteraction
from spindynapy.core.types import MagneticMaterial


geometry = CartesianGeometry(
    [
        CartesianMoment(
            CartesianCoordinates(np.array([3, 23402130491, 1])), CartesianDirection(1, 1, 1), 222222
        ),
        CartesianMoment(
            CartesianCoordinates(np.array([3, 1, 1])), CartesianDirection(1, 1, 1), 1
        ),
    ]
)
numpy_geometry = CartesianGeometry(
    np.array([[1, 1, 1, 1, 1, 1, 1], [1, 1, 1234123421, 1, 1, 1, 1], [1, 1, 1, 1, 1, 1, 1]])
)
print(str(numpy_geometry))

mat_reg = MaterialRegistry({1: MagneticMaterial(1.1)})
inter_reg = InteractionRegistry({1: ExchangeInteraction()})

sim = Simulation(geometry, CartesianLLGSolver(), mat_reg, inter_reg)
sim.simulate_many_steps(1)

print(sim.__repr__())

# geometry undefined
