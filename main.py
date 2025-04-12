from spindynapy.core.types.cartesian import CartesianDirection, CartesianSpin, CartesianCoordinates  # type: ignore

import numpy as np

from spindynapy.core.geometries.cartesian import CartesianGeometry
from spindynapy.core.simulation import Simulation


print(
    CartesianSpin(
        CartesianCoordinates(np.array([3, 1, 1])),
        CartesianDirection(1, 1, 1)
    )
)

geometry = CartesianGeometry(1)

sim = Simulation(geometry)

geometry.xxx = 100

del geometry

print(
    sim
)

# geometry undefined
