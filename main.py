from spindynapy.core.types.cartesian import CartesianDirection, CartesianSpin, CartesianCoordinates  # type: ignore

import numpy as np

print(
    CartesianSpin(
        CartesianCoordinates(np.array([3, 1, 1])),
        CartesianDirection(2, 2, 2)
    )
)
