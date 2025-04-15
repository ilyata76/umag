from spindynapy.core.types import CartesianDirection, CartesianMoment, CartesianCoordinates  # type: ignore

import numpy as np
import sys

from spindynapy.core.geometries import CartesianGeometry
from spindynapy.core.simulation import Simulation
from spindynapy.core.solvers import CartesianLLGSolver, SphericalLLGSolver
from spindynapy.core.registries import MaterialRegistry, InteractionRegistry
from spindynapy.core.interactions import ExchangeInteraction
from spindynapy.core.types import Material

mat_reg = MaterialRegistry({1: Material(1, 1.1), 2: Material(2, 2.2)})
inter_reg = InteractionRegistry({1: ExchangeInteraction()})

numpy_geometry = CartesianGeometry(
    np.array(
        [[1, 1, 1, 1, 1, 1, 2], [1, 1, 1234123421, 1, 1, 1, 1], [1, 1, 1, 1, 1, 1, 1]]
    ),
    material_registry=mat_reg,
)
print(numpy_geometry)
print(
    x := CartesianGeometry(
        moments=[
            CartesianMoment(
                coordinates=CartesianCoordinates(x=1.000, y=1.000, z=1.000),
                direction=CartesianDirection(x=1.000, y=1.000, z=1.000),
                material=Material(1, 1),
            ),
            CartesianMoment(
                coordinates=CartesianCoordinates(x=1.000, y=1.000, z=1234123421.000),
                direction=CartesianDirection(x=1.000, y=1.000, z=1.000),
                material=mat_reg.get_element(1),
            ),
            CartesianMoment(
                coordinates=CartesianCoordinates(x=1.000, y=1.000, z=1.000),
                direction=CartesianDirection(x=1.000, y=1.000, z=1.000),
                material=mat_reg.get_element(1),
            ),
        ]
    )
)

print(x[1].__repr__())

Simulation(numpy_geometry, CartesianLLGSolver(), mat_reg, inter_reg)

# 1. Линейная цепочка вдоль оси X
# 5 моментов на линии x = 0, 1, 2, 3, 4; y = z = 0
# Направление момента (sx, sy, sz) = (1, 0, 0), материал = 1
linear_chain = np.array([
    [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1],
    [1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1],
    [2.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1],
    [3.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1],
    [4.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1],
])
# Ожидания для getNeighbors:
# - cutoff_radius = 1.1: для index=0 → [1], для index=2 → [1, 3]
# - cutoff_radius = 2.1: для index=0 → [1, 2], для index=2 → [0, 1, 3, 4]
# - cutoff_radius = 0.5: для всех → []

# 2. Кубическая решётка 2x2x2
# 8 моментов в узлах куба с шагом 1.0: (0,0,0), (0,0,1), ..., (1,1,1)
# Направление момента (sx, sy, sz) = (0, 0, 1), материал = 2
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
# Ожидания для getNeighbors:
# - cutoff_radius = 1.1: для index=0 → [1, 2, 4]
# - cutoff_radius = 1.5: для index=0 → [1, 2, 4], для index=7 → [3, 5, 6]
# - cutoff_radius = 2.0: для index=0 → [1, 2, 3, 4, 5, 6, 7]

# 3. Случайное расположение
# 6 моментов с произвольными координатами
# Направление момента (sx, sy, sz) = (1, 1, 1), материал = 1
random_points = np.array([
    [0.1, 0.2, 0.3, 1.0, 1.0, 1.0, 1],
    [0.5, 0.6, 0.7, 1.0, 1.0, 1.0, 1],
    [2.0, 2.0, 2.0, 1.0, 1.0, 1.0, 1],
    [2.1, 2.2, 2.3, 1.0, 1.0, 1.0, 1],
    [5.0, 5.0, 5.0, 1.0, 1.0, 1.0, 1],
    [10.0, 10.0, 10.0, 1.0, 1.0, 1.0, 1],
])
# Ожидания для getNeighbors:
# - cutoff_radius = 0.5: для index=0 → [1], для index=2 → [3]
# - cutoff_radius = 2.0: для index=0 → [1], для index=2 → [3]
# - cutoff_radius = 10.0: для index=0 → [1, 2, 3]

# Тестирование линейной цепочки
geom_linear = CartesianGeometry(linear_chain, mat_reg)
print("Linear chain neighbors (index=2, cutoff=1.1):", geom_linear.get_neighbors(2, 1.1))
print("Linear chain neighbors (index=2, cutoff=2.1):", geom_linear.get_neighbors(2, 2.1))

# Тестирование кубической решётки
geom_cubic = CartesianGeometry(cubic_lattice, mat_reg)
print("Cubic lattice neighbors (index=0, cutoff=1.1):", geom_cubic.get_neighbors(0, 1.1))
print("Cubic lattice neighbors (index=7, cutoff=1.5):", geom_cubic.get_neighbors(7, 1.5))

# Тестирование случайного расположения
geom_random = CartesianGeometry(random_points, mat_reg)
print("Random points neighbors (index=0, cutoff=0.5):", geom_random.get_neighbors(0, 0.5))
print("Random points neighbors (index=2, cutoff=10.0):", geom_random.get_neighbors(2, 10.0))

cubic_3x3x3 = np.zeros((27, 7))  # 3x3x3 = 27 точек
index = 0
for x in range(3):
    for y in range(3):
        for z in range(3):
            cubic_3x3x3[index] = [x, y, z, 0.0, 0.0, 1.0, 1]
            index += 1

geom = CartesianGeometry(cubic_3x3x3, mat_reg)
print("Geometry:\n", geom)

# Тестирование getNeighbors
print("Neighbors for index=13 (center [1,1,1], cutoff=1.1):", geom.get_neighbors(13, 1.1))
print("Neighbors for index=13 (center [1,1,1], cutoff=1.5):", geom.get_neighbors(13, 1.5))
print("Neighbors for index=0 (corner [0,0,0], cutoff=1.1):", geom.get_neighbors(0, 1.1))
print("Neighbors for index=26 (corner [2,2,2], cutoff=2.0):", geom.get_neighbors(26, 2.0))
