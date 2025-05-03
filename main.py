"""
Параметры геометрии:
    - размер ячейки:                            <...>
    - размер:                                   <...>
    - материал:                                 <...>
    - начальное направление:                    <...>

Параметры решателя:
    - стратегия:                                <...>
    - шаг:                                      <...>

Параметры взаимодействия:
    - обмен:                                    <...>
    - демагнетизация:                           <...>
    - анизотропия:                              <...>
    - внешнее поле:                             <...>

Параметры симуляции:
    - временной шаг:                            <...>
    - шаги:                                     <...>
    - сохранять каждый шаг:                     <...>
    - обновлять макроячейки каждый шаг:         <...>
"""

from enum import Enum
from os import makedirs

from spindynapy.cartesian import (  # type: ignore # noqa
    MaterialRegistry,
    SolverStrategy,
    Geometry,
    LLGSolver,
    InteractionRegistry,
    ExchangeInteraction,
    ExternalInteraction,
    DemagnetizationInteraction,
    AnisotropyInteraction,
    Simulation,
)
from spindynapy.unit import XYZ, nano, femto  # type: ignore # noqa
from spindynapy.matlib import MaterialEnum, mat_lib  # type: ignore # noqa
from spindynapy.geometry import NumpyGeometryManager  # type: ignore # noqa


class InteractionEnum(Enum):
    EXCHANGE = 0
    DEMAGNETIZATION = 1
    ANISOTROPY = 2
    EXTERNAL = 3


path_dir = "./temp/___COBALT_100x100_RELAX___"
makedirs(path_dir, exist_ok=True)

numpy_geometry = NumpyGeometryManager.load_geometry(f"{path_dir}/INITIAL")
if numpy_geometry is None or not numpy_geometry.any():  # type: ignore
    numpy_geometry = NumpyGeometryManager.generate_parallelepiped_monomaterial_geometry(
        lattice_constant=XYZ(nano(0.2507), nano(0.2507), nano(0.2507)),
        size=XYZ(nano(100), nano(100), 0),
        material_number=mat_lib["Co"].get_number(),
        initial_direction=None,
    )
    NumpyGeometryManager.save_geometry(f"{path_dir}/INITIAL", numpy_geometry)

material_registry = MaterialRegistry({MaterialEnum.COBALT.value: mat_lib["Co"]})

geometry = Geometry(numpy_geometry, material_registry, macrocell_size=nano(1))  # type:ignore

solver = LLGSolver(strategy=SolverStrategy.HEUN)

interaction_registry = InteractionRegistry(
    {
        InteractionEnum.EXCHANGE.value: ExchangeInteraction(cutoff_radius=nano(0.51)),
        InteractionEnum.DEMAGNETIZATION.value: DemagnetizationInteraction(
            cutoff_radius=nano(20), strategy="macrocells"
        ),
        InteractionEnum.ANISOTROPY.value: AnisotropyInteraction(),
        # InteractionEnum.EXTERNAL.value: ExternalInteraction(0.5, 0.05, 0.0),
    }
)

simulation = Simulation(
    geometry=geometry,
    solver=solver,
    material_registry=material_registry,
    interaction_registry=interaction_registry,
    dt=femto(1),
)

steps, save_every_step, update_macrocells_every_step = 1_000_000, 10_000, 10

for i in range(steps):
    simulation.simulate_one_step(
        save_step=(i % save_every_step == 0), update_macrocells=(i % update_macrocells_every_step == 0)
    )
    if i % save_every_step == 0:
        step_data = simulation.get_steps()[-1]
        with open(f"{path_dir}/sconfiguration-{i:08d}.vvis", mode="w") as file:
            file.writelines(step_data.vvisString())
        NumpyGeometryManager.save_geometry(f"{path_dir}/geometry-{i:08d}", step_data.geometry.as_numpy())
