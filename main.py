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
    save_data,
    get_vvis_data,
    get_shot_data
)
from spindynapy.unit import XYZ, nano, femto  # type: ignore # noqa
from spindynapy.matlib import MaterialEnum, mat_lib  # type: ignore # noqa
from spindynapy.geometry import NumpyGeometryManager  # type: ignore # noqa


class InteractionEnum(Enum):
    EXCHANGE = 0
    DEMAGNETIZATION = 1
    ANISOTROPY = 2
    EXTERNAL = 3


path_dir = "./temp/_____TSTS_____"
makedirs(path_dir, exist_ok=True)

numpy_geometry = NumpyGeometryManager.load_geometry(f"{path_dir}/INITIAL")
if numpy_geometry is None or not numpy_geometry.any():  # type: ignore
    numpy_geometry = NumpyGeometryManager.generate_hcp_monomaterial_parallelepiped(
        lattice_constant=XYZ(nano(0.2507), nano(0.2507), nano(0.408)),
        size=XYZ(nano(10), nano(10), nano(0.1)),  # монослой
        material_number=mat_lib["Co"].get_number(),
        initial_direction=None,
        base_shift=None,
    )
    NumpyGeometryManager.save_geometry(f"{path_dir}/INITIAL", numpy_geometry)

material_registry = MaterialRegistry({MaterialEnum.COBALT.value: mat_lib["Co"]})

geometry = Geometry(numpy_geometry, material_registry, macrocell_size=nano(1))  # type:ignore

solver = LLGSolver(strategy=SolverStrategy.HEUN)

interaction_registry = InteractionRegistry(
    {
        InteractionEnum.EXCHANGE.value: ExchangeInteraction(cutoff_radius=nano(0.45)),
        InteractionEnum.DEMAGNETIZATION.value: DemagnetizationInteraction(
            cutoff_radius=nano(4), strategy="cutoff"
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

steps, save_every_step, update_macrocells_every_step = 500_000, 2500, 10

for i in range(steps):
    simulation.simulate_one_step(
        save_step=(i % save_every_step == 0), update_macrocells=(i % update_macrocells_every_step == 0)
    )
    if i % save_every_step == 0:
        step_data = simulation.get_steps()[-1]
        save_data(
            get_vvis_data(step_data, material_registry, interaction_registry),
            f"{path_dir}/sconfiguration-{i:08d}.vvis"
        )
        NumpyGeometryManager.save_geometry(f"{path_dir}/geometry-{i:08d}", step_data.geometry.as_numpy())
        save_data(get_shot_data(step_data, material_registry, interaction_registry),
                  f"{path_dir}/stepdata-{i:08d}.shot")
