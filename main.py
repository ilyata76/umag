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
import sys
from enum import Enum
from os import makedirs
import numpy as np

from spindynapy.cartesian import (  # type: ignore # noqa
    MaterialRegistry,
    SolverStrategy,
    Geometry,
    LLGSolver,
    InteractionRegistry,
    ExchangeInteraction,
    ExternalInteraction,
    DemagnetizationInteraction,
    DipoleDipoleInteraction,
    AnisotropyInteraction,
    ThermalInteraction,
    Simulation,
    save_data,
    get_vvis_data,
    get_shot_data
)
from spindynapy.unit import XYZ, nano, femto  # type: ignore # noqa
from spindynapy.matlib import MaterialEnum, mat_lib, lattice_lib  # type: ignore # noqa
from spindynapy.geometry import NumpyGeometryManager  # type: ignore # noqa
from spindynapy.logger import ScopedTimer, logger


class InteractionEnum(Enum):
    EXCHANGE = 0
    DEMAGNETIZATION = 1
    ANISOTROPY = 2
    EXTERNAL = 3
    THERMAL = 4


path_dir = "./temp/___TEST___"
makedirs(path_dir, exist_ok=True)

logger.set_stream(sys.stdout)

with ScopedTimer("Настройка окружения", always_flush=True, logger=logger):
    mat_lib[MaterialEnum.COBALT].exchange_constant_J = 5.5e-21
    numpy_geometry = NumpyGeometryManager.generate_hcp_monomaterial_parallelepiped(
        lattice_constant=lattice_lib[MaterialEnum.COBALT],
        size=XYZ(nano(5.0), nano(5.0), nano(5.0)),
        material_number=mat_lib[MaterialEnum.COBALT].get_number(),
        initial_direction=None,
        base_shift=None,
    )
    material_registry = MaterialRegistry({MaterialEnum.COBALT.value: mat_lib[MaterialEnum.COBALT]})
    geometry = Geometry(moments=numpy_geometry, material_registry=material_registry)  # type:ignore
    solver = LLGSolver(strategy=SolverStrategy.HEUN)
    interaction_registry = InteractionRegistry({
        InteractionEnum.EXCHANGE.value: ExchangeInteraction(cutoff_radius=nano(0.4)),
        InteractionEnum.DEMAGNETIZATION.value: DipoleDipoleInteraction(cutoff_radius=nano(10), strategy="cutoff"),
        InteractionEnum.EXTERNAL.value: ExternalInteraction(1, 0, 0),
        InteractionEnum.ANISOTROPY.value: AnisotropyInteraction()
    })
    simulation = Simulation(
        geometry=geometry,
        solver=solver,
        material_registry=material_registry,
        interaction_registry=interaction_registry,
        dt=femto(1.0),
    )

with ScopedTimer("ПОЛНАЯ СИМУЛЯЦИЯ", always_flush=True, logger=logger):
    steps, save_every_step, update_macrocells_every_step = 15, 100, 1

    for i in range(1, steps + 1):
        simulation.simulate_one_step(
            save_step=(i % save_every_step == 0),
            update_macrocells=(i % update_macrocells_every_step == 0)
        )

    # simulation.simulate_many_steps(
    #     steps=steps,
    #     save_every_step=save_every_step,
    #     update_macrocells_every_step=update_macrocells_every_step
    # )
