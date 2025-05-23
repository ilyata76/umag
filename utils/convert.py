#!/usr/bin/env python3
"""Собрать координаты + спины (метры!) и сохранить atoms_spins.npy"""

import numpy as np
from pathlib import Path

Å = 1e-10                       # 1 ангстрем = 10⁻¹⁰ м
MAT = 1                         # номер материала (Co)

coords = np.loadtxt(Path("atoms-coords.data"), skiprows=1)
spins  = np.loadtxt(Path("spins-00000000.data"), skiprows=1)

assert coords.shape[0] == spins.shape[0], "Разное число строк!"

xyz_m = coords[:, [2, 3, 4]] * Å            # z x y → x y z, перевод в метры
spins  = spins.astype(np.float64)

mat = np.full((xyz_m.shape[0], 1), MAT, dtype=np.int32)

data = np.hstack([xyz_m, spins, mat])
data = np.ascontiguousarray(data, dtype=np.float64)
data[:, 6] = data[:, 6].astype(np.int32)     # столбец material — int

np.save("atoms_spins.npy", data)
print("✔ atoms_spins.npy:", data.shape,
      "| max coord %.2f nm" % (data[:, :3].max()*1e9))