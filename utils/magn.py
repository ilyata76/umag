import glob, re, numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

root = Path("./__COBALT_10x10__DIRECT___")   # <-- ваш path_dir
shot_files = sorted(root.rglob("*.shot"))       # rglob = рекурсивный поиск

if not shot_files:
    raise SystemExit(f"файлы .shot не найдены в {root}")

t_ps, mx, mz = [], [], []
time_pat = re.compile(r"TIME \(s\):\s*([+\-.\deE]+)")
mag_pat  = re.compile(
    r"M =\s*\(\s*([+\-.\deE]+)\s*,\s*([+\-.\deE]+)\s*,\s*([+\-.\deE]+)\s*\)"
)

for fn in shot_files:
    txt = Path(fn).read_text()

    # время
    t_sec = float(time_pat.search(txt).group(1))
    t_ps.append(t_sec * 1e12)           # -> пс

    # компоненты магнитизации
    mx_, my_, mz_ = map(float, mag_pat.search(txt).groups())
    mx.append(mx_);  mz.append(mz_)

t_ps, mx, mz = map(np.asarray, (t_ps, mx, mz))

plt.plot(t_ps, mx, "o-", label=r"$S_x$")
plt.plot(t_ps, mz, "^--", label=r"$S_z$")
plt.xlabel("Time (ps)")
plt.ylabel("Magnetisation ($m/m_s$)")
plt.xlim(0, 150)
plt.ylim(-1.1, 1.1)
plt.legend()
plt.tight_layout()
plt.show()
