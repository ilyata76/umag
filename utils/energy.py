import re
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

folder = Path(".")
energy_pattern = re.compile(r"ENERGY\[\d+\|(\w+)\]:\s+([-+eE0-9.]+)")
time_pattern = re.compile(r"TIME \(s\):\s+([-+eE0-9.]+)")
step_pattern = re.compile(r"STEP:\s+(\d+)")

data = []

for file in sorted(folder.glob("stepdata-*.shot")):
    with open(file) as f:
        content = f.read()

    energies = dict(re.findall(energy_pattern, content))
    time_match = time_pattern.search(content)
    step_match = step_pattern.search(content)

    if time_match and step_match:
        time = float(time_match.group(1))
        step = int(step_match.group(1))
        data.append({
            "time_ps": time * 1e12,
            "step": step,
            "exchange": float(energies.get("EXCHANGE", 0)),
            "demag": float(energies.get("DEMAGNETIZATION", 0)) * 1e3,
            "anis": float(energies.get("ANISOTROPY", 0)) * 1e3,
        })

df = pd.DataFrame(data).sort_values("step")

# Добавляем суммарную энергию
# df["total"] = df["exchange"] + df["demag"] / 1e3 + df["anis"] / 1e3

# Построение графика
plt.figure(figsize=(10, 6))

markers = {"exchange": "o", "demag": "s", "anis": "D"}
labels = {
    "exchange": "Exchange",
    "demag": "Demagnetization *1e3",
    "anis": "Anisotropy *1e3"
}
colors = {"exchange": "C0", "demag": "C1", "anis": "C2"}

for key in ["exchange", "demag", "anis"]:
    plt.plot(df["step"], df[key], label=labels[key], marker=markers[key],
             markersize=4, linewidth=1, color=colors[key])

# plt.plot(df["step"], df["total"], label="Total Energy", color="black", linestyle="--", linewidth=1.2)

plt.xlabel("Simulation step")
plt.ylabel("Energy (J)")
plt.title("Energy evolution per step")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()
