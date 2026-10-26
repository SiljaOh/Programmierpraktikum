"""Overview figure for the README: Q3 of the GOR variants and an example prediction track.

    python3 tools/plot_overview.py docs/example_154l.txt docs/gor_overview.png
"""
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import ListedColormap  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402
import numpy as np  # noqa: E402

example, out = sys.argv[1:3]

# mean Q3 from our 5-fold cross-validation on CB513 (see README)
q3 = {"GOR III": 58.9, "GOR I": 61.4, "GOR V + III\n(MSA)": 66.2, "GOR V + I\n(MSA)": 67.1}

rows = {}
for line in open(example):
    key, _, seq = line.rstrip("\n").partition(" ")
    if key in ("SS", "gor1", "gor3", "gor4"):
        rows[key] = seq
true = rows["SS"]
tracks = [("DSSP (true)", true), ("GOR I", rows["gor1"]), ("GOR III", rows["gor3"]), ("GOR IV", rows["gor4"])]
code = {"H": 0, "E": 1, "C": 2, "-": 3}
cmap = ListedColormap(["#d62728", "#1f77b4", "#d9d9d9", "white"])

fig, (a, b) = plt.subplots(1, 2, figsize=(14, 3.8), gridspec_kw={"width_ratios": [1, 2.3]})
names = list(q3)
a.barh(names, list(q3.values()), color=["#9ecae1", "#6baed6", "#2171b5", "#08306b"])
for i, v in enumerate(q3.values()):
    a.text(v + 0.4, i, f"{v:.1f} %", va="center")
a.set(xlim=(50, 72), xlabel="Q3 accuracy (%)", title="5-fold cross-validation on CB513")
for s in ("top", "right"):
    a.spines[s].set_visible(False)

img = np.array([[code.get(c, 3) for c in s] for _, s in tracks])
b.imshow(img, cmap=cmap, aspect="auto", vmin=0, vmax=3, interpolation="nearest")
labels = []
for name, s in tracks:
    if name.startswith("DSSP"):
        labels.append(name)
        continue
    idx = [i for i, c in enumerate(s) if c != "-"]
    labels.append(f"{name}  (Q3 {100 * sum(s[i] == true[i] for i in idx) / len(idx):.0f} %)")
b.set_yticks(range(len(tracks)), labels)
b.set_xlabel("residue position")
b.set_title("Example: lysozyme (PDB 154L), predicted vs. DSSP secondary structure")
b.legend(handles=[Patch(fc="#d62728", label="helix (H)"), Patch(fc="#1f77b4", label="strand (E)"),
                  Patch(fc="#d9d9d9", label="coil (C)")], loc="upper center", bbox_to_anchor=(0.5, -0.28), ncol=3,
         frameon=False)
fig.tight_layout()
fig.savefig(out, dpi=150, bbox_inches="tight")
print("saved", out)
