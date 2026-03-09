import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
import numpy as np
from matplotlib.colors import LinearSegmentedColormap

rows = [
    ("gor1","n/a",3,2,61.5,60.3),("gor1","n/a",3,3,62.0,60.2),("gor1","n/a",3,4,61.9,58.3),
    ("gor1","n/a",4,2,61.5,60.1),("gor1","n/a",4,3,62.1,59.8),("gor1","n/a",4,4,61.9,57.7),
    ("gor1","n/a",5,2,61.1,59.5),("gor1","n/a",5,3,61.7,59.2),("gor1","n/a",5,4,61.5,56.9),
    ("gor1","n/a",6,2,60.7,59.2),("gor1","n/a",6,3,61.2,58.9),("gor1","n/a",6,4,61.1,56.4),
    ("gor3",0.01,3,2,59.9,56.5),("gor3",0.01,3,3,59.8,52.9),("gor3",0.01,3,4,59.7,50.2),
    ("gor3",0.01,4,2,58.9,54.3),("gor3",0.01,4,3,58.8,50.0),("gor3",0.01,4,4,58.7,46.8),
    ("gor3",0.01,5,2,58.4,52.5),("gor3",0.01,5,3,58.3,48.0),("gor3",0.01,5,4,58.2,44.7),
    ("gor3",0.01,6,2,57.2,51.2),("gor3",0.01,6,3,57.1,46.4),("gor3",0.01,6,4,57.0,43.0),
    ("gor3",0.1,3,2,59.9,56.5),("gor3",0.1,3,3,59.7,52.8),("gor3",0.1,3,4,59.7,50.1),
    ("gor3",0.1,4,2,58.8,54.2),("gor3",0.1,4,3,58.6,49.8),("gor3",0.1,4,4,58.6,46.7),
    ("gor3",0.1,5,2,58.2,52.4),("gor3",0.1,5,3,58.0,47.8),("gor3",0.1,5,4,58.0,44.6),
    ("gor3",0.1,6,2,56.9,50.9),("gor3",0.1,6,3,56.8,46.1),("gor3",0.1,6,4,56.7,42.8),
    ("gor3",0.5,3,2,58.7,55.7),("gor3",0.5,3,3,59.0,52.1),("gor3",0.5,3,4,58.7,49.1),
    ("gor3",0.5,4,2,57.9,53.9),("gor3",0.5,4,3,58.2,49.8),("gor3",0.5,4,4,57.9,46.4),
    ("gor3",0.5,5,2,57.2,52.0),("gor3",0.5,5,3,57.5,47.7),("gor3",0.5,5,4,57.2,44.1),
    ("gor3",0.5,6,2,56.0,50.5),("gor3",0.5,6,3,56.3,46.0),("gor3",0.5,6,4,56.0,42.3),
    ("gor3",0.8,3,2,58.1,55.2),("gor3",0.8,3,3,58.6,51.9),("gor3",0.8,3,4,58.6,49.2),
    ("gor3",0.8,4,2,57.2,53.1),("gor3",0.8,4,3,57.7,49.2),("gor3",0.8,4,4,57.7,46.1),
    ("gor3",0.8,5,2,56.6,51.3),("gor3",0.8,5,3,57.1,47.1),("gor3",0.8,5,4,57.1,43.8),
    ("gor3",0.8,6,2,55.4,49.8),("gor3",0.8,6,3,55.9,45.4),("gor3",0.8,6,4,55.9,42.1),
    ("gor3",1.0,3,2,57.8,55.0),("gor3",1.0,3,3,58.5,51.8),("gor3",1.0,3,4,58.5,49.0),
    ("gor3",1.0,4,2,56.8,52.9),("gor3",1.0,4,3,57.5,49.1),("gor3",1.0,4,4,57.6,45.9),
    ("gor3",1.0,5,2,56.3,51.3),("gor3",1.0,5,3,57.0,47.1),("gor3",1.0,5,4,57.1,43.9),
    ("gor3",1.0,6,2,55.2,49.9),("gor3",1.0,6,3,55.9,45.6),("gor3",1.0,6,4,55.9,42.3),
    ("gor3",2.0,3,2,55.6,53.0),("gor3",2.0,3,3,56.8,51.1),("gor3",2.0,3,4,57.1,48.5),
    ("gor3",2.0,4,2,54.5,50.8),("gor3",2.0,4,3,55.7,48.2),("gor3",2.0,4,4,56.0,45.2),
    ("gor3",2.0,5,2,53.8,49.1),("gor3",2.0,5,3,55.0,46.1),("gor3",2.0,5,4,55.4,43.0),
    ("gor3",2.0,6,2,52.8,48.0),("gor3",2.0,6,3,54.1,44.9),("gor3",2.0,6,4,54.4,41.7),
    ("gor4",0.01,3,2,56.4,43.1),("gor4",0.01,3,3,55.6,38.5),("gor4",0.01,3,4,54.7,35.6),
    ("gor4",0.01,4,2,54.7,39.7),("gor4",0.01,4,3,53.9,34.6),("gor4",0.01,4,4,53.0,31.5),
    ("gor4",0.01,5,2,53.5,36.8),("gor4",0.01,5,3,52.7,31.5),("gor4",0.01,5,4,51.8,28.2),
    ("gor4",0.01,6,2,51.8,33.6),("gor4",0.01,6,3,51.0,28.1),("gor4",0.01,6,4,50.1,24.8),
    ("gor4",0.1,3,2,57.8,48.3),("gor4",0.1,3,3,57.1,43.3),("gor4",0.1,3,4,56.1,38.9),
    ("gor4",0.1,4,2,56.8,46.2),("gor4",0.1,4,3,56.1,40.7),("gor4",0.1,4,4,55.1,36.1),
    ("gor4",0.1,5,2,55.4,43.4),("gor4",0.1,5,3,54.7,37.7),("gor4",0.1,5,4,53.7,32.8),
    ("gor4",0.1,6,2,53.9,40.3),("gor4",0.1,6,3,53.2,34.4),("gor4",0.1,6,4,52.2,29.5),
    ("gor4",0.5,3,2,51.8,49.4),("gor4",0.5,3,3,53.4,47.9),("gor4",0.5,3,4,53.7,44.2),
    ("gor4",0.5,4,2,51.3,48.2),("gor4",0.5,4,3,52.9,45.6),("gor4",0.5,4,4,53.1,41.3),
    ("gor4",0.5,5,2,50.3,46.5),("gor4",0.5,5,3,51.9,43.6),("gor4",0.5,5,4,52.1,39.1),
    ("gor4",0.5,6,2,49.1,44.8),("gor4",0.5,6,3,50.7,41.7),("gor4",0.5,6,4,51.0,37.1),
    ("gor4",0.8,3,2,40.8,38.1),("gor4",0.8,3,3,44.1,40.7),("gor4",0.8,3,4,45.8,38.6),
    ("gor4",0.8,4,2,40.4,37.5),("gor4",0.8,4,3,43.7,39.5),("gor4",0.8,4,4,45.5,37.0),
    ("gor4",0.8,5,2,38.6,35.0),("gor4",0.8,5,3,41.8,36.7),("gor4",0.8,5,4,43.6,34.0),
    ("gor4",0.8,6,2,37.2,33.1),("gor4",0.8,6,3,40.5,34.5),("gor4",0.8,6,4,42.3,31.6),
    ("gor4",1.0,3,2,35.9,32.5),("gor4",1.0,3,3,38.8,34.8),("gor4",1.0,3,4,40.7,33.7),
    ("gor4",1.0,4,2,35.2,31.3),("gor4",1.0,4,3,38.1,33.5),("gor4",1.0,4,4,40.0,31.8),
    ("gor4",1.0,5,2,33.9,29.3),("gor4",1.0,5,3,36.8,31.4),("gor4",1.0,5,4,38.7,29.5),
    ("gor4",1.0,6,2,33.1,28.2),("gor4",1.0,6,3,36.0,30.1),("gor4",1.0,6,4,38.0,28.0),
    ("gor4",2.0,3,2,26.0,20.4),("gor4",2.0,3,3,28.4,23.4),("gor4",2.0,3,4,30.6,23.3),
    ("gor4",2.0,4,2,26.1,20.2),("gor4",2.0,4,3,28.5,23.1),("gor4",2.0,4,4,30.8,23.3),
    ("gor4",2.0,5,2,26.0,20.1),("gor4",2.0,5,3,28.5,22.9),("gor4",2.0,5,4,30.7,23.1),
    ("gor4",2.0,6,2,26.0,20.1),("gor4",2.0,6,3,28.5,22.9),("gor4",2.0,6,4,30.7,23.1),
]
df = pd.DataFrame(rows, columns=["method","pseudo","minH","minE","Q3","SOV"])
df['pseudo_num'] = pd.to_numeric(df['pseudo'], errors='coerce')

BEST = {
    'GOR1': dict(pseudo=None,  minH=3, minE=2, Q3=61.5, SOV=60.3),
    'GOR3': dict(pseudo=0.01,  minH=3, minE=2, Q3=59.9, SOV=56.5),
    'GOR4': dict(pseudo=0.1,   minH=3, minE=2, Q3=57.8, SOV=48.3),
}

C1   = "#2d6a4f"
C3   = "#52b788"
C4   = "#1a6896"
STAR = "#d62828"
DARK = "#1a1a2e"
LIGHT_BG = "#fafbfc"

cmap_custom = LinearSegmentedColormap.from_list(
    "bio", ["#e8f5f0", "#52b788", "#1b4332"], N=256)

plt.rcParams.update({
    'font.family': 'DejaVu Sans',
    'figure.facecolor': LIGHT_BG,
    'axes.facecolor': 'white',
    'axes.spines.top': False,
    'axes.spines.right': False,
    'axes.spines.left': True,
    'axes.spines.bottom': True,
    'axes.labelsize': 11,
    'axes.titlesize': 13,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'axes.titlepad': 12,
})

minH_vals = [3,4,5,6]
minE_vals = [2,3,4]


fig1 = plt.figure(figsize=(15, 5.5), facecolor=LIGHT_BG)
fig1.suptitle("Einfluss von minH und minE auf den Q3-Score",
              fontsize=16, fontweight='bold', color=DARK, y=1.01)

gs = gridspec.GridSpec(1, 3, figure=fig1, wspace=0.35)

configs = [
    ("gor1", None, "GOR1",  "(kein Pseudocount)"),
    ("gor3", 0.01, "GOR3",  "(pseudo = 0.01)"),
    ("gor4", 0.1,  "GOR4",  "(pseudo = 0.1)"),
]

for idx, (method, pseudo, mtitle, psub) in enumerate(configs):
    ax = fig1.add_subplot(gs[idx])
    if pseudo is None:
        sub = df[df['method']==method]
    else:
        sub = df[(df['method']==method)&(df['pseudo_num']==pseudo)]
    pivot = sub.pivot_table(index='minH', columns='minE', values='Q3')
    pivot = pivot.reindex(index=minH_vals, columns=minE_vals)

    # Gemeinsame Skala über alle 3 (55–63)
    im = ax.imshow(pivot.values, cmap=cmap_custom, vmin=55, vmax=63,
                   aspect='auto', origin='upper')

    ax.set_xticks(range(len(minE_vals)))
    ax.set_xticklabels([f"minE={v}" for v in minE_vals], fontsize=9)
    ax.set_yticks(range(len(minH_vals)))
    ax.set_yticklabels([f"minH={v}" for v in minH_vals], fontsize=9)

    # Spines dezent
    for spine in ax.spines.values():
        spine.set_linewidth(0.8)
        spine.set_color('#ccc')

    bH = BEST[method.upper()]['minH']
    bE = BEST[method.upper()]['minE']
    bQ = BEST[method.upper()]['Q3']
    bS = BEST[method.upper()]['SOV']

    for i, mH in enumerate(minH_vals):
        for j, mE in enumerate(minE_vals):
            try:
                val = pivot.loc[mH, mE]
            except:
                continue
            if np.isnan(val): continue
            is_best = (mH==bH and mE==bE)
            ax.text(j, i, f"{val:.1f}",
                    ha='center', va='center', fontsize=10,
                    fontweight='bold' if is_best else 'normal',
                    color='white' if val > 59 else DARK)
            if is_best:
                ax.add_patch(plt.Rectangle(
                    (j-.48, i-.48), .96, .96,
                    fill=False, edgecolor=STAR, linewidth=2.8, zorder=10,
                    linestyle='-'))

    cbar = fig1.colorbar(im, ax=ax, shrink=0.82, pad=0.03)
    cbar.set_label('Q3 (%)', fontsize=9)
    cbar.ax.tick_params(labelsize=8)

    ax.set_title(f"{mtitle}", fontsize=14, fontweight='bold', color=DARK, pad=4)
    ax.text(0.5, -0.18,
            f"★  Best: minH={bH}, minE={bE}   Q3 = {bQ:.1f}%   SOV = {bS:.1f}%",
            transform=ax.transAxes, ha='center', fontsize=9,
            color=STAR, fontweight='bold',
            bbox=dict(boxstyle='round,pad=0.35', fc='#fff0f0',
                      ec=STAR, lw=1.0, alpha=0.9))
    ax.text(0.5, -0.28, psub, transform=ax.transAxes,
            ha='center', fontsize=8.5, color='#777', style='italic')

plt.savefig("plot1_minHE.png", dpi=160, bbox_inches='tight')
print("Saved: plot1_minHE.png")


fig2, ax = plt.subplots(figsize=(9, 5.5), facecolor=LIGHT_BG)
ax.set_facecolor('white')

PSEUDO_METHODS = [
    ("gor3", C3,  "o", "-",  "GOR3"),
    ("gor4", C4,  "s", "--", "GOR4"),
]

for method, color, marker, ls, label in PSEUDO_METHODS:
    sub = df[(df['method']==method)&(df['minH']==3)&(df['minE']==2)].dropna(subset=['pseudo_num'])
    sub = sub.sort_values('pseudo_num')
    # Q3 Linie
    ax.plot(sub['pseudo_num'], sub['Q3'], color=color, marker=marker,
            linewidth=2.5, markersize=8, label=f"{label}  Q3", linestyle=ls, zorder=4)
    # SOV gestrichelt
    ax.plot(sub['pseudo_num'], sub['SOV'], color=color, marker=marker,
            linewidth=1.4, markersize=5, linestyle=':', alpha=0.55,
            label=f"{label}  SOV", zorder=3)
    # Schattierung Q3/SOV Lücke
    ax.fill_between(sub['pseudo_num'], sub['SOV'], sub['Q3'],
                    color=color, alpha=0.08)

    # Besten Punkt
    bp = BEST[method.upper()]['pseudo']
    br = sub[sub['pseudo_num']==bp]
    if not br.empty:
        bq = br['Q3'].values[0]
        bs = br['SOV'].values[0]
        ax.scatter([bp], [bq], s=160, color=STAR, zorder=7,
                   edgecolors='white', linewidths=1.5)
        offset_x = 0.12 if method=='gor3' else 0.12
        offset_y = 2.0  if method=='gor3' else -4.5
        ax.annotate(
            f"★  {label} best\nPseudo = {bp}\nQ3 = {bq:.1f}%   SOV = {bs:.1f}%",
            xy=(bp, bq), xytext=(bp+offset_x, bq+offset_y),
            fontsize=8.5, color=STAR, fontweight='bold',
            bbox=dict(boxstyle='round,pad=0.4', fc='#fff0f0', ec=STAR, lw=1, alpha=0.92),
            arrowprops=dict(arrowstyle='->', color=STAR, lw=1.4))

# GOR1 Referenz
ax.axhline(61.5, color=C1, lw=2.0, ls='-.', zorder=2, label='GOR1  Q3 = 61.5%')
ax.axhline(60.3, color=C1, lw=1.2, ls=':', alpha=0.55, zorder=2, label='GOR1  SOV = 60.3%')
ax.text(1.85, 61.9, "GOR1 Q3",  color=C1, fontsize=8.5, fontweight='bold')
ax.text(1.85, 60.7, "GOR1 SOV", color=C1, fontsize=8,   alpha=0.7)

ax.set_xlabel("Pseudocount", fontweight='bold', fontsize=11)
ax.set_ylabel("Score (%)", fontweight='bold', fontsize=11)
ax.set_title("Einfluss des Pseudocounts auf Q3 und SOV\n(minH = 3, minE = 2)",
             fontweight='bold', color=DARK)
ax.set_xlim(-0.05, 2.25)
ax.set_ylim(15, 68)
ax.yaxis.grid(True, alpha=0.2, linestyle='--')
ax.set_axisbelow(True)

leg = ax.legend(fontsize=9, frameon=True, framealpha=0.95,
                ncol=2, loc='lower left',
                edgecolor='#ddd')

ax.text(0.5, -0.13,
        "Durchgezogene Linien = Q3  ·  Gestrichelt = SOV  ·  Fläche = Q3–SOV-Lücke  ·  ★ = bester Pseudocount",
        transform=ax.transAxes, ha='center', fontsize=8, color='#888', style='italic')

plt.tight_layout(rect=[0,0.05,1,1])
plt.savefig("plot2_pseudocount.png", dpi=160, bbox_inches='tight')
print("Saved: plot2_pseudocount.png")


fig3, ax = plt.subplots(figsize=(9, 6), facecolor=LIGHT_BG)
ax.set_facecolor('white')

methods  = ['GOR1', 'GOR3', 'GOR4']
q3s      = [BEST[m]['Q3']  for m in methods]
sovs     = [BEST[m]['SOV'] for m in methods]
cols     = [C1, C3, C4]
x = np.arange(len(methods))
w = 0.32

# Balken
b_q3  = ax.bar(x - w/2, q3s,  w, color=cols, alpha=1.0, zorder=3,
               edgecolor='white', linewidth=0.5)
b_sov = ax.bar(x + w/2, sovs, w, color=cols, alpha=0.42, zorder=3,
               edgecolor=cols, linewidth=1.8)

# Werte-Labels
for bar, v in zip(b_q3, q3s):
    ax.text(bar.get_x()+bar.get_width()/2, v+0.5,
            f"{v:.1f}%", ha='center', va='bottom',
            fontsize=12, fontweight='bold', color=DARK)
for bar, v in zip(b_sov, sovs):
    ax.text(bar.get_x()+bar.get_width()/2, v+0.5,
            f"{v:.1f}%", ha='center', va='bottom',
            fontsize=11, fontweight='bold', color=DARK, alpha=0.75)

# Parameter-Boxen direkt unter den Balken
param_info = [
    "minH = 3  ·  minE = 2\n(kein Pseudocount)",
    "pseudo = 0.01\nminH = 3  ·  minE = 2",
    "pseudo = 0.1\nminH = 3  ·  minE = 2",
]
box_colors = ['#e8f4f1', '#e8f4f1', '#e8f0f7']
for xi, txt, bc, c in zip(x, param_info, box_colors, cols):
    ax.text(xi, 3.5, txt, ha='center', va='bottom',
            fontsize=8.5, color='#444', linespacing=1.5,
            bbox=dict(boxstyle='round,pad=0.5', fc=bc, ec=c, lw=1.2, alpha=0.95))



ax.set_xticks(x)
ax.set_xticklabels(methods, fontsize=14, fontweight='bold', color=DARK)
ax.set_ylabel("Score (%)", fontweight='bold', fontsize=11)
ax.set_ylim(0, 70)
ax.set_title("Beste Parameterkombination pro GOR-Methode",
             fontweight='bold', color=DARK, fontsize=14)

ax.yaxis.grid(True, alpha=0.2, linestyle='--')
ax.set_axisbelow(True)

q3p  = mpatches.Patch(color='#555', alpha=1.0,  label='Q3  (voller Balken)')
sovp = mpatches.Patch(color='#555', alpha=0.42, label='SOV  (transparenter Balken)')
ax.legend(handles=[q3p, sovp], fontsize=10, frameon=True,
          framealpha=0.95, loc='upper right', edgecolor='#ddd')

for spine in ax.spines.values():
    spine.set_linewidth(0.8)
    spine.set_color('#ccc')

plt.tight_layout()
plt.savefig("plot3_best_overview.png", dpi=160, bbox_inches='tight')
print("Saved: plot3_best_overview.png")

plt.show()
print("\nDone. 3 plots saved:")
print("  plot1_minHE.png")
print("  plot2_pseudocount.png")
print("  plot3_best_overview.png")