#!/usr/bin/env python3
# Regenerate ALL manuscript figures in BLACK & WHITE (grayscale, print-safe).
# Sign/significance encoded by line style / hatch (not color) so figures are
# unambiguous in monochrome print. 300 dpi PNG + vector PDF.
import json, os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Ellipse, FancyArrowPatch, Patch
from matplotlib.lines import Line2D

plt.rcParams.update({
    "font.family": "DejaVu Serif",   # serif, matches Times-style manuscript
    "axes.edgecolor": "black", "axes.linewidth": 0.8,
    "savefig.facecolor": "white", "figure.facecolor": "white",
})

import os
# Portable: run from repo root (python3 make_figures.py)
PROJ = os.path.dirname(os.path.abspath(__file__))
RES = json.load(open(os.path.join(PROJ, "results/r_revision_results.json")))
DATA = os.path.join(PROJ, "data", "korean_happiness_qol_2021.csv")
OUT = [os.path.join(PROJ, "figures")]
for o in OUT: os.makedirs(o, exist_ok=True)

# grayscale palette
BLACK, DK, MID, LT, VLT, WHITE = "#000000", "#333333", "#777777", "#AAAAAA", "#DDDDDD", "#FFFFFF"

def stars(p):
    return "***" if p < .001 else "**" if p < .01 else "*" if p < .05 else "n.s."

def save(fig, name):
    for o in OUT:
        fig.savefig(os.path.join(o, name + ".png"), dpi=300, bbox_inches="tight")
        fig.savefig(os.path.join(o, name + ".pdf"), bbox_inches="tight")
    plt.close(fig)

# ---- pull M3 coefficients ----
coefs = {c["variable"]: c for c in RES["full_model_coefs"]}
order = ["subjective_class","log_income","health_mental","health_physical","autonomy",
         "female","age_c","age_sq_c","married","employed","yesterday_depress_r",
         "financial_worry_r","covid_economy"]
labels = {"subjective_class":"Subjective Class","log_income":"Log Income","health_mental":"Mental Health",
          "health_physical":"Physical Health","autonomy":"Autonomy","female":"Female","age_c":"Age",
          "age_sq_c":"Age²","married":"Married","employed":"Employed","yesterday_depress_r":"Low Depression (R)",
          "financial_worry_r":"Low Financial Worry (R)","covid_economy":"COVID Econ. Change"}
loadings = dict(zip(["Cantril Ladder","Life Satisfaction","Yesterday Happy","Eudaimonia"],
                    RES["measurement"]["std_loadings"]))
R2 = RES["competing"]["M3"]["r2"]

# ============================================================
# FIG 1: research-design schematic (B&W flowchart; uniform box width, no bold)
# ============================================================
fig, ax = plt.subplots(figsize=(9.8, 6.6)); ax.set_xlim(0,10); ax.set_ylim(0,7.4); ax.axis("off")
BW, CX = 9.4, 5.0                 # uniform full-width boxes, centered
LEFT, RIGHT = CX - BW/2, CX + BW/2
def fbox(cx, y, w, h, text, fs=8.6):
    ax.add_patch(FancyBboxPatch((cx-w/2, y-h/2), w, h,
        boxstyle="round,pad=0.02,rounding_size=0.05", fc=WHITE, ec=BLACK, lw=1.1))
    ax.text(cx, y, text, ha="center", va="center", fontsize=fs)
def varrow(x, y1, y2):
    ax.add_patch(FancyArrowPatch((x, y1), (x, y2), arrowstyle="-|>",
        mutation_scale=12, lw=1.1, color=BLACK))
fbox(CX, 6.85, BW, 0.6, "2021 Korean Happiness and Quality of Life Survey (NAFI / KOSSDA; N = 2,089)")
varrow(CX, 6.55, 6.15)
fbox(CX, 5.85, BW, 0.6, "Preprocessing: income 9,999 → missing; listwise deletion → N = 2,060; reverse-code C3, C8")
varrow(CX, 5.55, 5.15)
LWD = 4.55; LX = LEFT + LWD/2; RX = RIGHT - LWD/2   # two parallel boxes, outer ends aligned with full width
fbox(LX, 4.75, LWD, 0.72, "Measurement model\nCFA: latent happiness (4 indicators)")
fbox(RX, 4.75, LWD, 0.72, "ICC by region (3 groups)\n→ single-level SEM justified")
varrow(LX, 4.39, 3.97); varrow(RX, 4.39, 3.97)
fbox(CX, 3.65, BW, 0.6, "Competing SEMs:  M1 income-only   |   M2 SSC-only   |   M3 both (+ controls)")
varrow(CX, 3.35, 2.95)
fbox(CX, 2.55, BW, 0.74, "Robustness:  VIF | Harman & 1-factor CFA | ULMC (narrow / wide) | objective-SES composite |\noverlap-excluded | hierarchical R² | FIML | survey-weighted | mediation | cluster-robust | outliers", 8.0)
varrow(CX, 2.18, 1.74)
fbox(CX, 1.30, BW, 0.72, "Key finding: perceived social standing and psychological resources are more strongly\nassociated with happiness than household income", 8.6)
save(fig, "paper_overview")
print("Fig1 overview done")

# ============================================================
# FIG 2: PATH DIAGRAM (B&W: solid=+ , dashed=- , dotted=n.s.)
# Arrowheads land on DISTINCT points along the ellipse edge so each
# arrow's direction and thickness reads clearly (no convergence pile-up).
# ============================================================
fig, ax = plt.subplots(figsize=(13.5, 10.0)); ax.set_xlim(0,13.5); ax.set_ylim(-1.1,14); ax.axis("off")
n = len(order); ys = np.linspace(13.5, 0.6, n)
bx, bw, bh = 0.3, 3.5, 0.9
ex, ey = 7.5, 7.0
ea, eb = 2.0, 2.0           # ellipse (circle) radius
ix, iw, ih = 10.9, 2.4, 1.0
def edge_x(ty, sign):       # x on ellipse boundary at height ty (sign -1 left / +1 right)
    d = max(0.0, 1 - ((ty-ey)/eb)**2)
    return ex + sign*ea*np.sqrt(d)
for y, v in zip(ys, order):
    c = coefs[v]; b = c["std"]; sig = c["significant"]; st = stars(c["p"])
    box = FancyBboxPatch((bx, y-bh/2), bw, bh, boxstyle="round,pad=0.02,rounding_size=0.07",
                         fc=WHITE, ec=BLACK, lw=1.3 if sig else 0.8); ax.add_patch(box)
    ax.text(bx+bw/2, y+0.18, labels[v], ha="center", va="center", fontsize=11.5, fontweight="bold")
    ax.text(bx+bw/2, y-0.22, f"β = {b:+.3f} {st}", ha="center", va="center", fontsize=10.5,
            color=BLACK if sig else MID)
    ty = ey + (y-ey)*0.22                       # spread landing heights along the circle edge
    ty = max(ey-eb*0.85, min(ey+eb*0.85, ty))
    tx = edge_x(ty, -1) - 0.02
    if sig:
        acol = BLACK; lw = 1.0 + abs(b)*8.5; lsty = "-" if b > 0 else (0,(6,3))
    else:
        acol = MID; lw = 1.0; lsty = (0,(1,3))
    arr = FancyArrowPatch((bx+bw+0.06, y), (tx, ty),
                          arrowstyle="-|>", mutation_scale=13, lw=lw, color=acol,
                          linestyle=lsty, connectionstyle="arc3,rad=0.02",
                          shrinkA=0, shrinkB=2)
    ax.add_patch(arr)
el = Ellipse((ex, ey), 2*ea, 2*eb, fc=VLT, ec=BLACK, lw=1.6); ax.add_patch(el)
ax.text(ex, ey+0.34, "Latent\nHappiness", ha="center", va="center", fontsize=15.5, fontweight="bold")
ax.text(ex, ey-0.62, f"R² = {R2:.3f}", ha="center", va="center", fontsize=13, style="italic")
iys = np.linspace(10.7, 3.3, 4)
for y, (name, lam) in zip(iys, loadings.items()):
    box = FancyBboxPatch((ix, y-ih/2), iw, ih, boxstyle="round,pad=0.02,rounding_size=0.07",
                         fc=WHITE, ec=BLACK, lw=1.0); ax.add_patch(box)
    ax.text(ix+iw/2, y, name, ha="center", va="center", fontsize=11.5)
    ry = ey + (y-ey)*0.34
    rx = edge_x(ry, +1) + 0.02
    arr = FancyArrowPatch((rx, ry), (ix-0.06, y),
                          arrowstyle="-|>", mutation_scale=13, lw=0.9+lam*2.2, color=BLACK,
                          connectionstyle="arc3,rad=-0.03", shrinkA=2, shrinkB=0); ax.add_patch(arr)
    mx = rx + 0.55*(ix-0.06 - rx); my = ry + 0.55*(y - ry)
    ax.text(mx, my+0.26, f"{lam:.3f}", ha="center", va="center", fontsize=10.5, fontweight="bold",
            bbox=dict(boxstyle="square,pad=0.06", fc=WHITE, ec="none"))
leg = [Line2D([0],[0],color=BLACK,lw=4,ls="-",label="Positive (p < .05)"),
       Line2D([0],[0],color=BLACK,lw=4,ls=(0,(6,3)),label="Negative (p < .05)"),
       Line2D([0],[0],color=MID,lw=2.0,ls=(0,(1,3)),label="Non-significant"),
       Line2D([0],[0],color=BLACK,lw=3.5,ls="-",label="Factor loading")]
ax.legend(handles=leg, loc="upper center", bbox_to_anchor=(0.5,-0.02), ncol=2, frameon=False,
          fontsize=16, handlelength=3.2, columnspacing=3.0, handletextpad=0.8, labelspacing=1.0)
save(fig, "sem_path_diagram")
print("Fig2 path diagram done")

# ============================================================
# FIG 6: coefficient bar plot (B&W: black=+sig, hatched=-sig, gray=n.s.)
# ============================================================
sb = sorted(order, key=lambda v: coefs[v]["std"])
vals = [coefs[v]["std"] for v in sb]
fig, ax = plt.subplots(figsize=(8.5, 7))
yy = np.arange(len(sb))
for i, v in enumerate(sb):
    c = coefs[v]; b = c["std"]
    if c["significant"] and b > 0:
        ax.barh(i, b, color=BLACK, edgecolor=BLACK)
    elif c["significant"] and b < 0:
        ax.barh(i, b, color=WHITE, edgecolor=BLACK, hatch="////")
    else:
        ax.barh(i, b, color=LT, edgecolor=MID)
ax.set_yticks(yy); ax.set_yticklabels([labels[v] for v in sb], fontsize=10)
for sp in ["top", "right"]: ax.spines[sp].set_visible(False)   # remove top & right border
ax.axvline(0, ls="--", color=MID, lw=0.9)
for i, v in enumerate(sb):
    c = coefs[v]; st = stars(c["p"]); b = c["std"]
    ax.text(b + (0.006 if b>=0 else -0.006), i, st, va="center",
            ha="left" if b>=0 else "right", fontsize=8.5,
            fontweight="bold" if c["significant"] else "normal", color=BLACK)
ax.set_xlabel("Standardized coefficient (β)", fontsize=11)
ax.set_xlim(min(vals)-0.06, max(vals)+0.07)
leg = [Patch(fc=BLACK, ec=BLACK, label="Positive (p < .05)"),
       Patch(fc=WHITE, ec=BLACK, hatch="////", label="Negative (p < .05)"),
       Patch(fc=LT, ec=MID, label="Non-significant")]
ax.legend(handles=leg, loc="lower right", frameon=False, fontsize=9)
save(fig, "sem_coefficients")
print("Fig6 coefficient bars done")

# ============================================================
# FIG 5: correlation heatmap (B&W: shade = |r|, printed value carries sign)
# ============================================================
cm = np.array(RES["correlations"])
clab = ["Cantril","Satisfaction","Yest. Happy","Eudaimonia","Subj. Class","Log Income",
        "Mental Health","Phys. Health","Autonomy","Low Depr.(R)","Low Fin.Worry(R)","COVID Econ",
        "Female","Age","Age²","Married","Employed"]
mask = np.triu(np.ones_like(cm, dtype=bool), k=1)
absm = np.ma.array(np.abs(cm), mask=mask)
fig, ax = plt.subplots(figsize=(10, 9))
im = ax.imshow(absm, cmap="Greys", vmin=0, vmax=1)
for sp in ["top", "right"]: ax.spines[sp].set_visible(False)   # remove top & right border
ax.tick_params(top=False, right=False)
ax.set_xticks(range(len(clab))); ax.set_yticks(range(len(clab)))
ax.set_xticklabels(clab, rotation=45, ha="right", fontsize=8.3)
ax.set_yticklabels(clab, fontsize=8.3)
for i in range(len(clab)):
    for j in range(len(clab)):
        if not mask[i,j]:
            val = cm[i,j]
            ax.text(j, i, f"{val:.2f}".replace("0.",".").replace("-0.","-."),
                    ha="center", va="center", fontsize=6.3,
                    color="white" if abs(val)>0.6 else "black")
cb = fig.colorbar(im, ax=ax, shrink=0.7); cb.set_label("|Pearson r|  (sign shown in cell)")
save(fig, "correlation_heatmap")
print("Fig5 heatmap done")

# ============================================================
# data-driven figures (load raw survey)
# ============================================================
df = pd.read_csv(DATA, encoding="EUC-KR")
H = {"Cantril Ladder":"q1","Life Satisfaction":"C1","Yesterday's Happiness":"C2","Eudaimonic Well-Being":"C4"}

# FIG 3: distributions of the four happiness indicators (B&W bars)
fig, axes = plt.subplots(2, 2, figsize=(9.5, 6.6))
for ax, (name, col) in zip(axes.ravel(), H.items()):
    s = pd.to_numeric(df[col], errors="coerce").dropna()
    s = s[(s>=0)&(s<=10)]
    counts = s.value_counts(normalize=True).reindex(range(0,11), fill_value=0)*100
    ax.bar(range(0,11), counts.values, color=MID, edgecolor=BLACK, lw=0.6)
    ax.set_title(name, fontsize=10, fontweight="bold")
    ax.set_xticks(range(0,11)); ax.tick_params(labelsize=7.5)
    ax.set_xlabel("Response (0–10)", fontsize=8); ax.set_ylabel("%", fontsize=8)
    ax.text(0.03,0.95,f"M = {s.mean():.2f}\nSD = {s.std():.2f}", transform=ax.transAxes,
            va="top", ha="left", fontsize=8,
            bbox=dict(boxstyle="round,pad=0.3", fc=WHITE, ec=MID, lw=0.6))
    for sp in ["top","right"]: ax.spines[sp].set_visible(False)
fig.tight_layout()
save(fig, "happiness_distributions")
print("Fig3 distributions done")

# FIG 4: happiness indicators by region (B&W boxplots)
reg_map = {1:"Capital",2:"Non-Capital",3:"Jeju"}
df["_reg"] = pd.to_numeric(df["reg3"], errors="coerce").map(reg_map)
fig, axes = plt.subplots(2, 2, figsize=(9.5, 6.6))
for ax, (name, col) in zip(axes.ravel(), H.items()):
    vals = pd.to_numeric(df[col], errors="coerce")
    groups = [vals[(df["_reg"]==r) & vals.between(0,10)].dropna().values for r in ["Capital","Non-Capital","Jeju"]]
    bp = ax.boxplot(groups, labels=["Capital","Non-Capital","Jeju"], showmeans=True, widths=0.55,
                    patch_artist=True,
                    boxprops=dict(facecolor=VLT, edgecolor=BLACK, lw=0.9),
                    medianprops=dict(color=BLACK, lw=1.4),
                    whiskerprops=dict(color=BLACK, lw=0.8),
                    capprops=dict(color=BLACK, lw=0.8),
                    flierprops=dict(marker="o", markerfacecolor="none", markeredgecolor=MID, markersize=3, alpha=0.5),
                    meanprops=dict(marker="D", markerfacecolor=BLACK, markeredgecolor=BLACK, markersize=5))
    ax.set_title(name, fontsize=10, fontweight="bold")
    ax.tick_params(labelsize=8); ax.set_ylim(-0.5,10.5); ax.set_ylabel("Score", fontsize=8)
    for sp in ["top","right"]: ax.spines[sp].set_visible(False)
fig.tight_layout()
save(fig, "region_comparison")
print("Fig4 region boxplots done")

print("ALL B&W FIGURES SAVED")
