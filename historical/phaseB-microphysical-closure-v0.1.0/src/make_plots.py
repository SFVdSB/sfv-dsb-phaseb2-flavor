#!/usr/bin/env python3
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
OUT = ROOT / "plots"
OUT.mkdir(exist_ok=True)

c = pd.read_csv(RES / "baseline_microphysical_control_map.csv")
fig, ax = plt.subplots(figsize=(8, 4.5))
x = np.arange(len(c))
w = 0.36
ax.bar(x-w/2, c.canonical_fitted, width=w, label="Seven-control fit")
ax.bar(x+w/2, c.phaseA_formula, width=w, label="Phase-A formula")
ax.set_xticks(x, c.control, rotation=35, ha="right")
ax.set_ylabel("Control value")
ax.set_title("Benchmark controls: fitted versus Phase-A formula")
ax.legend()
fig.tight_layout()
fig.savefig(OUT / "control_map_comparison.png", dpi=180)
plt.close(fig)

r = pd.read_csv(RES / "refit_controls_all_walls.csv")
fig, ax = plt.subplots(figsize=(5.4, 5.0))
ax.scatter(r.pred_a_d1, r.fit_a_d1, s=26)
lo = min(r.pred_a_d1.min(), r.fit_a_d1.min())
hi = max(r.pred_a_d1.max(), r.fit_a_d1.max())
ax.plot([lo, hi], [lo, hi], linestyle="--")
ax.set_xlabel(r"Phase-A locking interval $\alpha(R_{\rm mix}-R_{\rm grad})$")
ax.set_ylabel(r"Refitted $a_{d1}$")
ax.set_title("51-wall audit of the locking relation")
fig.tight_layout()
fig.savefig(OUT / "locking_relation_51_walls.png", dpi=180)
plt.close(fig)

m = pd.read_csv(RES / "constrained_model_results.csv")
fig, ax = plt.subplots(figsize=(8.5, 4.8))
ax.bar(np.arange(len(m)), m.max_abs_percent_error)
ax.axhline(1.0, linestyle="--", label="1% threshold")
ax.set_xticks(np.arange(len(m)), m.model, rotation=35, ha="right")
ax.set_ylabel("Maximum absolute error (%)")
ax.set_title("Phase-B constrained closure tests")
ax.set_yscale("log")
ax.legend()
fig.tight_layout()
fig.savefig(OUT / "constrained_model_errors.png", dpi=180)
plt.close(fig)

f = pd.read_csv(RES / "formula_closure_all_walls.csv")
fig, ax = plt.subplots(figsize=(7.2, 4.7))
sc = ax.scatter(100*f.y, f.max_abs_percent_error, c=100*f.x, s=35)
ax.set_xlabel("Stiffness-coordinate shift y (%)")
ax.set_ylabel("Maximum target error (%)")
ax.set_title("Counterfactual wall variation of the frozen formula")
cb = fig.colorbar(sc, ax=ax)
cb.set_label("Driving-coordinate shift x (%)")
fig.tight_layout()
fig.savefig(OUT / "formula_counterfactual_wall_variation.png", dpi=180)
plt.close(fig)
