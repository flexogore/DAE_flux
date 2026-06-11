import os
import pickle
import random

# from cobra.io import read_sbml_model
import timeit

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from freeflux import Metabolite, Model, Reaction
from matplotlib.backends.backend_pdf import PdfPages

MODEL_FILE = "modreactions.tsv"
MEASURED_MDV_FILE = "mid_original.tsv"
MEASURED_FLUX_FILE = "measured_fluxes_final.tsv"

DILUTION_FROM = ["adpglcHu", "r5pHu", "n2pglycHu"]

modelo = Model("chlamys")
modelo.read_from_file(MODEL_FILE)

ifit = modelo.fitter("inst")

ifit.set_labeling_strategy(
    "co2E.ex", labeling_pattern=["1"], percentage=[0.50], purity=[0.95]
)

ifit.set_measured_MDVs_from_file(MEASURED_MDV_FILE)
ifit.set_measured_fluxes_from_file(MEASURED_FLUX_FILE)

# Set bounds: reversible can be ±1000, irreversible must be [0, 1000]
ifit.set_flux_bounds("all", bounds=[-1, 1])  # default
ifit.set_concentration_bounds("all", bounds=[0.001, 1])

ifit.prepare(n_jobs=11, dilution_from=DILUTION_FROM)
# ifit.prepare(n_jobs=1)
"""
res = ifit.solve(solver="ralg", ini_fluxes=None, fit_measured_fluxes=True)

print(res.optimization_successful)
res.chi2_test()

pd.Series(res.opt_net_fluxes).to_excel("estimated_net_fluxes_imputed.xlsx")
pd.Series(res.opt_total_fluxes).to_excel("estimated_total_fluxes_imputed.xlsx")
pd.Series(res.opt_concentrations).to_excel("estimated_concentrations_imputed.xlsx")

res.plot_normal_probability(show_fig=False, output_dir="Fitting_figures_imputed")
res.plot_simulated_vs_measured_MDVs(
    show_fig=False, output_dir="Fitting_figures_imputed"
)
res.plot_simulated_vs_measured_fluxes(
    show_fig=False, output_dir="Fitting_figures_imputed"
)
"""

res3 = ifit.solve_with_confidence_intervals(
    solver="ralg",
    ini_fluxes="estimated_net_fluxes_loose2.xlsx",
    fit_measured_fluxes=True,
    n_runs=90,
    n_jobs=11,
)
for which, name in [
    ("net", "netflux"),
    ("total", "totalflux"),
    ("conc", "concentration"),
]:
    try:
        cis = res3.estimate_confidence_intervals(which=which, confidence_level=0.95)
        df = pd.DataFrame(cis, index=["LB", "UB"]).T
        out = f"{name}_mc_CIs_original.xlsx"
        df.to_excel(out)
        print(f"[{which}] wrote {out} with {len(df)} rows")
    except Exception as e:
        print(f"[{which}] failed: {type(e).__name__}: {e}")
