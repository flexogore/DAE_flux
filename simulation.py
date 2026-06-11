import os
import random
import timeit

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from cobra.io import read_sbml_model
from freeflux import Model
from matplotlib.backends.backend_pdf import PdfPages

kokel = read_sbml_model("photo_model_last.xml")

concs = pd.read_excel("essentmets_new.xlsx")
concs = concs[["Abbreviation", "Concentration"]]

conc = concs.set_index("Abbreviation").to_dict()["Concentration"]


def convert_key(key):
    base, comp = key.split("[")
    comp = comp.rstrip("]").upper()
    base = base.replace("_", "")
    if base[0].isdigit():
        base = "n" + base
    return base + comp


conc = {convert_key(k): v for k, v in conc.items()}

# Activated to simulate randomized concentrations
concs_buf = dict()

for i in range(1):
    concs_buf = conc.copy()  # ← creates a new dict each iteration
    for key, value in concs_buf.items():
        coef = random.uniform(0.01, 100)
        concs_buf[key] = value * coef


fluxes = pd.read_csv("fluxes.csv").set_index("Var1")

fluxes.index.name = "Reaction"
fluxes = fluxes.mul(1000 / 3600)

emus = {
    # CBB cycle 
    "n2pglycH": "12",
    "glycltH": "12",  
    "rb15bpH": "12345",
    "ru5pDH": "12345",
    "r5pH": "12345",
    "xu5pDH": "12345",
    "e4pH": "1234",
    "s7pH": "1234567",
    "s17bpH": "1234567",  
    "dhapH": "123",
    "fdpBH": "123456",
    "f6pBH": "123456",
    "g6pAH": "123456",
    "g1pH": "123456",
    "adpglcH": "123456",
    "n3pgH": "123",
    "pepH": "123",  
    "pyrH": "123",  
    "glxH": "12", 
    "serLH": "123", 
    # Cytosol
    "n3pgC": "123",  
    "dhapC": "123",  
    "f6pBC": "123456",  
    "fdpBC": "123456", 
    "g1pC": "123456",  
    "g6pAC": "123456",  
    "pepC": "123",
    "pyrC": "123",
    "akgC": "12345",  
    "alaLC": "123",
    "aspLC": "1234",  
    "gluLC": "12345",
    "serLC": "123",
    "glxC": "12",
    "malLC": "1234",
    "oaaC": "1234",
    "udpgC": "123456",  
    # TCA 
    "pepM": "123",
    "pyrM": "123",
    "alaLM": "123",  
    "aspLM": "1234", 
    "gluLM": "12345", 
    "akgM": "12345",
    "citM": "123456",
    "icitM": "123456",
    "succM": "1234",
    "fumM": "1234",
    "malLM": "1234",
    "oaaM": "1234",
    "glxM": "12",
}

MODEL_FILE = "modreactions.tsv"

model_prod = Model("prod")
model_prod.read_from_file(MODEL_FILE)

isimucho = model_prod.simulator("inst")

isimucho.set_target_EMUs(emus)

isimucho.set_labeling_strategy(
    "co2E.ex", labeling_pattern=["1"], percentage=[0.5], purity=[0.95]
)

timepoints = [0, 5, 10, 20, 40]
isimucho.set_timepoints(timepoints)

counter = 0

for concid, value in concs.items():
    v = value
    isimucho.set_concentration(concid, v)

for fluxd in fluxes:
    
    # for a randomized MID subset
    """
    if counter > 800:
        break
    """
    
    flux_wb = fluxes[fluxd].to_dict()

    for fluxid, value in flux_wb.items():
        isimucho.set_flux(fluxid, value)

    isimucho.prepare(n_jobs=11)
    resprod = isimucho.simulate()

    # print(f'Simulation {counter} is done!')
    counter += 1

    frames = []
    for emu in resprod.simulated_EMUs:
        name = emu.split("_")[0]
        mdvs = [v.value for v in resprod.simulated_MDV(emu).values()]
        if len(mdvs) == 0:
            continue
        n_states = len(mdvs[0])
        cols = [f"{name}_{x}" for x in range(n_states)]
        frames.append(pd.DataFrame(mdvs, columns=cols))

    if frames:
        df = pd.concat(frames, axis=1)
        df.to_csv(os.path.join("MID_simulated/", f"{fluxd}.csv"), index=False)

directory = r'MID_simulated'

MID_8000 = pd.DataFrame()
rand_list = []

for filename in os.listdir(directory):
    filepath = os.path.join(directory, filename)
    file = pd.read_csv(filepath)
    rand_list.append(file)

MID_8000 = pd.concat(rand_list, axis=0)

MID_8000.to_csv('MID_8000.csv')
