import numpy as np
import pandas as pd

FRAGMENT_MAP = {
    "adpglcH_123456": "ADPG",
    "dhapH_123": "DHAP",
    "dhapC_123": "DHAP",
    "fdpBH_123456": "FBP",
    "fdpBC_123456": "FBP",
    "f6pBH_123456": "f6p",
    "g1pH_123456": "g1p",
    "g6pAH_123456": "G6P",
    "rb15bpH_12345": "RuBP",
    "ru5pDH_12345": "R5P",
    "s7pH_1234567": "S7P",
    "s17bpH_1234567": "SBP",
    "malLC_1234": "malate",
    "malLM_1234": "malate",
    "gluLC_12345": "glut",
    "gluLM_12345": "glut",
    "aspLC_1234": "asp",
    "aspLM_1234": "asp",
    "udpgC_123456": "UDPG",
    "n2pglycH_12": "2PG",
    "akgM_12345": "2-OG",
    "akgC_12345": "2-OG",
    "r5pH_12345": "PP",
    "succM_1234": "succ",
    "n3pgC_123": "3PGA",
    "n3pgH_123": "3PGA",
    "pepC_123": "PEP",
    "pepM_123": "PEP",
    "pepH_123": "PEP",
}

META = ["Light", "Alga", "Labeling time (sec)"]
TIME = "Labeling time (sec)"

# Load and filter
df = pd.read_csv("mixerIsotopomers.csv", sep=";")
df = df[df["Alga"] == "Chlamy"].reset_index(drop=True)
df = df.loc[:, ~df.columns.str.startswith("Unnamed")]

data_cols = [c for c in df.columns if c not in META]


# Get sorted isotopomer columns for a given prefix
def iso_cols(prefix):
    cols = [c for c in data_cols if c.startswith(prefix + "-")]
    cols.sort(key=lambda x: int(x.split("-")[-1]))
    return cols


prefixes = set(c.rsplit("-", 1)[0] for c in data_cols if "-" in c)
for prefix in prefixes:
    cols = iso_cols(prefix)
    block = df[cols].values.astype(float)
    for i in range(len(block)):
        row_sum = block[i].sum()
        if row_sum > 0:
            block[i] = block[i] / row_sum
        else:
            block[i] = np.nan
    df[cols] = block

time_points = sorted(df[TIME].unique())
mean_rows, std_rows = [], []
for t in time_points:
    subset = df[df[TIME] == t][data_cols]
    mean_rows.append(subset.mean())
    std_rows.append(subset.std(ddof=1))

df_mean = pd.DataFrame(mean_rows, index=time_points)
df_std = pd.DataFrame(std_rows, index=time_points)

arr = df_std.values.copy()
for j in range(arr.shape[1]):
    for i in range(arr.shape[0]):
        if arr[i, j] == 0 or np.isnan(arr[i, j]):
            # find max non-zero in this column
            col_vals = arr[:, j]
            nz = [v for v in col_vals if v > 0 and not np.isnan(v)]
            if not nz and j > 0:
                prev = arr[:, j - 1]
                nz = [v for v in prev if v > 0 and not np.isnan(v)]
            arr[i, j] = max(nz) if nz else np.nan

std_max = arr.max(axis=0)

# Build output rows
rows = []
for frag, prefix in FRAGMENT_MAP.items():
    cols = iso_cols(prefix)
    if not cols:
        continue
    col_indices = [data_cols.index(c) for c in cols]
    for t in time_points:
        mean_vals = df_mean.loc[t, cols].values
        mean_vals = np.nan_to_num(mean_vals, nan=0.0)
        sd_vals = std_max[col_indices]
        mean_str = ",".join(f"{v:.5f}" for v in mean_vals)
        sd_str = ",".join(f"{v:.5f}" for v in sd_vals)
        rows.append(
            {"#fragment_ID": frag, "time (s)": int(t), "mean": mean_str, "sd": sd_str}
        )

out = pd.DataFrame(rows, columns=["#fragment_ID", "time (s)", "mean", "sd"])
out.to_csv("mid_original.tsv", sep="\t", index=False)
print(f"Saved {len(out)} rows to 'mid_original.tsv'")
