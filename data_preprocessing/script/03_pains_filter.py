# === PAINS Filter on final JAK1 set ===
import os
import pandas as pd
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import FilterCatalog

# === CONFIG ===
IN_PATH   = r"C:\Users\AQUILJM\OneDrive - AbbVie Inc (O365)\Desktop\JAK1\data_preprocessing\jak1_selective_druglike_core.csv"
OUT_DIR   = r"C:\Users\AQUILJM\OneDrive - AbbVie Inc (O365)\Desktop\JAK1\data_preprocessing"
FIG_DIR   = r"C:\Users\AQUILJM\OneDrive - AbbVie Inc (O365)\Desktop\JAK1\reports\figures"
os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs(FIG_DIR, exist_ok=True)

# === Load data ===
df = pd.read_csv(IN_PATH)
print(f"Loaded final core dataset: {len(df)} molecules")

# === RDKit PAINS filters ===
params = FilterCatalog.FilterCatalogParams()
params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_A)
params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_B)
params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS_C)
catalog = FilterCatalog.FilterCatalog(params)

def pains_hits(smiles):
    m = Chem.MolFromSmiles(smiles)
    if not m:
        return []
    return [entry.GetDescription() for entry in catalog.GetMatches(m)]

# Apply filter
df["PAINS_hits"] = df["smiles"].astype(str).map(pains_hits)
df["PAINS_flag"] = df["PAINS_hits"].map(lambda x: len(x) > 0)

# Split sets
df_removed = df[df["PAINS_flag"]].copy()
df_kept    = df[~df["PAINS_flag"]].copy()

print(f"Before PAINS filter: {len(df)}")
print(f"Removed (PAINS+):   {len(df_removed)}")
print(f"Kept (PAINS–):      {len(df_kept)}")

# Save outputs
kept_csv    = os.path.join(OUT_DIR, "jak1_core_noPAINS.csv")
removed_csv = os.path.join(OUT_DIR, "jak1_PAINS_removed.csv")
df_kept.to_csv(kept_csv, index=False)
df_removed.to_csv(removed_csv, index=False)
print(f"Saved filtered set → {kept_csv}")
print(f"Saved removed set  → {removed_csv}")

# Summary of top PAINS alerts
from collections import Counter
all_hits = [h for sub in df_removed["PAINS_hits"] for h in sub]
hit_counts = Counter(all_hits).most_common(15)
print("\nTop PAINS alerts in removed compounds:")
for k, v in hit_counts:
    print(f"{k:<40} {v}")

# Plot pActivity distribution: kept vs removed
plt.figure(figsize=(7,4))
df_kept["pActivity"].hist(alpha=0.7, bins=30, label="Kept (PAINS–)")
df_removed["pActivity"].hist(alpha=0.7, bins=30, label="Removed (PAINS+)")
plt.xlabel("pActivity")
plt.ylabel("Count")
plt.title("pActivity Distribution: PAINS– vs PAINS+")
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, "pActivity_PAINS_filter.png"), dpi=150)
plt.show()

print("✅ PAINS filter complete.")

# Loaded final core dataset: 2609 molecules
# Before PAINS filter: 2609
# Removed (PAINS+):   37
# Kept (PAINS–):      2572
# Saved filtered set → C:\Users\AQUILJM\OneDrive - AbbVie Inc (O365)\Desktop\JAK1\data_preprocessing\jak1_core_noPAINS.csv
# Saved removed set  → C:\Users\AQUILJM\OneDrive - AbbVie Inc (O365)\Desktop\JAK1\data_preprocessing\jak1_PAINS_removed.csv

# Top PAINS alerts in removed compounds:
# anil_di_alk_A(478)                       29
# anil_di_alk_C(246)                       4
# azo_A(324)                               3
# anil_di_alk_E(186)                       1
