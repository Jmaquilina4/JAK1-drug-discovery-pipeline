# Selective, drug-like JAK1 filtering (logs + plots + saves)
import os, math, warnings
import pandas as pd, numpy as np
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, Lipinski

# === CONFIG ===
IN_PATH   = r"C:\Users\AQUILJM\OneDrive - AbbVie Inc (O365)\Desktop\JAK1\data_preprocessing\chembl_jak1_curated.csv"
# Optional: isoform activities (JAK2/JAK3/TYK2) from ChEMBL; columns needed: molecule_chembl_id, target_name, pActivity
ISO_PATH  = None  # e.g., r"C:\...\jak_all_isoforms_curated.csv" or leave None
OUT_DIR   = r"C:\Users\AQUILJM\OneDrive - AbbVie Inc (O365)\Desktop\JAK1\data_preprocessing"
FIG_DIR   = r"C:\Users\AQUILJM\OneDrive - AbbVie Inc (O365)\Desktop\JAK1\reports\figures"
os.makedirs(OUT_DIR, exist_ok=True); os.makedirs(FIG_DIR, exist_ok=True)

# Rules (tune as needed)
PACT_CORE_MIN = 7.0     # core keep
PACT_EXP_MIN  = 6.0     # expansion tier lower bound
MW_RANGE      = (250, 500)  # (strict alt: (300, 450))
LOGP_RANGE    = (0, 5)
HBD_MAX, HBA_MAX = 5, 10
ROT_MAX       = 10
TPSA_MAX      = 140
RO5_MAX_VIOL  = 1
SELECTIVITY_DELTA = 1.0  # pAct(JAK1) - max(other isoforms) >= 1.0

warnings.filterwarnings("ignore")

# === Helpers ===
def keep_big_frag(s):
    if pd.isna(s): return s
    parts = [p for p in str(s).split('.') if p]
    return max(parts, key=len) if parts else s

def mol_from_smiles(s):
    try: return Chem.MolFromSmiles(s)
    except: return None

def compute_descriptors(m):
    # Returns dict with MW, LogP, TPSA, HBD, HBA, RotB, Rings
    return dict(
        MW   = Descriptors.MolWt(m),
        LogP = Crippen.MolLogP(m),
        TPSA = Descriptors.TPSA(m),
        HBD  = Lipinski.NumHDonors(m),
        HBA  = Lipinski.NumHAcceptors(m),
        RotB = Lipinski.NumRotatableBonds(m),
        Rings= Lipinski.RingCount(m),
    )

def count_ro5_viol(d):
    v = 0
    if d["MW"]   > 500: v += 1
    if d["LogP"] > 5:   v += 1
    if d["HBD"]  > 5:   v += 1
    if d["HBA"]  > 10:  v += 1
    return v

def within_ranges(d):
    mw_ok   = MW_RANGE[0] <= d["MW"]   <= MW_RANGE[1]
    logp_ok = LOGP_RANGE[0] <= d["LogP"] <= LOGP_RANGE[1]
    hbd_ok  = d["HBD"] <= HBD_MAX
    hba_ok  = d["HBA"] <= HBA_MAX
    rot_ok  = d["RotB"] <= ROT_MAX
    tpsa_ok = d["TPSA"] <= TPSA_MAX
    ro5_ok  = count_ro5_viol(d) <= RO5_MAX_VIOL
    return mw_ok and logp_ok and hbd_ok and hba_ok and rot_ok and tpsa_ok and ro5_ok

def log_step(name, df):
    print(f"{name:<35} n={len(df):5d}")

# === Load base (JAK1) ===
df = pd.read_csv(IN_PATH)
need = {"smiles","pActivity"}
assert need.issubset(df.columns), f"Missing columns: {need - set(df.columns)}"
if "molecule_chembl_id" not in df.columns:
    df["molecule_chembl_id"] = pd.RangeIndex(len(df)).astype(str)

df["smiles"] = df["smiles"].astype(str).map(keep_big_frag)
df["pActivity"] = pd.to_numeric(df["pActivity"], errors="coerce")
df = df.dropna(subset=["smiles","pActivity"]).copy()
log_step("Loaded curated JAK1", df)

# === Tag biochemical assays preferred (kept earlier; ensure again)
if "assay_type" in df.columns:
    before = len(df); df = df[df["assay_type"].astype(str).str.upper().eq("B")]
    log_step("Assay type = B (binding)", df)

# === Optional: merge in isoform activities for selectivity
if ISO_PATH and os.path.exists(ISO_PATH):
    iso = pd.read_csv(ISO_PATH)
    # Expect: molecule_chembl_id, target_name, pActivity
    for c in ("molecule_chembl_id","target_name","pActivity"):
        assert c in iso.columns, f"{c} missing in isoform file"
    iso["target_name_norm"] = iso["target_name"].astype(str).str.lower()
    # map to buckets
    def iso_bucket(t):
        t = t.lower()
        if "jak2" in t: return "JAK2"
        if "jak3" in t: return "JAK3"
        if "tyk2" in t: return "TYK2"
        return None
    iso["iso"] = iso["target_name_norm"].map(iso_bucket)
    iso = iso[iso["iso"].notna()].copy()
    iso_pvt = (iso.pivot_table(index="molecule_chembl_id", columns="iso", values="pActivity", aggfunc="max")
                   .reset_index())
    df = df.merge(iso_pvt, on="molecule_chembl_id", how="left")
else:
    # No isoform file → leave cols absent; mark later as JAK1_only
    pass

# === Compute RDKit descriptors
mols, descs = [], []
for s in df["smiles"]:
    m = mol_from_smiles(s); mols.append(m)
    if m is None:
        descs.append(dict(MW=np.nan,LogP=np.nan,TPSA=np.nan,HBD=np.nan,HBA=np.nan,RotB=np.nan,Rings=np.nan))
    else:
        descs.append(compute_descriptors(m))
df = pd.concat([df, pd.DataFrame(descs)], axis=1)
df = df.dropna(subset=["MW","LogP","TPSA"])  # drop RDKit parse failures
log_step("RDKit-parsable", df)

# === Selectivity margin
other_cols = [c for c in ["JAK2","JAK3","TYK2"] if c in df.columns]
if other_cols:
    df["max_other"] = df[other_cols].max(axis=1)
    df["sel_margin"] = df["pActivity"] - df["max_other"]
    df["selective"] = df["sel_margin"] >= SELECTIVITY_DELTA
    df["selectivity_status"] = np.where(df["selective"], "Selective(≥1 log)", "Not_selective")
else:
    df["sel_margin"] = np.nan
    df["selective"] = False
    df["selectivity_status"] = "JAK1_only"

# === Filter by pActivity tiers
df_core = df[df["pActivity"] >= PACT_CORE_MIN].copy()
log_step("pActivity ≥ 7.0 (core)", df_core)

df_exp  = df[(df["pActivity"] >= PACT_EXP_MIN)].copy()
log_step("pActivity ≥ 6.0 (expansion)", df_exp)

# === Drug-likeness gates (Lipinski/Veber-ish + ≤1 Ro5 violation)
def passes_drug_like(row):
    d=dict(MW=row.MW, LogP=row.LogP, TPSA=row.TPSA, HBD=row.HBD, HBA=row.HBA, RotB=row.RotB)
    return within_ranges(d)

before = len(df_core)
df_core = df_core[df_core.apply(passes_drug_like, axis=1)]
log_step("Drug-likeness gates (core)", df_core)

before_exp = len(df_exp)
df_exp = df_exp[df_exp.apply(passes_drug_like, axis=1)]
log_step("Drug-likeness gates (expansion)", df_exp)

# === Remove duplicates: keep highest pAct per (mol, canonical SMILES string from RDKit)
def canon_smiles(s):
    m = mol_from_smiles(s)
    return Chem.MolToSmiles(m, isomericSmiles=True, canonical=True) if m else np.nan
df_core["canonical_smiles"] = df_core["smiles"].map(canon_smiles)
df_exp["canonical_smiles"]  = df_exp["smiles"].map(canon_smiles)

def dedupe_best(d):
    order = ["molecule_chembl_id","canonical_smiles","pActivity"]
    d = d.dropna(subset=["canonical_smiles"])
    d = d.sort_values(order, ascending=[True, True, False])
    return d.drop_duplicates(subset=["molecule_chembl_id","canonical_smiles"], keep="first")

df_core = dedupe_best(df_core)
df_exp  = dedupe_best(df_exp)
log_step("Deduped (core)", df_core)
log_step("Deduped (expansion)", df_exp)

# === Save outputs
core_csv = os.path.join(OUT_DIR, "jak1_selective_druglike_core.csv")
exp_csv  = os.path.join(OUT_DIR, "jak1_selective_druglike_expansion.csv")
df_core.to_csv(core_csv, index=False)
df_exp.to_csv(exp_csv, index=False)
print(f"Saved:\n  CORE → {core_csv}\n  EXP  → {exp_csv}")

# === Flag top-selective subset (if iso data present)
if other_cols:
    top_sel = df_core[df_core["selective"]].sort_values("sel_margin", ascending=False)
    top_sel_path = os.path.join(OUT_DIR, "jak1_top_selective_core.csv")
    top_sel.to_csv(top_sel_path, index=False)
    print(f"Saved top-selective CORE → {top_sel_path}")

# === Plots: MW, LogP, HBD, HBA, TPSA, pActivity (hist + box)
props = ["MW","LogP","HBD","HBA","TPSA","pActivity"]
for p in props:
    plt.figure(figsize=(7,4))
    df_core[p].dropna().hist(bins=40)
    plt.xlabel(p); plt.ylabel("Count"); plt.title(f"CORE {p} distribution")
    plt.tight_layout()
    plt.savefig(os.path.join(FIG_DIR, f"core_hist_{p}.png"), dpi=150)
    plt.show()

for p in props:
    plt.figure(figsize=(6,4))
    plt.boxplot(df_core[p].dropna().values)
    plt.ylabel(p); plt.title(f"CORE {p} boxplot")
    plt.tight_layout()
    plt.savefig(os.path.join(FIG_DIR, f"core_box_{p}.png"), dpi=150)
    plt.show()

print("✅ Done.")

# Loaded curated JAK1                 n= 4503
# Assay type = B (binding)            n= 4503
# RDKit-parsable                      n= 4503
# pActivity ≥ 7.0 (core)              n= 3546
# pActivity ≥ 6.0 (expansion)         n= 4172
# Drug-likeness gates (core)          n= 2609
# Drug-likeness gates (expansion)     n= 3134
# Deduped (core)                      n= 2609
# Deduped (expansion)                 n= 3134
# Saved:
  # CORE → C:\Users\AQUILJM\OneDrive - AbbVie Inc (O365)\Desktop\JAK1\data_preprocessing\jak1_selective_druglike_core.csv
  # EXP  → C:\Users\AQUILJM\OneDrive - AbbVie Inc (O365)\Desktop\JAK1\data_preprocessing\jak1_selective_druglike_expansion.csv
