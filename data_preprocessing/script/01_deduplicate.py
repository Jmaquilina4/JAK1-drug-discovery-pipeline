# JAK1 curation v2 (uses pActivity + SMILES)
import pandas as pd, os

IN_PATH  = r"C:\Users\AQUILJM\OneDrive - AbbVie Inc (O365)\Desktop\JAK1\data_preprocessing\chembl_jak1_clean.csv"
OUT_PATH = r"C:\Users\AQUILJM\OneDrive - AbbVie Inc (O365)\Desktop\JAK1\data_preprocessing\chembl_jak1_curated.csv"
os.makedirs(os.path.dirname(OUT_PATH), exist_ok=True)

df = pd.read_csv(IN_PATH)

# Ensure columns exist
assert {'smiles','pActivity'}.issubset(df.columns), f"Missing cols. Got: {list(df.columns)[:20]}"

# Clean SMILES (largest fragment)
def keep_big_frag(s):
    if pd.isna(s): return s
    parts=[p for p in str(s).split('.') if p]
    return max(parts,key=len) if parts else s
df['smiles'] = df['smiles'].astype(str).map(keep_big_frag)

# Make pActivity numeric + filter sane range
df['pActivity'] = pd.to_numeric(df['pActivity'], errors='coerce')
df = df.dropna(subset=['smiles','pActivity'])
df = df[df['pActivity'].between(4,11, inclusive='both')]

# Optional filters if present
if 'standard_type' in df:
    df = df[df['standard_type'].str.upper().isin({'IC50','KI','KD','EC50'})]
if 'standard_units' in df:
    df = df[df['standard_units'].astype(str).str.lower().eq('nm')]

# Dedupe: keep highest pActivity per (mol, smiles)
keys = [c for c in ['molecule_chembl_id','smiles'] if c in df.columns]
if keys:
    df = df.sort_values(keys+['pActivity'], ascending=[True]*len(keys)+[False]) \
           .drop_duplicates(subset=keys, keep='first')

df.to_csv(OUT_PATH, index=False)
print(f"✅ Curated dataset: {len(df)} | Saved → {OUT_PATH}")

# ✅ Curated dataset: 4503 | Saved → C:\Users\AQUILJM\OneDrive - AbbVie Inc (O365)\Desktop\JAK1\data_preprocessing\chembl_jak1_curated.csv
