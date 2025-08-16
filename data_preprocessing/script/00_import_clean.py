# JAK1 ChEMBL cleaner (semicolon CSV → pActivity → dedupe)
import pandas as pd, numpy as np, math, os, csv

# === CONFIG ===
FILE_PATH = r'C:\Users\AQUILJM\OneDrive - AbbVie Inc (O365)\Desktop\JAK1\JAK1_Chembl_data.csv'
OUT_DIR   = r'C:\Users\AQUILJM\OneDrive - AbbVie Inc (O365)\Desktop\JAK1\data_preprocessing'
OUT_CSV   = os.path.join(OUT_DIR, 'chembl_jak1_clean.csv')
os.makedirs(OUT_DIR, exist_ok=True)

# === LOAD (semicolon + tolerate bad quotes) ===
df = pd.read_csv(
    FILE_PATH, sep=';', engine='python',
    quoting=csv.QUOTE_NONE, on_bad_lines='skip',
    encoding='utf-8', encoding_errors='ignore'
)

# === NORMALIZE HEADERS ===
df.columns = (df.columns.astype(str)
              .str.strip().str.replace('"','', regex=False)
              .str.replace(r'\s+', '_', regex=True)
              .str.lower())

# === STRIP EMBEDDED QUOTES IN CELLS (no applymap warning) ===
def strip_col(s):
    return s.str.replace('"','', regex=False).str.strip() if s.dtype=='object' else s
df = df.apply(strip_col)

# === CAST NUMERIC FIELDS (if present) ===
for c in ['standard_value','standard_text_value','pchembl_value']:
    if c in df: df[c] = pd.to_numeric(df[c], errors='coerce')

# === FILTERS ===
allowed = {'IC50','KI','KD','EC50'}
df = df[df['standard_type'].str.upper().isin(allowed)]
df = df[df['standard_relation'].isin(['=','~',"'='"]) | df['standard_relation'].isna()]
df = df[df['standard_units'].str.lower().eq('nm')]
df = df[df['assay_type'].astype(str).str.upper().eq('B')]

# === pActivity (prefer pChEMBL / text; else convert nM) ===
def to_pact(r):
    for c in ('pchembl_value','standard_text_value'):
        v = r.get(c, np.nan)
        if pd.notna(v) and 0 < float(v) <= 12:
            return float(v)
    v = r.get('standard_value', np.nan)
    return 9.0 - math.log10(float(v)) if pd.notna(v) and v > 0 else np.nan

df['pActivity'] = df.apply(to_pact, axis=1)
df = df[df['pActivity'].between(0, 12, inclusive='both')]

# === DE-SALT SMILES (largest fragment) ===
def keep_big_frag(s):
    if pd.isna(s): return s
    parts = [p for p in str(s).split('.') if p]
    return max(parts, key=len) if parts else s
df['smiles'] = df['smiles'].astype(str).map(keep_big_frag)

# === DEDUPE: median pActivity per (molecule, smiles) ===
subset = ['molecule_chembl_id', 'smiles']
agg = {c: 'first' for c in df.columns if c not in subset + ['pActivity','standard_value','standard_text_value','pchembl_value']}
clean = (df.groupby(subset, dropna=False)
           .agg({**agg, 'pActivity':'median'})
           .reset_index())

# === SAVE ===
clean.to_csv(OUT_CSV, index=False)
print(f'✅ Cleaned molecules: {len(clean)} | Saved: {OUT_CSV}')

# ✅ Cleaned molecules: 4508 | Saved: C:\Users\AQUILJM\OneDrive - AbbVie Inc (O365)\Desktop\JAK1\data_preprocessing\chembl_jak1_clean.csv
