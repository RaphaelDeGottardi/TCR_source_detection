import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import scirpy as ir
import argparse
import os

AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")

DB_FILE = "/home/dego/.cache/iggytop_airr/mcpas_latest/McPAS-TCR.csv"

CDR3a = "CDR3.alpha.aa"
CDR3b = "CDR3.beta.aa"
PEPTIDE = "Epitope.peptide"

TARGET = None
TARGET_POS = "1"

"""
DB_FILE = "/home/dego/.cache/iggytop_airr/cedar_latest/receptor_full_v3.zip.unzip/tcr_full_v3.csv"

CDR3a = "Chain 1 CDR3 Calculated"
CDR3b = "Chain 2 CDR3 Calculated"
PEPTIDE = "Epitope Name"

TARGET = None
TARGET_POS = "1"


DB_FILE = "/home/dego/.cache/iggytop_airr/tcr3d_latest/tcr_complexes_data.tsv"

CDR3a = "CDR3_alpha"
CDR3b = "CDR3_beta"
PEPTIDE = "Epitope"

TARGET = None
TARGET_POS = "1"

DB_FILE = "/home/dego/iggytop/t2pmhc_train.tsv"

CDR3a = "cdr3a"
CDR3b = "cdr3b"
PEPTIDE = "peptide"

TARGET = "binder"
TARGET_POS = 1


DB_FILE = "/home/dego/.cache/iggytop_airr/trait_latest/Interactive_TCR-pMHC_Pairs.zip_20250312.zip.unzip/20250312-TRAIT_search_download.xlsx"

CDR3a = "CDR3α"
CDR3b =  "CDR3β"
PEPTIDE = "Epitope"

TARGET = None
TARGET_POS = 1
"""
        
def is_valid_peptide_sequence(seq: str) -> bool:
    """Checks if a given sequence is a valid peptide sequence."""
    if isinstance(seq, str) and len(seq) > 2:
        return all([aa in AMINO_ACIDS for aa in seq])
    else:
        return False

def process_cdr3_sequence(seq: str, is_igh: bool = False) -> str | None:
    if seq is None or pd.isna(seq):
        return None

    # Clean and normalize the sequence
    seq = str(seq).upper().strip().replace(" ", "").replace("\n", "")

    # Validate that the sequence contains only valid amino acids (optional: define valid AAs if needed)
    if not is_valid_peptide_sequence(seq):
        return None

    # Check if sequence has a valid CDR3 format
    starts_with_c = seq.startswith("C")
    ends_with_fw = seq.endswith("F") or (is_igh and seq.endswith("W"))

    if starts_with_c and ends_with_fw:
        return seq

    # Pad the sequence appropriately
    seq = seq.lstrip("C")  # remove leading C if already present
    if is_igh:
        seq = seq.rstrip("FW")  # remove existing F or W if present
        return f"C{seq}W"
    else:
        seq = seq.rstrip("F")
        return f"C{seq}F"

def run_coverage_report(dataset_path, alpha_col, beta_col, epitope_col, cache_dir="/home/dego/.cache/iggytop_airr"):
    merged_path = os.path.join(cache_dir, "merged_anndata.h5ad")
    
    print(f"Loading Iggytop dataset from {merged_path}...")
    try:
        adata_merged = sc.read_h5ad(merged_path)
    except Exception as e:
        print(f"Error loading Iggytop dataset: {e}")
        return

    print(f"Loading target dataset from {dataset_path}...")
    if dataset_path.endswith(".tsv") or dataset_path.endswith(".txt"):
        target_df = pd.read_csv(dataset_path, sep="\t")
    elif dataset_path.endswith(".csv"):
        if CDR3a == "Chain 1 CDR3 Calculated":
            target_df = pd.read_csv(dataset_path, header=[0, 1], dtype=str)
            # Combine multi-index header into single strings
            target_df.columns = [
                ' '.join([str(c).strip() for c in col if str(c).strip() and not str(c).lower().startswith('unnamed')])
                for col in target_df.columns.values
            ]
        else:
            target_df = pd.read_csv(dataset_path, dtype=str)
    elif dataset_path.endswith(".xlsx") or dataset_path.endswith(".xls"):
        target_df = pd.read_excel(dataset_path)
    else:
        print(f"Unsupported file format for target dataset: {dataset_path}")
        return
    if TARGET:
        target_df = target_df[target_df[TARGET] == TARGET_POS].copy()
    
    # Handle the specific case where epitope_col is not in target_df
    # but might be derivable from peptide_HLA
    if epitope_col not in target_df.columns:
        if "peptide_HLA" in target_df.columns:
            print("Extracting peptide from 'peptide_HLA'...")
            target_df[["peptide", "HLA"]] = target_df["peptide_HLA"].str.split(r"\s+", n=1, expand=True)
            epitope_col = "peptide"
        else:
            print(f"Error: Column '{epitope_col}' not found in dataset and no fallback available.")
            return

    # Preprocessing
    print("Preprocessing sequences...")
    try:
        match_df = target_df[[alpha_col, beta_col, epitope_col, "PMID"]].copy()        
        match_df.columns = ["CDR3a", "CDR3b", "peptide", "PMID"]
    except KeyError as e:
        print(f"Error: PMID not found in target dataset: {e}")
        match_df = target_df[[alpha_col, beta_col, epitope_col]].copy()        
        match_df.columns = ["CDR3a", "CDR3b", "peptide"]

    
    for col in ["CDR3a", "CDR3b", "peptide"]:
        match_df[col] = match_df[col].astype("string").str.strip().str.upper()
    
    for cdr_col in ["CDR3a", "CDR3b"]:
        match_df[cdr_col] = match_df[cdr_col].apply(lambda x: process_cdr3_sequence(x, is_igh=False))

    # Drop rows where we can't build a key
    initial_len = len(match_df)
    match_df = match_df.dropna(subset=["peptide"])
    match_df = match_df.reset_index(drop=True)
    if len(match_df) < initial_len:
        print(f"Dropped {initial_len - len(match_df)} rows with missing epitope sequence.")

    # Iggytop keys
    print("Extracting keys from Iggytop...")
    iggy_match_df = pd.DataFrame(index=adata_merged.obs_names)
    iggy_match_df["CDR3a"] = ir.get.airr(adata_merged, "junction_aa", chain="VJ_1")
    iggy_match_df["CDR3b"] = ir.get.airr(adata_merged, "junction_aa", chain="VDJ_1")
    iggy_match_df["peptide"] = adata_merged.obs["epitope_sequence"].astype("string")

    def make_key_set(df, cols):
        return set(
            df[list(cols)]
            .dropna(subset=list(cols))
            .drop_duplicates()
            .itertuples(index=False, name=None)
        )

    full_cols = ["CDR3a", "CDR3b", "peptide"]
    iggy_full_set = make_key_set(iggy_match_df, full_cols)
    iggy_a_pep = make_key_set(iggy_match_df, ["CDR3a", "peptide"])
    iggy_b_pep = make_key_set(iggy_match_df, ["CDR3b", "peptide"])

    # Match target dataset
    print("Matching records...")
    
    # We will categorize each record into exactly one category
    # Order of priority:
    # 1. Missing Epitope
    # 2. Multi-epitope (peptide contains "+")
    # 3. No CDR3 info (both null)
    # 4. Exact Match (Full)
    # 5. Partial Match (Alpha or Beta)
    # 6. Non-matched but from 10x
    # 7. Rest (No match)

    results = []
    for i, row in match_df.iterrows():
        a = row["CDR3a"]
        b = row["CDR3b"]
        pep = row["peptide"]
        pmid = str(row["PMID"]) if "PMID" in row and not pd.isna(row["PMID"]) else ""
        
        category = "Rest"
        
        if pd.isna(pep) or pep == "":
            category = "No Epitope"
        elif "+" in str(pep):
            category = "Multi-epitope (+)"
        elif (pd.isna(a) or a == "") and (pd.isna(b) or b == ""):
            category = "No CDR3 info"
        else:
            # Check for matches
            full_match = (a, b, pep) in iggy_full_set
            if full_match:
                category = "Exact Match"
            else:
                a_pep_match = (a, pep) in iggy_a_pep if not pd.isna(a) else False
                b_pep_match = (b, pep) in iggy_b_pep if not pd.isna(b) else False
                
                if a_pep_match or b_pep_match:
                    category = "Partial Match"
                elif pmid.startswith("https://www.10x"):
                    category = "Unmatched (10x)"
                else:
                    category = "Rest"
            
        results.append(category)
    
    match_df["category"] = results
    total_records = len(match_df)
    
    counts = match_df["category"].value_counts()
    
    def print_stat(label, category_name):
        count = counts.get(category_name, 0)
        pct = (count / total_records * 100) if total_records > 0 else 0
        print(f"{label}: {count:,} ({pct:.2f}%)")

    print("\n=== Iggytop Coverage Report ===")
    print(f"Target Dataset: {dataset_path}")
    print(f"Total Records: {total_records:,}")
    print("-" * 30)
    
    print_stat("Exact matches (Full)", "Exact Match")
    print_stat("Partial matches (Alpha or Beta)", "Partial Match")
    print_stat("No CDR3 info (both null)", "No CDR3 info")
    print_stat("No epitope", "No Epitope")
    print_stat("Multi-epitope (+)", "Multi-epitope (+)")
    print_stat("Unmatched (from 10x)", "Unmatched (10x)")
    print_stat("Rest (No match found)", "Rest")
    
    print("=" * 30)
    
    rest_df = match_df[match_df["category"] == "Rest"]
    if not rest_df.empty:
        print("\nSample of 'Rest' (unmatched and not 10x) records:")
        print(rest_df.head(10))
    
    tenx_unmatched = match_df[match_df["category"] == "Unmatched (10x)"]
    if not tenx_unmatched.empty:
        print("\nSample of Unmatched (10x) records:")
        print(tenx_unmatched.head(10))



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Report Iggytop coverage for a dataset.")
    parser.add_argument("--dataset", type=str, default= DB_FILE, help="Path to the CSV dataset")
    parser.add_argument("--alpha_col", type=str, default= CDR3a, help="Column name for alpha CDR3")
    parser.add_argument("--beta_col", type=str, default= CDR3b, help="Column name for beta CDR3")
    parser.add_argument("--epitope_col", type=str, default= PEPTIDE, help="Column name for epitope/peptide")
    parser.add_argument("--cache_dir", type=str, default="/home/dego/.cache/iggytop_airr", help="Cache directory for Iggytop")
    parser.add_argument("--target_col", type=str, default= TARGET, help="Column name for target column")
    args = parser.parse_args()
    run_coverage_report(args.dataset, args.alpha_col, args.beta_col, args.epitope_col, args.cache_dir)
