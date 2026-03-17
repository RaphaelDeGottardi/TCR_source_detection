import os
import pandas as pd
import scanpy as sc
import anndata as ad
from typing import Optional

def load_target_dataset(dataset_path: str, alpha_col: str, beta_col: str, epitope_col: str, target_col: Optional[str] = None, target_pos: str = "1") -> pd.DataFrame:
    """Loads and basic-preprocesses the target dataset from various formats."""
    print(f"Loading target dataset from {dataset_path}...")
    if dataset_path.endswith(".tsv") or dataset_path.endswith(".txt"):
        df = pd.read_csv(dataset_path, sep="\t")
    elif dataset_path.endswith(".csv"):
        # Specific handling for multi-index headers sometimes seen in CDR3 datasets
        try:
            df = pd.read_csv(dataset_path, header=[0, 1], dtype=str)
            df.columns = [
                ' '.join([str(c).strip() for c in col if str(c).strip() and not str(c).lower().startswith('unnamed')])
                for col in df.columns.values
            ]
        except:
            df = pd.read_csv(dataset_path, dtype=str)
    elif dataset_path.endswith(".xlsx") or dataset_path.endswith(".xls"):
        df = pd.read_excel(dataset_path)
    else:
        raise ValueError(f"Unsupported file format: {dataset_path}")

    if target_col and target_col in df.columns:
        df = df[df[target_col].astype(str) == str(target_pos)].copy()
    
    # Extract peptide from peptide_HLA if needed
    if epitope_col not in df.columns and "peptide_HLA" in df.columns:
        print("Extracting peptide from 'peptide_HLA'...")
        extracted = df["peptide_HLA"].str.split(r"\s+", n=1, expand=True)
        df[epitope_col] = extracted[0]
        if "HLA" not in df.columns:
            df["HLA"] = extracted[1]

    # Ensure required columns exist
    missing = [col for col in [alpha_col, beta_col, epitope_col] if col not in df.columns]
    if missing:
        raise KeyError(f"Missing columns in dataset: {missing}")

    return df

def load_iggytop_anndata(merged_path: str) -> ad.AnnData:
    """Loads the Iggytop merged AnnData file."""
    if not os.path.exists(merged_path):
        raise FileNotFoundError(f"Iggytop dataset not found at {merged_path}")
    print(f"Loading Iggytop dataset from {merged_path}...")
    return sc.read_h5ad(merged_path)

def load_tenx_reference(tenx_path: str = "10x_raw.csv") -> pd.DataFrame:
    """Loads the 10X reference dataset."""
    if not os.path.exists(tenx_path):
        print(f"Warning: 10X reference not found at {tenx_path}")
        return pd.DataFrame()
    print(f"Loading 10X reference from {tenx_path}...")
    df = pd.read_csv(tenx_path)
    
    # Pre-extract peptide from peptide_HLA if present
    if "peptide" not in df.columns and "peptide_HLA" in df.columns:
        df["peptide"] = df["peptide_HLA"].str.split(r"\s+", n=1, expand=True)[0]
    
    return df
