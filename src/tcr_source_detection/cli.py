import argparse
import os
import pandas as pd
import json
from .loader import load_target_dataset, load_iggytop_anndata, load_tenx_reference
from .normalize import process_cdr3_sequence
from .matcher import TCRMatcher
from .reporter import generate_coverage_report
from .config import save_config, get_config, list_configs

DEFAULT_IGGYTOP_PATH = "/home/dego/.cache/iggytop_airr/merged_anndata.h5ad"
DEFAULT_TENX_PATH = "10x_raw.csv"

def run_pipeline():
    parser = argparse.ArgumentParser(description="TCR Source Detection and Coverage Report")
    
    # Configuration selection and management
    config_group = parser.add_argument_group("Configuration Management")
    config_group.add_argument("--config", help="Name of a saved configuration to use")
    config_group.add_argument("--save-config", help="Save the current CLI arguments as a named configuration")
    config_group.add_argument("--list-configs", action="store_true", help="List all saved configurations and exit")
    
    # Target dataset parameters
    target_group = parser.add_argument_group("Target Dataset Parameters")
    target_group.add_argument("--dataset", help="Path to the target TCR dataset (CSV, TSV, XLSX)")
    target_group.add_argument("--alpha_col", help="Column name for alpha CDR3")
    target_group.add_argument("--beta_col", help="Column name for beta CDR3")
    target_group.add_argument("--epitope_col", help="Column name for epitope/peptide")
    target_group.add_argument("--target_col", help="Column name for filtering (optional)")
    target_group.add_argument("--target_pos", default=None, help="Positive value for target_col filtering")
    
    # External configuration
    ref_group = parser.add_argument_group("External Configuration")
    ref_group.add_argument("--iggytop_path", default=DEFAULT_IGGYTOP_PATH, help="Path to the Iggytop h5ad file")
    ref_group.add_argument("--tenx_path", default=DEFAULT_TENX_PATH, help="Path to 10X reference CSV")
    ref_group.add_argument("--output_log", default="coverage_report.log", help="Path to log file")
    
    args = parser.parse_args()

    if args.list_configs:
        list_configs()
        return

    # 0. Load Configuration if provided
    current_config = {}
    if args.config:
        loaded = get_config(args.config)
        if loaded:
            print(f"Loading configuration '{args.config}'...")
            current_config.update(loaded)
        else:
            print(f"Error: Configuration '{args.config}' not found.")
            return

    # CLI arguments override loaded config
    # We collect all current parameters to potentially save them
    params = {
        "dataset": args.dataset or current_config.get("dataset"),
        "alpha_col": args.alpha_col or current_config.get("alpha_col"),
        "beta_col": args.beta_col or current_config.get("beta_col"),
        "epitope_col": args.epitope_col or current_config.get("epitope_col"),
        "target_col": args.target_col or current_config.get("target_col"),
        "target_pos": args.target_pos or current_config.get("target_pos"),
        "iggytop_path": args.iggytop_path if args.iggytop_path != DEFAULT_IGGYTOP_PATH else (current_config.get("iggytop_path") or args.iggytop_path),
        "tenx_path": args.tenx_path if args.tenx_path != DEFAULT_TENX_PATH else (current_config.get("tenx_path") or args.tenx_path),
        "output_log": args.output_log if args.output_log != "coverage_report.log" else (current_config.get("output_log") or args.output_log)
    }

    if args.save_config:
        # Don't save empty/None values for required fields
        save_config(args.save_config, {k: v for k, v in params.items() if v is not None})
        if not params["dataset"]: # If we are only saving, we might exit or continue
            print("Config saved. To run, ensure --dataset is provided.")
            return

    # Validate required parameters
    required = ["dataset", "alpha_col", "beta_col", "epitope_col"]
    missing = [r for r in required if not params.get(r)]
    if missing:
        parser.error(f"The following arguments are required: {', '.join('--'+m for m in missing)}")
    
    # 1. Load data
    try:
        target_df = load_target_dataset(
            params["dataset"], 
            params["alpha_col"], 
            params["beta_col"], 
            params["epitope_col"], 
            params["target_col"], 
            params["target_pos"]
        )
    except Exception as e:
        print(f"Error loading target dataset: {e}")
        return

    # 2. Load references
    try:
        iggy_adata = load_iggytop_anndata(params["iggytop_path"])
        tenx_df = load_tenx_reference(params["tenx_path"])
    except Exception as e:
        print(f"Error loading reference datasets: {e}")
        return

    # 3. Match and Categorize
    print("Preprocessing sequences and matching...")
    matcher = TCRMatcher(iggy_adata=iggy_adata, tenx_df=tenx_df)
    
    # Column extraction and normalization
    pmid_col = "PMID" if "PMID" in target_df.columns else None
    
    results = []
    source_infos = []
    
    # Pre-calculate normalized columns for efficiency and debugging
    target_df["CDR3a_norm"] = target_df[params["alpha_col"]].apply(lambda x: process_cdr3_sequence(x))
    target_df["CDR3b_norm"] = target_df[params["beta_col"]].apply(lambda x: process_cdr3_sequence(x))
    target_df["Peptide_norm"] = target_df[params["epitope_col"]].astype("string").str.strip().str.upper()
    
    for i, row in target_df.iterrows():
        a = row["CDR3a_norm"]
        b = row["CDR3b_norm"]
        pep = row["Peptide_norm"]
        pmid = row[pmid_col] if pmid_col else ""
        
        category, source_info = matcher.match_record(a, b, pep, pmid)
        results.append(category)
        source_infos.append(source_info)
    
    target_df["Category"] = results
    target_df["SourceInfo"] = source_infos
    
    # 4. Generate Reports
    generate_coverage_report(target_df, params["iggytop_path"], params["dataset"], params["output_log"])

if __name__ == "__main__":
    run_pipeline()
