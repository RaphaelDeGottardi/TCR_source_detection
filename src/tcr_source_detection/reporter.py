import pandas as pd
import os
import json
from datetime import datetime

def generate_coverage_report(df: pd.DataFrame, source_path: str, dataset_path: str, output_log: str = "coverage_report.log"):
    """
    Prints summary stats and writes detailed logs/results to files.
    """
    total_records = len(df)
    counts = df["Category"].value_counts()
    
    # Save detailed CSV
    output_csv = "detailed_matching_results.csv"
    df.to_csv(output_csv, index=False)
    
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    
    def log_and_print(msg, file_handle):
        print(msg)
        file_handle.write(msg + "\n")

    with open(output_log, "a") as f:
        log_and_print(f"\n{'='*50}", f)
        log_and_print(f"REPORT GENERATED: {timestamp}", f)
        log_and_print(f"Source Dataset: {source_path}", f)
        log_and_print(f"Target Dataset: {dataset_path}", f)
        log_and_print(f"Total Records: {total_records:,}", f)
        log_and_print(f"{'-'*30}", f)
        
        categories = [
            "Exact Match (Full)", 
            "Partial Match (A/B)", 
            "Unmatched (10X Exclusion)", 
            "No CDR3 info", 
            "No Epitope", 
            "Multi-epitope (+)", 
            "Rest"
        ]
        
        for cat in categories:
            count = counts.get(cat, 0)
            pct = (count / total_records * 100) if total_records > 0 else 0
            log_and_print(f"{cat:25}: {count:,} ({pct:.2f}%)", f)
        
        log_and_print(f"{'='*50}", f)
        log_and_print(f"Detailed results saved to: {output_csv}", f)
        
    print(f"\nDetailed report appended to: {output_log}")

    # Save unexplained (Rest) records to a separate file
    rest_df = df[df["Category"] == "Rest"]
    if not rest_df.empty:
        rest_csv = "unexplained_records.csv"
        rest_df.to_csv(rest_csv, index=False)
        print(f"Unexplained records (Category: 'Rest') saved to: {rest_csv}")
