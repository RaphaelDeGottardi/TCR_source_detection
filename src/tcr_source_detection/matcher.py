import pandas as pd
import numpy as np
import scirpy as ir
from typing import Set, Tuple, List, Dict
from .normalize import process_cdr3_sequence

TCR_TRIPLET = Tuple[str, str, str]  # (CDR3a, CDR3b, Peptide)
CDR3_PEPTIDE = Tuple[str, str]       # (CDR3, Peptide)

class TCRMatcher:
    """Core logic for matching TCR sequences against reference databases."""
    
    def __init__(self, iggy_adata=None, tenx_df=None):
        self.iggy_full_set: Set[TCR_TRIPLET] = set()
        self.iggy_a_pep: Set[CDR3_PEPTIDE] = set()
        self.iggy_b_pep: Set[CDR3_PEPTIDE] = set()
        self.iggy_source_map: Dict[TCR_TRIPLET, str] = {}
        self.tenx_full_set: Set[TCR_TRIPLET] = set()
        self.tenx_a_pep: Set[CDR3_PEPTIDE] = set()
        self.tenx_b_pep: Set[CDR3_PEPTIDE] = set()
        
        if iggy_adata is not None:
            self._prepare_iggy_sets(iggy_adata)
        if tenx_df is not None:
            self._prepare_tenx_sets(tenx_df)

    def _prepare_iggy_sets(self, adata):
        """Prepare sets and source mappings from Iggytop AnnData."""
        print("Preprocessing Iggytop references...")
        
        subset_df = pd.DataFrame(index=adata.obs_names)
        subset_df["CDR3a"] = ir.get.airr(adata, "junction_aa", chain="VJ_1")
        subset_df["CDR3b"] = ir.get.airr(adata, "junction_aa", chain="VDJ_1")
        subset_df["peptide"] = adata.obs["epitope_sequence"].astype("string")
        
        # Add source information columns to subset_df if they exist
        for col in ["PMID", "source"]:
            if col in adata.obs.columns:
                subset_df[col] = adata.obs[col]
        
        # We also need source information
        source_col = "PMID" if "PMID" in subset_df.columns else "source"
        if source_col not in subset_df.columns:
            source_col = None
            
        def build_set(df, cols, source_col=None):
            records = df[list(cols)].dropna(subset=list(cols)).drop_duplicates()
            if source_col:
                subset_with_source = df[list(cols) + [source_col]].dropna(subset=list(cols))
                # For multiple sources for the same triplet, join them
                grouped = subset_with_source.groupby(list(cols))[source_col].unique().apply(lambda x: ", ".join(map(str, x)))
                source_dict = grouped.to_dict()
                return set(records.itertuples(index=False, name=None)), source_dict
            return set(records.itertuples(index=False, name=None)), {}

        self.iggy_full_set, self.iggy_source_map = build_set(subset_df, ["CDR3a", "CDR3b", "peptide"], source_col)
        self.iggy_a_pep, _ = build_set(subset_df, ["CDR3a", "peptide"])
        self.iggy_b_pep, _ = build_set(subset_df, ["CDR3b", "peptide"])

    def _prepare_tenx_sets(self, tenx_df):
        """Prepare sets from 10X reference CSV."""
        print("Preprocessing 10X reference...")
        
        # 10X columns are different: cdr3_TRA, cdr3_TRB, peptide (extracted)
        norm_df = tenx_df.copy()
        norm_df["CDR3a"] = norm_df["cdr3_TRA"].apply(lambda x: process_cdr3_sequence(x, is_igh=False))
        norm_df["CDR3b"] = norm_df["cdr3_TRB"].apply(lambda x: process_cdr3_sequence(x, is_igh=False))
        norm_df["peptide"] = norm_df["peptide"].astype("string").str.strip().str.upper()
        
        def build_set(df, cols):
            return set(df[list(cols)].dropna(subset=list(cols)).drop_duplicates().itertuples(index=False, name=None))
            
        self.tenx_full_set = build_set(norm_df, ["CDR3a", "CDR3b", "peptide"])
        self.tenx_a_pep = build_set(norm_df, ["CDR3a", "peptide"])
        self.tenx_b_pep = build_set(norm_df, ["CDR3b", "peptide"])

    def match_record(self, a_norm: str, b_norm: str, pep_norm: str, pmid: str = "") -> Tuple[str, str]:
        """
        Categorizes a single TCR record.
        Returns: (category, source_info)
        """
        source_info = ""
        
        if pd.isna(pep_norm) or pep_norm == "":
            return "No Epitope", ""
        
        if "+" in str(pep_norm):
            return "Multi-epitope (+)", ""
            
        if (pd.isna(a_norm) or a_norm == "") and (pd.isna(b_norm) or b_norm == ""):
            return "No CDR3 info", ""

        # Pre-check for 10X exclusion from input pmid or sequence match in local 10X ref
        triplet = (a_norm, b_norm, pep_norm)
        is_tenx_by_seq = triplet in self.tenx_full_set or \
                        ((a_norm, pep_norm) in self.tenx_a_pep if not pd.isna(a_norm) else False) or \
                        ((b_norm, pep_norm) in self.tenx_b_pep if not pd.isna(b_norm) else False)
        
        if str(pmid).startswith("https://www.10x") or is_tenx_by_seq:
            return "Unmatched (10X Exclusion)", "10X"

        # Priority 1: Exact Match in Iggytop
        if triplet in self.iggy_full_set:
            source_info = self.iggy_source_map.get(triplet, "")
            
            # Sub-check: If matching part of Iggytop that is actually 10X
            if "https://www.10xgenomics.com" in source_info or "no_pmid_1036521" in source_info:
                return "Unmatched (10X Exclusion)", f"Iggytop/10X ({source_info})"
                
            return "Exact Match (Full)", source_info

        # Priority 2: Partial Match in Iggytop
        a_pep_match = (a_norm, pep_norm) in self.iggy_a_pep if not pd.isna(a_norm) else False
        b_pep_match = (b_norm, pep_norm) in self.iggy_b_pep if not pd.isna(b_norm) else False
        if a_pep_match or b_pep_match:
            return "Partial Match (A/B)", ""

        # Priority 4: Everything Else
        return "Rest", ""
