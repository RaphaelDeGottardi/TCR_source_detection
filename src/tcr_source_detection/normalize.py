import pandas as pd

AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")

def is_valid_peptide_sequence(seq: str) -> bool:
    """Checks if a given sequence is a valid peptide sequence."""
    if isinstance(seq, str) and len(seq) > 2:
        return all([aa in AMINO_ACIDS for aa in seq])
    else:
        return False

def process_cdr3_sequence(seq: str, is_igh: bool = False) -> str | None:
    """
    Normalizes a CDR3 sequence:
    - Upper case, stripped of whitespace.
    - Validates amino acids.
    - Ensures it starts with 'C' and ends with 'F' (or 'W' if IGH).
    """
    if seq is None or pd.isna(seq):
        return None

    # Clean and normalize the sequence
    seq = str(seq).upper().strip().replace(" ", "").replace("\n", "")

    # Validate that the sequence contains only valid amino acids
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
