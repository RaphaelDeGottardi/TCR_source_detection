# TCR Source Detection

A pipeline for checking TCR sequence provenance (CDR3a, CDR3b, and Peptide) against the [Iggytop](https://github.com/biocypher/iggytop) and 10X (from [ITRAP(https://github.com/mnielLab/ITRAP_benchmark/blob/main/raw.csv)])databases. It helps identify whether a dataset contains previously known records or excludes them based on their presence in the 10X dataset.

## Features

- **Normalization**: Standardizes CDR3 and peptide sequences (ensures proper start/end characters for TCR chains), analogous to IggyTop.
- **Matching**: Performs exact and partial (single-chain) matching against the Iggytop database.
- **10X Exclusion**: Robustly checks for presence in the `10x_raw.csv` dataset (from [ITRAP(https://github.com/mnielLab/ITRAP_benchmark/blob/main/raw.csv)]), using sequence information to catch records even when metadata is missing.
- **Detailed Reporting**: Generates summary statistics and detailed provenance logs (CSV/Log format).

## Installation

This project uses `uv` for dependency management.

```bash
# Install dependencies
uv sync
```

Also go to [Iggytop](https://github.com/biocypher/iggytop) and follow the instructions to create the merged_adata.h5ad. Add the path to this file to cli.py > DEFAULT_IGGYTOP_PATH. you need to also update the path in the tcr_configs.json if you are using it.
Pro-tip (outcomment the filtering for 10X records in base_adapter to obtain all the records from the [10Xdataset](https://www.10xgenomics.com/library/a14cde) found in the databases)

## Usage

You can run the pipeline using `main.py` or the provided CLI script.

### Example

using a dataset from t2pmhc(https://github.com/qbic-pipelines/t2pmhc/tree/main/data)
(please download the file to this repo to try it out)
```bash
uv run python main.py \
    --dataset t2pmhc_train.tsv \
    --alpha_col cdr3a \
    --beta_col cdr3b \
    --epitope_col peptide \
    --target_col binder \
    --target_pos 1 \
    --save-config t2pmhc_project
```

### Configuration Management

The pipeline supports saving and loading configurations from a `tcr_configs.json` file.

- `--save-config <name>`: Save current CLI arguments as a named configuration.
- `--config <name>`: Load parameters from a saved configuration.
- `--list-configs`: List all saved configurations and exit.

**Example: Running with a saved config**
```bash
uv run python main.py --config t2pmhc_project
```

### Arguments

#### Configuration Management
- `--config`: Name of a saved configuration to use.
- `--save-config`: Save the current CLI arguments as a named configuration.
- `--list-configs`: List all saved configurations and exit.

#### Target Dataset Parameters
- `--dataset`: Path to the input file (CSV, TSV, or XLSX).
- `--alpha_col`: Column name for the Alpha CDR3 sequence.
- `--beta_col`: Column name for the Beta CDR3 sequence.
- `--epitope_col`: Column name for the Epitope/Peptide sequence.
- `--target_col`: (Optional) Column used for filtering (e.g., 'binder').
- `--target_pos`: (Optional) The value to filter by in the target column (default: "1").

#### External Configuration
- `--iggytop_path`: (Optional) Path to the Iggytop `.h5ad` file.
- `--tenx_path`: (Optional) Path to the `10x_raw.csv` reference.
- `--output_log`: (Optional) Path to the log file (default: `coverage_report.log`).

## Output

1. **Console Summary**: Percentages of Exact, Partial, and 10X matches.
2. **`coverage_report.log`**: Appends the summary statistics and configuration for each run.
3. **`detailed_matching_results.csv`**: A copy of the input dataset with added categorization and source provenance columns.
4. **`unexplained_records.csv`**: A separate file containing all records categorized as "Rest" (records that matched neither Iggytop nor 10X databases).

