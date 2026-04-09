# Evo2 MVP Variant Scoring

A pipeline for scoring genetic variants from the Million Veteran Program (MVP) fine-mapping dataset using the [Evo2](https://github.com/evo-design/evo) genomic foundation model. The pipeline computes Evo2 delta scores (reference vs. variant log-likelihood differences) and benchmarks them against established conservation metrics (phastCons, phyloP, GERP).

## Repository Structure

```
Evo2-MVP-Variant-Scoring/
├── MVP_Low-High_Variant_Test/       # Discriminative performance experiment
│   ├── mvp_variants_test.py         # Main Evo2 scoring script (prepare/process modes)
│   ├── filter_and_match_variants.py # Filter MVP data and match high/low PIP variantstop
│   ├── get_conservation_scores.py   # Fetch phastCons, phyloP, GERP scores via REST APIs
│   ├── mvp_plot_script.R            # Faceted boxplots of Evo2 scores by model/context
│   ├── conservation_regression_plots.R  # Logistic regression and forest plots
│   ├── run_MVP_evo2_test.sh         # Single-job SLURM submission (for testing)
│   ├── run_MVP_evo2_worker.sh       # SLURM worker script (chunked processing)
│   └── launch_MVP_jobs.sh           # Batch launcher for multiple model/context configs
│
└── MVP_Timing_Test/                 # Computational efficiency experiment
    ├── timing_variants_test.py      # Evo2 scoring with timing instrumentation
    ├── sample_variants.py           # Stratified sampling (50 high + 50 low PIP)
    ├── timing_plot_script.R         # Timing analysis plots
    ├── run_timing_evo2_worker.sh    # SLURM worker with timing/status tracking
    └── launch_timing_jobs.sh        # Launcher for full timing experiment matrix
```

## Experimental Tracks

### MVP Low-High Variant Test
Evaluates whether Evo2 delta scores can differentiate high-confidence causal variants (PIP ≥ 0.95) from low-confidence variants (PIP < 0.05). Variants are matched on minor allele frequency (MAF), p-value, and VEP annotation using k-nearest neighbors to control for confounders. Results are compared against conservation scores from phastCons, phyloP, and GERP.

### MVP Timing Test
Benchmarks the computational efficiency of Evo2 inference across different configurations: model size (7B vs. 40B), context window length (8 Kbp–512 Kbp), batch/chunk size (10–100 variants), and GPU parallelism strategy (tensor parallelism × context parallelism combinations).

## Data Inputs
| File | Description | Source |
|------|-------------|--------|
| `Data_S1.xlsx` | MVP fine-mapping summary statistics (EUR population) | [Dryad](https://doi.org/10.5061/dryad.zgmsbcck4) |
| GRCh38 FASTA + index | Reference genome (`.fa` + `.fai`) | NCBI/Ensembl |

## Environment and Dependencies

### Python
```
pandas, numpy, scipy, scikit-learn
torch
pyfaidx
requests
biopython
bionemo  # NVIDIA BioNeMo (provides Evo2 model loading and inference)
```

### R
```r
tidyverse, ggplot2, ggpubr, scales, patchwork, car
```

### Hardware and HPC
- NVIDIA GPUs with CUDA support (8 GPUs per job recommended)
- FP8 inference requires compute capability 8.9+ (Ada Lovelace / Hopper architecture)
- SLURM job scheduler
- Apptainer/Singularity container runtime

### Required Environment Variables
To use the code in this repository, you will likely need to update the file following locations embedded in the scripts to ones suited to your application:
```bash
HF_HOME              # Hugging Face model cache directory
NEMO_CACHE_DIR       # NVIDIA NeMo cache directory
BIONEMO_CACHE_DIR    # BioNeMo cache directory
NGC_CLI_API_KEY      # NVIDIA GPU Cloud API key (for model downloads)
PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True
```

## Workflow: Low-High Variant Test

### Step 1 — Filter and match variants
Filters the MVP dataset for EUR-population variants and matches high PIP (≥ 0.95) to low PIP (< 0.05) variants by MAF, p-value, and VEP annotation.

```bash
python filter_and_match_variants.py -i Data_S1.xlsx -o matched_variants.csv
```

### Step 2 — Fetch conservation scores
Retrieves phastCons100way, phyloP100way, and GERP RS scores for each variant. Uses myVariant.info as the primary source with UCSC bigWig tracks as a fallback.

```bash
python get_conservation_scores.py -i matched_variants.csv -o conservation_scores.csv
```

### Step 3 — Score variants with Evo2
Submit SLURM jobs to run Evo2 inference. Use the test script for a quick single job, or the batch launcher to sweep over multiple model and context length configurations.

```bash
# Quick test (single job)
sbatch run_MVP_evo2_test.sh

# Full experiment (multiple model/context configurations)
bash launch_MVP_jobs.sh
```

Each job runs in two phases automatically:
1. **prepare** — Extracts reference and alternate sequences from the FASTA and writes batched input files
2. **process** — Loads PyTorch prediction outputs, computes delta scores (`var_log_probs - ref_log_probs`), and writes a results CSV

### Step 4 — Visualize results
```bash
Rscript mvp_plot_script.R           # Faceted boxplots by model size × context length
Rscript conservation_regression_plots.R  # Logistic regression, forest plot, pairwise correlations
```

## Workflow: Timing Test

### Step 1 — Sample variants
Creates a stratified sample of 100 variants (50 high PIP, 50 low PIP) for use in the timing experiments.

```bash
python sample_variants.py
```

### Step 2 — Submit timing experiments
Launches a matrix of 36 jobs covering:
- 32 jobs: 2 model sizes × 4 chunk sizes × 4 TP/CP parallelism combinations (fixed 8 Kbp context)
- 4 jobs: Context length series (8 Kbp–512 Kbp, fixed TP=4, CP=2, chunk=100)

```bash
bash launch_timing_jobs.sh
```

### Step 3 — Visualize timing results
```bash
Rscript timing_plot_script.R
```
Produces three plots: parallelism strategy comparison, run time vs. context length, and chunk size vs. run time.

## Outputs

| Script | Output File | Format | Description |
|--------|------------|--------|-------------|
| `filter_and_match_variants.py` | `matched_variants_output.csv` | CSV | Balanced high/low PIP variant set |
| `get_conservation_scores.py` | `MVP_conservation_scores.csv` | CSV | phastCons, phyloP, GERP scores per variant |
| `mvp_variants_test.py` | `MVP_variant_scores.*.csv` | CSV | All input columns + `ref_log_probs`, `var_log_probs`, `evo2_delta_score`, `class` |
| `mvp_plot_script.R` | Multiple PNGs + `collated_results.csv` | PNG/CSV | Faceted boxplots; collated scores across all model/context runs |
| `conservation_regression_plots.R` | 5 PNG files | PNG | Regression panels, forest plot, pairwise correlation plots, VIF bar chart, MAF vs Evo2 correlation |
| `timing_variants_test.py` | `timing_scores.*.csv` | CSV | Scored variants with timing metadata |
| `timing_plot_script.R` | 3 PNG files | PNG | Parallelism, context length, chunk size timing plots |
