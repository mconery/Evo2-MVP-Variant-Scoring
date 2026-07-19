################################################################################
# finemap_CARMA_t2d.R
#
# CARMA fine-mapping for MVP T2D EUR loci.
# Adapted from Betanfer/paper_analyses/.../CARMA/finemap_CARMA.R.
#
# Usage:
#   Rscript finemap_CARMA_t2d.R \
#       <sumstats_tsv> <vcor1_file> <out_file> \
#       [--use-priors] [--prior-file <file>] [--input-alpha <float>] \
#       [--outlier-check TRUE|FALSE] [--lambda <float>] [--rho <float>] \
#       [--seed <int>]
#
# Required positional arguments:
#   1: summary statistics TSV  (columns: SNP_ID, BETA, SE, N, ...)
#   2: plink2 .vcor1 LD matrix file
#   3: output TSV file path
#
# Optional keyword arguments (any order after position 3):
#   --use-priors             enable Evo2-informed priors (requires --prior-file)
#   --prior-file <file>      per-variant prior weight TSV (SNP_ID, prior_weight)
#   --input-alpha <float>    weight of functional prior in CARMA (default: 0.1)
#   --outlier-check <bool>   enable outlier detection (default: TRUE)
#   --lambda <float>         Poisson prior on causal count (default: 1)
#   --rho <float>            credible set confidence (default: 0.99)
#   --seed <int>             random seed (default: 5)
#
# Output columns:
#   SNP_ID, Z_SCORE, PIP, CS_ID (0 = not in any CS), PRIOR_WEIGHT (NA if uniform)
################################################################################

library(CARMA)

# ---------------------------------------------------------------------------
# 1) Parse arguments
# ---------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  cat("Usage: Rscript finemap_CARMA_t2d.R sumstats.tsv matrix.vcor1 out.tsv [options]\n")
  quit(save = "no")
}

sumstats_file <- args[1]
ld_file       <- args[2]
out_file      <- args[3]

# Defaults
use_priors      <- FALSE
prior_file      <- NULL
input_alpha     <- 0.1
outlier_check   <- TRUE
lambda          <- 1.0
rho             <- 0.99
random_seed     <- 5L

# Parse keyword args
i <- 4
while (i <= length(args)) {
  switch(args[i],
    "--use-priors"     = { use_priors <- TRUE },
    "--prior-file"     = { i <- i + 1; prior_file <- args[i] },
    "--input-alpha"    = { i <- i + 1; input_alpha <- as.numeric(args[i]) },
    "--outlier-check"  = { i <- i + 1
                           first <- tolower(substr(args[i], 1, 1))
                           outlier_check <- (first == "t" || first == "y") },
    "--lambda"         = { i <- i + 1; lambda <- as.numeric(args[i]) },
    "--rho"            = { i <- i + 1; rho <- as.numeric(args[i]) },
    "--seed"           = { i <- i + 1; random_seed <- as.integer(args[i]) }
  )
  i <- i + 1
}

# ---------------------------------------------------------------------------
# 2) Validate inputs
# ---------------------------------------------------------------------------

error_ind <- FALSE

if (!file.exists(sumstats_file)) {
  cat(paste0("ERROR: sumstats file not found: ", sumstats_file, "\n"))
  error_ind <- TRUE
}

if (!file.exists(ld_file)) {
  cat(paste0("ERROR: LD matrix not found: ", ld_file, "\n"))
  error_ind <- TRUE
}

if (!dir.exists(dirname(out_file))) {
  cat(paste0("ERROR: output directory does not exist: ", dirname(out_file), "\n"))
  error_ind <- TRUE
}

if (use_priors) {
  if (is.null(prior_file) || !file.exists(prior_file)) {
    cat(paste0("ERROR: --use-priors set but prior file not found: ", prior_file, "\n"))
    error_ind <- TRUE
  }
}

if (error_ind) quit(save = "no")

setwd(dirname(out_file))

# ---------------------------------------------------------------------------
# 3) Load LD matrix and variant list
# ---------------------------------------------------------------------------

vars_file <- paste0(ld_file, ".vars")
if (!file.exists(vars_file)) {
  cat(paste0("ERROR: variant list not found: ", vars_file, "\n"))
  quit(save = "no")
}

var_ids <- read.csv(vars_file, sep = "\t", header = FALSE)[, 1]
R <- as.matrix(read.csv(ld_file, sep = "\t", header = FALSE))
rownames(R) <- var_ids
colnames(R) <- var_ids

# ---------------------------------------------------------------------------
# 4) Load and align summary statistics
# ---------------------------------------------------------------------------

ss <- read.table(sumstats_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

required <- c("SNP_ID", "BETA", "SE")
missing_cols <- required[!required %in% colnames(ss)]
if (length(missing_cols) > 0) {
  cat(paste0("ERROR: missing columns in sumstats: ", paste(missing_cols, collapse = ", "), "\n"))
  quit(save = "no")
}

ss <- ss[ss$SNP_ID %in% var_ids, ]
if (nrow(ss) == 0) {
  cat("ERROR: no variants from LD matrix found in summary statistics.\n")
  quit(save = "no")
}

ss <- ss[match(var_ids, ss$SNP_ID), ]
missing_mask <- is.na(ss$SNP_ID)
if (any(missing_mask)) {
  cat(paste0("WARNING: ", sum(missing_mask), " LD matrix variants not in sumstats; dropping.\n"))
  keep <- !missing_mask
  ss    <- ss[keep, ]
  R     <- R[keep, keep]
  var_ids <- var_ids[keep]
}

z_scores <- ss$BETA / ss$SE
z.list     <- list(z_scores)
ld.list    <- list(as.matrix(R))
lambda.list <- list(lambda)

# ---------------------------------------------------------------------------
# 5) Load priors (optional)
# ---------------------------------------------------------------------------

w.list <- NULL

if (use_priors) {
  pw <- read.table(prior_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  # Align prior weights to the order of variants in the LD matrix / sumstats
  weights_aligned <- pw$prior_weight[match(var_ids, pw$SNP_ID)]
  # Variants with no Evo2 score get a neutral weight of 1.0
  weights_aligned[is.na(weights_aligned)] <- 1.0
  w.list <- list(weights_aligned)
  cat(paste0("Priors loaded: ", sum(!is.na(pw$prior_weight[match(var_ids, pw$SNP_ID)])),
             "/", length(var_ids), " variants have Evo2 weights (input.alpha=", input_alpha, ")\n"))
}

# ---------------------------------------------------------------------------
# 6) Define CARMA caller
# ---------------------------------------------------------------------------

call_carma <- function(z.list, ld.list, lambda.list, rho, outlier_check, w.list, input_alpha) {
  suppressWarnings(
    CARMA(
      z.list         = z.list,
      ld.list        = ld.list,
      w.list         = w.list,
      lambda.list    = lambda.list,
      input.alpha    = if (!is.null(w.list)) input_alpha else 0,
      outlier.switch = outlier_check,
      rho.index      = rho,
      BF.index       = 1000000,
      Max.Model.Dim  = 10000000
    )
  )
}

# ---------------------------------------------------------------------------
# 7) Execute CARMA (retry loop)
# ---------------------------------------------------------------------------

attempts <- 0
max_attempts <- 10
success <- FALSE
carma_time_seconds <- NA

while (!success && attempts < max_attempts) {
  tryCatch({
    set.seed(random_seed)
    t_start <- Sys.time()
    carma_results <- call_carma(z.list, ld.list, lambda.list, rho, outlier_check, w.list, input_alpha)
    t_end <- Sys.time()
    carma_time_seconds <- as.numeric(difftime(t_end, t_start, units = "secs"))
    success <- TRUE
  }, error = function(e) {
    attempts    <<- attempts + 1
    random_seed <<- random_seed * 2L
    cat(paste0("WARNING: attempt ", attempts, " failed; retrying with seed ", random_seed, "\n"))
  })
}

# ---------------------------------------------------------------------------
# 8) Build output data frame
# ---------------------------------------------------------------------------

prior_weight_col <- if (!is.null(w.list)) w.list[[1]] else rep(NA_real_, length(var_ids))

if (!success) {
  cat(paste0("ERROR: CARMA failed after ", max_attempts, " attempts.\n"))
  out_data <- data.frame(
    SNP_ID       = var_ids,
    Z_SCORE      = z_scores,
    PIP          = NA_real_,
    CS_ID        = 0L,
    PRIOR_WEIGHT = prior_weight_col
  )
} else {
  cat(paste0("SUCCESS: CARMA completed in ", round(carma_time_seconds, 2), "s after ",
             attempts, " attempt(s).\n"))

  pips <- carma_results[[1]]$PIPs

  cs_ids <- integer(length(var_ids))
  cs_list <- carma_results[[1]]$"Credible set"[[2]]
  if (length(cs_list) > 0) {
    for (s in seq_along(cs_list)) {
      cs_ids[cs_list[[s]]] <- as.integer(s)
    }
  }

  out_data <- data.frame(
    SNP_ID       = var_ids,
    Z_SCORE      = z_scores,
    PIP          = pips,
    CS_ID        = cs_ids,
    PRIOR_WEIGHT = prior_weight_col
  )
}

write.table(out_data, file = out_file,
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

cat(paste0("Output written to: ", out_file, "\n"))
cat(paste0("Variants in any CS: ", sum(out_data$CS_ID > 0), "\n"))
