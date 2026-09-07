################################################################################
# compare_finemapping.R
#
# Collates CARMA fine-mapping results (with and without Evo2 priors) and
# compares them against each other and against the original MVP phenome-wide
# mapping results from Supplementary Table S11.
#
# Metrics:
#   - Per-locus credible set (CS) sizes and counts
#   - Top-PIP variant per locus
#   - Jaccard overlap between uniform-prior CS and Evo2-prior CS per locus
#   - Overlap of each new CS with the original S11 EUR credible set
#   - Paired Wilcoxon signed-rank test: CS size, uniform vs Evo2 prior
#   - Paired Wilcoxon signed-rank test: Jaccard overlap with original S11 CS,
#     uniform vs Evo2 prior
#   - Capture (presence in CS) and assigned PIP of original high-confidence
#     (PIP > 0.95) S11 variants under each new approach, with a paired
#     Wilcoxon signed-rank test on assigned PIP
#
# Outputs (written to $BASE/collation/):
#   per_locus_summary.tsv       -- one row per locus
#   cs_comparison.tsv           -- one row per locus × CS signal
#   s11_highpip_capture.tsv     -- one row per original high-PIP (>0.95) signal
#   aggregate_metrics.txt       -- printed summary statistics + Wilcoxon tests
#   plots/cs_size_violin.pdf    -- includes paired Wilcoxon p-value
#   plots/pip_scatter.pdf
#   plots/jaccard_histogram.pdf
#   plots/venn_overlap.pdf      -- pooled CS-variant overlap: S11 vs uniform vs Evo2
#   plots/locuszoom/{locus_id}_locuszoom.pdf
#                               -- one per locus, 3 stacked panels: GWAS -log10(p),
#                                  uniform-prior PIP, Evo2-prior PIP, sharing a
#                                  genomic-position x-axis
################################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(readr)
  library(purrr)
  library(ggvenn)
  library(patchwork)
})

BASE       <- "/grand/GeomicVar/mconery/evo2_variant_scoring_mapping"
S11_FILE   <- paste0(BASE, "/mapping_loci/Supplementary_Table-S11.txt")
NO_PRIOR   <- paste0(BASE, "/carma_results/without_priors")
WITH_PRIOR <- paste0(BASE, "/carma_results/with_priors")
LOCI_FILE  <- paste0(BASE, "/loci_definition/t2d_eur_loci.tsv")
VARIANT_LIST_DIR <- paste0(BASE, "/variant_lists")
OUT_DIR    <- paste0(BASE, "/collation")

dir.create(paste0(OUT_DIR, "/plots"), recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# 1) Load loci list
# ---------------------------------------------------------------------------

loci <- read_tsv(LOCI_FILE, col_types = cols(.default = "c"), show_col_types = FALSE)
cat(sprintf("Loaded %d loci\n", nrow(loci)))

# ---------------------------------------------------------------------------
# 2) Load S11 original fine-mapping results for T2D EUR
# ---------------------------------------------------------------------------

s11_raw <- read_tsv(S11_FILE, col_types = cols(.default = "c"), show_col_types = FALSE)
s11_raw <- s11_raw %>% filter(phenotype != "phenotype")  # drop duplicate header

s11_t2d <- s11_raw %>%
  filter(phenotype == "Phe_250_2",
         !is.na(EUR.best_variant) & EUR.best_variant != "") %>%
  select(locus_id = locus, signal,
         EUR.best_variant, EUR.max_overall_pip, EUR.variant_ids) %>%
  mutate(signal = as.integer(signal),
         EUR.max_overall_pip = as.numeric(EUR.max_overall_pip))

# Parse EUR.variant_ids: comma-separated list of variants in the original CS
# (NOTE: fixed from semicolon -- the S11 file delimits multi-variant fields
# with commas, e.g. "10:114758349:C:T,10:114754071:T:C")
s11_t2d <- s11_t2d %>%
  mutate(orig_cs_variants = str_split(EUR.variant_ids, ","))

cat(sprintf("S11: %d T2D EUR signals across %d unique loci\n",
            nrow(s11_t2d), n_distinct(s11_t2d$locus_id)))

# High-confidence original signals: those whose best variant reached
# PIP > 0.95 in the original SuSiE (in-sample LD) analysis. S11 only reports
# a PIP for each signal's best variant, so "high-PIP variant" here means
# a signal-level best variant clearing that threshold.
s11_highpip <- s11_t2d %>%
  filter(!is.na(EUR.max_overall_pip), EUR.max_overall_pip > 0.95) %>%
  select(locus_id, signal, orig_variant = EUR.best_variant, orig_pip = EUR.max_overall_pip)

cat(sprintf("S11: %d high-PIP (>0.95) original signals across %d loci\n",
            nrow(s11_highpip), n_distinct(s11_highpip$locus_id)))

# ---------------------------------------------------------------------------
# Helper: load CARMA output for all loci from a directory
# ---------------------------------------------------------------------------

load_carma_dir <- function(dir_path, label) {
  files <- list.files(dir_path, pattern = "\\.carma\\.tsv$", full.names = TRUE)
  if (length(files) == 0) {
    warning(sprintf("No CARMA outputs found in %s", dir_path))
    return(tibble())
  }
  purrr::map_dfr(files, function(f) {
    locus_id <- str_replace(basename(f), "\\.carma\\.tsv$", "")
    tryCatch(
      read_tsv(f, col_types = cols(.default = "c"), show_col_types = FALSE) %>%
        mutate(locus_id = locus_id,
               approach = label,
               PIP     = as.numeric(PIP),
               Z_SCORE = as.numeric(Z_SCORE),
               CS_ID   = as.integer(CS_ID),
               PRIOR_WEIGHT = suppressWarnings(as.numeric(PRIOR_WEIGHT))),
      error = function(e) {
        warning(sprintf("Failed to load %s: %s", f, e$message))
        tibble()
      }
    )
  })
}

cat("Loading uniform-prior results ...\n")
res_no  <- load_carma_dir(NO_PRIOR,   "uniform")
cat(sprintf("  %d loci loaded\n", n_distinct(res_no$locus_id)))

cat("Loading Evo2-prior results ...\n")
res_evo <- load_carma_dir(WITH_PRIOR, "evo2_prior")
cat(sprintf("  %d loci loaded\n", n_distinct(res_evo$locus_id)))

# ---------------------------------------------------------------------------
# 2b) Diagnostics: do S11 and CARMA actually key on the same IDs?
# ---------------------------------------------------------------------------
# If every locus shows zero S11 overlap, that is almost always a join/ID-format
# problem rather than a real biological result -- check both possible failure
# points before trusting jaccard_orig_vs_* downstream.

cat("\n--- Diagnostic: S11 vs CARMA identifier formats ---\n")
cat("S11 locus_id examples:     ",
    paste(head(unique(s11_t2d$locus_id), 5), collapse = " | "), "\n")
cat("CARMA locus_id examples:   ",
    paste(head(unique(res_no$locus_id), 5), collapse = " | "), "\n")

n_locus_matches <- length(intersect(s11_t2d$locus_id, res_no$locus_id))
cat(sprintf("Locus IDs shared between S11 and CARMA: %d / %d S11 loci, %d / %d CARMA loci\n",
            n_locus_matches, n_distinct(s11_t2d$locus_id),
            n_locus_matches, n_distinct(res_no$locus_id)))
if (n_locus_matches == 0) {
  warning(paste(
    "No locus_id overlap at all between S11 ('locus' column, e.g.",
    "'chr10.114250001.115250000') and CARMA output (derived from filenames).",
    "Every S11 overlap metric (jaccard_orig_vs_uniform/evo2, high-PIP capture)",
    "will be degenerate (NA/0) until locus_id formats are reconciled -- check",
    "delimiter (dot vs underscore vs colon/dash), chr prefix, and coordinate",
    "convention (start/end vs signal-specific ranges) on both sides."
  ))
}

s11_variant_examples  <- head(unique(unlist(s11_t2d$orig_cs_variants)), 5)
carma_variant_examples <- head(unique(res_no$SNP_ID), 5)
cat("S11 variant ID examples:   ", paste(s11_variant_examples, collapse = " | "), "\n")
cat("CARMA variant ID examples: ", paste(carma_variant_examples, collapse = " | "), "\n")

n_variant_matches <- length(intersect(unique(unlist(s11_t2d$orig_cs_variants)),
                                       unique(res_no$SNP_ID)))
cat(sprintf("Variant IDs shared between S11 and CARMA, genome-wide: %d\n", n_variant_matches))
if (n_variant_matches == 0) {
  warning(paste(
    "No variant ID overlap at all between S11 and CARMA output, even before",
    "restricting to matching loci -- CARMA output appears to use rsIDs while",
    "S11 uses chr:pos:ref:alt positional IDs. Resolved below by translating",
    "CARMA's rsIDs via the per-locus variant list files."
  ))
}
cat("--- End diagnostic ---\n\n")

# ---------------------------------------------------------------------------
# 2c) Translate CARMA rsIDs to chr:pos:ref:alt IDs via variant list files
# ---------------------------------------------------------------------------
# CARMA output identifies variants by rsID (SNP_ID); S11 identifies them by
# chr:pos:ref:alt. Each locus has a variant list file (named the same way as
# the locus itself, e.g. "chr10.114250001.115250000.variants.tsv") with
# columns SNP_ID/CHR/POS/REF/ALT that let us build that mapping per locus.

cat("Loading variant ID maps (rsID -> chr:pos:ref:alt) from variant list files ...\n")

all_res_loci <- union(unique(res_no$locus_id), unique(res_evo$locus_id))

variant_map <- purrr::map_dfr(all_res_loci, function(lid) {
  f <- paste0(VARIANT_LIST_DIR, "/", lid, ".variants.tsv")
  if (!file.exists(f)) {
    warning(sprintf("No variant list file found for locus %s (expected %s)", lid, f))
    return(tibble())
  }
  tryCatch(
    read_tsv(f, col_types = cols(.default = "c"), show_col_types = FALSE) %>%
      transmute(locus_id = lid,
                SNP_ID   = SNP_ID,
                CHR      = CHR,
                POS      = POS,
                PVAL     = PVAL,
                pos_id   = paste(CHR, POS, REF, ALT, sep = ":")),
    error = function(e) {
      warning(sprintf("Failed to load variant list %s: %s", f, e$message))
      tibble()
    }
  )
})

cat(sprintf("  Loaded variant ID maps for %d loci (%d variant rows total)\n",
            n_distinct(variant_map$locus_id), nrow(variant_map)))

res_no <- res_no %>%
  left_join(variant_map %>% select(locus_id, SNP_ID, pos_id), by = c("locus_id", "SNP_ID"))
res_evo <- res_evo %>%
  left_join(variant_map %>% select(locus_id, SNP_ID, pos_id), by = c("locus_id", "SNP_ID"))

n_unmapped_no  <- sum(is.na(res_no$pos_id))
n_unmapped_evo <- sum(is.na(res_evo$pos_id))
if (n_unmapped_no > 0 || n_unmapped_evo > 0) {
  warning(sprintf(
    "%d uniform-prior and %d Evo2-prior CARMA variants had no matching rsID in their locus's variant list file; these will not contribute to S11 overlap metrics.",
    n_unmapped_no, n_unmapped_evo
  ))
}

# Post-translation sanity check: overlap should now be non-trivial if the
# translation worked.
n_variant_matches_post <- length(intersect(unique(unlist(s11_t2d$orig_cs_variants)),
                                            unique(na.omit(res_no$pos_id))))
cat(sprintf("Post-translation: %d variant IDs now shared between S11 and CARMA (uniform), genome-wide\n\n",
            n_variant_matches_post))

# ---------------------------------------------------------------------------
# 3) Per-locus CS summary
# ---------------------------------------------------------------------------

summarize_locus <- function(df, approach_label) {
  df %>%
    group_by(locus_id) %>%
    summarise(
      approach          = approach_label,
      n_variants        = n(),
      n_cs              = n_distinct(CS_ID[CS_ID > 0]),
      total_cs_size     = sum(CS_ID > 0),
      median_cs_size    = if (n_distinct(CS_ID[CS_ID > 0]) == 0) NA_real_ else {
        cs_sizes <- table(CS_ID[CS_ID > 0])
        median(as.integer(cs_sizes))
      },
      top_pip           = max(PIP, na.rm = TRUE),
      top_variant       = SNP_ID[which.max(PIP)],
      n_singleton_cs    = {
        if (n_distinct(CS_ID[CS_ID > 0]) == 0) 0L
        else sum(table(CS_ID[CS_ID > 0]) == 1L)
      },
      .groups = "drop"
    )
}

sum_no  <- summarize_locus(res_no,  "uniform")
sum_evo <- summarize_locus(res_evo, "evo2_prior")

# ---------------------------------------------------------------------------
# 4) Jaccard overlap between uniform and Evo2 CS per locus
# ---------------------------------------------------------------------------

# NOTE: returns NA if EITHER set is empty, not just if both are. A 0 here
# should mean "we compared two non-trivial sets and they truly don't overlap" --
# not "one side had no data" (e.g. because of a join/parsing mismatch). Treating
# an empty-vs-nonempty case as a hard 0 would silently mask exactly that failure
# mode, which is otherwise indistinguishable from genuine non-overlap.
jaccard <- function(a, b) {
  if (length(a) == 0 || length(b) == 0) return(NA_real_)
  length(intersect(a, b)) / length(union(a, b))
}

loci_with_both <- intersect(sum_no$locus_id, sum_evo$locus_id)

jaccard_df <- purrr::map_dfr(loci_with_both, function(lid) {
  cs_no  <- res_no  %>% filter(locus_id == lid, CS_ID > 0) %>% pull(SNP_ID)
  cs_evo <- res_evo %>% filter(locus_id == lid, CS_ID > 0) %>% pull(SNP_ID)
  tibble(
    locus_id = lid,
    jaccard_uniform_vs_evo2 = jaccard(cs_no, cs_evo),
    cs_size_uniform         = length(cs_no),
    cs_size_evo2            = length(cs_evo),
    cs_size_delta           = length(cs_evo) - length(cs_no)
  )
})

# ---------------------------------------------------------------------------
# 5) Overlap with original S11 credible sets
# ---------------------------------------------------------------------------

overlap_s11 <- purrr::map_dfr(loci_with_both, function(lid) {
  orig_variants <- s11_t2d %>%
    filter(locus_id == lid) %>%
    pull(orig_cs_variants) %>%
    unlist() %>%
    unique()

  # Use pos_id (chr:pos:ref:alt), not the raw rsID SNP_ID, since S11 identifies
  # variants positionally -- see section 2c for the rsID -> pos_id translation.
  cs_no  <- res_no  %>% filter(locus_id == lid, CS_ID > 0) %>% pull(pos_id) %>% na.omit()
  cs_evo <- res_evo %>% filter(locus_id == lid, CS_ID > 0) %>% pull(pos_id) %>% na.omit()

  s11_best_variants <- s11_t2d$EUR.best_variant[s11_t2d$locus_id == lid]

  tibble(
    locus_id                     = lid,
    jaccard_orig_vs_uniform      = jaccard(orig_variants, cs_no),
    jaccard_orig_vs_evo2         = jaccard(orig_variants, cs_evo),
    n_orig_variants              = length(orig_variants),
    orig_top_in_uniform_cs       = any(s11_best_variants %in% cs_no),
    orig_top_in_evo2_cs          = any(s11_best_variants %in% cs_evo)
  )
})

# ---------------------------------------------------------------------------
# 5b) Capture of original high-PIP (>0.95) variants by each new approach
# ---------------------------------------------------------------------------
# For each S11 signal with a best variant PIP > 0.95, check whether uniform-
# and Evo2-prior CARMA put that variant in a credible set, and what PIP each
# approach assigned it. A variant absent from an approach's output entirely
# is treated as PIP = 0 for that approach (it was not captured at all).

lookup_pip <- function(df, lid, variant) {
  hit <- df %>% filter(locus_id == lid, pos_id == variant)
  if (nrow(hit) == 0) return(tibble(pip = NA_real_, in_cs = FALSE))
  tibble(pip = max(hit$PIP, na.rm = TRUE), in_cs = any(hit$CS_ID > 0, na.rm = TRUE))
}

highpip_capture <- purrr::pmap_dfr(
  list(s11_highpip$locus_id, s11_highpip$signal, s11_highpip$orig_variant, s11_highpip$orig_pip),
  function(lid, sig, variant, orig_pip) {
    u <- lookup_pip(res_no,  lid, variant)
    e <- lookup_pip(res_evo, lid, variant)
    tibble(
      locus_id      = lid,
      signal        = sig,
      orig_variant  = variant,
      orig_pip      = orig_pip,
      uniform_pip   = u$pip,
      uniform_in_cs = u$in_cs,
      evo2_pip      = e$pip,
      evo2_in_cs    = e$in_cs
    )
  }
)

highpip_capture <- highpip_capture %>%
  mutate(uniform_pip_filled = ifelse(is.na(uniform_pip), 0, uniform_pip),
         evo2_pip_filled    = ifelse(is.na(evo2_pip),    0, evo2_pip))

write_tsv(highpip_capture, paste0(OUT_DIR, "/s11_highpip_capture.tsv"))
cat(sprintf("Written: %s/s11_highpip_capture.tsv\n", OUT_DIR))

# ---------------------------------------------------------------------------
# 6) Assemble per-locus summary table
# ---------------------------------------------------------------------------

per_locus <- loci %>%
  left_join(sum_no  %>% select(locus_id, n_cs_uniform = n_cs,
                                top_pip_uniform = top_pip,
                                top_var_uniform = top_variant,
                                n_singleton_uniform = n_singleton_cs,
                                total_cs_size_uniform = total_cs_size),
            by = "locus_id") %>%
  left_join(sum_evo %>% select(locus_id, n_cs_evo2 = n_cs,
                                top_pip_evo2 = top_pip,
                                top_var_evo2 = top_variant,
                                n_singleton_evo2 = n_singleton_cs,
                                total_cs_size_evo2 = total_cs_size),
            by = "locus_id") %>%
  left_join(jaccard_df,  by = "locus_id") %>%
  left_join(overlap_s11, by = "locus_id")

write_tsv(per_locus, paste0(OUT_DIR, "/per_locus_summary.tsv"))
cat(sprintf("Written: %s/per_locus_summary.tsv\n", OUT_DIR))

# ---------------------------------------------------------------------------
# 6b) Statistical tests (paired by locus, since each locus is mapped once
#     under each approach)
# ---------------------------------------------------------------------------

# CS size: uniform vs Evo2 prior
cs_size_paired <- per_locus %>%
  filter(!is.na(total_cs_size_uniform), !is.na(total_cs_size_evo2))

wt_cs <- if (nrow(cs_size_paired) >= 2) {
  wilcox.test(cs_size_paired$total_cs_size_evo2,
              cs_size_paired$total_cs_size_uniform,
              paired = TRUE)
} else NULL

# Overlap with original S11 CS: uniform vs Evo2 prior
jaccard_paired <- overlap_s11 %>%
  filter(!is.na(jaccard_orig_vs_uniform), !is.na(jaccard_orig_vs_evo2))

wt_jaccard <- if (nrow(jaccard_paired) >= 2) {
  wilcox.test(jaccard_paired$jaccard_orig_vs_evo2,
              jaccard_paired$jaccard_orig_vs_uniform,
              paired = TRUE)
} else NULL

# PIP assigned to original high-PIP variants: uniform vs Evo2 prior
wt_pip <- if (nrow(highpip_capture) >= 2) {
  wilcox.test(highpip_capture$evo2_pip_filled,
              highpip_capture$uniform_pip_filled,
              paired = TRUE)
} else NULL

# ---------------------------------------------------------------------------
# 7) Per-CS comparison table
# ---------------------------------------------------------------------------

get_cs_rows <- function(df, approach_label) {
  df %>%
    filter(CS_ID > 0) %>%
    group_by(locus_id, CS_ID) %>%
    summarise(
      approach       = approach_label,
      cs_size        = n(),
      top_pip        = max(PIP, na.rm = TRUE),
      top_variant    = SNP_ID[which.max(PIP)],
      top_variant_pos_id = pos_id[which.max(PIP)],
      variants       = paste(SNP_ID, collapse = ";"),
      variants_pos_id = paste(pos_id, collapse = ";"),
      .groups = "drop"
    )
}

cs_no  <- get_cs_rows(res_no,  "uniform")
cs_evo <- get_cs_rows(res_evo, "evo2_prior")
cs_all <- bind_rows(cs_no, cs_evo)
write_tsv(cs_all, paste0(OUT_DIR, "/cs_comparison.tsv"))
cat(sprintf("Written: %s/cs_comparison.tsv\n", OUT_DIR))

# ---------------------------------------------------------------------------
# 8) Aggregate metrics
# ---------------------------------------------------------------------------

agg_metrics <- function(label, sumdf, jdf) {
  cat(sprintf("\n=== %s ===\n", label))
  cat(sprintf("  Loci with ≥1 CS:          %d / %d\n",
              sum(sumdf$n_cs > 0, na.rm = TRUE), nrow(sumdf)))
  cat(sprintf("  Loci with ≥1 singleton CS: %d\n",
              sum(sumdf$n_singleton_cs > 0, na.rm = TRUE)))
  cat(sprintf("  Median total CS size:       %.1f\n",
              median(sumdf$total_cs_size, na.rm = TRUE)))
  cat(sprintf("  Mean total CS size:         %.1f\n",
              mean(sumdf$total_cs_size, na.rm = TRUE)))
  cat(sprintf("  Median top PIP:             %.3f\n",
              median(sumdf$top_pip, na.rm = TRUE)))
}

sink(paste0(OUT_DIR, "/aggregate_metrics.txt"))
cat("=== CARMA Fine-Mapping Comparison: T2D EUR (MVP) ===\n")
agg_metrics("Uniform prior", sum_no, jaccard_df)
agg_metrics("Evo2 prior",    sum_evo, jaccard_df)

if (nrow(jaccard_df) > 0) {
  cat("\n=== Uniform vs Evo2 CS Comparison ===\n")
  cat(sprintf("  Loci where Evo2 CS is smaller:  %d\n", sum(jaccard_df$cs_size_delta < 0, na.rm = TRUE)))
  cat(sprintf("  Loci where Evo2 CS is same:     %d\n", sum(jaccard_df$cs_size_delta == 0, na.rm = TRUE)))
  cat(sprintf("  Loci where Evo2 CS is larger:   %d\n", sum(jaccard_df$cs_size_delta > 0, na.rm = TRUE)))
  cat(sprintf("  Median Jaccard (uniform vs Evo2): %.3f\n",
              median(jaccard_df$jaccard_uniform_vs_evo2, na.rm = TRUE)))
}

if (nrow(overlap_s11) > 0) {
  cat("\n=== Overlap with Original S11 Mapping ===\n")
  cat(sprintf("  Loci where S11 top var in uniform CS: %d\n",
              sum(overlap_s11$orig_top_in_uniform_cs, na.rm = TRUE)))
  cat(sprintf("  Loci where S11 top var in Evo2 CS:   %d\n",
              sum(overlap_s11$orig_top_in_evo2_cs, na.rm = TRUE)))
  cat(sprintf("  Median Jaccard (orig vs uniform):     %.3f\n",
              median(overlap_s11$jaccard_orig_vs_uniform, na.rm = TRUE)))
  cat(sprintf("  Median Jaccard (orig vs Evo2):        %.3f\n",
              median(overlap_s11$jaccard_orig_vs_evo2, na.rm = TRUE)))
}

cat("\n=== Wilcoxon signed-rank test: total CS size (Evo2 vs uniform), paired by locus ===\n")
if (!is.null(wt_cs)) {
  cat(sprintf("  n loci paired: %d\n", nrow(cs_size_paired)))
  cat(sprintf("  V = %.1f, p = %.4g\n", wt_cs$statistic, wt_cs$p.value))
} else {
  cat("  Skipped: fewer than 2 loci with CS results under both approaches.\n")
}

cat("\n=== Wilcoxon signed-rank test: overlap with original S11 CS (Evo2 vs uniform), paired by locus ===\n")
if (!is.null(wt_jaccard)) {
  cat(sprintf("  n loci paired: %d\n", nrow(jaccard_paired)))
  cat(sprintf("  Loci better overlapped by uniform: %d\n",
              sum(jaccard_paired$jaccard_orig_vs_uniform > jaccard_paired$jaccard_orig_vs_evo2)))
  cat(sprintf("  Loci better overlapped by Evo2:     %d\n",
              sum(jaccard_paired$jaccard_orig_vs_evo2 > jaccard_paired$jaccard_orig_vs_uniform)))
  cat(sprintf("  Loci tied:                          %d\n",
              sum(jaccard_paired$jaccard_orig_vs_uniform == jaccard_paired$jaccard_orig_vs_evo2)))
  cat(sprintf("  V = %.1f, p = %.4g\n", wt_jaccard$statistic, wt_jaccard$p.value))
} else {
  cat("  Skipped: fewer than 2 loci with Jaccard overlap under both approaches.\n")
}

cat("\n=== Capture of original high-PIP (>0.95) variants ===\n")
if (nrow(highpip_capture) > 0) {
  cat(sprintf("  Total high-PIP original signals: %d\n", nrow(highpip_capture)))
  cat(sprintf("  Captured in a credible set, uniform: %d (%.1f%%)\n",
              sum(highpip_capture$uniform_in_cs), 100 * mean(highpip_capture$uniform_in_cs)))
  cat(sprintf("  Captured in a credible set, Evo2:    %d (%.1f%%)\n",
              sum(highpip_capture$evo2_in_cs), 100 * mean(highpip_capture$evo2_in_cs)))
  cat(sprintf("  Median PIP assigned, uniform (0 if not found): %.3f\n",
              median(highpip_capture$uniform_pip_filled)))
  cat(sprintf("  Median PIP assigned, Evo2 (0 if not found):    %.3f\n",
              median(highpip_capture$evo2_pip_filled)))
  if (!is.null(wt_pip)) {
    cat(sprintf("  Wilcoxon signed-rank (Evo2 vs uniform PIP), paired by variant: V = %.1f, p = %.4g\n",
                wt_pip$statistic, wt_pip$p.value))
  } else {
    cat("  Wilcoxon test skipped: fewer than 2 high-PIP variants.\n")
  }
} else {
  cat("  No original signals exceeded PIP > 0.95.\n")
}
sink()
cat(sprintf("Written: %s/aggregate_metrics.txt\n", OUT_DIR))

# ---------------------------------------------------------------------------
# 9) Plots
# ---------------------------------------------------------------------------

# Display labels for the two approaches, used on plots only -- the underlying
# "approach" values ("uniform" / "evo2_prior") are left untouched everywhere
# else in the script (joins, grouping, file outputs).
approach_display_levels <- c("Uniform Prior", "Evo2 Prior")
recode_approach_display <- function(x) {
  factor(recode(x, "uniform" = "Uniform Prior", "evo2_prior" = "Evo2 Prior"),
         levels = approach_display_levels)
}

plot_data <- bind_rows(
  sum_no  %>% select(locus_id, approach, total_cs_size, top_pip),
  sum_evo %>% select(locus_id, approach, total_cs_size, top_pip)
) %>%
  filter(!is.na(total_cs_size)) %>%
  mutate(approach_label = recode_approach_display(approach))

# CS size violin
# Paired Wilcoxon signed-rank test (wt_cs, computed in section 6b) compares
# total CS size within-locus between the uniform- and Evo2-prior analyses.
# Drawn as a manual comparison bracket (two ticks + a bar + a label) using
# base ggplot2's annotate() rather than the ggsignif package, to avoid an
# extra dependency.
cs_signif_label <- if (!is.null(wt_cs)) {
  sprintf("p = %s", format.pval(wt_cs$p.value, digits = 3, eps = 1e-4))
} else {
  "p = NA"
}

cs_y_max  <- if (nrow(plot_data) > 0) max(plot_data$total_cs_size, na.rm = TRUE) else 1
bracket_y <- cs_y_max * 1.08
tip_len   <- cs_y_max * 0.03
label_y   <- bracket_y + cs_y_max * 0.05

p1 <- ggplot(plot_data, aes(x = approach_label, y = total_cs_size, fill = approach_label)) +
  geom_violin(alpha = 0.6, trim = FALSE) +
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.8) +
  annotate("segment", x = 1, xend = 1, y = bracket_y - tip_len, yend = bracket_y) +
  annotate("segment", x = 2, xend = 2, y = bracket_y - tip_len, yend = bracket_y) +
  annotate("segment", x = 1, xend = 2, y = bracket_y, yend = bracket_y) +
  annotate("text", x = 1.5, y = label_y, label = cs_signif_label, vjust = 0, size = 4) +
  scale_fill_manual(values = c("Uniform Prior" = "#4393C3", "Evo2 Prior" = "#D6604D")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
  labs(x = "Approach", y = "Total CS size (variants)") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")
ggsave(paste0(OUT_DIR, "/plots/cs_size_violin.pdf"), p1, width = 5, height = 5)

# PIP scatter: uniform vs Evo2
if (nrow(per_locus) > 0 && all(c("top_pip_uniform", "top_pip_evo2") %in% colnames(per_locus))) {
  p2 <- ggplot(per_locus %>% filter(!is.na(top_pip_uniform), !is.na(top_pip_evo2)),
               aes(x = top_pip_uniform, y = top_pip_evo2)) +
    geom_point(alpha = 0.6, size = 1.5) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
    labs(x = "Top PIP (Uniform Prior)", y = "Top PIP (Evo2 Prior)") +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    theme_bw(base_size = 12)
  ggsave(paste0(OUT_DIR, "/plots/pip_scatter.pdf"), p2, width = 5, height = 5)
}

# Jaccard histogram
if (nrow(jaccard_df) > 0) {
  p3 <- ggplot(jaccard_df %>% filter(!is.na(jaccard_uniform_vs_evo2)),
               aes(x = jaccard_uniform_vs_evo2)) +
    geom_histogram(bins = 20, fill = "#4393C3", colour = "white") +
    labs(x = "Jaccard index", y = "Locus count") +
    theme_bw(base_size = 12)
  ggsave(paste0(OUT_DIR, "/plots/jaccard_histogram.pdf"), p3, width = 5, height = 4)
}

# Venn diagram: overlap between original S11 CS variants and each new
# approach's CS variants, pooled across all loci. Uses pos_id (chr:pos:ref:alt)
# throughout since that's the ID space shared with S11 -- pooling is safe
# because loci are non-overlapping genomic windows, so a given pos_id cannot
# appear in more than one locus.
orig_all_variants <- unique(unlist(s11_t2d$orig_cs_variants))
uniform_all_variants <- res_no  %>% filter(CS_ID > 0) %>% pull(pos_id) %>% na.omit() %>% unique()
evo2_all_variants    <- res_evo %>% filter(CS_ID > 0) %>% pull(pos_id) %>% na.omit() %>% unique()

if (length(orig_all_variants) > 0 && length(uniform_all_variants) > 0 && length(evo2_all_variants) > 0) {
  venn_sets <- list(
    "Original (S11)" = orig_all_variants,
    "Uniform Prior"  = uniform_all_variants,
    "Evo2 Prior"     = evo2_all_variants
  )
  p4 <- ggvenn(venn_sets,
               fill_color = c("#66C2A5", "#4393C3", "#D6604D"),
               stroke_size = 0.6,
               set_name_size = 4,
               text_size = 4) +
    theme(legend.position = "none")
  ggsave(paste0(OUT_DIR, "/plots/venn_overlap.pdf"), p4, width = 6, height = 6)
} else {
  warning("Skipped Venn diagram: at least one of the three variant sets (S11, uniform, Evo2) is empty.")
}

cat("\nAll outputs written to:", OUT_DIR, "\n")

# ---------------------------------------------------------------------------
# 10) Per-locus 3-panel plots: GWAS -log10(p), uniform-prior PIP, Evo2-prior PIP
# ---------------------------------------------------------------------------
# One "locuszoom-style" figure per locus, three panels stacked on a shared
# genomic-position x-axis: GWAS -log10(p) on top, uniform-prior PIP in the
# middle, Evo2-prior PIP on the bottom. Base points are light grey; a
# variant's point gets a colored FILL if it's a uniform-prior CS member, and
# a colored RING (border) if it's an Evo2-prior CS member. Points in neither
# CS render as plain light grey (fill and border both grey, so no
# distinguishable ring). The two encodings use two different hue families
# (shades of red for uniform-prior CS fills, shades of blue for Evo2-prior CS
# rings) so it's easy to tell at a glance which analysis a highlight belongs
# to; different shades within each family distinguish different CS indices.
# Palettes are re-built per locus since CS numbering is locus-local, not
# shared genome-wide.

red_palette  <- function(n) colorRampPalette(c("#FFB3B3", "#800000"))(n)
blue_palette <- function(n) colorRampPalette(c("#A6D0FF", "#00204D"))(n)

# Mechanical label rule per your spec: title case, underscores -> spaces.
# Note this will render acronyms as "Pval"/"Pip" rather than "PVAL"/"PIP" --
# let me know if you'd rather I hand-curate those instead.
axis_label <- function(x) str_to_title(gsub("_", " ", x))

locuszoom_dir <- paste0(OUT_DIR, "/plots/locuszoom")
dir.create(locuszoom_dir, showWarnings = FALSE, recursive = TRUE)

build_panel <- function(data, yvar, x_range, fill_values, colour_values,
                         is_bottom, x_title = NULL) {
  ggplot(data, aes(x = POS_num, y = .data[[yvar]])) +
    geom_point(aes(fill = fill_grp, colour = ring_grp),
               shape = 21, size = 2.8, stroke = 1.1) +
    scale_fill_manual(values = fill_values, guide = "none") +
    scale_colour_manual(values = colour_values, guide = "none") +
    scale_x_continuous(limits = x_range, labels = scales::comma) +
    labs(y = axis_label(yvar), x = x_title) +
    theme_bw(base_size = 18) +
    theme(
      axis.title.y     = element_text(size = 18, face = "bold"),
      axis.text.y      = element_text(size = 14),
      axis.title.x     = element_text(size = 18, face = "bold"),
      axis.text.x      = if (is_bottom) element_text(size = 14) else element_blank(),
      axis.ticks.x     = if (is_bottom) element_line() else element_blank(),
      panel.grid.minor = element_blank()
    )
}

cat(sprintf("Generating %d per-locus 3-panel plots ...\n", length(loci_with_both)))

n_locuszoom_written <- 0
for (lid in loci_with_both) {

  vm <- variant_map %>% filter(locus_id == lid)
  if (nrow(vm) == 0) {
    warning(sprintf("Skipping locuszoom plot for %s: no variant list data found.", lid))
    next
  }

  locus_df <- vm %>%
    mutate(POS_num  = suppressWarnings(as.numeric(POS)),
           PVAL_num = suppressWarnings(as.numeric(PVAL))) %>%
    left_join(res_no  %>% select(locus_id, SNP_ID, uniform_pip = PIP, uniform_cs = CS_ID),
              by = c("locus_id", "SNP_ID")) %>%
    left_join(res_evo %>% select(locus_id, SNP_ID, evo2_pip = PIP, evo2_cs = CS_ID),
              by = c("locus_id", "SNP_ID")) %>%
    mutate(
      # A GWAS p-value of exactly 0 reflects numerical underflow (MVP is a
      # very large cohort), not a true p = 0 -- capped at -log10(p) = 320
      # rather than plotted as Inf.
      neg_log10_pval = ifelse(is.na(PVAL_num), NA_real_,
                               ifelse(PVAL_num <= 0, 320, -log10(PVAL_num))),
      fill_grp = ifelse(!is.na(uniform_cs) & uniform_cs > 0, as.character(uniform_cs), "None"),
      ring_grp = ifelse(!is.na(evo2_cs)    & evo2_cs    > 0, as.character(evo2_cs),    "None")
    )

  # Draw highlighted (in-CS) points last so they render on top of the grey background
  locus_df <- locus_df %>%
    mutate(is_highlighted = fill_grp != "None" | ring_grp != "None") %>%
    arrange(is_highlighted)

  fill_levels <- sort(setdiff(unique(locus_df$fill_grp), "None"))
  ring_levels <- sort(setdiff(unique(locus_df$ring_grp), "None"))

  pal_fill <- if (length(fill_levels) > 0) setNames(red_palette(length(fill_levels)),  fill_levels) else character(0)
  pal_ring <- if (length(ring_levels) > 0) setNames(blue_palette(length(ring_levels)), ring_levels) else character(0)

  fill_values <- c(pal_fill, "None" = "grey80")
  # NOTE: "None" must NOT map to NA here. geom_point() treats colour as a
  # required aesthetic even with shape 21 (fill + border), so an NA colour
  # value causes ggplot to silently drop that row entirely via
  # remove_missing() -- not just render an invisible border. Using the same
  # grey as the fill keeps the row visible while showing no distinguishable
  # ring (fill and border are indistinguishable in that colour).
  colour_values <- c(pal_ring, "None" = "grey80")

  x_range <- range(locus_df$POS_num, na.rm = TRUE)
  chr_num <- unique(locus_df$CHR)[1]

  p_top <- build_panel(locus_df, "neg_log10_pval", x_range, fill_values, colour_values,
                        is_bottom = FALSE)
  p_mid <- build_panel(locus_df, "uniform_pip", x_range, fill_values, colour_values,
                        is_bottom = FALSE)
  p_bot <- build_panel(locus_df, "evo2_pip", x_range, fill_values, colour_values,
                        is_bottom = TRUE, x_title = paste0("Chromosome ", chr_num))

  combined <- p_top / p_mid / p_bot

  ggsave(paste0(locuszoom_dir, "/", lid, "_locuszoom.pdf"), combined, width = 8, height = 10)
  n_locuszoom_written <- n_locuszoom_written + 1
}

cat(sprintf("Written %d locuszoom-style plots to %s\n", n_locuszoom_written, locuszoom_dir))