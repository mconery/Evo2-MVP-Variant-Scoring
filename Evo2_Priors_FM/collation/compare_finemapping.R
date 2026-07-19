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
#
# Outputs (written to $BASE/collation/):
#   per_locus_summary.tsv       -- one row per locus
#   cs_comparison.tsv           -- one row per locus × CS signal
#   aggregate_metrics.txt       -- printed summary statistics
#   plots/cs_size_violin.pdf
#   plots/pip_scatter.pdf
#   plots/jaccard_histogram.pdf
################################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(readr)
  library(purrr)
})

BASE       <- "/grand/GeomicVar/mconery/evo2_variant_scoring_mapping"
S11_FILE   <- paste0(BASE, "/mapping_loci/Supplementary_Table-S11.txt")
NO_PRIOR   <- paste0(BASE, "/carma_results/without_priors")
WITH_PRIOR <- paste0(BASE, "/carma_results/with_priors")
LOCI_FILE  <- paste0(BASE, "/loci_definition/t2d_eur_loci.tsv")
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

# Parse EUR.variant_ids: semicolon-separated list of variants in the original CS
s11_t2d <- s11_t2d %>%
  mutate(orig_cs_variants = str_split(EUR.variant_ids, ";"))

cat(sprintf("S11: %d T2D EUR signals across %d unique loci\n",
            nrow(s11_t2d), n_distinct(s11_t2d$locus_id)))

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

jaccard <- function(a, b) {
  if (length(a) == 0 && length(b) == 0) return(NA_real_)
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

  cs_no  <- res_no  %>% filter(locus_id == lid, CS_ID > 0) %>% pull(SNP_ID)
  cs_evo <- res_evo %>% filter(locus_id == lid, CS_ID > 0) %>% pull(SNP_ID)

  tibble(
    locus_id                     = lid,
    jaccard_orig_vs_uniform      = jaccard(orig_variants, cs_no),
    jaccard_orig_vs_evo2         = jaccard(orig_variants, cs_evo),
    n_orig_variants              = length(orig_variants),
    orig_top_in_uniform_cs       = any(s11_t2d$EUR.best_variant[s11_t2d$locus_id == lid] %in% cs_no),
    orig_top_in_evo2_cs          = any(s11_t2d$EUR.best_variant[s11_t2d$locus_id == lid] %in% cs_evo)
  )
})

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
# 7) Per-CS comparison table
# ---------------------------------------------------------------------------

get_cs_rows <- function(df, approach_label) {
  df %>%
    filter(CS_ID > 0) %>%
    group_by(locus_id, CS_ID) %>%
    summarise(
      approach    = approach_label,
      cs_size     = n(),
      top_pip     = max(PIP, na.rm = TRUE),
      top_variant = SNP_ID[which.max(PIP)],
      variants    = paste(SNP_ID, collapse = ";"),
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
sink()
cat(sprintf("Written: %s/aggregate_metrics.txt\n", OUT_DIR))

# ---------------------------------------------------------------------------
# 9) Plots
# ---------------------------------------------------------------------------

plot_data <- bind_rows(
  sum_no  %>% select(locus_id, approach, total_cs_size, top_pip),
  sum_evo %>% select(locus_id, approach, total_cs_size, top_pip)
) %>% filter(!is.na(total_cs_size))

# CS size violin
p1 <- ggplot(plot_data, aes(x = approach, y = total_cs_size, fill = approach)) +
  geom_violin(alpha = 0.6, trim = FALSE) +
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.8) +
  scale_fill_manual(values = c("uniform" = "#4393C3", "evo2_prior" = "#D6604D")) +
  labs(title = "Total credible set size per locus",
       x = "Approach", y = "Total CS size (variants)") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")
ggsave(paste0(OUT_DIR, "/plots/cs_size_violin.pdf"), p1, width = 5, height = 5)

# PIP scatter: uniform vs Evo2
if (nrow(per_locus) > 0 && all(c("top_pip_uniform", "top_pip_evo2") %in% colnames(per_locus))) {
  p2 <- ggplot(per_locus %>% filter(!is.na(top_pip_uniform), !is.na(top_pip_evo2)),
               aes(x = top_pip_uniform, y = top_pip_evo2)) +
    geom_point(alpha = 0.6, size = 1.5) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
    labs(title = "Top-variant PIP: uniform vs Evo2 prior",
         x = "Top PIP (uniform prior)", y = "Top PIP (Evo2 prior)") +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    theme_bw(base_size = 12)
  ggsave(paste0(OUT_DIR, "/plots/pip_scatter.pdf"), p2, width = 5, height = 5)
}

# Jaccard histogram
if (nrow(jaccard_df) > 0) {
  p3 <- ggplot(jaccard_df %>% filter(!is.na(jaccard_uniform_vs_evo2)),
               aes(x = jaccard_uniform_vs_evo2)) +
    geom_histogram(bins = 20, fill = "#4393C3", colour = "white") +
    labs(title = "Jaccard similarity: uniform CS vs Evo2-prior CS",
         x = "Jaccard index", y = "Locus count") +
    theme_bw(base_size = 12)
  ggsave(paste0(OUT_DIR, "/plots/jaccard_histogram.pdf"), p3, width = 5, height = 4)
}

cat("\nAll outputs written to:", OUT_DIR, "\n")
