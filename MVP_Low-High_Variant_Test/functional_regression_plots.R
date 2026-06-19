# ============================================================================
# Logistic Regression: High vs Low PIP ~ Functional Annotations + Evo2 Scores
# ============================================================================
# Central question: Does Evo2's delta score predict High vs Low PIP status
# independently of established functional annotations (CADD, RegulomeDB, etc.)?
#
# Three annotation tiers (parallel to conservation_regression_plots.R):
#   1. All-variant: CADD_phred + RegulomeDB_prob + Evo2 (n~7,000+)
#   2. Near-transcript: DANN/FitCons/GenoCanyon/FATHMM_XF + Evo2 (n~600)
#   3. Missense-only: 6 missense scores + Evo2 (n~450-500)
# Plus categorical annotations (ENCODE cCRE, RegulomeDB rank) and MAF check.
# ============================================================================

library(tidyverse)
library(ggplot2)
library(scales)
library(patchwork)
library(car)

# ============================================================================
# CONFIGURATION
# ============================================================================

func_file  <- "C:/Users/mitch/Documents/Argonne/Variant_Scoring/MVP_functional_annotations.csv"
evo2_file  <- "C:/Users/mitch/Documents/Argonne/Variant_Scoring/scored_MVP_results/MVP_variant_scores.7b_arc_longcontext_model.16384bp_context.csv"
output_dir <- "C:/Users/mitch/Documents/Argonne/Variant_Scoring/scored_MVP_figures"

col_low  <- "#6baed6"
col_high <- "#08519c"

# ============================================================================
# SHARED HELPERS (mirrors conservation_regression_plots.R)
# ============================================================================

theme_mvp <- function() {
  theme_bw() +
    theme(
      plot.title        = element_text(hjust = 0.5, face = "bold", size = 18),
      plot.margin       = margin(5, 20, 5, 5, "mm"),
      axis.text         = element_text(size = 12),
      axis.text.x       = element_text(angle = 45, hjust = 1),
      axis.title        = element_text(size = 13),
      strip.text.x      = element_text(face = "bold", size = 13),
      strip.text.y      = element_text(face = "bold", size = 12, angle = 0),
      legend.position   = "none",
      panel.grid.major  = element_blank(),
      panel.grid.minor  = element_blank()
    )
}

format_pval <- function(p) {
  paste0("p=", formatC(p, format = "e", digits = 1))
}

sig_stars <- function(p) {
  case_when(p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", TRUE ~ "ns")
}

# Violin panel for a single predictor
make_violin_panel <- function(pred_name, pred_label, p_label_str, beta_val, data, beta_prefix = "β") {
  annotation_str <- paste0(beta_prefix, " = ", sprintf("%+.2f", beta_val), "\n", p_label_str)
  ggplot(data, aes(x = class, y = .data[[pred_name]], fill = class)) +
    geom_violin(trim = FALSE) +
    stat_summary(
      fun = median, fun.min = median, fun.max = median,
      geom = "crossbar", width = 0.3, color = "black", linewidth = 0.8
    ) +
    annotate("text",
             x = 1.5, y = Inf,
             label = annotation_str,
             vjust = 1.5, size = 5.5, fontface = "bold",
             lineheight = 0.9) +
    scale_fill_manual(values = c("Low PIP" = col_low, "High PIP" = col_high)) +
    labs(x = NULL, y = pred_label) +
    theme_mvp() +
    theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16))
}

# Forest plot with side text panel (base + patchwork composition)
make_forest_plot <- function(results_df, title_str,
                             col_values  = c("Functional" = "#6baed6", "Evo2" = "#08519c"),
                             col_col     = "var_type",
                             legend_labs = c("Functional annotation", "Evo2 model")) {
  fdf <- results_df %>%
    mutate(label = factor(label, levels = rev(label)))

  plot_base <- ggplot(fdf, aes(x = or, y = label, color = .data[[col_col]])) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.8) +
    geom_errorbar(aes(xmin = ci_lo, xmax = ci_hi),
                  width = 0.2, linewidth = 0.8, orientation = "y") +
    geom_point(size = 3.5) +
    scale_x_log10() +
    scale_color_manual(
      values = col_values,
      name   = NULL,
      labels = legend_labs
    ) +
    labs(
      title = title_str,
      x     = "Odds Ratio (95% CI, standardised predictors)",
      y     = NULL
    ) +
    theme_bw() +
    theme(
      plot.title       = element_text(hjust = 0.5, face = "bold", size = 18),
      plot.margin      = margin(5, 2, 5, 5, "mm"),
      axis.text        = element_text(size = 12),
      axis.text.x      = element_text(angle = 45, hjust = 1),
      axis.text.y      = element_text(size = 12),
      axis.title       = element_text(size = 13),
      legend.position  = "bottom",
      legend.text      = element_text(size = 12),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )

  sdf <- fdf %>%
    mutate(side_label = paste0(
      "β = ", sprintf("%+.2f", beta), "\n",
      p_label, "  ", stars
    ))

  plot_side <- ggplot(sdf, aes(x = 0, y = label)) +
    geom_text(aes(label = side_label),
              hjust = 0, size = 4.5, fontface = "bold", lineheight = 0.9) +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    labs(title = "", x = NULL, y = NULL) +
    theme_void() +
    theme(plot.margin = margin(5, 5, 5, 2, "mm"),
          axis.text = element_blank(), axis.ticks = element_blank())

  plot_base + plot_side + plot_layout(widths = c(3, 1))
}

# ============================================================================
# STEP 1: LOAD AND MERGE DATA
# ============================================================================

message("Loading functional annotations...")
df_func <- read_csv(func_file, show_col_types = FALSE)

message("Loading Evo2 scores (7b-arc-longcontext, 16384 bp)...")
df_evo2 <- read_csv(evo2_file, show_col_types = FALSE) %>%
  select(`MVP ID`, evo2_delta_score)

df_raw <- df_func %>%
  inner_join(df_evo2, by = "MVP ID")

message("Merged rows: ", nrow(df_raw))

# Coding VEP consequences (Ensembl standard)
coding_consequences <- c(
  "missense_variant", "synonymous_variant", "stop_gained", "stop_lost",
  "start_lost", "frameshift_variant", "inframe_insertion", "inframe_deletion",
  "splice_donor_variant", "splice_acceptor_variant"
)

# RegulomeDB ordinal encoding: 1a (strongest evidence) -> 7 (no evidence)
rank_order <- c("1a"=1,"1b"=2,"1c"=3,"1d"=4,"1e"=5,"1f"=6,
                "2a"=7,"2b"=8,"3a"=9,"3b"=10,"4"=11,"5"=12,"6"=13,"7"=14)

df <- df_raw %>%
  rename(EAF = `EAF Population`, VEP = `VEP Annotation`) %>%
  mutate(
    class        = str_replace(Set, "_", " "),
    pip_high     = as.integer(Set == "High_PIP"),
    maf          = pmin(EAF, 1 - EAF),
    is_missense  = grepl("missense_variant", VEP),
    variant_type = if_else(VEP %in% coding_consequences, "Coding", "Non-coding"),
    rdb_rank_num = rank_order[RegulomeDB_rank],
    has_ccre     = as.integer(!is.na(ENCODE_cCRE))
  )

message("Total variants: ", nrow(df))
message("Class: High PIP=", sum(df$pip_high), ", Low PIP=", sum(1 - df$pip_high))
message("Evo2 coverage: ", sum(!is.na(df$evo2_delta_score)), "/", nrow(df))
message("Missense variants: ", sum(df$is_missense))

# ============================================================================
# STEP 2: ALL-VARIANT LOGISTIC REGRESSIONS (CADD + RegulomeDB + Evo2)
# ============================================================================

df_allvar <- df %>% drop_na(CADD_phred, RegulomeDB_prob, evo2_delta_score)
message("\nAll-variant complete cases (CADD_phred + RegulomeDB_prob + Evo2): ", nrow(df_allvar))

allvar_preds <- list(
  CADD_phred       = list(col = "CADD_phred",       label = "CADD Phred Score"),
  RegulomeDB_prob  = list(col = "RegulomeDB_prob",   label = "RegulomeDB Probability Score"),
  evo2_delta_score = list(col = "evo2_delta_score",  label = "Evo2 Delta Score\n(7b-arc-longcontext, 16384 bp)")
)

message("\n--- All-variant individual regressions ---")

allvar_individual <- map_dfr(names(allvar_preds), function(name) {
  meta      <- allvar_preds[[name]]
  col       <- meta$col
  label     <- meta$label
  use_scale <- col == "evo2_delta_score"
  fml_str   <- if (use_scale) paste0("pip_high ~ scale(`", col, "`)") else paste0("pip_high ~ `", col, "`")
  fit  <- glm(as.formula(fml_str), data = df_allvar, family = binomial)
  beta <- coef(fit)[2]
  or   <- exp(beta)
  ci   <- exp(confint.default(fit)[2, ])
  pval <- summary(fit)$coefficients[2, "Pr(>|z|)"]
  message(sprintf("  %-25s  beta=%+.3f  OR=%.3f [%.3f, %.3f]  %s  %s",
                  name, beta, or, ci[1], ci[2], format_pval(pval), sig_stars(pval)))
  tibble(predictor = name, label = label, beta = beta, or = or,
         ci_lo = ci[1], ci_hi = ci[2], p_value = pval,
         p_label = format_pval(pval), stars = sig_stars(pval), model = "Individual",
         var_type = if_else(name == "evo2_delta_score", "Evo2", "Functional"),
         beta_prefix = if_else(use_scale, "β/SD", "β"))
})

message("\n--- All-variant joint regression (standardised predictors) ---")

fit_joint_allvar <- glm(
  pip_high ~ scale(CADD_phred) + scale(RegulomeDB_prob) + scale(evo2_delta_score),
  data = df_allvar, family = binomial
)

coef_joint_allvar <- summary(fit_joint_allvar)$coefficients[-1, , drop = FALSE]
ci_joint_allvar   <- confint.default(fit_joint_allvar)[-1, , drop = FALSE]

allvar_joint <- tibble(
  predictor = names(allvar_preds),
  label     = map_chr(allvar_preds, "label"),
  beta      = coef_joint_allvar[, "Estimate"],
  or        = exp(coef_joint_allvar[, "Estimate"]),
  ci_lo     = exp(ci_joint_allvar[, 1]),
  ci_hi     = exp(ci_joint_allvar[, 2]),
  p_value   = coef_joint_allvar[, "Pr(>|z|)"],
  p_label   = format_pval(p_value),
  stars     = sig_stars(p_value),
  model     = "Joint",
  var_type  = if_else(predictor == "evo2_delta_score", "Evo2", "Functional")
)

walk(seq_len(nrow(allvar_joint)), function(i) {
  r <- allvar_joint[i, ]
  message(sprintf("  %-25s  OR=%.3f [%.3f, %.3f]  %s  %s",
                  r$predictor, r$or, r$ci_lo, r$ci_hi, r$p_label, r$stars))
})

message("\n--- VIF (all-variant joint model) ---")
vif_allvar <- car::vif(fit_joint_allvar)
vif_allvar_df <- tibble(
  predictor = names(allvar_preds),
  label     = map_chr(allvar_preds, "label"),
  vif       = as.numeric(vif_allvar),
  var_type  = if_else(predictor == "evo2_delta_score", "Evo2", "Functional")
) %>%
  mutate(
    disp_label = str_remove(label, "\n.*"),
    disp_label = factor(disp_label, levels = disp_label)
  )

walk(seq_len(nrow(vif_allvar_df)), function(i) {
  message(sprintf("  %-25s  VIF = %.3f", vif_allvar_df$predictor[i], vif_allvar_df$vif[i]))
})

message("\n=== All-Variant Summary Table ===")
bind_rows(allvar_individual, allvar_joint) %>%
  select(model, predictor, or, ci_lo, ci_hi, p_value, stars) %>%
  mutate(across(c(or, ci_lo, ci_hi), ~round(.x, 3))) %>%
  print(n = Inf)

# --- Build all-variant plots ---
message("\nBuilding all-variant violin panels...")

violin_allvar <- pmap(
  list(pred_name   = allvar_individual$predictor,
       pred_label  = allvar_individual$label,
       p_label_str = allvar_individual$p_label,
       beta_val    = allvar_individual$beta,
       beta_prefix = allvar_individual$beta_prefix),
  make_violin_panel, data = df_allvar
)

plot_violins_allvar <- wrap_plots(violin_allvar, ncol = 3) +
  plot_annotation(
    title = "High vs Low PIP: Functional Annotation & Evo2 Score Distributions\n(β and p-values from individual logistic regressions)",
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18))
  )

ggsave(file.path(output_dir, "functional_individual_regressions_allvar.png"),
       plot_violins_allvar, width = 15, height = 5, dpi = 300)
message("Saved: functional_individual_regressions_allvar.png")

message("Building all-variant joint forest plot...")

plot_forest_allvar <- make_forest_plot(
  allvar_joint %>% mutate(label = str_remove(label, "\n.*")),
  "Joint Logistic Regression: High vs Low PIP\n(CADD + RegulomeDB + Evo2, Standardised ORs, 95% CI)"
)

ggsave(file.path(output_dir, "functional_joint_regression_forest_allvar.png"),
       plot_forest_allvar, width = 12, height = 5, dpi = 300)
message("Saved: functional_joint_regression_forest_allvar.png")

plot_vif_allvar <- ggplot(vif_allvar_df, aes(x = disp_label, y = vif, fill = var_type)) +
  geom_col(width = 0.6, color = "white") +
  geom_hline(yintercept = 5,  linetype = "dashed", color = "orange",    linewidth = 0.8) +
  geom_hline(yintercept = 10, linetype = "dashed", color = "firebrick", linewidth = 0.8) +
  geom_text(aes(label = sprintf("%.3f", vif)), vjust = -0.4, size = 4, fontface = "bold") +
  annotate("text", x = Inf, y = 5,  label = "VIF = 5",  hjust = 1.1, vjust = -0.4,
           color = "orange",    size = 3.5) +
  annotate("text", x = Inf, y = 10, label = "VIF = 10", hjust = 1.1, vjust = -0.4,
           color = "firebrick", size = 3.5) +
  scale_fill_manual(values = c("Functional" = col_low, "Evo2" = col_high), name = NULL,
                    labels = c("Functional annotation", "Evo2 model")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(title = "Variance Inflation Factors — All-Variant Joint Model",
       x = NULL, y = "VIF") +
  theme_mvp() +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 0, hjust = 0.5))

ggsave(file.path(output_dir, "functional_vif_allvar.png"),
       plot_vif_allvar, width = 8, height = 5, dpi = 300)
message("Saved: functional_vif_allvar.png")

# ============================================================================
# STEP 3: NEAR-TRANSCRIPT INDIVIDUAL REGRESSIONS (Evo2 included)
# ============================================================================

df_nt <- df %>%
  drop_na(DANN_score, FitCons_score, GenoCanyon_score, FATHMM_XF, evo2_delta_score)

message("\nNear-transcript complete cases: ", nrow(df_nt),
        " (dbNSFP near-transcript only; non-random missingness)")

nt_preds <- list(
  DANN_score       = list(col = "DANN_score",       label = "DANN Score"),
  FitCons_score    = list(col = "FitCons_score",     label = "FitCons Score"),
  GenoCanyon_score = list(col = "GenoCanyon_score",  label = "GenoCanyon Score"),
  FATHMM_XF        = list(col = "FATHMM_XF",         label = "FATHMM-XF Score"),
  evo2_delta_score = list(col = "evo2_delta_score",  label = "Evo2 Delta Score\n(7b-arc-longcontext, 16384 bp)")
)

message("\n--- Near-transcript individual regressions ---")

nt_individual <- map_dfr(names(nt_preds), function(name) {
  meta      <- nt_preds[[name]]
  col       <- meta$col
  label     <- meta$label
  use_scale <- col == "evo2_delta_score"
  fml_str   <- if (use_scale) paste0("pip_high ~ scale(`", col, "`)") else paste0("pip_high ~ `", col, "`")
  fit  <- glm(as.formula(fml_str), data = df_nt, family = binomial)
  beta <- coef(fit)[2]
  or   <- exp(beta)
  ci   <- exp(confint.default(fit)[2, ])
  pval <- summary(fit)$coefficients[2, "Pr(>|z|)"]
  message(sprintf("  %-25s  beta=%+.3f  OR=%.3f [%.3f, %.3f]  %s  %s",
                  name, beta, or, ci[1], ci[2], format_pval(pval), sig_stars(pval)))
  tibble(predictor = name, label = label, beta = beta, or = or,
         ci_lo = ci[1], ci_hi = ci[2], p_value = pval,
         p_label = format_pval(pval), stars = sig_stars(pval), model = "Individual",
         var_type = if_else(name == "evo2_delta_score", "Evo2", "Functional"),
         beta_prefix = if_else(use_scale, "β/SD", "β"))
})

message("Building near-transcript violin panels...")

violin_nt <- pmap(
  list(pred_name   = nt_individual$predictor,
       pred_label  = nt_individual$label,
       p_label_str = nt_individual$p_label,
       beta_val    = nt_individual$beta,
       beta_prefix = nt_individual$beta_prefix),
  make_violin_panel, data = df_nt
)

# 5 panels in ncol=3: row 1 has DANN, FitCons, GenoCanyon; row 2 has FATHMM-XF, Evo2
plot_violins_nt <- wrap_plots(violin_nt, ncol = 3) +
  plot_annotation(
    title = "Near-Transcript Variants: Functional & Evo2 Score Distributions\n(β and p-values from individual logistic regressions; no joint model due to non-random missingness)",
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))
  )

# Individual forest plot (no joint model for near-transcript)
plot_nt_forest <- make_forest_plot(
  nt_individual %>% mutate(label = str_remove(label, "\n.*")),
  "Near-Transcript Individual Regressions: High vs Low PIP\n(Raw ORs, 95% CI — individual models only; n≈600)"
)

plot_nt_combined <- plot_violins_nt / plot_nt_forest +
  plot_layout(heights = c(2, 1))

ggsave(file.path(output_dir, "functional_neartranscript_regressions.png"),
       plot_nt_combined, width = 15, height = 14, dpi = 300)
message("Saved: functional_neartranscript_regressions.png")

# ============================================================================
# STEP 4: MISSENSE-ONLY ANALYSIS (all 6 scores + Evo2)
# ============================================================================

df_missense <- df %>%
  filter(is_missense) %>%
  drop_na(SIFT_score, PolyPhen2_HVAR, REVEL_score, AlphaMissense,
          BayesDel_score, PrimateAI_score, evo2_delta_score)

message("\nMissense complete cases: ", nrow(df_missense))
message("Missense class: High PIP=", sum(df_missense$pip_high),
        ", Low PIP=", sum(1 - df_missense$pip_high))

# ============================================================================
# REGRESSION TIER COVERAGE SUMMARY
# ============================================================================

message("\n=== Regression Tier Coverage Summary ===")
coverage_tbl <- tibble(
  Tier        = c("All-variant", "Near-transcript", "Missense"),
  Predictors  = c(
    "CADD_phred, RegulomeDB_prob, Evo2",
    "DANN, FitCons, GenoCanyon, FATHMM-XF, Evo2",
    "SIFT, PolyPhen2, REVEL, AlphaMissense, BayesDel, PrimateAI, Evo2"
  ),
  N_total    = c(nrow(df_allvar),            nrow(df_nt),            nrow(df_missense)),
  N_high_pip = c(sum(df_allvar$pip_high),    sum(df_nt$pip_high),    sum(df_missense$pip_high)),
  N_low_pip  = c(sum(1-df_allvar$pip_high),  sum(1-df_nt$pip_high),  sum(1-df_missense$pip_high)),
  Pct_of_all = round(c(nrow(df_allvar), nrow(df_nt), nrow(df_missense)) / nrow(df) * 100, 1)
)
print(coverage_tbl, width = Inf)

missense_preds <- list(
  SIFT_score       = list(col = "SIFT_score",       label = "SIFT Score\n(lower = more damaging)"),
  PolyPhen2_HVAR   = list(col = "PolyPhen2_HVAR",   label = "PolyPhen2 HVAR Score"),
  REVEL_score      = list(col = "REVEL_score",       label = "REVEL Score"),
  AlphaMissense    = list(col = "AlphaMissense",     label = "AlphaMissense Score"),
  BayesDel_score   = list(col = "BayesDel_score",    label = "BayesDel Score"),
  PrimateAI_score  = list(col = "PrimateAI_score",   label = "PrimateAI Score"),
  evo2_delta_score = list(col = "evo2_delta_score",  label = "Evo2 Delta Score\n(7b-arc-longcontext, 16384 bp)")
)

message("\n--- Missense individual regressions ---")

missense_individual <- map_dfr(names(missense_preds), function(name) {
  meta      <- missense_preds[[name]]
  col       <- meta$col
  label     <- meta$label
  use_scale <- col == "evo2_delta_score"
  fml_str   <- if (use_scale) paste0("pip_high ~ scale(`", col, "`)") else paste0("pip_high ~ `", col, "`")
  fit  <- glm(as.formula(fml_str), data = df_missense, family = binomial)
  beta <- coef(fit)[2]
  or   <- exp(beta)
  ci   <- exp(confint.default(fit)[2, ])
  pval <- summary(fit)$coefficients[2, "Pr(>|z|)"]
  message(sprintf("  %-25s  beta=%+.3f  OR=%.3f [%.3f, %.3f]  %s  %s",
                  name, beta, or, ci[1], ci[2], format_pval(pval), sig_stars(pval)))
  tibble(predictor = name, label = label, beta = beta, or = or,
         ci_lo = ci[1], ci_hi = ci[2], p_value = pval,
         p_label = format_pval(pval), stars = sig_stars(pval), model = "Individual",
         var_type = if_else(name == "evo2_delta_score", "Evo2", "Functional"),
         beta_prefix = if_else(use_scale, "β/SD", "β"))
})

message("\n--- Missense joint regression (standardised predictors) ---")

fit_joint_missense <- glm(
  pip_high ~ scale(SIFT_score) + scale(PolyPhen2_HVAR) + scale(REVEL_score) +
             scale(AlphaMissense) + scale(BayesDel_score) + scale(PrimateAI_score) +
             scale(evo2_delta_score),
  data = df_missense, family = binomial
)

coef_joint_missense <- summary(fit_joint_missense)$coefficients[-1, , drop = FALSE]
ci_joint_missense   <- confint.default(fit_joint_missense)[-1, , drop = FALSE]

missense_joint <- tibble(
  predictor = names(missense_preds),
  label     = map_chr(missense_preds, "label"),
  beta      = coef_joint_missense[, "Estimate"],
  or        = exp(coef_joint_missense[, "Estimate"]),
  ci_lo     = exp(ci_joint_missense[, 1]),
  ci_hi     = exp(ci_joint_missense[, 2]),
  p_value   = coef_joint_missense[, "Pr(>|z|)"],
  p_label   = format_pval(p_value),
  stars     = sig_stars(p_value),
  model     = "Joint",
  var_type  = if_else(predictor == "evo2_delta_score", "Evo2", "Functional")
)

walk(seq_len(nrow(missense_joint)), function(i) {
  r <- missense_joint[i, ]
  message(sprintf("  %-25s  OR=%.3f [%.3f, %.3f]  %s  %s",
                  r$predictor, r$or, r$ci_lo, r$ci_hi, r$p_label, r$stars))
})

message("\n--- VIF (missense joint model) ---")
vif_missense <- car::vif(fit_joint_missense)
vif_missense_df <- tibble(
  predictor  = names(missense_preds),
  label      = map_chr(missense_preds, "label"),
  vif        = as.numeric(vif_missense),
  var_type   = if_else(predictor == "evo2_delta_score", "Evo2", "Functional")
) %>%
  mutate(
    disp_label = str_remove(label, "\n.*"),
    disp_label = factor(disp_label, levels = disp_label)
  )

walk(seq_len(nrow(vif_missense_df)), function(i) {
  message(sprintf("  %-25s  VIF = %.3f", vif_missense_df$predictor[i], vif_missense_df$vif[i]))
})

message("\n=== Missense Summary Table ===")
bind_rows(missense_individual, missense_joint) %>%
  select(model, predictor, or, ci_lo, ci_hi, p_value, stars) %>%
  mutate(across(c(or, ci_lo, ci_hi), ~round(.x, 3))) %>%
  print(n = Inf)

# --- Build missense plots ---
message("Building missense violin panels...")

violin_missense <- pmap(
  list(pred_name   = missense_individual$predictor,
       pred_label  = missense_individual$label,
       p_label_str = missense_individual$p_label,
       beta_val    = missense_individual$beta,
       beta_prefix = missense_individual$beta_prefix),
  make_violin_panel, data = df_missense
)

# 7 panels in ncol=4: row 1 = 4, row 2 = 3
plot_violins_missense <- wrap_plots(violin_missense, ncol = 4) +
  plot_annotation(
    title = "Missense Variants: Functional & Evo2 Score Distributions\n(β and p-values from individual logistic regressions)",
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18))
  )

ggsave(file.path(output_dir, "functional_missense_violins.png"),
       plot_violins_missense, width = 16, height = 10, dpi = 300)
message("Saved: functional_missense_violins.png")

message("Building missense joint forest plot...")

plot_forest_missense <- make_forest_plot(
  missense_joint %>% mutate(label = str_remove(label, "\n.*")),
  "Missense Joint Logistic Regression: High vs Low PIP\n(6 Missense Scores + Evo2, Standardised ORs, 95% CI)"
)

ggsave(file.path(output_dir, "functional_missense_forest.png"),
       plot_forest_missense, width = 14, height = 5, dpi = 300)
message("Saved: functional_missense_forest.png")

plot_vif_missense <- ggplot(vif_missense_df, aes(x = disp_label, y = vif, fill = var_type)) +
  geom_col(width = 0.6, color = "white") +
  geom_hline(yintercept = 5,  linetype = "dashed", color = "orange",    linewidth = 0.8) +
  geom_hline(yintercept = 10, linetype = "dashed", color = "firebrick", linewidth = 0.8) +
  geom_text(aes(label = sprintf("%.3f", vif)), vjust = -0.4, size = 4, fontface = "bold") +
  annotate("text", x = Inf, y = 5,  label = "VIF = 5",  hjust = 1.1, vjust = -0.4,
           color = "orange",    size = 3.5) +
  annotate("text", x = Inf, y = 10, label = "VIF = 10", hjust = 1.1, vjust = -0.4,
           color = "firebrick", size = 3.5) +
  scale_fill_manual(values = c("Functional" = col_low, "Evo2" = col_high), name = NULL,
                    labels = c("Functional annotation", "Evo2 model")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(title = "Variance Inflation Factors — Missense Joint Model",
       x = NULL, y = "VIF") +
  theme_mvp() +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 30, hjust = 1))

ggsave(file.path(output_dir, "functional_vif_missense.png"),
       plot_vif_missense, width = 9, height = 5, dpi = 300)
message("Saved: functional_vif_missense.png")

# ============================================================================
# STEP 5: CATEGORICAL ANNOTATION PLOTS (ENCODE cCRE + RegulomeDB rank)
# ============================================================================

message("\n--- Categorical annotation plots ---")

# ENCODE cCRE binary logistic
fit_ccre <- glm(pip_high ~ has_ccre, data = df, family = binomial)
ccre_beta <- coef(fit_ccre)[2]
ccre_or   <- exp(ccre_beta)
ccre_ci   <- exp(confint.default(fit_ccre)[2, ])
ccre_pval <- summary(fit_ccre)$coefficients[2, "Pr(>|z|)"]
chisq_ccre <- chisq.test(table(df$class, df$has_ccre))
message(sprintf("  ENCODE cCRE (binary): OR=%.3f [%.3f, %.3f]  %s  %s",
                ccre_or, ccre_ci[1], ccre_ci[2], format_pval(ccre_pval), sig_stars(ccre_pval)))
message(sprintf("  Chi-square: X2=%.2f  df=%d  %s",
                chisq_ccre$statistic, chisq_ccre$parameter, format_pval(chisq_ccre$p.value)))

# ENCODE cCRE bar chart
ccre_types <- c("PLS", "pELS", "dELS", "CTCF-only", "DNase-H3K4me3")

df_ccre <- df %>%
  filter(!is.na(ENCODE_cCRE)) %>%
  mutate(
    ccre_type = str_extract(ENCODE_cCRE, paste(ccre_types, collapse = "|")),
    ccre_type = if_else(is.na(ccre_type), "Other", ccre_type),
    ccre_type = factor(ccre_type, levels = c(ccre_types, "Other"))
  )

class_totals <- df %>% count(class, name = "total")

ccre_props <- df_ccre %>%
  count(class, ccre_type) %>%
  left_join(class_totals, by = "class") %>%
  mutate(prop = n / total)

plot_ccre_bar <- ggplot(ccre_props, aes(x = ccre_type, y = prop, fill = class)) +
  geom_col(position = "dodge", width = 0.7) +
  scale_fill_manual(values = c("Low PIP" = col_low, "High PIP" = col_high)) +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.12))) +
  annotate("text", x = Inf, y = Inf,
           label = paste0("Binary logistic OR = ", sprintf("%.2f", ccre_or),
                          "\n", format_pval(ccre_pval), "  ", sig_stars(ccre_pval),
                          "\nChi-sq p", format_pval(chisq_ccre$p.value)),
           hjust = 1.05, vjust = 1.3, size = 4.5, fontface = "bold", lineheight = 0.9) +
  labs(title = "ENCODE cCRE Overlap by Variant Class",
       x = "cCRE Category", y = "Proportion of Variants",
       fill = NULL) +
  theme_mvp() +
  theme(legend.position = "bottom", legend.text = element_text(size = 12),
        axis.text.x = element_text(angle = 30, hjust = 1))

# RegulomeDB rank distribution
message("RegulomeDB non-NA variants: ", sum(!is.na(df$RegulomeDB_rank)))

rdb_raw_counts <- df %>%
  filter(!is.na(RegulomeDB_rank)) %>%
  count(class, RegulomeDB_rank) %>%
  group_by(class) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

top_ranks <- rdb_raw_counts %>%
  group_by(RegulomeDB_rank) %>%
  summarise(total = sum(n), .groups = "drop") %>%
  arrange(desc(total)) %>%
  slice_head(n = 8) %>%
  pull(RegulomeDB_rank)

plot_rdb_bar <- rdb_raw_counts %>%
  filter(RegulomeDB_rank %in% top_ranks) %>%
  mutate(RegulomeDB_rank = factor(RegulomeDB_rank, levels = top_ranks)) %>%
  ggplot(aes(x = RegulomeDB_rank, y = prop, fill = class)) +
  geom_col(position = "dodge", width = 0.7) +
  scale_fill_manual(values = c("Low PIP" = col_low, "High PIP" = col_high)) +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.12))) +
  labs(title = "RegulomeDB Rank Distribution by Variant Class\n(Top 8 Ranks)",
       x = "RegulomeDB Rank", y = "Proportion of Variants (within class)",
       fill = NULL) +
  theme_mvp() +
  theme(legend.position = "bottom", legend.text = element_text(size = 12),
        axis.text.x = element_text(angle = 0, hjust = 0.5))

# RegulomeDB ordinal logistic
df_rdb_num <- df %>% drop_na(rdb_rank_num)
fit_rdb <- glm(pip_high ~ scale(rdb_rank_num), data = df_rdb_num, family = binomial)
rdb_beta <- coef(fit_rdb)[2]
rdb_or   <- exp(rdb_beta)
rdb_ci   <- exp(confint.default(fit_rdb)[2, ])
rdb_pval <- summary(fit_rdb)$coefficients[2, "Pr(>|z|)"]
message(sprintf("  RegulomeDB rank (ordinal): OR=%.3f [%.3f, %.3f]  %s  %s",
                rdb_or, rdb_ci[1], rdb_ci[2], format_pval(rdb_pval), sig_stars(rdb_pval)))

plot_categorical <- plot_ccre_bar + plot_rdb_bar +
  plot_layout(ncol = 2) +
  plot_annotation(
    title = "Categorical Functional Annotations: ENCODE cCRE & RegulomeDB Rank",
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18))
  )

ggsave(file.path(output_dir, "functional_categorical_plots.png"),
       plot_categorical, width = 14, height = 6, dpi = 300)
message("Saved: functional_categorical_plots.png")

# ============================================================================
# STEP 6: PAIRWISE CORRELATION PLOTS
# ============================================================================

message("\nBuilding pairwise correlation plots...")

make_corr_panel <- function(x_col, y_col, x_lab, y_lab, data) {
  ct    <- cor.test(data[[x_col]], data[[y_col]], method = "pearson")
  annot <- paste0("r = ", sprintf("%.2f", ct$estimate), "\n", format_pval(ct$p.value))
  ggplot(data, aes(x = .data[[x_col]], y = .data[[y_col]], color = class)) +
    geom_point(alpha = 0.4, size = 1.2) +
    geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.8) +
    annotate("text", x = -Inf, y = Inf, label = annot,
             hjust = -0.1, vjust = 1.4, size = 5, fontface = "bold", lineheight = 0.9) +
    scale_color_manual(values = c("Low PIP" = col_low, "High PIP" = col_high), name = NULL) +
    labs(x = x_lab, y = y_lab) +
    theme_mvp() +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
          legend.position = "bottom", legend.text = element_text(size = 12))
}

# 3 all-variant pairs (CADD, RegulomeDB, Evo2)
corr_panels_allvar <- list(
  make_corr_panel("CADD_phred",      "RegulomeDB_prob",  "CADD Phred Score",           "RegulomeDB Probability Score", df_allvar),
  make_corr_panel("CADD_phred",      "evo2_delta_score", "CADD Phred Score",            "Evo2 Delta Score",             df_allvar),
  make_corr_panel("RegulomeDB_prob", "evo2_delta_score", "RegulomeDB Probability Score","Evo2 Delta Score",             df_allvar)
)

plot_corr_allvar <- wrap_plots(corr_panels_allvar, ncol = 3, guides = "collect") +
  plot_annotation(
    title = "Pairwise Correlations: CADD, RegulomeDB Probability, and Evo2 Score",
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
                  legend.position = "bottom")
  ) &
  theme(legend.position = "bottom")

ggsave(file.path(output_dir, "functional_pairwise_correlations_allvar.png"),
       plot_corr_allvar, width = 15, height = 5, dpi = 300)
message("Saved: functional_pairwise_correlations_allvar.png")

# Missense 7x7 correlation heatmap (all 6 scores + Evo2)
message("Building missense correlation heatmap...")

missense_score_cols <- c("SIFT_score", "PolyPhen2_HVAR", "REVEL_score",
                         "AlphaMissense", "BayesDel_score", "PrimateAI_score",
                         "evo2_delta_score")
missense_labels <- c("SIFT", "PolyPhen2", "REVEL", "AlphaMissense",
                     "BayesDel", "PrimateAI", "Evo2")

cor_mat <- cor(df_missense[, missense_score_cols], use = "complete.obs")

cor_long <- as.data.frame(as.table(cor_mat)) %>%
  rename(x = Var1, y = Var2, r = Freq) %>%
  mutate(
    x = factor(missense_labels[match(as.character(x), missense_score_cols)], levels = missense_labels),
    y = factor(missense_labels[match(as.character(y), missense_score_cols)], levels = rev(missense_labels))
  )

plot_heatmap <- ggplot(cor_long, aes(x = x, y = y, fill = r)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.2f", r)), size = 3.5, fontface = "bold") +
  scale_fill_gradient2(low = col_low, mid = "white", high = col_high,
                       midpoint = 0, limits = c(-1, 1), name = "Pearson r") +
  labs(title = "Missense Score Correlation Matrix\n(including Evo2 Delta Score)",
       x = NULL, y = NULL) +
  theme_bw() +
  theme(
    plot.title      = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.text.x     = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y     = element_text(size = 12),
    panel.grid      = element_blank(),
    legend.position = "right"
  )

ggsave(file.path(output_dir, "functional_missense_corr_heatmap.png"),
       plot_heatmap, width = 9, height = 8, dpi = 300)
message("Saved: functional_missense_corr_heatmap.png")

# ============================================================================
# STEP 7: MAF vs EVO2 CONFOUNDING CHECK (mirrors conservation script)
# ============================================================================

message("\nBuilding MAF vs Evo2 correlation plots...")

df_maf <- df %>%
  drop_na(maf, evo2_delta_score)

message("MAF analysis dataset: ", nrow(df_maf))

run_cor <- function(data, label) {
  ct <- cor.test(data$maf, data$evo2_delta_score, method = "pearson")
  tibble(subset = label, r = ct$estimate, p = ct$p.value, n = nrow(data))
}

cor_maf_evo2 <- bind_rows(
  run_cor(df_maf,                                        "All variants"),
  run_cor(filter(df_maf, class        == "High PIP"),    "High PIP"),
  run_cor(filter(df_maf, class        == "Low PIP"),     "Low PIP"),
  run_cor(filter(df_maf, variant_type == "Coding"),      "Coding"),
  run_cor(filter(df_maf, variant_type == "Non-coding"),  "Non-coding")
)

message("\n--- MAF vs Evo2 Pearson correlations ---")
walk(seq_len(nrow(cor_maf_evo2)), function(i) {
  r <- cor_maf_evo2[i, ]
  message(sprintf("  %-20s  r = %+.2f  %s  n = %d",
                  r$subset, r$r, format_pval(r$p), r$n))
})

# Panel A: all variants
cor_all_label <- paste0("r = ", sprintf("%.2f", cor_maf_evo2$r[1]),
                        "     ", format_pval(cor_maf_evo2$p[1]),
                        "     n = ", cor_maf_evo2$n[1])

plot_maf_all <- ggplot(df_maf, aes(x = maf, y = evo2_delta_score)) +
  geom_point(alpha = 0.35, size = 1.2, color = col_high) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.9) +
  annotate("text", x = -Inf, y = -Inf, label = cor_all_label,
           hjust = -0.05, vjust = -0.5, size = 4, fontface = "bold") +
  labs(title = "All Variants",
       x = "Minor Allele Frequency", y = "Evo2 Delta Score") +
  theme_mvp() +
  theme(legend.position = "none")

# Panel B: 2x2 facet
facet_order <- c("High PIP", "Low PIP", "Coding", "Non-coding")

df_facets <- bind_rows(
  df_maf %>% filter(class        == "High PIP")   %>% mutate(subset = "High PIP"),
  df_maf %>% filter(class        == "Low PIP")    %>% mutate(subset = "Low PIP"),
  df_maf %>% filter(variant_type == "Coding")     %>% mutate(subset = "Coding"),
  df_maf %>% filter(variant_type == "Non-coding") %>% mutate(subset = "Non-coding")
) %>%
  mutate(subset = factor(subset, levels = facet_order))

cor_facet_df <- cor_maf_evo2 %>%
  filter(subset %in% facet_order) %>%
  mutate(
    subset = factor(subset, levels = facet_order),
    label  = paste0("r = ", sprintf("%.2f", r), "     ",
                    format_pval(p), "     n = ", n)
  )

plot_maf_facets <- ggplot(df_facets, aes(x = maf, y = evo2_delta_score)) +
  geom_point(alpha = 0.35, size = 1.0, color = col_high) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.9) +
  geom_text(data = cor_facet_df,
            aes(x = -Inf, y = -Inf, label = label),
            hjust = -0.05, vjust = -0.5,
            size = 3.5, color = "black", fontface = "bold",
            inherit.aes = FALSE) +
  facet_wrap(~subset, ncol = 2) +
  labs(title = "Variant Subsets",
       x = "Minor Allele Frequency", y = "Evo2 Delta Score") +
  theme_mvp() +
  theme(legend.position  = "none",
        strip.text       = element_text(face = "bold", size = 12),
        strip.background = element_blank())

plot_maf_evo2 <- plot_maf_all + plot_maf_facets +
  plot_layout(ncol = 2, widths = c(1, 2)) +
  plot_annotation(
    title = "MAF vs Evo2 Delta Score Correlation",
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18))
  )

ggsave(file.path(output_dir, "functional_maf_evo2_correlation.png"),
       plot_maf_evo2, width = 16, height = 6, dpi = 300)
message("Saved: functional_maf_evo2_correlation.png")

# --- Evo2 logistic: 3-context comparison forest plot ---
message("\n--- Evo2 logistic regression: context comparison ---")

fit_evo2_only <- glm(pip_high ~ scale(evo2_delta_score),
                     data = df_maf, family = binomial)
fit_evo2_maf  <- glm(pip_high ~ scale(evo2_delta_score) + scale(maf),
                     data = df_maf, family = binomial)

extract_row <- function(fit, raw_name, display_label, model_label, is_maf = FALSE) {
  coefs <- summary(fit)$coefficients
  cis   <- confint.default(fit)
  tibble(
    raw_name = raw_name,
    label    = display_label,
    beta     = coefs[raw_name, "Estimate"],
    or       = exp(coefs[raw_name, "Estimate"]),
    ci_lo    = exp(cis[raw_name, 1]),
    ci_hi    = exp(cis[raw_name, 2]),
    p_value  = coefs[raw_name, "Pr(>|z|)"],
    p_label  = format_pval(p_value),
    stars    = sig_stars(p_value),
    model    = model_label,
    var_type = if_else(is_maf, "MAF", "Evo2")
  )
}

# Evo2 from the all-variant joint model (functional adjustment) — re-extract
coef_jv <- summary(fit_joint_allvar)$coefficients
ci_jv   <- confint.default(fit_joint_allvar)
evo2_joint_row <- tibble(
  raw_name = "scale(evo2_delta_score)",
  label    = "Evo2 — adj. for CADD + RegulomeDB",
  beta     = coef_jv["scale(evo2_delta_score)", "Estimate"],
  or       = exp(coef_jv["scale(evo2_delta_score)", "Estimate"]),
  ci_lo    = exp(ci_jv["scale(evo2_delta_score)", 1]),
  ci_hi    = exp(ci_jv["scale(evo2_delta_score)", 2]),
  p_value  = coef_jv["scale(evo2_delta_score)", "Pr(>|z|)"],
  p_label  = format_pval(coef_jv["scale(evo2_delta_score)", "Pr(>|z|)"]),
  stars    = sig_stars(coef_jv["scale(evo2_delta_score)", "Pr(>|z|)"]),
  model    = "Functional joint",
  var_type = "Evo2"
)

maf_cmp_df <- bind_rows(
  extract_row(fit_evo2_only, "scale(evo2_delta_score)", "Evo2 — unadjusted",      "Evo2 only"),
  extract_row(fit_evo2_maf,  "scale(evo2_delta_score)", "Evo2 — adj. for MAF",    "Evo2 + MAF"),
  extract_row(fit_evo2_maf,  "scale(maf)",              "MAF — adj. for Evo2",     "Evo2 + MAF", is_maf = TRUE),
  evo2_joint_row
) %>%
  mutate(
    label = factor(label, levels = c(
      "Evo2 — unadjusted",
      "Evo2 — adj. for MAF",
      "MAF — adj. for Evo2",
      "Evo2 — adj. for CADD + RegulomeDB"
    ))
  )

walk(seq_len(nrow(maf_cmp_df)), function(i) {
  r <- maf_cmp_df[i, ]
  message(sprintf("  %-45s  OR=%.3f [%.3f, %.3f]  %s  %s",
                  as.character(r$label), r$or, r$ci_lo, r$ci_hi, r$p_label, r$stars))
})

# Forest plot
plot_maf_forest_base <- ggplot(maf_cmp_df, aes(x = or, y = label, color = var_type)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.8) +
  geom_errorbar(aes(xmin = ci_lo, xmax = ci_hi),
                width = 0.2, linewidth = 0.8, orientation = "y") +
  geom_point(size = 3.5) +
  scale_x_log10() +
  scale_color_manual(
    values = c("Evo2" = col_high, "MAF" = col_low),
    name   = NULL,
    labels = c("Evo2 coefficient", "MAF coefficient")
  ) +
  labs(
    title = "Evo2 Score: Effect of Adjusting for MAF and Functional Annotations\n(Standardised ORs, 95% CI)",
    x     = "Odds Ratio (95% CI, standardised predictors)",
    y     = NULL
  ) +
  theme_bw() +
  theme(
    plot.title       = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.margin      = margin(5, 2, 5, 5, "mm"),
    axis.text        = element_text(size = 12),
    axis.text.x      = element_text(angle = 45, hjust = 1),
    axis.text.y      = element_text(size = 12),
    axis.title       = element_text(size = 13),
    legend.position  = "bottom",
    legend.text      = element_text(size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

maf_side_df <- maf_cmp_df %>%
  mutate(side_label = paste0("β = ", sprintf("%+.2f", beta), "\n",
                             p_label, "  ", stars))

plot_maf_side <- ggplot(maf_side_df, aes(x = 0, y = label)) +
  geom_text(aes(label = side_label),
            hjust = 0, size = 4.5, fontface = "bold", lineheight = 0.9) +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  labs(title = "", x = NULL, y = NULL) +
  theme_void() +
  theme(plot.margin = margin(5, 5, 5, 2, "mm"),
        axis.text = element_blank(), axis.ticks = element_blank())

plot_maf_forest <- plot_maf_forest_base + plot_maf_side +
  plot_layout(widths = c(3, 1))

ggsave(file.path(output_dir, "functional_maf_regression_forest.png"),
       plot_maf_forest, width = 13, height = 5, dpi = 300)
message("Saved: functional_maf_regression_forest.png")

message("\n=== All Done ===")
