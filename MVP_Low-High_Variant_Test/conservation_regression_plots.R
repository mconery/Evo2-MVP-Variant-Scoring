# ============================================================================
# Logistic Regression: High vs Low PIP ~ Conservation + Evo2 Scores
# ============================================================================
# Fits individual and joint logistic regressions predicting High vs Low PIP
# status from three conservation metrics (phastCons, phyloP, GERP) and the
# Evo2 7b-arc-longcontext delta score (16384 bp context window).
#
# Produces two plots in the same style as mvp_plot_script.R:
#   1. conservation_individual_regressions.png  — 2x2 violin panel (individual models)
#   2. conservation_joint_regression_forest.png — forest plot (joint model)
# ============================================================================

library(tidyverse)
library(ggplot2)
library(scales)
library(patchwork)
library(car)

# ============================================================================
# CONFIGURATION
# ============================================================================

conservation_file <- "C:/Users/mitch/Documents/Argonne/Variant_Scoring/MVP_conservation_scores.csv"
evo2_file         <- "C:/Users/mitch/Documents/Argonne/Variant_Scoring/scored_MVP_results/MVP_variant_scores.7b_arc_longcontext_model.16384bp_context.csv"
output_dir        <- "C:/Users/mitch/Documents/Argonne/Variant_Scoring/scored_MVP_figures"

# Colors matching mvp_plot_script.R
col_low  <- "#6baed6"
col_high <- "#08519c"

# ============================================================================
# SHARED THEME (mirrors mvp_plot_script.R)
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

# ============================================================================
# STEP 1: LOAD AND MERGE DATA
# ============================================================================

message("Loading conservation scores...")
conservation <- read_csv(conservation_file, show_col_types = FALSE)

message("Loading Evo2 scores (7b-arc-longcontext, 16384 bp)...")
evo2 <- read_csv(evo2_file, show_col_types = FALSE) %>%
  select(`MVP ID`, evo2_delta_score, class)

# Harmonise class labels: conservation uses "High_PIP" / "Low_PIP" (underscores)
conservation <- conservation %>%
  mutate(class = str_replace(Set, "_", " "))   # "High PIP" / "Low PIP"

# Join on MVP ID
df <- conservation %>%
  inner_join(evo2, by = "MVP ID", suffix = c("", ".evo2")) %>%
  mutate(
    class     = coalesce(class, class.evo2),   # both should match; keep one
    pip_high  = as.integer(class == "High PIP")
  ) %>%
  select(`MVP ID`, class, pip_high,
         phastCons100way, phyloP100way, GERP_RS, evo2_delta_score,
         `EAF Population`, `VEP Annotation`)

message("Merged rows: ", nrow(df))
message("Class distribution: ", paste(names(table(df$class)), table(df$class), sep = "=", collapse = ", "))

# Complete cases across all four predictors
df_complete <- df %>% drop_na(phastCons100way, phyloP100way, GERP_RS, evo2_delta_score)
message("Complete cases (all 4 predictors present): ", nrow(df_complete))

# Coding VEP consequences (standard Ensembl/VEP exonic terms)
coding_consequences <- c(
  "missense_variant", "synonymous_variant", "stop_gained", "stop_lost",
  "start_lost", "frameshift_variant", "inframe_insertion", "inframe_deletion",
  "splice_donor_variant", "splice_acceptor_variant"
)

df_maf <- df %>%
  mutate(
    maf          = pmin(`EAF Population`, 1 - `EAF Population`),
    variant_type = if_else(`VEP Annotation` %in% coding_consequences,
                           "Coding", "Non-coding"),
    log10_evo2   = if_else(evo2_delta_score != 0,
                           log10(abs(evo2_delta_score)), NA_real_)
  ) %>%
  drop_na(maf, evo2_delta_score)

message("MAF analysis dataset: ", nrow(df_maf), " variants")
message("VEP annotation breakdown:\n",
        paste(capture.output(print(table(df_maf$variant_type))), collapse = "\n"))

# ============================================================================
# STEP 2: INDIVIDUAL LOGISTIC REGRESSIONS
# ============================================================================

predictors <- list(
  phastCons100way  = list(col = "phastCons100way",  label = "phastCons (100-way vertebrate)"),
  phyloP100way     = list(col = "phyloP100way",     label = "phyloP (100-way vertebrate)"),
  GERP_RS          = list(col = "GERP_RS",          label = "GERP++ RS"),
  evo2_delta_score = list(col = "evo2_delta_score", label = "Evo2 Delta Score\n(7b-arc-longcontext, 16384 bp)")
)

message("\n--- Individual logistic regressions ---")

individual_results <- map_dfr(names(predictors), function(name) {
  meta  <- predictors[[name]]
  col   <- meta$col
  label <- meta$label

  formula <- as.formula(paste0("pip_high ~ `", col, "`"))
  fit     <- glm(formula, data = df_complete, family = binomial)
  coef_df <- summary(fit)$coefficients

  # Second row = predictor (skip intercept)
  beta  <- coef(fit)[2]
  or    <- exp(beta)
  ci    <- exp(confint.default(fit)[2, ])
  pval  <- coef_df[2, "Pr(>|z|)"]

  message(sprintf("  %-30s  beta=%+.3f  OR=%.3f [%.3f, %.3f]  %s  %s",
                  name, beta, or, ci[1], ci[2], format_pval(pval), sig_stars(pval)))

  tibble(
    predictor = name,
    label     = label,
    beta      = beta,
    or        = or,
    ci_lo     = ci[1],
    ci_hi     = ci[2],
    p_value   = pval,
    p_label   = format_pval(pval),
    stars     = sig_stars(pval),
    model     = "Individual"
  )
})

# ============================================================================
# STEP 3: JOINT LOGISTIC REGRESSION (STANDARDISED PREDICTORS)
# ============================================================================

message("\n--- Joint logistic regression (standardised predictors) ---")

fit_joint <- glm(
  pip_high ~ scale(phastCons100way) + scale(phyloP100way) +
             scale(GERP_RS) + scale(evo2_delta_score),
  data   = df_complete,
  family = binomial
)

coef_joint <- summary(fit_joint)$coefficients[-1, , drop = FALSE]   # drop intercept
ci_joint   <- confint.default(fit_joint)[-1, , drop = FALSE]

joint_results <- tibble(
  predictor = names(predictors),
  label     = map_chr(predictors, "label"),
  beta      = coef_joint[, "Estimate"],
  or        = exp(coef_joint[, "Estimate"]),
  ci_lo     = exp(ci_joint[, 1]),
  ci_hi     = exp(ci_joint[, 2]),
  p_value   = coef_joint[, "Pr(>|z|)"],
  p_label   = format_pval(p_value),
  stars     = sig_stars(p_value),
  model     = "Joint"
)

walk(seq_len(nrow(joint_results)), function(i) {
  r <- joint_results[i, ]
  message(sprintf("  %-30s  OR=%.3f [%.3f, %.3f]  %s  %s",
                  r$predictor, r$or, r$ci_lo, r$ci_hi, r$p_label, r$stars))
})

# ============================================================================
# STEP 3b: VARIANCE INFLATION FACTORS (Joint Model)
# ============================================================================

message("\n--- Variance Inflation Factors (joint model) ---")

vif_vals <- car::vif(fit_joint)   # named numeric vector

# Build tidy VIF table — strip scale() wrapper from names, map to labels
vif_df <- tibble(
  predictor = names(predictors),
  label     = map_chr(predictors, "label"),
  vif       = as.numeric(vif_vals),
  var_type  = if_else(predictor == "evo2_delta_score", "Evo2", "Conservation")
)

walk(seq_len(nrow(vif_df)), function(i) {
  message(sprintf("  %-30s  VIF = %.3f", vif_df$predictor[i], vif_df$vif[i]))
})

# --- VIF bar chart ---
vif_df <- vif_df %>%
  mutate(
    label = str_replace(label, " \\(", "\n("),   # wrap at opening parenthesis,
    label = str_replace(label, "7b-arc-longcontext", "7b long-context"),
    label = factor(label, levels = label)
  )

plot_vif <- ggplot(vif_df, aes(x = label, y = vif, fill = var_type)) +
  geom_col(width = 0.6, color = "white") +
  geom_hline(yintercept = 5,  linetype = "dashed", color = "orange",    linewidth = 0.8) +
  geom_hline(yintercept = 10, linetype = "dashed", color = "firebrick", linewidth = 0.8) +
  geom_text(aes(label = sprintf("%.2f", vif)),
            vjust = -0.4, size = 4, fontface = "bold") +
  annotate("text", x = Inf, y = 5,  label = "VIF = 5",  hjust = 1.1, vjust = -0.4,
           color = "orange",    size = 3.5) +
  annotate("text", x = Inf, y = 10, label = "VIF = 10", hjust = 1.1, vjust = -0.4,
           color = "firebrick", size = 3.5) +
  scale_fill_manual(
    values = c("Conservation" = col_low, "Evo2" = col_high),
    name   = NULL,
    labels = c("Conservation metric", "Evo2 model")
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    title = "Variance Inflation Factors \u2014 Joint Logistic Regression",
    x     = NULL,
    y     = "VIF"
  ) +
  theme_mvp() +
  theme(
    legend.position = "bottom",
    axis.text.x     = element_text(angle = 0, hjust = 0.5)
  )

# ============================================================================
# STEP 4: SUMMARY TABLE
# ============================================================================

message("\n=== Summary Table ===")
bind_rows(individual_results, joint_results) %>%
  select(model, predictor, or, ci_lo, ci_hi, p_value, stars) %>%
  mutate(across(c(or, ci_lo, ci_hi), ~round(.x, 3))) %>%
  print(n = Inf)

# ============================================================================
# STEP 5: PLOT A — INDIVIDUAL VIOLIN PANELS (2×2)
# ============================================================================

message("\nBuilding violin panel (individual regressions)...")

make_violin_panel <- function(pred_name, pred_label, p_label_str, beta_val) {
  annotation_str <- paste0(
    "\u03b2 = ", sprintf("%+.3f", beta_val), "\n",
    p_label_str
  )
  ggplot(df_complete,
         aes(x = class,
             y = .data[[pred_name]],
             fill = class)) +
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
    theme(
      axis.text  = element_text(size = 14),
      axis.title = element_text(size = 16)
    )
}

violin_panels <- pmap(
  list(individual_results$predictor,
       individual_results$label,
       individual_results$p_label,
       individual_results$beta),
  make_violin_panel
)

plot_violins <- wrap_plots(violin_panels, ncol = 2) +
  plot_annotation(
    title = "High vs Low PIP: Conservation & Evo2 Score Distributions\n(\u03b2 and p-values from individual logistic regressions)",
    theme = theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 18)
    )
  )

# ============================================================================
# STEP 6: PLOT B — JOINT MODEL FOREST PLOT + SIDE TABLE
# ============================================================================

message("Building forest plot (joint regression)...")

# Predictor display order (bottom to top on y-axis)
forest_df <- joint_results %>%
  mutate(
    label    = factor(label, levels = rev(label)),
    var_type = if_else(predictor == "evo2_delta_score", "Evo2", "Conservation")
  )

# --- 6a: forest plot (no in-plot text) ---
plot_forest_base <- ggplot(forest_df,
                           aes(x = or, y = label, color = var_type)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.8) +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi),
                 height = 0.2, linewidth = 0.8) +
  geom_point(size = 3.5) +
  scale_x_log10() +
  scale_color_manual(
    values  = c("Conservation" = col_low, "Evo2" = col_high),
    name    = NULL,
    labels  = c("Conservation metric", "Evo2 model")
  ) +
  labs(
    title = "Joint Logistic Regression: High vs Low PIP\n(Standardised Odds Ratios, 95% CI)",
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

# --- 6b: side text panel showing β and p-value per row ---
side_df <- forest_df %>%
  mutate(
    side_label = paste0(
      "\u03b2 = ", sprintf("%+.3f", beta), "\n",
      p_label, "  ", stars
    )
  )

plot_side <- ggplot(side_df, aes(x = 0, y = label)) +
  geom_text(aes(label = side_label),
            hjust = 0, size = 4.5, fontface = "bold",
            lineheight = 0.9) +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  labs(title = "", x = NULL, y = NULL) +
  theme_void() +
  theme(
    plot.title    = element_text(hjust = 0, face = "bold", size = 13,
                                 margin = margin(b = 4)),
    plot.margin   = margin(5, 5, 5, 2, "mm"),
    # add bottom margin to align with forest plot x-axis space
    axis.text     = element_blank(),
    axis.ticks    = element_blank()
  )

# --- 6c: compose with patchwork ---
plot_forest <- plot_forest_base + plot_side +
  plot_layout(widths = c(3, 1))

# ============================================================================
# STEP 7: SAVE PLOTS
# ============================================================================

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  message("Created output directory: ", output_dir)
}

violin_path <- file.path(output_dir, "conservation_individual_regressions.png")
ggsave(violin_path, plot_violins, width = 12, height = 10, dpi = 300)
message("Saved: ", violin_path)

forest_path <- file.path(output_dir, "conservation_joint_regression_forest.png")
ggsave(forest_path, plot_forest, width = 12, height = 5, dpi = 300)
message("Saved: ", forest_path)

vif_path <- file.path(output_dir, "conservation_vif_barplot.png")
ggsave(vif_path, plot_vif, width = 8, height = 5, dpi = 300)
message("Saved: ", vif_path)

message("\n=== Done ===")

# ============================================================================
# STEP 8: PAIRWISE CORRELATION PLOTS
# ============================================================================

message("\nBuilding pairwise correlation plots...")

pairs_def <- list(
  list(x = "phastCons100way",  y = "phyloP100way",     x_lab = "phastCons (100-way vertebrate)", y_lab = "phyloP (100-way vertebrate)",             pair_type = "Conservation vs Conservation"),
  list(x = "phastCons100way",  y = "GERP_RS",           x_lab = "phastCons (100-way vertebrate)", y_lab = "GERP++ RS",                               pair_type = "Conservation vs Conservation"),
  list(x = "phyloP100way",     y = "GERP_RS",           x_lab = "phyloP (100-way vertebrate)",    y_lab = "GERP++ RS",                               pair_type = "Conservation vs Conservation"),
  list(x = "phastCons100way",  y = "evo2_delta_score",  x_lab = "phastCons (100-way vertebrate)", y_lab = "Evo2 Delta Score\n(7b-arc-longcontext, 16384 bp)", pair_type = "Conservation vs Evo2"),
  list(x = "phyloP100way",     y = "evo2_delta_score",  x_lab = "phyloP (100-way vertebrate)",    y_lab = "Evo2 Delta Score\n(7b-arc-longcontext, 16384 bp)", pair_type = "Conservation vs Evo2"),
  list(x = "GERP_RS",          y = "evo2_delta_score",  x_lab = "GERP++ RS",                      y_lab = "Evo2 Delta Score\n(7b-arc-longcontext, 16384 bp)", pair_type = "Conservation vs Evo2")
)

make_corr_panel <- function(pair) {
  x_col <- pair$x
  y_col <- pair$y

  ct    <- cor.test(df_complete[[x_col]], df_complete[[y_col]], method = "pearson")
  r_val <- sprintf("r = %.3f", ct$estimate)
  p_val <- format_pval(ct$p.value)
  annot <- paste0(r_val, "\n", p_val)

  ggplot(df_complete,
         aes(x = .data[[x_col]], y = .data[[y_col]], color = class)) +
    geom_point(alpha = 0.4, size = 1.5) +
    geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.8) +
    annotate("text",
             x = -Inf, y = Inf,
             label = annot,
             hjust = -0.1, vjust = 1.4,
             size = 5, fontface = "bold",
             lineheight = 0.9) +
    scale_color_manual(
      values = c("Low PIP" = col_low, "High PIP" = col_high),
      name   = NULL
    ) +
    labs(x = pair$x_lab, y = pair$y_lab) +
    theme_mvp() +
    theme(
      axis.text.x     = element_text(angle = 0, hjust = 0.5),
      legend.position = "bottom",
      legend.text     = element_text(size = 12)
    )
}

corr_panels <- map(pairs_def, make_corr_panel)

plot_corr <- wrap_plots(corr_panels, ncol = 3, nrow = 2, guides = "collect") +
  plot_annotation(
    title    = "Pairwise Correlations: Conservation Metrics & Evo2 Score",
    subtitle = "Row 1: Conservation vs Conservation    |    Row 2: Conservation vs Evo2",
    theme    = theme(
      plot.title    = element_text(hjust = 0.5, face = "bold", size = 18),
      plot.subtitle = element_text(hjust = 0.5, size = 13),
      legend.position = "bottom"
    )
  ) &
  theme(legend.position = "bottom")

corr_path <- file.path(output_dir, "conservation_pairwise_correlations.png")
ggsave(corr_path, plot_corr, width = 14, height = 10, dpi = 300)
message("Saved: ", corr_path)

# ============================================================================
# STEP 9: MAF vs EVO2 SCORE CORRELATION
# ============================================================================

message("\nBuilding MAF vs Evo2 correlation plot...")

# --- 9a: Pearson correlations for each subset ---
run_cor <- function(data, label) {
  ct <- cor.test(data$maf, data$evo2_delta_score, method = "pearson")
  tibble(
    subset = label,
    r      = ct$estimate,
    p      = ct$p.value,
    n      = nrow(data)
  )
}

cor_results <- bind_rows(
  run_cor(df_maf,                                       "All variants"),
  run_cor(filter(df_maf, class        == "High PIP"),   "High PIP"),
  run_cor(filter(df_maf, class        == "Low PIP"),    "Low PIP"),
  run_cor(filter(df_maf, variant_type == "Coding"),     "Coding"),
  run_cor(filter(df_maf, variant_type == "Non-coding"), "Non-coding")
)

message("\n--- MAF vs Evo2 Pearson correlations ---")
walk(seq_len(nrow(cor_results)), function(i) {
  r <- cor_results[i, ]
  message(sprintf("  %-20s  r = %+.3f  %s  n = %d",
                  r$subset, r$r, format_pval(r$p), r$n))
})

# --- 9b: Left panel — All variants ---
cor_all_label <- paste0("r = ", sprintf("%.3f", cor_results$r[1]),
                        "     ", format_pval(cor_results$p[1]),
                        "     n = ", cor_results$n[1])

plot_maf_all <- ggplot(df_maf, aes(x = maf, y = evo2_delta_score)) +
  geom_point(alpha = 0.35, size = 1.2, color = col_high) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.9) +
  annotate("text", x = -Inf, y = -Inf, label = cor_all_label,
           hjust = -0.05, vjust = -0.5, size = 4, fontface = "bold") +
  labs(
    title = "All Variants",
    x     = "Minor Allele Frequency",
    y     = "Evo2 Delta Score"
  ) +
  theme_mvp() +
  theme(legend.position = "none")

# --- 9c: Right panel — 2x2 faceted scatter ---
facet_order <- c("High PIP", "Low PIP", "Coding", "Non-coding")

df_facets <- bind_rows(
  df_maf %>% filter(class        == "High PIP")    %>% mutate(subset = "High PIP"),
  df_maf %>% filter(class        == "Low PIP")     %>% mutate(subset = "Low PIP"),
  df_maf %>% filter(variant_type == "Coding")      %>% mutate(subset = "Coding"),
  df_maf %>% filter(variant_type == "Non-coding")  %>% mutate(subset = "Non-coding")
) %>%
  mutate(subset = factor(subset, levels = facet_order))

cor_facet_df <- cor_results %>%
  filter(subset %in% facet_order) %>%
  mutate(
    subset = factor(subset, levels = facet_order),
    label  = paste0("r = ", sprintf("%.3f", r), "     ", format_pval(p),
                    "     n = ", n)
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
  labs(
    title = "Variant Subsets",
    x     = "Minor Allele Frequency",
    y     = "Evo2 Delta Score"
  ) +
  theme_mvp() +
  theme(
    legend.position  = "none",
    strip.text       = element_text(face = "bold", size = 12),
    strip.background = element_blank()
  )

# --- 9d: Combine panels with patchwork ---
plot_maf <- plot_maf_all + plot_maf_facets +
  plot_layout(ncol = 2, widths = c(1, 2)) +
  plot_annotation(
    title = "MAF vs Evo2 Delta Score Correlation",
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18))
  )

maf_path <- file.path(output_dir, "conservation_maf_evo2_correlation.png")
ggsave(maf_path, plot_maf, width = 16, height = 6, dpi = 300)
message("Saved: ", maf_path)

# ============================================================================
# STEP 10: MAF vs LOG10(|EVO2 DELTA SCORE|) CORRELATION
# ============================================================================

message("\nBuilding MAF vs log10(|Evo2 delta score|) correlation plot...")

run_cor_y <- function(data, label, y_col) {
  y_col <- rlang::ensym(y_col)
  ct <- cor.test(data$maf, dplyr::pull(data, !!y_col), method = "pearson")
  tibble(subset = label, r = ct$estimate, p = ct$p.value, n = nrow(data))
}

df_maf_log <- df_maf %>% drop_na(log10_evo2)
message("log10 analysis dataset: ", nrow(df_maf_log),
        " variants (", nrow(df_maf) - nrow(df_maf_log), " zeros dropped)")

# --- 10a: Pearson correlations for each subset ---
cor_results_log <- bind_rows(
  run_cor_y(df_maf_log,                                       "All variants", log10_evo2),
  run_cor_y(filter(df_maf_log, class        == "High PIP"),   "High PIP",     log10_evo2),
  run_cor_y(filter(df_maf_log, class        == "Low PIP"),    "Low PIP",      log10_evo2),
  run_cor_y(filter(df_maf_log, variant_type == "Coding"),     "Coding",       log10_evo2),
  run_cor_y(filter(df_maf_log, variant_type == "Non-coding"), "Non-coding",   log10_evo2)
)

message("\n--- MAF vs log10(|Evo2|) Pearson correlations ---")
walk(seq_len(nrow(cor_results_log)), function(i) {
  r <- cor_results_log[i, ]
  message(sprintf("  %-20s  r = %+.3f  %s  n = %d",
                  r$subset, r$r, format_pval(r$p), r$n))
})

# --- 10b: Left panel — All variants ---
cor_log_all_label <- paste0("r = ", sprintf("%.3f", cor_results_log$r[1]),
                            "     ", format_pval(cor_results_log$p[1]),
                            "     n = ", cor_results_log$n[1])

plot_maf_log_all <- ggplot(df_maf_log, aes(x = maf, y = log10_evo2)) +
  geom_point(alpha = 0.35, size = 1.2, color = col_high) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.9) +
  scale_x_log10() +
  annotate("text", x = -Inf, y = -Inf, label = cor_log_all_label,
           hjust = -0.05, vjust = -0.5, size = 4, fontface = "bold") +
  labs(
    title = "All Variants",
    x     = expression(log[10]("MAF")),
    y     = expression(log[10]("|Evo2 Delta Score|"))
  ) +
  theme_mvp() +
  theme(legend.position = "none")

# --- 10c: Right panel — 2x2 faceted scatter ---
df_facets_log <- bind_rows(
  df_maf_log %>% filter(class        == "High PIP")    %>% mutate(subset = "High PIP"),
  df_maf_log %>% filter(class        == "Low PIP")     %>% mutate(subset = "Low PIP"),
  df_maf_log %>% filter(variant_type == "Coding")      %>% mutate(subset = "Coding"),
  df_maf_log %>% filter(variant_type == "Non-coding")  %>% mutate(subset = "Non-coding")
) %>%
  mutate(subset = factor(subset, levels = facet_order))

cor_facet_log_df <- cor_results_log %>%
  filter(subset %in% facet_order) %>%
  mutate(
    subset = factor(subset, levels = facet_order),
    label  = paste0("r = ", sprintf("%.3f", r), "     ", format_pval(p),
                    "     n = ", n)
  )

plot_maf_log_facets <- ggplot(df_facets_log, aes(x = maf, y = log10_evo2)) +
  geom_point(alpha = 0.35, size = 1.0, color = col_high) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.9) +
  scale_x_log10() +
  geom_text(data = cor_facet_log_df,
            aes(x = -Inf, y = -Inf, label = label),
            hjust = -0.05, vjust = -0.5,
            size = 3.5, color = "black", fontface = "bold",
            inherit.aes = FALSE) +
  facet_wrap(~subset, ncol = 2) +
  labs(
    title = "Variant Subsets",
    x     = expression(log[10]("MAF")),
    y     = expression(log[10]("|Evo2 Delta Score|"))
  ) +
  theme_mvp() +
  theme(
    legend.position  = "none",
    strip.text       = element_text(face = "bold", size = 12),
    strip.background = element_blank()
  )

# --- 10d: Combine panels with patchwork ---
plot_maf_log <- plot_maf_log_all + plot_maf_log_facets +
  plot_layout(ncol = 2, widths = c(1, 2)) +
  plot_annotation(
    title = expression("MAF vs " * log[10]("|Evo2 Delta Score|") * " Correlation"),
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18))
  )

maf_log_path <- file.path(output_dir, "conservation_maf_log10evo2_correlation.png")
ggsave(maf_log_path, plot_maf_log, width = 16, height = 6, dpi = 300)
message("Saved: ", maf_log_path)

# ============================================================================
# STEP 11: EVO2 LOGISTIC REGRESSION — WITH AND WITHOUT MAF COVARIATE
# ============================================================================

message("\n--- Evo2 logistic regression: with and without MAF covariate ---")

# --- 11a: Fit models ---
fit_evo2_only <- glm(pip_high ~ scale(evo2_delta_score),
                     data = df_maf, family = binomial)
fit_evo2_maf  <- glm(pip_high ~ scale(evo2_delta_score) + scale(maf),
                     data = df_maf, family = binomial)

# --- 11b: Extract results ---
extract_coefs <- function(fit, model_label, predictors_keep) {
  coefs <- summary(fit)$coefficients
  cis   <- confint.default(fit)
  rows  <- rownames(coefs)[rownames(coefs) %in% predictors_keep]
  tibble(
    raw_name = rows,
    beta     = coefs[rows, "Estimate"],
    or       = exp(coefs[rows, "Estimate"]),
    ci_lo    = exp(cis[rows, 1]),
    ci_hi    = exp(cis[rows, 2]),
    p_value  = coefs[rows, "Pr(>|z|)"],
    p_label  = format_pval(p_value),
    stars    = sig_stars(p_value),
    model    = model_label
  )
}

maf_forest_df <- bind_rows(
  extract_coefs(fit_evo2_only,
                "Evo2 only",
                "scale(evo2_delta_score)"),
  extract_coefs(fit_evo2_maf,
                "Evo2 + MAF",
                c("scale(evo2_delta_score)", "scale(maf)"))
) %>%
  mutate(
    label = case_when(
      raw_name == "scale(evo2_delta_score)" & model == "Evo2 only"  ~ "Evo2 \u2014 without MAF",
      raw_name == "scale(evo2_delta_score)" & model == "Evo2 + MAF" ~ "Evo2 \u2014 with MAF",
      raw_name == "scale(maf)"                                       ~ "MAF \u2014 with Evo2"
    ),
    label = factor(label, levels = c("Evo2 \u2014 without MAF",
                                     "Evo2 \u2014 with MAF",
                                     "MAF \u2014 with Evo2"))
  )

walk(seq_len(nrow(maf_forest_df)), function(i) {
  r <- maf_forest_df[i, ]
  message(sprintf("  %-35s  OR=%.3f [%.3f, %.3f]  %s  %s",
                  r$label, r$or, r$ci_lo, r$ci_hi, r$p_label, r$stars))
})

# --- 11c: Forest plot base ---
plot_maf_forest_base <- ggplot(maf_forest_df,
                               aes(x = or, y = label, color = model)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.8) +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi),
                 height = 0.2, linewidth = 0.8) +
  geom_point(size = 3.5) +
  scale_x_log10() +
  scale_color_manual(
    values = c("Evo2 only" = col_high, "Evo2 + MAF" = col_low),
    name   = NULL
  ) +
  labs(
    title = "Evo2 Score: High vs Low PIP\n(Standardised ORs, 95% CI \u2014 with and without MAF covariate)",
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

# --- 11d: Side text panel ---
maf_side_df <- maf_forest_df %>%
  mutate(side_label = paste0("\u03b2 = ", sprintf("%+.3f", beta), "\n",
                             p_label, "  ", stars))

plot_maf_side <- ggplot(maf_side_df, aes(x = 0, y = label)) +
  geom_text(aes(label = side_label),
            hjust = 0, size = 4.5, fontface = "bold", lineheight = 0.9) +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  labs(title = "", x = NULL, y = NULL) +
  theme_void() +
  theme(
    plot.margin = margin(5, 5, 5, 2, "mm"),
    axis.text   = element_blank(),
    axis.ticks  = element_blank()
  )

# --- 11e: Compose and save ---
plot_maf_forest <- plot_maf_forest_base + plot_maf_side +
  plot_layout(widths = c(3, 1))

maf_forest_path <- file.path(output_dir, "conservation_maf_regression_forest.png")
ggsave(maf_forest_path, plot_maf_forest, width = 12, height = 4, dpi = 300)
message("Saved: ", maf_forest_path)

message("\n=== All Done ===")
