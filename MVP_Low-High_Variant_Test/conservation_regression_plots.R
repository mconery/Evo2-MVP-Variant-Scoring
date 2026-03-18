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
         phastCons100way, phyloP100way, GERP_RS, evo2_delta_score)

message("Merged rows: ", nrow(df))
message("Class distribution: ", paste(names(table(df$class)), table(df$class), sep = "=", collapse = ", "))

# Complete cases across all four predictors
df_complete <- df %>% drop_na(phastCons100way, phyloP100way, GERP_RS, evo2_delta_score)
message("Complete cases (all 4 predictors present): ", nrow(df_complete))

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

message("\n=== Done ===")
