# ============================================================================
# R Script to Plot Evo2 Timing Results
# ============================================================================
# Produces three scatter plots (with connected lines) from timing experiments
# measuring Evo2 run time under varying parallelism, context length, and chunk
# size conditions.  Visual theme matches mvp_plot_script.R.
# ============================================================================

library(tidyverse)
library(ggplot2)

# ============================================================================
# CONFIGURATION
# ============================================================================

input_file <- "C:/Users/mitch/Documents/Argonne/Variant_Scoring/timing_MVP_results/timing_results.csv"
output_dir <- "C:/Users/mitch/Documents/Argonne/Variant_Scoring/timing_MVP_results"

plot_width  <- 8
plot_height <- 6
output_format <- "png"

# Colors: match Low PIP / High PIP palette from mvp_plot_script.R
color_7b  <- "#6baed6"   # same color used for all single-line plots
color_40b <- "#08519c"

# ============================================================================
# SHARED THEME
# ============================================================================

timing_theme <- theme_bw() +
  theme(
    plot.title       = element_text(hjust = 0.5, face = "bold", size = 18),
    plot.margin      = margin(5, 20, 5, 5, "mm"),
    axis.text        = element_text(size = 12),
    axis.title       = element_text(size = 13),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# ============================================================================
# HELPER
# ============================================================================

save_plot <- function(plot, filename, width = plot_width, height = plot_height) {
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    message("Created output directory: ", output_dir)
  }
  out_path <- file.path(output_dir, paste0(filename, ".", output_format))
  ggsave(filename = out_path, plot = plot, width = width, height = height, dpi = 300)
  message("Saved: ", out_path)
}

# ============================================================================
# LOAD DATA
# ============================================================================

timing <- read_csv(input_file, show_col_types = FALSE) %>%
  mutate(time_minutes = duration_seconds / 60)

# ============================================================================
# PLOT 1: Tensor vs. Context Parallelism
# Filter: 7b model, context = 8192 bp, chunk_size = 10
# X: tp/cp combination  |  Y: time in minutes
# ============================================================================

plot1_data <- timing %>%
  filter(
    model      == "7b_arc_longcontext",
    context_bp == 8192,
    chunk_size == 10
  ) %>%
  arrange(tp_size) %>%
  mutate(
    parallelism_label = factor(
      paste0("TP=", tp_size, "\nCP=", cp_size),
      levels = paste0("TP=", tp_size, "\nCP=", cp_size)
    )
  )

plot1 <- ggplot(plot1_data, aes(x = parallelism_label, y = time_minutes, group = 1)) +
  geom_line(color = color_7b, linewidth = 1.2) +
  geom_point(color = color_7b, size = 3.5) +
  labs(
    title = "Tensor vs. Context Parallelism\n(7b Arc Long-Context, 8,192 bp Context, Chunk Size = 10)",
    x     = "Tensor / Context Parallelism",
    y     = "Run Time (minutes)"
  ) +
  timing_theme

save_plot(plot1, "timing_parallelism_vs_time")

# ============================================================================
# PLOT 2: Context Length vs. Time
# Filter: 7b model, tp=4, cp=2, chunk_size=100
# X: context_bp  |  Y: time in minutes
# ============================================================================

context_bp_to_label <- function(bp) {
  dplyr::case_when(
    bp == 8192   ~ "8 Kbp",
    bp == 16384  ~ "16 Kbp",
    bp == 65536  ~ "64 Kbp",
    bp == 131072 ~ "128 Kbp",
    bp == 524288 ~ "512 Kbp",
    TRUE         ~ paste0(bp, " bp")
  )
}

plot2_data <- timing %>%
  filter(
    model      == "7b_arc_longcontext",
    tp_size    == 4,
    cp_size    == 2,
    chunk_size == 100
  ) %>%
  arrange(context_bp)

plot2_breaks <- sort(unique(plot2_data$context_bp))

plot2 <- ggplot(plot2_data, aes(x = context_bp, y = time_minutes, group = 1)) +
  geom_line(color = color_7b, linewidth = 1.2) +
  geom_point(color = color_7b, size = 3.5) +
  scale_x_log10(
    breaks = plot2_breaks,
    labels = context_bp_to_label(plot2_breaks)
  ) +
  labs(
    title = "Run Time vs. Context Length\n(7b Arc Long-Context, TP=4, CP=2, Chunk Size = 100)",
    x     = "Sequence Length",
    y     = "Run Time (minutes)"
  ) +
  timing_theme

save_plot(plot2, "timing_context_length_vs_time")

plot2_linear <- ggplot(plot2_data, aes(x = context_bp, y = time_minutes, group = 1)) +
  geom_line(color = color_7b, linewidth = 1.2) +
  geom_point(color = color_7b, size = 3.5) +
  scale_x_continuous(
    breaks = plot2_breaks,
    labels = context_bp_to_label(plot2_breaks)
  ) +
  labs(
    title = "Run Time vs. Context Length\n(7b Arc Long-Context, TP=4, CP=2, Chunk Size = 100)",
    x     = "Sequence Length",
    y     = "Run Time (minutes)"
  ) +
  timing_theme

save_plot(plot2_linear, "timing_context_length_vs_time_linear")

# ============================================================================
# PLOT 3: Chunk Size vs. Time (both models)
# Filter: tp=8, cp=1, context=8192 bp
# X: chunk_size  |  Y: time in minutes  |  Color: model
# ============================================================================

plot3_data <- timing %>%
  filter(
    tp_size    == 8,
    cp_size    == 1,
    context_bp == 8192
  ) %>%
  arrange(model, chunk_size) %>%
  mutate(chunk_size = factor(chunk_size, levels = c(10, 20, 50, 100)))

model_colors <- c("7b_arc_longcontext"  = color_7b,
                  "40b_arc_longcontext" = color_40b)
model_labels <- c("7b_arc_longcontext"  = "7b Arc Long-Context",
                  "40b_arc_longcontext" = "40b Arc Long-Context")

plot3 <- ggplot(plot3_data, aes(x = chunk_size, y = time_minutes,
                                color = model, group = model)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3.5) +
  scale_color_manual(values = model_colors, labels = model_labels) +
  labs(
    title = "Run Time vs. Chunk Size\n(TP=8, CP=1, 8,192 bp Context)",
    x     = "Chunk Size",
    y     = "Run Time (minutes)",
    color = NULL
  ) +
  timing_theme +
  theme(
    legend.position  = "bottom",
    legend.text      = element_text(size = 11)
  )

save_plot(plot3, "timing_chunk_size_vs_time")

message("\n=== All timing plots saved to: ", output_dir, " ===")
