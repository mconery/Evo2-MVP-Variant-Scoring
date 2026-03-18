# ============================================================================
# R Script to Plot Evo2 Variant Scoring Results
# ============================================================================
# This script reads multiple Evo2 variant scoring results files, combines them,
# maps VEP annotations to grouped categories, and creates faceted boxplots
# comparing score distributions between Low PIP and High PIP classes.
# Includes Wilcoxon rank-sum test p-values for each facet.
# Uses pseudo-log scale for better visualization of small differences.
# ============================================================================

# Load required libraries
library(tidyverse)
library(ggplot2)
library(ggpubr)  # For statistical annotations
library(scales)  # For pseudo-log transformation

# ============================================================================
# CONFIGURATION
# ============================================================================

# Directory containing the results files
results_dir <- "C:/Users/mitch/Documents/Argonne/Variant_Scoring/scored_MVP_results"

# Output directory for plots
output_dir <- "C:/Users/mitch/Documents/Argonne/Variant_Scoring/scored_MVP_figures"

# VEP annotations mapping file (located in output directory)
vep_annotations_file <- "C:/Users/mitch/Documents/Argonne/Variant_Scoring/vep_annotations.txt"

# Conservation scores file (phastCons, phyloP, GERP)
conservation_file <- "C:/Users/mitch/Documents/Argonne/Variant_Scoring/MVP_conservation_scores.csv"

# Output plot settings
plot_width <- 12
plot_height <- 8
output_format <- "png"  # Can be "png", "pdf", or "svg"

# Pseudo-log scale parameter (sigma)
# This determines the linear region around zero
# Smaller values = more logarithmic behavior near zero
# Typical range: 0.00001 to 0.0001 for your data
pseudolog_sigma <- 0.00005

# ============================================================================
# FUNCTION DEFINITIONS
# ============================================================================

#' Load and parse all results files from the directory
#'
#' @param dir_path Path to directory containing results files
#' @return A combined data frame with all results and metadata
load_results_files <- function(dir_path) {
  # Get all files matching the naming pattern
  file_pattern <- "MVP_variant_scores\\.(\\w+)_model\\.(\\d+)bp_context\\.csv$"
  files <- list.files(dir_path, pattern = file_pattern, full.names = TRUE)
  
  if (length(files) == 0) {
    stop("No files found matching the expected naming pattern in: ", dir_path)
  }
  
  message("Found ", length(files), " results files")
  
  # Read and combine all files
  combined_data <- map_dfr(files, function(file) {
    # Extract model size and context size from filename
    filename <- basename(file)
    matches <- str_match(filename, "MVP_variant_scores\\.(\\w+)_model\\.(\\d+)bp_context")
    
    model_size <- matches[1, 2]
    context_size <- as.numeric(matches[1, 3])
    
    message("Reading: ", filename, " (Model: ", model_size, ", Context: ", context_size, "bp)")
    
    # Read the file
    data <- read_csv(file, show_col_types = FALSE)
    
    # Add metadata columns
    data <- data %>%
      mutate(
        model_size = model_size,
        context_size = context_size
      )
    
    return(data)
  })
  
  return(combined_data)
}

#' Load VEP annotations mapping and merge with results
#'
#' @param data Results data frame
#' @param vep_file Path to VEP annotations mapping file
#' @return Data frame with grouped annotations added
add_grouped_annotations <- function(data, vep_file) {
  # Check if VEP file exists
  if (!file.exists(vep_file)) {
    stop("VEP annotations file not found: ", vep_file)
  }
  
  # Read VEP mapping file
  vep_mapping <- read_tsv(vep_file, show_col_types = FALSE)
  
  # Rename columns for clarity
  vep_mapping <- vep_mapping %>%
    select(`VEP Annotation`, `Grouped Annotation`, Priority, Coding)
  
  # Merge with results data
  data_with_groups <- data %>%
    left_join(vep_mapping, by = c("VEP Annotation" = "VEP Annotation"))
  
  # Check for any unmapped annotations
  unmapped <- data_with_groups %>%
    filter(is.na(`Grouped Annotation`)) %>%
    pull(`VEP Annotation`) %>%
    unique()
  
  if (length(unmapped) > 0) {
    warning("The following VEP annotations were not mapped: ", paste(unmapped, collapse = ", "))
  }
  
  return(data_with_groups)
}

#' Calculate Wilcoxon test p-values for each facet combination
#'
#' @param data Data frame with results
#' @return Data frame with p-values for each model_size/context_size combination
calculate_wilcoxon_pvalues <- function(data) {
  # Calculate p-values for each combination
  pvalues <- data %>%
    group_by(model_size, context_size) %>%
    summarise(
      p_value = tryCatch({
        # Check if we have both classes
        classes <- unique(class)
        if (length(classes) < 2) {
          return(NA)
        }
        
        # Get scores for each class
        low_pip <- evo2_delta_score[class == "Low PIP"]
        high_pip <- evo2_delta_score[class == "High PIP"]
        
        # Need at least 2 observations in each group
        if (length(low_pip) < 2 | length(high_pip) < 2) {
          return(NA)
        }
        
        # Perform Wilcoxon test
        test_result <- wilcox.test(low_pip, high_pip, exact = FALSE, alternative = "greater")
        test_result$p.value
      }, error = function(e) {
        return(NA)
      }),
      .groups = "drop"
    )
  
  # Format p-values for display
  pvalues <- pvalues %>%
    mutate(
      p_label = case_when(
        is.na(p_value) ~ "N/A",
        TRUE ~ formatC(p_value, format = "e", digits = 1)
      ),
      significance = case_when(
        is.na(p_value) ~ "",
        p_value < 0.001 ~ "***",
        p_value < 0.01 ~ "**",
        p_value < 0.05 ~ "*",
        TRUE ~ "ns"
      ),
      p_display = paste0("p=", p_label)
    )
  
  return(pvalues)
}

#' Create faceted boxplot of Evo2 delta scores with Wilcoxon test p-values
#' Uses pseudo-log scale for y-axis
#'
#' @param data Data frame with results
#' @param vep_filter Optional grouped annotation to filter by (NULL for all)
#' @param coding_filter Optional coding filter: 'coding', 'non-coding', or NULL for all
#' @param title Plot title
#' @param output_directory Directory to save plots
#' @param sigma Pseudo-log scale parameter (controls linear region near zero)
#' @return ggplot object
create_faceted_boxplot <- function(data, vep_filter = NULL, title = NULL,
                                   output_directory = NULL, coding_filter = NULL,
                                   long_context_only = FALSE, sigma = pseudolog_sigma,
                                   max_context_size = NULL) {
  # Filter by VEP annotation if specified
  if (!is.null(vep_filter)) {
    data <- data %>% filter(`Grouped Annotation` == vep_filter)
    
    if (nrow(data) == 0) {
      stop("No data found for grouped annotation: ", vep_filter)
    }
    
    if (is.null(title)) {
      title <- paste0("Evo2 Delta Scores by Model and Context Size\n(", vep_filter, ")")
    }
  } else if (!is.null(coding_filter)){
    data <- data %>% filter(`Coding` == ifelse(coding_filter == "coding", 1, 0))
    
    if (nrow(data) == 0) {
      stop("No data found for coding_filter: ", coding_filter)
    }
    
    if (is.null(title)) {
      title <- paste0("Evo2 Delta Scores by Model and Context Size\n(", str_to_title(coding_filter), " Variants)")
    }
  }else {
    if (is.null(title)) {
      title <- "Evo2 Delta Scores by Model and Context Size\n(All Variants)"
    }
  } 
  
  # Create factor for model size with proper ordering
  data <- data %>%
    mutate(
      model_size = str_replace_all(str_replace_all(str_replace_all(model_size, "_", " "), "arc ", ""), "longcontext", "Long-Context"),
      model_size_factor = factor(model_size, levels = c("1b", "7b", "7b Long-Context", "40b", "40b Long-Context")),
      context_size_factor = factor(ifelse(context_size!=1000000, paste0(context_size, " bp"), "1 Mb"), levels = ifelse(sort(unique(context_size))!=1000000, paste0(sort(unique(context_size)), " bp"), "1 Mb"))
    )

  # Optionally restrict to long-context models only
  if (long_context_only) {
    data <- data %>% filter(grepl("Long-Context", model_size))
    if (nrow(data) == 0) stop("No long-context model data found")
  }

  # Optionally cap context size
  if (!is.null(max_context_size)) {
    data <- data %>% filter(context_size <= max_context_size)
    if (nrow(data) == 0) stop("No data remaining after context size filter (<= ", max_context_size, ")")
  }

  # Calculate Wilcoxon test p-values
  message("Calculating Wilcoxon test p-values...")
  pvalues <- calculate_wilcoxon_pvalues(data)
  
  # Add factor columns to pvalues for merging
  pvalues <- pvalues %>%
    mutate(
      model_size_factor = factor(model_size, levels = c("1b", "7b", "7b Long-Context", "40b", "40b Long-Context")),
      context_size_factor = factor(ifelse(context_size!=1000000, paste0(context_size, " bp"), "1 Mb"), levels = ifelse(sort(unique(context_size))!=1000000, paste0(sort(unique(context_size)), " bp"), "1 Mb"))
    )
  
  # Print p-values
  message("\nWilcoxon test results:")
  print(pvalues %>% select(model_size, context_size, p_value, significance))
  
  # Count samples for each facet combination
  counts <- data %>%
    group_by(model_size_factor, context_size_factor, class) %>%
    summarise(n = n(), .groups = "drop")
  
  # Create the plot with pseudo-log scale
  p <- ggplot(data, aes(x = class, y = evo2_delta_score, fill = class)) +
    geom_violin() +
    stat_summary(
      fun = median, fun.min = median, fun.max = median,
      geom = "crossbar", width = 0.3, color = "black", linewidth = 0.8
    ) +
    # Add p-value labels pinned to top of each facet panel
    geom_text(
      data = pvalues,
      aes(x = 1.5, y = Inf, label = p_display),
      inherit.aes = FALSE,
      vjust = 1.5,
      size = 4.5,
      fontface = "bold"
    ) +
    # Apply pseudo-log scale transformation to y-axis, showing only 0
    scale_y_continuous(
      trans = pseudo_log_trans(sigma = sigma),
      breaks = 0,
      labels = "0"
    ) +
    facet_grid(context_size_factor ~ model_size_factor, scales = "free") +
    labs(
      title = title,
      x = NULL,
      y = "Evo2 Delta Score (pseudo-log scale)"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
      plot.margin = margin(5, 20, 5, 5, "mm"),
      axis.text = element_text(size = 12),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title = element_text(size = 13),
      strip.text.x = element_text(face = "bold", size = 13),
      strip.text.y = element_text(face = "bold", size = 12, angle = 0),
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    scale_fill_manual(values = c("Low PIP" = "#6baed6", "High PIP" = "#08519c"))
  
  # Print sample counts
  message("\nSample counts per facet:")
  print(counts)
  
  return(p)
}

#' Save plot to file in specified directory
#'
#' @param plot ggplot object
#' @param filename Output filename (without extension)
#' @param width Plot width in inches
#' @param height Plot height in inches
#' @param format Output format (png, pdf, svg)
#' @param output_directory Directory to save the plot
save_plot <- function(plot, filename, width = 12, height = 8, format = "png", output_directory = output_dir) {
  # Create output directory if it doesn't exist
  if (!dir.exists(output_directory)) {
    dir.create(output_directory, recursive = TRUE)
    message("Created output directory: ", output_directory)
  }
  
  # Create full output path
  output_file <- file.path(output_directory, paste0(filename, ".", format))
  
  ggsave(
    filename = output_file,
    plot = plot,
    width = width,
    height = height,
    dpi = 300
  )
  
  message("Plot saved to: ", output_file)
}

#' Build a wide-format collated results table from all model/context combinations
#'
#' @param data Combined data frame from load_results_files()
#' @return Wide data frame with 13 metadata columns + 3 score columns per combo
build_collated_table <- function(data) {
  base_cols <- c("MVP ID", "RSID", "BP", "BP38", "VEP Annotation", "CHR",
                 "EAF Population", "Beta Population", "P-Value Population",
                 "Overall PIP", "CS-Level Pip", "mu", "Set")

  # One copy of variant metadata (identical across all files)
  base <- data %>%
    filter(model_size == first(model_size), context_size == first(context_size)) %>%
    select(all_of(base_cols))

  score_cols <- c("ref_log_probs", "var_log_probs", "evo2_delta_score")

  combos <- data %>%
    distinct(model_size, context_size) %>%
    arrange(model_size, context_size) %>%
    filter(
      (grepl("longcontext", model_size) & context_size <= 524288) |
      (!grepl("longcontext", model_size) & context_size <= 131072)
    )

  message("Collating ", nrow(combos), " model/context combinations")

  for (i in seq_len(nrow(combos))) {
    ms <- combos$model_size[i]
    cs <- combos$context_size[i]
    suffix <- paste0("_", ms, "_", cs, "bp")

    scores <- data %>%
      filter(model_size == ms, context_size == cs) %>%
      select(`MVP ID`, all_of(score_cols)) %>%
      rename_with(~ paste0(.x, suffix), all_of(score_cols))

    base <- base %>% left_join(scores, by = "MVP ID")
  }

  return(base)
}

# ============================================================================
# MAIN SCRIPT
# ============================================================================

main <- function() {
  message("=== Evo2 Variant Scoring Analysis ===\n")
  
  # Step 1: Load all results files
  message("Step 1: Loading results files...")
  results_data <- load_results_files(results_dir)
  message("Total variants loaded: ", nrow(results_data))
  message("Unique model sizes: ", paste(unique(results_data$model_size), collapse = ", "))
  message("Unique context sizes: ", paste(unique(results_data$context_size), collapse = ", "))
  message("Unique classes: ", paste(unique(results_data$class), collapse = ", "))
  
  # Step 2: Add grouped VEP annotations
  message("\nStep 2: Adding grouped VEP annotations...")
  message("VEP annotations file: ", vep_annotations_file)
  results_data <- add_grouped_annotations(results_data, vep_annotations_file)
  message("Grouped annotations available: ", 
          paste(sort(unique(results_data$`Grouped Annotation`[!is.na(results_data$`Grouped Annotation`)])), 
                collapse = ", "))
  
  # Step 3: Create plots
  message("\nStep 3: Creating plots...\n")
  message("Output directory: ", output_dir)
  message("Using pseudo-log scale with sigma = ", pseudolog_sigma)

  for (lc_only in c(FALSE, TRUE)) {
    lc_suffix <- if (lc_only) "_long_context" else ""
    lc_label  <- if (lc_only) " (long-context models only)" else ""
    max_context <- if (lc_only) 524288 else 131072

    # All variants
    message("Creating plot for all variants", lc_label, "...")
    plot_all <- create_faceted_boxplot(
      results_data, vep_filter = NULL,
      long_context_only = lc_only, output_directory = output_dir,
      max_context_size = max_context
    )
    save_plot(plot_all, paste0("evo2_scores_all_variants", lc_suffix),
              width = plot_width, height = plot_height,
              format = output_format, output_directory = output_dir)

    # Coding / non-coding variants
    message("Creating coding/non-coding plots", lc_label, "...")
    for (coding_option in c("coding", "non-coding")) {
      safe_filename <- paste0(
        "evo2_scores_",
        tolower(gsub("[^A-Za-z0-9_]", "_", coding_option)),
        lc_suffix
      )
      tryCatch({
        plot_group <- create_faceted_boxplot(
          results_data, coding_filter = coding_option,
          long_context_only = lc_only, output_directory = output_dir,
          max_context_size = max_context
        )
        save_plot(plot_group, safe_filename,
                  width = plot_width, height = plot_height,
                  format = output_format, output_directory = output_dir)
      }, error = function(e) {
        message("Error creating plot for ", coding_option, ": ", e$message)
      })
    }
  }
  
  # VEP-level plots for long-context models
  message("\nCreating VEP-level plots for long-context models...")
  vep_subdir <- file.path(output_dir, "vep-level_predictions")
  vep_annotations <- sort(unique(results_data$`Grouped Annotation`[!is.na(results_data$`Grouped Annotation`)]))
  message("VEP annotations to plot: ", paste(vep_annotations, collapse = ", "))

  for (vep_ann in vep_annotations) {
    safe_name <- tolower(gsub("[^A-Za-z0-9_]", "_", vep_ann))
    filename <- paste0("evo2_scores_vep_", safe_name, "_long_context")
    message("Creating VEP-level plot for: ", vep_ann)
    tryCatch({
      plot_vep <- create_faceted_boxplot(
        results_data,
        vep_filter        = vep_ann,
        long_context_only = TRUE,
        max_context_size  = 524288
      )
      save_plot(plot_vep, filename,
                width = plot_width, height = plot_height,
                format = output_format, output_directory = vep_subdir)
    }, error = function(e) {
      message("Skipping ", vep_ann, ": ", e$message)
    })
  }

  # Step 4: Build and save collated results table
  message("\nStep 4: Building collated results table...")
  collated <- build_collated_table(results_data)

  # Join conservation scores (phastCons, phyloP, GERP) by MVP ID
  message("Adding conservation scores from: ", conservation_file)
  conservation <- read_csv(conservation_file, show_col_types = FALSE) %>%
    select(`MVP ID`, phastCons100way, phyloP100way, GERP_RS)
  collated <- collated %>%
    left_join(conservation, by = "MVP ID")
  n_with_conservation <- sum(!is.na(collated$phastCons100way))
  message("Conservation scores joined: ", n_with_conservation, "/", nrow(collated),
          " variants with phastCons100way")

  collated_file <- file.path(output_dir, "evo2_scores_collated.csv")
  write_csv(collated, collated_file)
  message("Collated table saved to: ", collated_file)
  message("Dimensions: ", nrow(collated), " rows x ", ncol(collated), " columns")

  message("\n=== Analysis Complete ===")
  message("All plots have been saved to: ", output_dir)
  
  # Return the data for further analysis if needed
  invisible(results_data)
}

# ============================================================================
# EXECUTION
# ============================================================================

# Run the main function
results <- main()
