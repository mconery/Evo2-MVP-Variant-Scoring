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
output_dir <- "C:/Users/mitch/Documents/Argonne/Variant_Scoring"

# VEP annotations mapping file (located in output directory)
vep_annotations_file <- paste(output_dir, "vep_annotations.txt", sep = "/")

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
        p_value < 0.001 ~ formatC(p_value, format = "e", digits = 2),
        p_value < 0.01 ~ formatC(p_value, format = "e", digits = 2),
        p_value < 0.05 ~ formatC(p_value, format = "e", digits = 2),
        TRUE ~ formatC(p_value, format = "e", digits = 2)
      ),
      significance = case_when(
        is.na(p_value) ~ "",
        p_value < 0.001 ~ "***",
        p_value < 0.01 ~ "**",
        p_value < 0.05 ~ "*",
        TRUE ~ "ns"
      ),
      p_display = paste0(p_label, " ", significance)
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
                                   output_directory = NULL, coding_filter = NULL, sigma = pseudolog_sigma) {
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
      title <- paste0("Evo2 Delta Scores by Model and Context Size\n(", capitalize(coding_filter), " Variants)")
    }
  }else {
    if (is.null(title)) {
      title <- "Evo2 Delta Scores by Model and Context Size\n(All Variants)"
    }
  } 
  
  # Create factor for model size with proper ordering
  data <- data %>%
    mutate(
      model_size_factor = factor(model_size, levels = c("1b", "7b", "40b")),
      context_size_factor = factor(paste0(context_size, " bp"), levels = paste0(sort(unique(context_size)), " bp"))
    )
  
  # Calculate Wilcoxon test p-values
  message("Calculating Wilcoxon test p-values...")
  pvalues <- calculate_wilcoxon_pvalues(data)
  
  # Add factor columns to pvalues for merging
  pvalues <- pvalues %>%
    mutate(
      model_size_factor = factor(model_size, levels = c("1b", "7b", "40b")),
      context_size_factor = factor(paste0(context_size, " bp"), levels = paste0(sort(unique(context_size)), " bp"))
    )
  
  # Print p-values
  message("\nWilcoxon test results:")
  print(pvalues %>% select(model_size, context_size, p_value, significance))
  
  # Count samples for each facet combination
  counts <- data %>%
    group_by(model_size_factor, context_size_factor, class) %>%
    summarise(n = n(), .groups = "drop")
  
  # Calculate y-position for p-value labels in pseudo-log space
  # Position labels at top of transformed data range
  y_max <- max(data$evo2_delta_score, na.rm = TRUE)
  y_min <- min(data$evo2_delta_score, na.rm = TRUE)
  
  # For pseudo-log transformed space, we need to be careful with label positioning
  # Place label slightly above the maximum value
  if (y_max > 0) {
    label_y <- y_max * 1.3
  } else {
    label_y <- y_max * 0.7
  }
  
  # Create the plot with pseudo-log scale
  p <- ggplot(data, aes(x = class, y = evo2_delta_score, fill = class)) +
    geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.3) +
    # Add p-value labels to each facet
    geom_text(
      data = pvalues,
      aes(x = 1.5, y = label_y, label = p_display),
      inherit.aes = FALSE,
      size = 3,
      fontface = "bold"
    ) +
    # Apply pseudo-log scale transformation to y-axis
    scale_y_continuous(
      trans = pseudo_log_trans(sigma = sigma),
      labels = label_number(accuracy = 0.0001)
    ) +
    facet_grid(context_size_factor ~ model_size_factor,
               labeller = labeller(
                 context_size_factor = function(x) paste("Context:", x),
                 model_size_factor = function(x) paste("Model:", x)
               )) +
    labs(
      title = title,
      x = "PIP Class",
      y = "Evo2 Delta Score (pseudo-log scale)",
      fill = "PIP Class"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(face = "bold", size = 10),
      legend.position = "bottom"
    ) +
    scale_fill_brewer(palette = "Set2")
  
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
  
  # Plot 1: Overall results (all variants)
  message("Creating plot for all variants...")
  plot_all <- create_faceted_boxplot(results_data, vep_filter = NULL, output_directory = output_dir)
  save_plot(plot_all, "evo2_scores_all_variants", 
            width = plot_width, height = plot_height, format = output_format, output_directory = output_dir)
  
  # Plot 2-3: Plots for coding/non-coding variants
  message("Creating plot for coding/non-coding variants...")
  for (coding_option in c("coding", "non-coding")) {
    # Create safe filename
    safe_filename <- gsub("[^A-Za-z0-9_]", "_", coding_option)
    safe_filename <- paste0("evo2_scores_", tolower(safe_filename))
    
    #Make plot
    tryCatch({
      plot_group <- create_faceted_boxplot(results_data, coding_filter = coding_option, output_directory = output_dir)
      save_plot(plot_group, safe_filename, 
                width = plot_width, height = plot_height, format = output_format, output_directory = output_dir)
    }, error = function(e) {
      message("Error creating plot for ", coding_option, ": ", e$message)
    })
  }
  
  # Plot 4-N: Individual plots for each grouped annotation
  grouped_annotations <- unique(results_data$`Grouped Annotation`[!is.na(results_data$`Grouped Annotation`)])
  
  for (annotation in sort(grouped_annotations)) {
    message("\nCreating plot for: ", annotation)
    
    # Create safe filename
    safe_filename <- gsub("[^A-Za-z0-9_]", "_", annotation)
    safe_filename <- paste0("evo2_scores_", tolower(safe_filename))
    
    tryCatch({
      plot_group <- create_faceted_boxplot(results_data, vep_filter = annotation, output_directory = output_dir)
      save_plot(plot_group, safe_filename, 
                width = plot_width, height = plot_height, format = output_format, output_directory = output_dir)
    }, error = function(e) {
      message("Error creating plot for ", annotation, ": ", e$message)
    })
  }
  
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
