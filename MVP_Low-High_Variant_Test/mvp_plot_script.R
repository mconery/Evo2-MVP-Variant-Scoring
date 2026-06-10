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

#' Calculate percentage difference in median scores between High and Low PIP
#'
#' @param data Data frame with results
#' @return Data frame with median_high, median_low, pct_label per model_size/context_size
calculate_median_pct_diff <- function(data) {
  data %>%
    group_by(model_size, context_size) %>%
    summarise(
      median_high = median(evo2_delta_score[class == "High PIP"], na.rm = TRUE),
      median_low  = median(evo2_delta_score[class == "Low PIP"],  na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      pct_label = case_when(
        is.na(median_high) | is.na(median_low) ~ "Δmedian: N/A",
        abs(median_low) < 1e-10               ~ "Δmedian: N/A",
        TRUE ~ sprintf("Δmedian: %+.0f%%",
                       (median_high - median_low) / abs(median_low) * 100)
      )
    )
}

#' Compute ROC curve data per model/context facet
#'
#' High PIP is the positive class; lower evo2_delta_score predicts positive.
#'
#' @param data Data frame with model_size_factor, context_size_factor, class, evo2_delta_score
#' @return Long-format data frame with fpr, tpr, auc_roc, auc_label per facet
compute_roc_data <- function(data) {
  groups <- data %>% distinct(model_size_factor, context_size_factor)
  result_list <- vector("list", nrow(groups))

  for (i in seq_len(nrow(groups))) {
    ms_f <- as.character(groups$model_size_factor[i])
    cs_f <- as.character(groups$context_size_factor[i])

    grp <- data %>%
      filter(as.character(model_size_factor) == ms_f,
             as.character(context_size_factor) == cs_f)

    n       <- nrow(grp)
    classes <- unique(grp$class)

    degenerate <- n < 5 || length(classes) < 2
    if (!degenerate) {
      ord   <- order(grp$evo2_delta_score)
      y     <- as.integer(grp$class[ord] == "High PIP")
      n_pos <- sum(y)
      n_neg <- n - n_pos
      degenerate <- n_pos == 0 || n_neg == 0
    }

    if (degenerate) {
      result_list[[i]] <- data.frame(
        model_size_factor = ms_f, context_size_factor = cs_f,
        fpr = NA_real_, tpr = NA_real_,
        auc_roc = NA_real_, auc_label = "AUC: N/A",
        stringsAsFactors = FALSE
      )
      next
    }

    tpr_vals <- cumsum(y) / n_pos
    fpr_vals <- cumsum(1L - y) / n_neg

    # Deduplicate at tied-score boundaries
    boundary <- c(diff(grp$evo2_delta_score[ord]) != 0, TRUE)
    tpr_vals <- tpr_vals[boundary]
    fpr_vals <- fpr_vals[boundary]

    fpr_vals <- c(0, fpr_vals)
    tpr_vals <- c(0, tpr_vals)

    auc_roc <- sum(diff(fpr_vals) *
                   (tpr_vals[-1] + tpr_vals[-length(tpr_vals)]) / 2)

    result_list[[i]] <- data.frame(
      model_size_factor = ms_f, context_size_factor = cs_f,
      fpr = fpr_vals, tpr = tpr_vals,
      auc_roc = auc_roc,
      auc_label = sprintf("AUC = %.2f", auc_roc),
      stringsAsFactors = FALSE
    )
  }

  out <- bind_rows(result_list)
  out$model_size_factor   <- factor(out$model_size_factor,
                                    levels = levels(data$model_size_factor))
  out$context_size_factor <- factor(out$context_size_factor,
                                    levels = levels(data$context_size_factor))
  out
}

#' Compute precision-recall curve data per model/context facet
#'
#' High PIP is the positive class; lower evo2_delta_score predicts positive.
#'
#' @param data Data frame with model_size_factor, context_size_factor, class, evo2_delta_score
#' @return Long-format data frame with recall, precision, baseline, auc_pr, auc_label per facet
compute_pr_data <- function(data) {
  groups <- data %>% distinct(model_size_factor, context_size_factor)
  result_list <- vector("list", nrow(groups))

  for (i in seq_len(nrow(groups))) {
    ms_f <- as.character(groups$model_size_factor[i])
    cs_f <- as.character(groups$context_size_factor[i])

    grp <- data %>%
      filter(as.character(model_size_factor) == ms_f,
             as.character(context_size_factor) == cs_f)

    n       <- nrow(grp)
    classes <- unique(grp$class)
    n_pos   <- sum(grp$class == "High PIP")

    if (n < 5 || length(classes) < 2 || n_pos == 0) {
      result_list[[i]] <- data.frame(
        model_size_factor = ms_f, context_size_factor = cs_f,
        recall = NA_real_, precision = NA_real_,
        baseline = NA_real_, auc_pr = NA_real_, auc_label = "AUC: N/A",
        stringsAsFactors = FALSE
      )
      next
    }

    ord            <- order(grp$evo2_delta_score)
    y              <- as.integer(grp$class[ord] == "High PIP")
    recall_vals    <- cumsum(y) / n_pos
    precision_vals <- cumsum(y) / seq_len(n)

    recall_vals    <- c(0, recall_vals)
    precision_vals <- c(1, precision_vals)

    baseline <- n_pos / n
    auc_pr   <- sum(diff(recall_vals) *
                    (precision_vals[-1] + precision_vals[-length(precision_vals)]) / 2)

    result_list[[i]] <- data.frame(
      model_size_factor = ms_f, context_size_factor = cs_f,
      recall = recall_vals, precision = precision_vals,
      baseline = baseline, auc_pr = auc_pr,
      auc_label = sprintf("AUC-PR = %.2f", auc_pr),
      stringsAsFactors = FALSE
    )
  }

  out <- bind_rows(result_list)
  out$model_size_factor   <- factor(out$model_size_factor,
                                    levels = levels(data$model_size_factor))
  out$context_size_factor <- factor(out$context_size_factor,
                                    levels = levels(data$context_size_factor))
  out
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
      context_size_factor = factor(ifelse(context_size!=1000000, paste0(prettyNum(context_size, big.mark=","), " bp"), "1 Mb"), levels = ifelse(sort(unique(context_size))!=1000000, paste0(prettyNum(sort(unique(context_size)), big.mark=","), " bp"), "1 Mb"))
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

  # Calculate Wilcoxon test p-values and median percentage differences
  message("Calculating Wilcoxon test p-values...")
  pvalues <- calculate_wilcoxon_pvalues(data)

  message("Calculating median percentage differences...")
  median_diffs <- calculate_median_pct_diff(data)

  # Build separate annotation dataframes for top (p-value) and bottom (median diff)
  annotations <- pvalues %>%
    mutate(
      model_size_factor = factor(model_size, levels = c("1b", "7b", "7b Long-Context", "40b", "40b Long-Context")),
      context_size_factor = factor(ifelse(context_size!=1000000, paste0(prettyNum(context_size, big.mark=","), " bp"), "1 Mb"), levels = ifelse(sort(unique(context_size))!=1000000, paste0(prettyNum(sort(unique(context_size)), big.mark=","), " bp"), "1 Mb"))
    )

  median_annot <- median_diffs %>%
    mutate(
      model_size_factor = factor(model_size, levels = c("1b", "7b", "7b Long-Context", "40b", "40b Long-Context")),
      context_size_factor = factor(ifelse(context_size!=1000000, paste0(prettyNum(context_size, big.mark=","), " bp"), "1 Mb"), levels = ifelse(sort(unique(context_size))!=1000000, paste0(prettyNum(sort(unique(context_size)), big.mark=","), " bp"), "1 Mb"))
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
    # P-value pinned to the top of each facet panel
    geom_text(
      data = annotations,
      aes(x = 1.5, y = Inf, label = p_display),
      inherit.aes = FALSE,
      vjust = 1.5,
      size = 4.5,
      fontface = "bold"
    ) +
    # Median % difference pinned to the bottom of each facet panel
    geom_text(
      data = median_annot,
      aes(x = 1.5, y = -Inf, label = pct_label),
      inherit.aes = FALSE,
      vjust = -0.5,
      size = 4,
      fontface = "italic"
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
      plot.margin = margin(5, 20, 10, 5, "mm"),
      axis.text = element_text(size = 12),
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

#' Create faceted ROC curve plot for High PIP vs Low PIP discrimination
#'
#' @param data Data frame with results
#' @param coding_filter Optional: 'coding', 'non-coding', or NULL for all variants
#' @param title Plot title (auto-generated if NULL)
#' @param long_context_only Restrict to long-context models
#' @param max_context_size Upper cap on context size
#' @return ggplot object
create_faceted_roc_plot <- function(data, coding_filter = NULL, title = NULL,
                                    long_context_only = FALSE,
                                    max_context_size = NULL) {
  if (!is.null(coding_filter)) {
    data <- data %>% filter(Coding == ifelse(coding_filter == "coding", 1, 0))
    if (nrow(data) == 0) stop("No data found for coding_filter: ", coding_filter)
    if (is.null(title)) {
      title <- paste0("ROC Curves — Evo2 Delta Score\n(",
                      str_to_title(coding_filter), " Variants)")
    }
  } else {
    if (is.null(title)) title <- "ROC Curves — Evo2 Delta Score\n(All Variants)"
  }

  data <- data %>%
    mutate(
      model_size = str_replace_all(str_replace_all(str_replace_all(str_replace_all(model_size, "_", " "), "arc ", ""), "longcontext", "Long-Context"), "Long-Context", "LC"),
      model_size_factor = factor(model_size, levels = c("1b", "7b", "7b LC", "40b", "40b LC")),
      context_size_factor = factor(ifelse(context_size!=1000000, paste0(prettyNum(context_size, big.mark=","), " bp"), "1 Mb"), levels = ifelse(sort(unique(context_size))!=1000000, paste0(prettyNum(sort(unique(context_size)), big.mark=","), " bp"), "1 Mb"))
    )

  if (long_context_only) {
    data <- data %>% filter(grepl("LC", model_size))
    if (nrow(data) == 0) stop("No long-context model data found")
  }

  if (!is.null(max_context_size)) {
    data <- data %>% filter(context_size <= max_context_size)
    if (nrow(data) == 0) stop("No data remaining after context size filter")
  }

  roc_df <- compute_roc_data(data)

  auc_df <- roc_df %>%
    group_by(model_size_factor, context_size_factor) %>%
    slice(1) %>%
    ungroup() %>%
    select(model_size_factor, context_size_factor, auc_label)

  ggplot(roc_df, aes(x = fpr, y = tpr)) +
    geom_line(color = "#08519c", linewidth = 0.9, na.rm = TRUE) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                color = "grey50", linewidth = 0.6) +
    geom_text(
      data = auc_df,
      aes(x = 0.65, y = 0.08, label = auc_label),
      inherit.aes = FALSE, size = 4, fontface = "bold"
    ) +
    facet_grid(context_size_factor ~ model_size_factor) +
    coord_fixed() +
    scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1),
                       labels = c("0", "0.5", "1")) +
    scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1),
                       labels = c("0", "0.5", "1")) +
    labs(title = title, x = "False Positive Rate", y = "True Positive Rate") +
    theme_bw() +
    theme(
      plot.title       = element_text(hjust = 0.5, face = "bold", size = 18),
      plot.margin      = margin(5, 20, 5, 5, "mm"),
      axis.text        = element_text(size = 12),
      axis.title       = element_text(size = 13),
      strip.text.x     = element_text(face = "bold", size = 13),
      strip.text.y     = element_text(face = "bold", size = 12, angle = 0),
      legend.position  = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
}

#' Create faceted precision-recall curve plot for High PIP vs Low PIP discrimination
#'
#' @param data Data frame with results
#' @param coding_filter Optional: 'coding', 'non-coding', or NULL for all variants
#' @param title Plot title (auto-generated if NULL)
#' @param long_context_only Restrict to long-context models
#' @param max_context_size Upper cap on context size
#' @return ggplot object
create_faceted_pr_plot <- function(data, coding_filter = NULL, title = NULL,
                                   long_context_only = FALSE,
                                   max_context_size = NULL) {
  if (!is.null(coding_filter)) {
    data <- data %>% filter(Coding == ifelse(coding_filter == "coding", 1, 0))
    if (nrow(data) == 0) stop("No data found for coding_filter: ", coding_filter)
    if (is.null(title)) {
      title <- paste0("Precision-Recall Curves — Evo2 Delta Score\n(",
                      str_to_title(coding_filter), " Variants)")
    }
  } else {
    if (is.null(title)) title <- "Precision-Recall Curves — Evo2 Delta Score\n(All Variants)"
  }

  data <- data %>%
    mutate(
      model_size = str_replace_all(str_replace_all(str_replace_all(model_size, "_", " "), "arc ", ""), "longcontext", "Long-Context"),
      model_size_factor = factor(model_size, levels = c("1b", "7b", "7b Long-Context", "40b", "40b Long-Context")),
      context_size_factor = factor(ifelse(context_size!=1000000, paste0(prettyNum(context_size, big.mark=","), " bp"), "1 Mb"), levels = ifelse(sort(unique(context_size))!=1000000, paste0(prettyNum(sort(unique(context_size)), big.mark=","), " bp"), "1 Mb"))
    )

  if (long_context_only) {
    data <- data %>% filter(grepl("Long-Context", model_size))
    if (nrow(data) == 0) stop("No long-context model data found")
  }

  if (!is.null(max_context_size)) {
    data <- data %>% filter(context_size <= max_context_size)
    if (nrow(data) == 0) stop("No data remaining after context size filter")
  }

  pr_df <- compute_pr_data(data)

  auc_df <- pr_df %>%
    group_by(model_size_factor, context_size_factor) %>%
    slice(1) %>%
    ungroup() %>%
    select(model_size_factor, context_size_factor, auc_label)

  baseline_df <- pr_df %>%
    group_by(model_size_factor, context_size_factor) %>%
    slice(1) %>%
    ungroup() %>%
    select(model_size_factor, context_size_factor, baseline)

  ggplot(pr_df, aes(x = recall, y = precision)) +
    geom_line(color = "#08519c", linewidth = 0.9, na.rm = TRUE) +
    geom_hline(
      data = baseline_df,
      aes(yintercept = baseline),
      linetype = "dashed", color = "firebrick", linewidth = 0.6,
      na.rm = TRUE
    ) +
    geom_text(
      data = auc_df,
      aes(x = 0.5, y = 0.05, label = auc_label),
      inherit.aes = FALSE, size = 4, fontface = "bold"
    ) +
    facet_grid(context_size_factor ~ model_size_factor) +
    scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1),
                       labels = c("0", "0.5", "1")) +
    scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1),
                       labels = c("0", "0.5", "1")) +
    labs(title = title, x = "Recall", y = "Precision") +
    theme_bw() +
    theme(
      plot.title       = element_text(hjust = 0.5, face = "bold", size = 18),
      plot.margin      = margin(5, 20, 5, 5, "mm"),
      axis.text        = element_text(size = 12),
      axis.title       = element_text(size = 13),
      strip.text.x     = element_text(face = "bold", size = 13),
      strip.text.y     = element_text(face = "bold", size = 12, angle = 0),
      legend.position  = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
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

    # ROC curves
    message("Creating ROC curve plots", lc_label, "...")
    tryCatch({
      save_plot(
        create_faceted_roc_plot(results_data,
                                long_context_only = lc_only,
                                max_context_size  = max_context),
        paste0("evo2_roc_all_variants", lc_suffix),
        width = plot_width, height = plot_height,
        format = output_format, output_directory = output_dir
      )
    }, error = function(e) {
      message("Error creating ROC (all variants): ", e$message)
    })

    for (coding_option in c("coding", "non-coding")) {
      tryCatch({
        save_plot(
          create_faceted_roc_plot(results_data, coding_filter = coding_option,
                                  long_context_only = lc_only,
                                  max_context_size  = max_context),
          paste0("evo2_roc_", tolower(gsub("[^A-Za-z0-9_]", "_", coding_option)), lc_suffix),
          width = plot_width, height = plot_height,
          format = output_format, output_directory = output_dir
        )
      }, error = function(e) {
        message("Error creating ROC (", coding_option, "): ", e$message)
      })
    }

    # Precision-recall curves
    message("Creating precision-recall curve plots", lc_label, "...")
    tryCatch({
      save_plot(
        create_faceted_pr_plot(results_data,
                               long_context_only = lc_only,
                               max_context_size  = max_context),
        paste0("evo2_pr_all_variants", lc_suffix),
        width = plot_width, height = plot_height,
        format = output_format, output_directory = output_dir
      )
    }, error = function(e) {
      message("Error creating PR (all variants): ", e$message)
    })

    for (coding_option in c("coding", "non-coding")) {
      tryCatch({
        save_plot(
          create_faceted_pr_plot(results_data, coding_filter = coding_option,
                                 long_context_only = lc_only,
                                 max_context_size  = max_context),
          paste0("evo2_pr_", tolower(gsub("[^A-Za-z0-9_]", "_", coding_option)), lc_suffix),
          width = plot_width, height = plot_height,
          format = output_format, output_directory = output_dir
        )
      }, error = function(e) {
        message("Error creating PR (", coding_option, "): ", e$message)
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
