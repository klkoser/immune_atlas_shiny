# Backend Rscript for generating plots

# ---- Load packages ----
library(tidyverse)
library(ggpubr)

# ---- Load datasets ----
data_obj <- read.csv("data/filtered_data.csv")

# ---- Tab 1: Generate "Volcano PLots" ----

# Color palette for cell type/population
cell_type_palette <- c(
  "B Cells" = "#E41A1C",       # Red
  "Basophils" = "#377EB8",     # Blue
  "CD11b-/CD16-" = "#4DAF4A",  # Green
  "CD11b+/CD16-" = "#984EA3",  # Purple
  "CD16+/CD11b-" = "#FF7F00",  # Orange
  "CD16+/CD11b+" = "#FFF010",  # Yellow
  "CD1c+ B cells" = "#A65628", # Brown
  "CD4+/CD8+" = "#F781BF",     # Pink
  "CD4+T cells" = "#999999",  # Gray
  "CD7+/HLA-DR-" = "#66C2A5",  # Teal
  "CD8+T cells" = "#FC8D62",  # Salmon
  "Neutrophils" = "#8DA0CB",   # Lavender
  "pDCs" = "#E78AC3"          # Magenta
)

#' Generate a Volcano Plot as a ggplot Object
#'
#' Calculates median difference and variance per population and generates a volcano plot.
#'
#' @param stimulus The stimulus/condition to filter by (default "TNFa").
#' @param data_obj A data frame containing the input data.
#' @param output_dir Directory where the plot will be saved if save_plot is TRUE.
#' @param save_plot Logical. If TRUE, saves the plot to disk; defaults to FALSE.
#' @return A ggplot object.
#' @export
plot_median_variance <- function(stimulus = "TNFa", data_obj = data_obj,
                                 output_dir = "./plots/", save_plot = FALSE) {

  p <- data_obj %>%
    group_by(population, reagent, Condition) %>%
    summarise(median_diff = median(value_diff),
              median_var = median(variance)) %>%
    ungroup() %>%
    filter(Condition == stimulus) %>%
    ggplot(aes(x = median_diff, y = median_var)) +
    geom_point(aes(color = population), size = 2) +
    theme_classic() +
    scale_color_manual(values = cell_type_palette, name = "Cell Type") +
    labs(x = "Median Difference\n(From Basal)",
         y = "Median Variance",
         title = paste0("Median stimulation by variance for ", stimulus))

  if (save_plot) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    plot_path <- file.path(output_dir, paste0("volcano_", stimulus, ".svg"))
    ggsave(plot_path, p, width = 6, height = 4, dpi = 300, device = "svg")
    message("Plot saved to: ", plot_path)
  }

  return(p)
}

# ---- Tab 2: Generate Boxplots ----



#' Generate a Boxplot as a ggplot Object
#'
#' Creates a boxplot for a given cell type and stimulus. Optionally splits the plot by gender.
#'
#' @param cell_type The target cell type (default "CD4+T cells").
#' @param stimulus The condition (default "TNFa").
#' @param gender Logical. If TRUE, splits the boxplot by gender; otherwise, plots both together.
#' @param output_dir Directory where the plot will be saved if save_plot is TRUE.
#' @param save_plot Logical. If TRUE, saves the plot to disk; defaults to FALSE.
#' @return A ggplot object.
#' @export
generate_boxplot <- function(cell_type = "CD4+T cells", stimulus = "TNFa", gender = TRUE,
                             output_dir = "./plots/", save_plot = FALSE) {

  # Subset data
  specific_data <- data_obj[data_obj$population == cell_type & data_obj$Condition == stimulus, ]

  if (gender) {
    # When splitting by gender, add significance using a Wilcoxon test.
    p <- ggplot(specific_data, aes(x = reorder(reagent, -value_diff), y = value_diff, fill = Gender)) +
      geom_hline(yintercept = 0, color = "red") +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(size = 0.2, position = position_jitterdodge(jitter.width = 0.2)) +
      stat_compare_means(method = "wilcox.test", label = "p.signif", color = "red") +  # significance annotations
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(x = "Reagent Name", y = "Median Difference from Basal")
  } else {
    p <- ggplot(specific_data, aes(x = reorder(reagent, -value_diff), y = value_diff)) +
      geom_hline(yintercept = 0, color = "red") +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(size = 0.2, width = 0.2) +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(x = "Reagent Name", y = "Median Difference from Basal")
  }

  if (save_plot) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    plot_path <- file.path(output_dir, paste0("boxplot_", stimulus, ".svg"))
    ggsave(plot_path, p, width = 6, height = 4, dpi = 300, device = "svg")
    message("Boxplot saved to: ", plot_path)
  }

  return(p)
}



# ---- Tab 3: Generate correlation plots ----


#' Generate a Correlation Plot as a ggplot Object
#'
#' Produces a correlation plot between two reagents for a specific cell type and stimulus.
#'
#' @param cell_type The cell type (default "Neutrophils").
#' @param stimulus The condition (default "TNFa").
#' @param read1 The first reagent (default "pP38").
#' @param read2 The second reagent (default "pErk1/2").
#' @param output_dir Directory where the plot will be saved if save_plot is TRUE.
#' @param save_plot Logical. If TRUE, saves the plot to disk; defaults to FALSE.
#' @return A ggplot object.
#' @export
correlation_plot <- function(cell_type = "Neutrophils", stimulus = "TNFa",
                             read1 = "pP38", read2 = "pErk1/2",
                             output_dir = "./plots/", save_plot = FALSE) {

  # Sanitize names for potential file saving
  safe_cell_type <- gsub("[+/]", "_", cell_type)
  safe_read_1 <- gsub("[+/]", "_", read1)
  safe_read_2 <- gsub("[+/]", "_", read2)

  # Input validation
  if (!cell_type %in% unique(data_obj$population)) {
    stop("Invalid cell type '", cell_type, "'. Available types:\n",
         paste(sort(unique(data_obj$population)), collapse = "\n"))
  }
  if (!stimulus %in% unique(data_obj$Condition)) {
    stop("Invalid stimulus '", stimulus, "'. Available types:\n",
         paste(sort(unique(data_obj$Condition)), collapse = "\n"))
  }
  if (!read1 %in% unique(data_obj$reagent) || !read2 %in% unique(data_obj$reagent)) {
    stop("Invalid reagent. Available types:\n",
         paste(sort(unique(data_obj$reagent)), collapse = "\n"))
  }

  # Filter data for the specified cell type and condition
  data_obj_specific <- subset(data_obj, population == cell_type & Condition == stimulus)

  # Separate data for each reagent
  read1_data <- subset(data_obj_specific, reagent == read1)
  read2_data <- subset(data_obj_specific, reagent == read2)

  cat(sprintf("Found %d donors for %s and %d donors for %s\n",
              length(unique(read1_data$Donor)), read1,
              length(unique(read2_data$Donor)), read2))

  # Merge data by Donor and Gender
  matched_data <- merge(
    read1_data[, c("Donor", "Gender", "value_diff")],
    read2_data[, c("Donor", "Gender", "value_diff")],
    by = c("Donor", "Gender"),
    suffixes = c("_read1", "_read2")
  )

  # Prepare data for plotting
  plot_data <- data.frame(
    read1_vals = matched_data$value_diff_read1,
    read2_vals = matched_data$value_diff_read2,
    Gender = matched_data$Gender
  )

  if (nrow(plot_data) < 3) {
    warning(sprintf("Only %d matching donors found. Need at least 3 for correlation.",
                    nrow(plot_data)))
    return(NULL)
  }

  plot_data <- na.omit(plot_data)

  # Fit a linear model and compute statistics
  corr_fit <- lm(read2_vals ~ read1_vals, data = plot_data)
  r_squared <- summary(corr_fit)$r.squared
  pearson_corr <- cor(plot_data$read1_vals, plot_data$read2_vals, method = "pearson")
  spearman_corr <- cor(plot_data$read1_vals, plot_data$read2_vals, method = "spearman")

  # Build the plot
  p <- ggplot(plot_data, aes(x = read1_vals, y = read2_vals)) +
    geom_point(size = 3, alpha = 0.6) +
    geom_smooth(method = "lm", se = TRUE, color = "red1") +
    labs(
      title = paste("Correlation Plot:", cell_type, "-", stimulus),
    #  subtitle = paste("n =", nrow(plot_data), "donors"),
      x = paste(read1, "response \n(Median Difference from Basal)"),
      y = paste(read2, "response \n(Median Difference from Basal)")
    ) +
    annotate("text",
             x = -Inf, y = Inf,
             hjust = -0.1, vjust = 1.5,
             label = sprintf("R² = %.3f\nPearson = %.3f\nSpearman = %.3f",
                             r_squared, pearson_corr, spearman_corr)) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(color = "gray50"),
      panel.grid.minor = element_blank()
    )

  if (save_plot) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    plot_path <- file.path(output_dir, paste0("correlation_", safe_cell_type, "_",
                                              stimulus, "_", safe_read_1, "_", safe_read_2, ".svg"))
    ggsave(plot_path, p, width = 6, height = 4, dpi = 300, device = "svg")
    message("Correlation plot saved to: ", plot_path)
  }

  return(p)
}


