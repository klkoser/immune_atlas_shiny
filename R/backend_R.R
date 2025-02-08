# Backend Rscript for generating plots

# ---- Load packages ----
library(tidyverse)

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

# Plot Function
plot_median_variance <- function(stimulus = "TNFa", data_obj = data_obj, output_dir = "./plots/"){

  # Directory to save
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  plot_path <- file.path(output_dir, paste0("volcano_", stimulus, ".svg"))

  p <- data_obj %>%
    group_by(population, reagent, Condition) %>%
    summarise(median_diff = median(value_diff),
              median_var = median(variance)) %>% ungroup() %>%
    filter(Condition == stimulus) %>%
    ggplot(aes(x = median_diff, y = median_var)) +
    geom_point(aes(color = population), size = 2) +
    theme_classic() +
    scale_color_manual(values = cell_type_palette, name ="Cell Type") +
    labs(x = "Median Difference\n(From Basal)",
         y = "Median Variance",
         title = paste0("Median stimulation by variance for ", stimulus))

  ggsave(plot_path, p, width = 6, height = 4, dpi = 300, device = "svg")

  return(plot_path) # Return file path

}

# Call Function
#plot_median_variance(stimulus = "TNFa", data_obj = data_obj)


# ---- Tab 2: Generate Boxplots ----


### cell_type is the target cell type, selected from dropdown (such as Neutrophil)
### stimulus is the condition, selected from dropdown (such as TNFa)
### gender is a boolean, TRUE or FALSE
# TRUE will split the boxplot by gender
# FALSE will simply plot both together

generate_boxplot <- function(cell_type="CD4+T cells", stimulus = "TNFa", gender =T, output_dir = "./plots/") {

  # Check if directory for plots exists and make if not
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  plot_path <- file.path(output_dir, paste0("boxplot_", stimulus, ".svg"))

  # Subset data
  specific_data <- data_obj[data_obj$population == cell_type & data_obj$Condition == stimulus,]

  # Generate plot
  if (gender) {

    p_2 <- ggplot(specific_data, aes(x = reorder(reagent, -value_diff), y = value_diff, fill = Gender)) +
      geom_hline(yintercept = 0, color = "red") +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(size = .2, position = position_jitterdodge(jitter.width = .2)) +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(x = "Reagent Name", y = "Median Difference from Basal")
    #  stat_compare_means(method = "wilcox", label = "..p..")

    ggsave(plot_path, p_2, width = 6, height = 4, dpi = 300, device = "svg")

    return(plot_path) # Return file path


  } else {

    p_2 <- ggplot(specific_data, aes(x = reorder(reagent, -value_diff), y = value_diff)) +
      geom_hline(yintercept = 0, color = "red") +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(size = .2, width = 0.2) +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(x = "Reagent Name", y = "Median Difference from Basal")

    ggsave(plot_path, p_2, width = 6, height = 4, dpi = 300, device = "svg")

    return(plot_path) # Return file path


  }
}

# Call function
#generate_boxplot()
#generate_boxplot(cell_type = celltype,
#                stimulus = condition,
#                 gender = gender)

# ---- Tab 3: Generate Correlation Plots ----

# Define the correlation plot function
correlation_plot <- function(cell_type = "Neutrophils", stimulus = "TNFa", read1 = "pP38", read2 = "pErk1/2", output_dir = "./plots/") {

  # Ensure output directory exists
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  # Sanitize filename: Remove special characters (`/`, `+`)
  safe_cell_type <- gsub("[+/]", "_", cell_type)  # Replace `/` and `+` with `_`
  safe_read_1 <- gsub("[+/]", "_", read1)  # Replace `/` and `+` with `_`
  safe_read_2 <- gsub("[+/]", "_", read2)  # Replace `/` and `+` with `_`
  plot_path <- file.path(output_dir, paste0("correlation_", safe_cell_type, "_", stimulus, "_", safe_read_1, "_", safe_read_2, ".svg"))


  # Check if directory for plots exists and make if not
  # if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  # plot_path <- file.path(output_dir, paste0("correlation_", gsub("\\+", "", cell_type), "_", stimulus, ".svg"))
  #
  # Input validation - check that all parameters acutally exist in the data
  if (!cell_type %in% unique(data_obj$population)) {
    # throws error and stops exeuction if we run into invalid cell types, also print out the valid types (unique)so we know what we can use
    stop("Invalid cell type '", cell_type, "'. Available types:\n",
         paste(sort(unique(data_obj$population)), collapse="\n"))
  }
  if (!stimulus %in% unique(data_obj$Condition)) {
    stop("Invalid stimulus '", stimulus, "'. Available types:\n",
         paste(sort(unique(data_obj$Condition)), collapse="\n"))
  }
  if (!read1 %in% unique(data_obj$reagent) || !read2 %in% unique(data_obj$reagent)) {
    stop("Invalid reagent. Available types:\n",
         paste(sort(unique(data_obj$reagent)), collapse="\n"))
  }

  # Filter data for the specific cell type and condition
  data_obj_specific <- subset(data_obj, population == cell_type & Condition == stimulus)

  # Create separate dataframes for each reagent
  read1_data <- subset(data_obj_specific, reagent == read1)
  read2_data <- subset(data_obj_specific, reagent == read2)

  # Print debugging info
  cat(sprintf("Found %d donors for %s and %d donors for %s\n",
              length(unique(read1_data$Donor)), read1,
              length(unique(read2_data$Donor)), read2))

  # Match the data by Donor and Gender
  matched_data <- merge(
    read1_data[, c("Donor", "Gender", "value_diff")],
    read2_data[, c("Donor", "Gender", "value_diff")],
    by = c("Donor", "Gender"),
    suffixes = c("_read1", "_read2")
  )

  # Create plot data frame
  plot_data <- data.frame(
    read1_vals = matched_data$value_diff_read1,
    read2_vals = matched_data$value_diff_read2,
    Gender = matched_data$Gender
  )


  # Check if we have matching data points
  if (nrow(plot_data) < 3) {
    warning(sprintf("Only %d matching donors found. Need at least 3 for correlation.",
                    nrow(plot_data)))
    return(list(r_squared = NA, pearson = NA, spearman = NA, n_donors = nrow(plot_data)))
  }

  # Remove any NA values
  plot_data <- na.omit(plot_data)

  # Fit a linear regression model: read2_vals as dependent var, read1_vals as independent var
  corr_fit <- lm(read2_vals ~ read1_vals, data = plot_data)

  # Extract R-squared value from the linear model summary
  r_squared <- summary(corr_fit)$r.squared

  # Pearson corrrelation coef - linear corr between variables
  pearson_corr <- cor(plot_data$read1_vals, plot_data$read2_vals, method = "pearson")

  # Spearman rank correlation - meassures monotonic relatonship, robust to outliers
  spearman_corr <- cor(plot_data$read1_vals, plot_data$read2_vals, method = "spearman")

  # Create a ggplot object with plot_data as the data source, mapping read1_vals to x-axis and read2_vals to y-axis
  p <- ggplot(plot_data, aes(x = read1_vals, y = read2_vals)) +
    # Add scatter points: size=3 makes them visible, alpha=0.6 makes them slightly transparent
    geom_point(size = 3, alpha = 0.6) +
    # Add a linear regression line (method="lm") in red color with standrd error bands (se=TRUE)
    geom_smooth(method = "lm", se = TRUE, color = "red1") +
    #geom_smooth(method = "lm", se = TRUE, aes(group = Gender, color = Gender)) +
   # stat_cor(method = "pearson")
    # Add labels to the plot
    labs(
      # Main title
      title = paste("Correlation Plot:", cell_type, "-", stimulus),
      # Subtitle
      subtitle = paste("n =", nrow(plot_data), "donors"),

      x = paste(read1, "response \n(Median Difference from Basal)"),

      y = paste(read2, "response \n(Median Difference from Basal)")
    ) +
    # annotations
    annotate("text",
             x = -Inf,          # Place text at the lefttmost edge
             y = Inf,           # Place text at the topmost edge
             hjust = -0.1,      # Slight horizontal offset to prevvent touching the edge
             vjust = 1.5,       # Vertical offset to positon below the top edge
             # Format stat values with 3 decimals
             label = sprintf("R² = %.3f\nPearson = %.3f\nSpearman = %.3f",
                             r_squared, pearson_corr, spearman_corr)) +
    # Use a minimal theme
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold"),     # Make title bold
      plot.subtitle = element_text(color = "gray50"), # Make subtitle grey
      panel.grid.minor = element_blank()            # Remove minor grid lines
    )

  # Display the plot
  ggsave(plot_path, p, width = 6, height = 4, dpi = 300, device = "svg")

  return(plot_path)

  # Return statistics (without printing)
  invisible(list(
    r_squared = r_squared,     # R-squared value from linear regresssion
    pearson = pearson_corr,    # Pearson correlation coeficient
    spearman = spearman_corr,  # Spearman corrrelation coefficient
    n_donors = nrow(plot_data) # Number of donors in the analysis
  ))
}

# Call function
#correlation_plot(cell_type = celltype, stimulus = condition, read1 = reagent_1, read2 = reagent_2)
