# 1.Packages loading and parameters parsing ----------------------------------------------------------------
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(yaml))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(htmltools))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(gridExtra))

option_list <- list(
  make_option(c("-r", "--rds"), type="character",
              help="read_assignment_result"),
  make_option(c("-0", "--outdir"), type="character",
              help="outdir"),
  make_option(c("-a","--annotation"), type="character",
              help="Path of barcode annotation file")
)
option_list <- parse_args(OptionParser(option_list=option_list))

barcode_annotation_path <- option_list$annotation
outpath <- option_list$outdir
tmp_rds_path <- option_list$rds

## Initialize log file
logfile <- paste0(outpath, "/doublets_separation_plotting.processing_log.txt")
logwrite <- function(msg) {
  write(msg, file = logfile, append = TRUE)
}




# 2.Load files ---------------------------------------------------------------------------------------------

## load RDS file
dgecounts <- readRDS(tmp_rds_path)
## load barcode annotation
barcode_annotation <- read.table(barcode_annotation_path, header = TRUE, sep = "\t")

logwrite("Barcode annotation file loaded successfully.\n")
logwrite("--------------------First few lines of Barcode annotation file--------------------:\n")
logwrite(paste(head(barcode_annotation, n=3), collapse = "\n"))




# 3.Statistical analysis -----------------------------------------------------------------------------------

# Initialize an empty read_assignment_result data.table
read_assignment_result <- data.table(RG = character(), Tsc_count = integer(), Nalm6_count = integer(), other_count = integer())

# Iterate through all data.tables
for (dt in dgecounts) {
  # Use data.table's .N to count Tsc and Nalm6 occurrences
  counts <- dt[, .(Tsc_count = sum(readcall == "Tsc"), 
                   Nalm6_count = sum(readcall == "Nalm6"), 
                   other_count = sum(readcall == "other")),
                by = RG]
  
# Bind the current data.table's results to the overall result
  read_assignment_result <- rbind(read_assignment_result, counts, use.names = TRUE, fill = TRUE)
}

# Group by RG and sum the counts
final_read_assignment_result <- read_assignment_result[, .(Tsc_count = sum(Tsc_count, na.rm = TRUE), 
                            Nalm6_count = sum(Nalm6_count, na.rm = TRUE),
                            other_count = sum(other_count, na.rm = TRUE)), 
                        by = RG]
barcode_annotation_subset <- barcode_annotation[, c("condition", "RG")]
final_read_assignment_result_with_barcodeinfo <- merge(final_read_assignment_result, barcode_annotation_subset, by = "RG")
# Sort by 'condition' first, and then by 'Tsc_count'
final_read_assignment_result_sorted <- final_read_assignment_result_with_barcodeinfo[order(final_read_assignment_result_with_barcodeinfo$condition, 
                                                           final_read_assignment_result_with_barcodeinfo$Tsc_count), ]


# Save the final result as a text file
fwrite(final_read_assignment_result_sorted, paste0(outpath,"/final_read_assignment_result_sorted.txt"), sep = "\t")

#colnames(final_read_assignment_result_sorted)

# Initialize counters
total_obs <- 0
total_Tsc <- 0
total_Nalm6 <- 0
total_other <- 0

# Iterate through all data.tables
for (dt in dgecounts) {
  # Add the number of observations in the current data.table to total_obs
  total_obs <- total_obs + nrow(dt)
  
  # Count Tsc and Nalm6 in the current data.table
  total_Tsc <- total_Tsc + sum(dt$readcall == "Tsc")
  total_Nalm6 <- total_Nalm6 + sum(dt$readcall == "Nalm6")
  total_other <- total_other + sum(dt$readcall == "other")
}

# Print the results
logwrite("--------------------Statistics of separated reads--------------------:\n")
logwrite(paste("Total number of observations:", total_obs, "\n"))
logwrite(paste("Total number of Tsc:", total_Tsc, "\n"))
logwrite(paste("Total number of Nalm6:", total_Nalm6, "\n"))
logwrite(paste("Total number of other:", total_other, "\n"))



# 4.Scatter plot -------------------------------------------------------------------------------------

#print(head(barcode_annotation))
#colnames(barcode_annotation)

barcode_annotation_subset <- barcode_annotation %>% select(RG, WellID)

final_read_assignment_result_sorted <- final_read_assignment_result_sorted %>%
  left_join(barcode_annotation_subset, by = "RG")

#colnames(final_read_assignment_result_sorted)
#print(head(final_read_assignment_result_sorted))

determine_celltype <- function(condition, WellID) {
  condition_suffix <- str_split(condition, "_", simplify = TRUE)[, 2]
  condition_prefix <- str_split(condition, "_", simplify = TRUE)[, 1]
  
  if (!is.na(condition_suffix) && condition_suffix == "mono") {
    return(condition)
  } else if (!is.na(WellID)) {
    if (any(str_detect(WellID, "[ABC]"))) {
      return(paste0(condition_prefix, "_T_sc"))
    } else if (any(str_detect(WellID, "[DEF]"))) {
      return(paste0(condition_prefix, "_Nalm6_sc"))
    } else if (any(str_detect(WellID, "[GHIJKLMNOP]"))) {
      return(paste0(condition_prefix, "_doublets"))
    }
  }
  return(NA)
}

final_read_assignment_result_sorted <- final_read_assignment_result_sorted %>%
  rowwise() %>%
  mutate(celltype = determine_celltype(condition, WellID)) %>%
  ungroup() %>%
  mutate(
    time = str_extract(condition, "\\d+min|\\d+hr"),
    time = replace_na(time, "0")
  ) %>%
  filter(!str_detect(condition, "Mart1")) %>%
  mutate(
    Nalm6_count = as.numeric(Nalm6_count),
    Tsc_count = as.numeric(Tsc_count)
  )

# Get the axis range
x_min <- min(final_read_assignment_result_sorted$Nalm6_count, na.rm = TRUE)
x_max <- max(final_read_assignment_result_sorted$Nalm6_count, na.rm = TRUE)
y_min <- min(final_read_assignment_result_sorted$Tsc_count, na.rm = TRUE)
y_max <- max(final_read_assignment_result_sorted$Tsc_count, na.rm = TRUE)

# Get all celltypes
celltypes <- unique(final_read_assignment_result_sorted$celltype)


html_list <- list()
# Traverse each celltype and generate images
for (celltype in celltypes) {
  subset <- final_read_assignment_result_sorted %>% filter(celltype == !!celltype)
  
  plot <- ggplot(subset, aes(x = Nalm6_count, y = Tsc_count, color = time)) +
    geom_point(size = 2, alpha = 0.6) +
    labs(title = paste("Celltype:", celltype), x = "Nalm6 Count", y = "Tsc Count") +
    xlim(x_min, x_max) +
    ylim(y_min, y_max) +
    theme_minimal()
  
  filename <- paste0(outpath, "/scatterplot_", celltype, ".png")
  ggsave(filename, plot, width = 5, height = 5)
  
  # Add images to HTML list
  html_list <- append(html_list, tags$div(
    tags$h3(paste("Celltype:", celltype)),
    tags$img(src = paste0("scatterplot_", celltype, ".png"), width = "500px", height = "500px")
  ))
}

# Generate HTML file
html_page <- tags$html(
  tags$head(tags$title("Scatterplots by Celltype")),
  tags$body(
    tags$h1("Scatterplots by Celltype"),
    html_list
  )
)

html_file <- paste0(outpath, "/scatterplots_summary.html")
write(as.character(html_page), html_file)

cat("Plots and HTML file saved to:", outpath, "\n")



# 5.Violin plot ---------------------------------------------------------------------

# Load and preprocess data
file_path_2 <- "/mnt/nvme/home/zhouhui/project_workspace/Tcell_NALM6_project/reads_counts_and_testing/read_counts_per_cell_total.csv"
sum_read_count_per_cell <- read.csv(file_path_2)

# Rename 'RG' to 'Barcode'
final_read_assignment_result_sorted <- final_read_assignment_result_sorted %>% rename(Barcode = RG)

# Merge datasets on 'Barcode'
final_read_assignment_result_sorted <- left_join(
  final_read_assignment_result_sorted,
  sum_read_count_per_cell %>% select(Barcode, Read_Count),
  by = "Barcode"
)

# Create new columns for further analysis
final_read_assignment_result_sorted <- final_read_assignment_result_sorted %>%
  mutate(
    Read_with_SNP_counts = Tsc_count + Nalm6_count + other_count,
    Read_with_SNP_proportion = Read_with_SNP_counts / Read_Count,
    Read_Tsc_Nalm6 = Tsc_count + Nalm6_count,
    Read_Tsc_Nalm6_proportion = Read_Tsc_Nalm6 / Read_Count
  )

# Filter data for each cell type
data_Tsc_sc <- filter(final_read_assignment_result_sorted, celltype == 'Tsc_sc')
data_Nalm6_mono <- filter(final_read_assignment_result_sorted, celltype == 'Nalm6_mono')
data_Tsc_doublets <- filter(final_read_assignment_result_sorted, celltype == 'Tsc_doublets')
data_Tsc_mono <- filter(final_read_assignment_result_sorted, celltype == 'Tsc_mono')

# Initialize an HTML list to store plots
html_list <- list()

# Define a function to create and save violin plots
plot_violin <- function(input_data, celltype) {
  # Define plot colors
  colors <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd")
  
  # Define columns to be plotted
  cols_to_plot <- c("Tsc_count", "Nalm6_count", "Read_Count", "Read_with_SNP_proportion","Read_Tsc_Nalm6_proportion")
  
  # Initialize a list to store individual plots
  plot_list <- list()
  
  # Create violin plots for each column
  for (col in cols_to_plot) {
    plot <- ggplot(input_data, aes(x = celltype, y = !!sym(col), fill = celltype)) +
      geom_violin(trim = TRUE, alpha = 0.7) +
      geom_boxplot(width = 0.1, alpha = 0.8, outlier.color = "black", outlier.size = 1.5) +
      labs(
        title = paste("Violin Plot for", col, "-", celltype),
        y = col,
        x = NULL
      ) +
      theme_minimal(base_size = 15) +
      scale_fill_manual(values = colors) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.text.x = element_blank(),
        legend.position = "none"
      )
    
    # Save individual plots
    filename <- paste0(outpath, "/violin_plot_", celltype, "_", col, ".png")
    ggsave(filename, plot, width = 6, height = 5)
    
    # Append individual plot to HTML list
    html_list <<- append(html_list, tags$div(
      tags$h3(paste("Celltype:", celltype, " - ", col)),
      tags$img(src = paste0("violin_plot_", celltype, "_", col, ".png"), width = "500px", height = "500px")
    ))
    
    # Add plot to the plot list
    plot_list[[col]] <- plot
  }
  
  # Combine individual plots into a grid
  combined_plot <- grid.arrange(grobs = plot_list, ncol = 2)
  combined_filename <- paste0(outpath, "/combined_violin_plot_", celltype, ".png")
  ggsave(combined_filename, combined_plot, width = 12, height = 10)
  
  # Append the combined plot to the HTML list
  html_list <<- append(html_list, tags$div(
    tags$h3(paste("Combined Violin Plot for Celltype:", celltype)),
    tags$img(src = paste0("combined_violin_plot_", celltype, ".png"), width = "1000px", height = "800px")
  ))
}

# Generate plots for each subset of data
plot_violin(data_Nalm6_mono, "Nalm6_mono")
plot_violin(data_Tsc_mono, "Tsc_mono")
plot_violin(data_Tsc_doublets, "Tsc_doublets")
plot_violin(data_Tsc_sc, "Tsc_sc")
plot_violin(final_read_assignment_result_sorted, "All_Data")

# Create the HTML page for visualization
html_page <- tags$html(
  tags$head(tags$title("Violin Plots by Celltype")),
  tags$body(
    tags$h1("Violin Plots by Celltype"),
    html_list
  )
)

# Save the HTML file
html_file <- paste0(outpath, "/violin_plots_summary.html")
write(as.character(html_page), html_file)

cat("Plots and HTML report saved to:", outpath, "\n")