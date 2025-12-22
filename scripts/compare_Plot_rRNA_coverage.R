
# Load necessary libraries
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(Hmisc))
library(jsonlite)

option_list = list(
  make_option(c("-i", "--input"), type = "character", help = "GSE id"),
  make_option(c("-o", "--output"), type = "character", help = "Output directory for plots and data"),
  make_option(c("-s", "--samples"), type = "character", help = "Comma-separated list of sample names"),
  make_option(c("-g", "--groups"), type = "character", help = "String representation of the groups dictionary")
)
opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

# Get input parameters
GSE_ID <- opt$input
BASIC_OUTPUT_DIR <- opt$output
SAMPLES <- strsplit(opt$samples, ",")[[1]]
GROUPS <- fromJSON(opt$groups)
case_samples <- GROUPS$case
control_samples <- GROUPS$control



# Define rRNA subtypes and their lengths
rRNA_lengths_df <- data.frame(
  Subtype = c("5S", "5.8S", "12S", "16S", "18S", "28S", "45S"),
  Length = c(121, 159, 954, 1559, 1869, 5070, 13357)
)

# --- Main Script ---
MATRIX_FILE_PATH <- file.path(BASIC_OUTPUT_DIR, "stats", paste0(GSE_ID, "_rRNA_expression_matrix.txt"))
df_matrix <- read.delim(MATRIX_FILE_PATH, header = TRUE, sep = "\t")

rRNA_subtypes_to_plot <- rRNA_lengths_df$Subtype

OUTPUT_DIR <- file.path(BASIC_OUTPUT_DIR, "plots")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

for (subtype in rRNA_subtypes_to_plot) {
  message(paste0("Processing ", subtype, "..."))
  rRNA_len <- rRNA_lengths_df[rRNA_lengths_df$Subtype == subtype, "Length"]
  subtype_data <- df_matrix[df_matrix$Subtype == subtype, ]

  sample_cols <- intersect(c(case_samples, control_samples), colnames(df_matrix))
  coverage_mat <- matrix(0, nrow = rRNA_len, ncol = length(sample_cols))
  colnames(coverage_mat) <- sample_cols

  for (i in seq_len(nrow(subtype_data))) {
    start_pos <- as.numeric(subtype_data$Start[i])
    frag_len <- as.numeric(subtype_data$Length[i])
    if (!is.finite(start_pos) || !is.finite(frag_len)) next
    covered <- seq(from = start_pos, to = start_pos + frag_len - 1)
    covered <- covered[covered >= 1 & covered <= rRNA_len]
    for (j in seq_along(sample_cols)) {
      rpm <- as.numeric(subtype_data[[sample_cols[j]]][i])
      if (!is.finite(rpm) || rpm == 0) next
      coverage_mat[covered, j] <- coverage_mat[covered, j] + rpm
    }
  }

  # Calculate mean coverage
  type_vec <- sapply(sample_cols, function(sample) ifelse(sample %in% case_samples, "case", "control"))
  mean_case <- rowMeans(coverage_mat[, type_vec == "case", drop = FALSE], na.rm = TRUE)
  mean_control <- rowMeans(coverage_mat[, type_vec == "control", drop = FALSE], na.rm = TRUE)

  # Output
  out_df <- data.frame(length = 1:rRNA_len, coverage_mat, mean.case = mean_case, mean.control = mean_control)
  out_path <- file.path(OUTPUT_DIR, paste0(GSE_ID, "_", subtype, "_Coverage_Data.txt"))
  write.table(out_df, file = out_path, sep = "\t", row.names = FALSE, quote = FALSE)
  message(paste0("Saved coverage data for ", subtype, " to ", out_path))

  # Plot
  plot_long <- tidyr::pivot_longer(out_df[, 1:(ncol(out_df) - 2)], cols = sample_cols, names_to = "sample", values_to = "RPM")
  plot_long$Type <- sapply(plot_long$sample, function(sample) ifelse(sample %in% case_samples, "case", "control"))
  p <- ggplot(plot_long, aes(x = length, y = RPM, color = Type, fill = Type)) +
    stat_summary(fun = "mean", geom = "line", size = 1) +
    labs(x = "Length", y = "RPM", title = paste0(GSE_ID, " - ", subtype, " rRNA Coverage")) +
    theme_bw() +
    theme(
      axis.line = element_line(),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      panel.background = element_rect(fill = "transparent", colour = NA),
      plot.background = element_rect(fill = "transparent", colour = NA),
      legend.title = element_blank(),
      legend.background = element_rect(fill = "transparent", colour = NA),
      legend.box.background = element_rect(fill = "transparent", colour = NA),
      axis.text = element_text(size = 10, color = "black"),
      axis.title = element_text(size = 10),
      plot.title = element_text(hjust = 0.5)
    ) +
    scale_color_manual(values = c("case" = "#2a3379", "control" = "#f57e2a")) +
    scale_fill_manual(values = c("case" = "#2a3379", "control" = "#f57e2a"))
  
  output_pdf_path <- file.path(OUTPUT_DIR, paste0(GSE_ID, "_", subtype, "_Coverage.pdf"))
  ggsave(output_pdf_path, plot = p, width = 8, height = 5, dpi = 300)
  message(paste0("Saved plot for ", subtype, " to ", output_pdf_path))

  output_png_path <- file.path(OUTPUT_DIR, paste0(GSE_ID, "_", subtype, "_Coverage.png"))
  ggsave(output_png_path, plot = p, width = 8, height = 5, dpi = 300)
}

message("Plots finished.")
