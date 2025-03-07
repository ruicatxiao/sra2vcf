#!/usr/bin/env Rscript

# Load necessary libraries
library(ggplot2)
library(reshape2)
library(dplyr)
library(tools)

# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Please provide the path to the TSV file, output directory, and output prefix as arguments.")
}
tsv_file <- args[1]
output_dir <- args[2]
output_prefix <- args[3]

# Read the TSV file
data <- read.table(tsv_file, header=TRUE, sep="\t", stringsAsFactors = FALSE)

# Initialize columns for counts
data$ref_count <- NA
data$alt_count <- NA

# Process each row
for (i in 1:nrow(data)) {
  AD <- data$AD[i]
  DP <- as.numeric(data$DP[i])

  if (!is.na(DP) && DP >= 20) {
    if (!is.na(AD) && AD != 'NA') {
      # SNP Variant
      ad_values <- as.numeric(unlist(strsplit(AD, ',')))
      if (length(ad_values) >= 2 && !any(is.na(ad_values))) {
        data$ref_count[i] <- ad_values[1]
        data$alt_count[i] <- sum(ad_values[-1])
      } else {
        data$ref_count[i] <- NA
        data$alt_count[i] <- NA
      }
    } else if (!is.na(data$IDV[i]) && data$IDV[i] != 'NA') {
      # INDEL Variant
      IDV <- as.numeric(data$IDV[i])
      data$alt_count[i] <- IDV
      data$ref_count[i] <- DP - IDV
    } else {
      # Sufficient coverage but no variant detected (e.g., INDEL with IDV=0 or SNP with AF=0)
      if (!is.na(data$IDV[i]) && data$IDV[i] == '0') {
        # INDEL with IDV=0
        data$alt_count[i] <- 0
        data$ref_count[i] <- DP
      } else if (!is.na(AD) && AD != 'NA') {
        # SNP with AF=0
        ad_values <- as.numeric(unlist(strsplit(AD, ',')))
        if (length(ad_values) >= 2 && !any(is.na(ad_values))) {
          data$ref_count[i] <- ad_values[1]
          data$alt_count[i] <- sum(ad_values[-1])
        } else {
          data$ref_count[i] <- DP
          data$alt_count[i] <- 0
        }
      } else {
        # No variant detected, report all reads as reference
        data$ref_count[i] <- DP
        data$alt_count[i] <- 0
      }
    }
  } else {
    # Insufficient coverage or missing DP
    data$ref_count[i] <- NA
    data$alt_count[i] <- NA
  }
}

# Remove samples with NA counts (insufficient coverage)
data_filtered <- data %>% filter(!is.na(ref_count) & !is.na(alt_count))

# Handle cases where DP is zero or counts are NA
data_filtered$ref_count[is.na(data_filtered$ref_count) | data_filtered$DP == 0] <- 0
data_filtered$alt_count[is.na(data_filtered$alt_count) | data_filtered$DP == 0] <- 0

# Calculate proportions
data_filtered$ref_prop <- ifelse(data_filtered$DP > 0, data_filtered$ref_count / data_filtered$DP, 0)
data_filtered$alt_prop <- ifelse(data_filtered$DP > 0, data_filtered$alt_count / data_filtered$DP, 0)

# Sort data by decreasing 'alt_prop' to highlight samples with higher alternative allele frequency
data_sorted <- data_filtered %>% arrange(desc(alt_prop))
sample_order <- data_sorted$sampleID

# For absolute stacked bar chart
# Melt the data for plotting
data_melt <- melt(data_sorted[, c("sampleID", "ref_count", "alt_count")], id.vars = "sampleID")

# Set factor levels for sampleID to maintain order
data_melt$sampleID <- factor(data_melt$sampleID, levels = sample_order)

# Plot absolute counts
p1 <- ggplot(data_melt, aes(x=sampleID, y=(value+1), fill=variable)) +
  geom_bar(stat="identity") +
  labs(x="Sample ID", y="Counts (log scale)", title=paste0("Absolute Counts of Reference and Alternative\nVariant: ", output_prefix)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("ref_count"="blue", "alt_count"="red"), labels=c("Reference", "Alternative")) +
  scale_y_log10() +
  theme(plot.title = element_text(hjust = 0.5))

# Save the plot to PDF
plot_file_abs <- file.path(output_dir, paste0(output_prefix, "_absolute_counts.pdf"))
ggsave(plot_file_abs, plot=p1, width=10, height=6)

# For proportion stacked bar chart
# Melt the data for plotting
data_melt_prop <- melt(data_sorted[, c("sampleID", "ref_prop", "alt_prop")], id.vars = "sampleID")

# Set factor levels for sampleID to maintain order
data_melt_prop$sampleID <- factor(data_melt_prop$sampleID, levels = sample_order)

# Plot proportions
p2 <- ggplot(data_melt_prop, aes(x=sampleID, y=value, fill=variable)) +
  geom_bar(stat="identity") +
  labs(x="Sample ID", y="Proportion", title=paste0("Proportion of Reference and Alternative\nVariant: ", output_prefix)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("ref_prop"="blue", "alt_prop"="red"), labels=c("Reference", "Alternative")) +
  theme(plot.title = element_text(hjust = 0.5))

# Save the plot to PDF
plot_file_prop <- file.path(output_dir, paste0(output_prefix, "_proportion.pdf"))
ggsave(plot_file_prop, plot=p2, width=10, height=6)

# List samples with insufficient coverage
insufficient_coverage_samples <- data %>% filter(is.na(ref_count) | is.na(alt_count)) %>% select(sampleID, DP)

if (nrow(insufficient_coverage_samples) > 0) {
  message("Samples with insufficient coverage (DP < 20) are not included in the plots:")
  print(insufficient_coverage_samples)
}

# Optionally, save the list of samples with insufficient coverage
insufficient_coverage_file <- file.path(output_dir, paste0(output_prefix, "_insufficient_coverage.txt"))
write.table(insufficient_coverage_samples, insufficient_coverage_file, sep="\t", row.names=FALSE, quote=FALSE)
