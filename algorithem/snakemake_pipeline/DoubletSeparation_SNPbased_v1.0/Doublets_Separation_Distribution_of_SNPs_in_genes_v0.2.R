# 1.Load necessary libraries -----------------------------------------------------------------
suppressPackageStartupMessages(library(data.table))    # For reading and handling data tables
suppressPackageStartupMessages(library(rtracklayer))    # For reading GTF files
suppressPackageStartupMessages(library(GenomicRanges))  # For genomic ranges operations
suppressPackageStartupMessages(library(ggplot2))        # For plotting
suppressPackageStartupMessages(library(optparse))       # For parsing command-line arguments




# 2.Define command-line options ----------------------------------------------------------------
option_list <- list(
  make_option(c("-G", "--gtf"), type = "character",
              help = "GTF file"),
  make_option(c("-v", "--vcf"), type = "character",
              help = "SNP list (VCF file) with variant annotation. Mandatory"),
  make_option(c("-o", "--outputDir"), type = "character",
              help = "Directory to save output files. Mandatory")
)
## Parse command-line arguments
option_list <- parse_args(OptionParser(option_list = option_list))

## Check if mandatory parameters are provided
if (any(is.null(option_list$vcf), is.null(option_list$gtf), is.null(option_list$outputDir))) {
  stop("All mandatory parameters must be provided. See script usage (--help)")
}

## Print out the parameters for verification
#print(option_list)

## Retrieve input file paths and other parameters from the command-line arguments
vcf_file <- option_list$vcf      # Path to VCF file
output_dir <- option_list$outputDir  # Output directory
gtf_file <- option_list$gtf

if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
}

## Initialize log file
logfile <- paste0(output_dir, "/SNP_distribution_per_gene.processing_log.txt")
logwrite <- function(msg) {
  write(msg, file = logfile, append = TRUE)
}
logwrite("Module1 SNP_distribution_per_gene\n")




# 3.Read and process GTF file using rtracklayer ---------------------------------------------
print("Reading GTF")
gtf <- import(gtf_file)

print("Processing GTF...")
## Extract gene information from GTF
genes <- gtf[gtf$type == "gene"]  # Keep only gene information

logwrite("GTF is imported and processed!\n")
logwrite("Gene annotation Information\n")
logwrite(paste(head(genes, n=5), collapse = "\n"))
#print("Gene Information")
#print(head(genes))  # Check the gene data





# 4.Read and process VCF file -----------------------------------------------------------------
print("Reading Variants...")
## Check if the VCF file is gzipped and read the data accordingly
if (grepl(vcf_file, pattern = ".gz$")) {
  vcf_data <- fread(cmd = paste("zcat", vcf_file, " | grep -v '^##'"), header = TRUE)
} else {
  vcf_data <- fread(cmd = paste("grep -v '^##'", vcf_file), header = TRUE)
}

print("Processing Variants...")
## Extract chromosome and position information from VCF
vcf_positions <- vcf_data[, c("#CHROM", "POS")]
## Rename columns for easier access
setnames(vcf_positions, old = "#CHROM", new = "CHROM")

## Filter out chromosomes containing "KI" or "GL"
vcf_positions_filtered <- vcf_positions[!grepl("^(KI|GL)", vcf_positions$CHROM), ]

## Print out filtered VCF positions
#print("Filtered VCF Positions")
#print(head(vcf_positions_filtered))  # Check the SNP information

## Check structure of filtered VCF positions
#print("Structure of VCF Positions")
#str(vcf_positions_filtered)
logwrite("VCF is imported and processed!\n")
logwrite("VCF Information\n")
logwrite(paste(head(vcf_positions_filtered, n=5), collapse = "\n"))

## Create SNP ranges as a GRanges object
snp_gr <- GRanges(seqnames = vcf_positions_filtered$CHROM,
                  ranges = IRanges(start = vcf_positions_filtered$POS, end = vcf_positions_filtered$POS))

#print("SNP GRanges")
#print(snp_gr)




# 5.Find overlaps between SNPs and genes --------------------------------------------------------
snp_in_genes <- findOverlaps(snp_gr, genes)
logwrite("Finding overlaps between SNPs and genes is done!\n")
logwrite("Overlaps between SNPs and genes\n")
logwrite(capture.output(str(snp_in_genes)))
#print("SNP-Gene Overlaps")
#str(snp_in_genes)

## Get indices of SNPs and genes that overlap
snp_idx <- queryHits(snp_in_genes)  # SNP indices
gene_idx <- subjectHits(snp_in_genes)  # Gene indices

## Extract gene names from the gene metadata
gene_names <- mcols(genes)$gene_name

## Count the number of SNPs per gene
snp_per_gene <- table(gene_idx)

## Create a data frame containing SNP counts and gene names
snp_counts_df <- data.frame(
  Gene_Index = as.numeric(names(snp_per_gene)),  # Gene indices
  Gene_Name = gene_names[as.numeric(names(snp_per_gene))],  # Gene names
  SNP_Count = as.vector(snp_per_gene)           # Corresponding SNP count
)

## Sort by SNP count in descending order and select top 20 genes
top_20_genes <- snp_counts_df[order(-snp_counts_df$SNP_Count), ][1:20, ]
logwrite("View the top 20 genes with the most SNPs\n")
logwrite(paste(head(top_20_genes, n=20), collapse = "\n"))


## Calculate the total number of SNPs
total_snp_count <- sum(snp_counts_df$SNP_Count)

## Output the total SNP count
logwrite("Total SNP Count:")
logwrite(total_snp_count)



# 6.Save results and plotting -------------------------------------------------------------------
# Save the SNP counts per gene to a CSV file in the specified output directory
output_file <- paste0(output_dir, "/", "snp_counts_per_gene.csv")
write.csv(snp_counts_df, output_file, row.names = FALSE)


# Perform frequency analysis on the SNP count distribution
freq_table <- table(snp_counts_df$SNP_Count)
freq_percentage <- prop.table(freq_table) * 100

# Convert frequency results to a data frame
freq_df <- data.frame(
  SNP_Count = as.numeric(names(freq_table)),  # Different SNP counts
  Frequency = as.vector(freq_table),         # Frequency of each count
  Percentage = as.vector(freq_percentage)    # Percentage of each count
)

# Plot the SNP count frequency distribution
ggplot(freq_df, aes(x = SNP_Count, y = Frequency)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  labs(
    title = "Distribution of SNP Counts per Gene",
    x = "SNP Count",
    y = "Number of Genes"
  ) +
  theme_minimal()

# Save the plot as a PDF file in the specified output directory
output_plot <- paste0(output_dir, "/","snp_count_frequency_distribution.pdf")
ggsave(output_plot, width = 10, height = 6)
