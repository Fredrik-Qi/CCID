# 1.Description
#   This R script is designed for separating the reads from the bam file based on the vcf file
#   (including different heterozygous and homozygous SNPs between two cell lines) 
#   so that the mixed reads in the doublets can be accurately assigned to the corresponding cells

# 2.New functionalities of read-separation algorithm(v0.4 version)
#   1) Read Pair Detection
#     - If one of read1 and read2 is assigned to a specific type of cell line, the other will also be included.

#   2) Separation of Reads Based on both heterozygous and homozygous SNPs
#     - The script also includes functionality to separate reads based on the presence of both 
#       heterozygous and homozygous SNPs rather than just homozygous SNPs

#   3) New parameter("-f","--subfolder") to specify the path of the subfolder

#   4) New parameter("-a","--annotation")

# 3.Usages
#   1)Inputs
#     - path_bam: Path to the input BAM file.
#     - vcf: VCF file containing the SNPs to be used for read separation.
#       - chroms_todo: A vector of chromosome names to filter the data.
#     - cellBCs: A vector of cell barcodes used to filter the reads.
#     - yaml: Other configuration information to ensure consistency with zUMIs pipeline.
#       - ncores: Number of cores to use for parallel processing.
#       - out_dir: Output directory.
#       - project: Project name.
#       - read_layout: Read layout (PE or SE).     
#      - annotation:

#   2)Outputs
#     - processing_log.txt: A log file containing information about the processing steps.
#     - "chr"_out_reads.txt: A file containing the results of read separation with readcall(TdT or Nalm6) and	nVars as annotations.
#     - "chr"_out_read_pair.txt: A file containing the results of read separation for pair reads 
#     - tmp.rds: A file containing the records of separated read(Sum of all "chr"_out_read_pair.txt files)
#     - separated_counts.rds: A file containing the count tables for the  reads.
#
#      

#   3)Example
#   Rscript DoubletSeparation_SNPbased_v0.4_processing_read_pair.R \
#   --yaml a.yaml \
#   --vcf b.vcf.gz \
#   --bam c.bam \
#   -c d.txt 

# 4.Rscript structure
# 1) Packages loading
# 2) log info
# 3) Function
# 4) Input variables
# 5) VCF loading
# 6) Bam loading and read separation
# 7）Results

# 5.Function structure
#   lapply(chroms_todo, read_separation){
#     load_reads_from_bam{
#       parse_flag
#     }
#     SNPdeconv{
#       variant_parsing
#       read_decision
#     }
  
#   }

#   makeWide

# 1.Packages loading ----------------------------------------------------------------
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(GenomicAlignments))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(Rsamtools))
suppressPackageStartupMessages(library(yaml))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(htmltools))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(gridExtra))


# 2.log info ----------------------------------------------------------------
Sys.time()
print("SNP based expression deconvolution of doublets")
print("Zhouhui & Christoph")
print("v 0.6.0")
print(timestamp())



# 3.Function -----------------------------------------------
load_reads_from_bam <- function(path_bam, cellBCs, ncores, chroms_todo){
  idxstats <- Rsamtools::idxstatsBam(path_bam)
  if("*" %in% idxstats$seqnames){
    idxstats <- idxstats[idxstats$seqnames != "*", ]
    idxstats$seqnames <- as.character(idxstats$seqnames)
  }
  idxstats <- idxstats[idxstats$seqnames %in% chroms_todo,]
  taglist <- c("BC", "UB","GE", "GI")
  whatlist <- c("qname", "pos","cigar","seq")
  
  rsamtools_reads <- parallel::mclapply(1:nrow(idxstats), function(x) {
    parms <- ScanBamParam(tag=taglist,
                          what=whatlist,
                          tagFilter = list(BC = cellBCs),
                          which = GRanges(seqnames = idxstats[x,"seqnames"], ranges = IRanges(1,idxstats[x,"seqlength"])))
    
    
    dat <- scanBam(file = path_bam, param = parms)
    
    dt <- data.table(readID = dat[[1]]$qname,chrom = idxstats[x,"seqnames"], pos = dat[[1]]$pos, cigar = dat[[1]]$cigar, seq = as.character(dat[[1]]$seq), RG = dat[[1]]$tag$BC, UB = dat[[1]]$tag$UB, GE = dat[[1]]$tag$GE, GEin = dat[[1]]$tag$GI)
    rm(dat)
    gc(verbose = F)
    dt <- dt[! (is.na(GE) & is.na(GEin)) ]
    
    return(dt)
  }, mc.cores = ncores)
  
  return(rsamtools_reads)
}



variant_parsing <- function(reads, variant_positions){
  #parse all cigars to reference seq
  ops <- c("M", "=", "X")
  ranges_on_ref <- cigarRangesAlongReferenceSpace(reads$cigar, pos=reads$pos, ops=ops)
  ranges_on_query <- cigarRangesAlongQuerySpace(reads$cigar, ops=ops)
  gc(verbose = F)
  range_group <- togroup(PartitioningByWidth(ranges_on_ref))
  ranges_on_ref <- unlist(ranges_on_ref, use.names=FALSE)
  ranges_on_query <- unlist(ranges_on_query, use.names=FALSE)
  query2ref_shift <- start(ranges_on_ref) - start(ranges_on_query)
  
  var_pos <- variant_positions
  hits <- findOverlaps(var_pos, ranges_on_ref)
  hits_at_in_x <- var_pos[queryHits(hits)] - query2ref_shift[subjectHits(hits)]
  hits_group <- range_group[subjectHits(hits)]
  fetched_bases <- subseq(reads[hits_group,]$seq, start=hits_at_in_x, width=1L)
  
  #now add everything together in the output data.table
  out_vars <- data.table(
    obs_base = fetched_bases,
    pos = var_pos[queryHits(hits)]
  )
  out_vars[, c("RG","UB","GE","GEin","readID") := reads[hits_group, c("RG","UB", "GE","GEin", "readID"), with = F] ]
  
  out_vars <- out_vars[obs_base %in% c("A","C","G","T") ]
  setnames(out_vars,"pos","POS")
  
  return(out_vars)
}

read_decision <- function(basecalls){
  if(length(basecalls) == 1){
    return(basecalls)
  }else{
    ux <- unique(basecalls)
    basecall_summary <- tabulate(match(basecalls, ux))
    names(basecall_summary) <- ux
    majority_basecall <- ux[which.max(basecall_summary)]
    if(basecall_summary[majority_basecall]/sum(basecall_summary) >= 0.66){
      return(majority_basecall)
    }else{
      return("other")
    }
  }
}

SNPdeconv <- function(vcf_chunk, reads_chunk, cellBCs, read_layout ){
  
  #get the bases that were observed in our reads according to the positions in VCF
  out_vars <- variant_parsing(reads = reads_chunk, variant_positions = vcf_chunk$POS)
  
  #combine VCF with observations
  out_vars <- merge(out_vars,vcf_chunk,by = "POS" )
  
  #categorize:
  out_vars[              , basecall := "other"][
          obs_base == REF, basecall := "REF"][
          obs_base == ALT, basecall := "ALT"]
  
  out_vars_filename <- paste0(log_dir, "/", vcf_chunk$CHROM[1], "_out_vars.txt")
  write.table(out_vars, file = out_vars_filename, sep = "\t", row.names = FALSE, quote = FALSE)

  #translate ref/alt into cell line:
  cell_lines <- tail(colnames(vcf_chunk),2)

  # The default setting is "other"
  out_vars[, basecall2 := "other"]

  for (i in seq_len(nrow(out_vars))) {
    # Get the basecall and cell_lines values ​​for the current row
    current_basecall <- out_vars[i, basecall]
    genotype_1 <- out_vars[i, get(cell_lines[1])]
    genotype_2 <- out_vars[i, get(cell_lines[2])]
    
    # Logic of Decision Tree 
    if (current_basecall == "REF") {
      if (genotype_1 %in% c("./.", ".|.", "0/0", "0|0") & genotype_2 %in% c("1/1", "1|1")) {
        out_vars[i, basecall2 := cell_lines[1]]
      } else if (genotype_1 %in% c("1/1", "1|1") & genotype_2 %in% c("./.", ".|.", "0/0", "0|0")) {
        out_vars[i, basecall2 := cell_lines[2]]
      } else if (genotype_1 %in% c("1/0", "0/1", "0|1", "1|0") & genotype_2 %in% c("1/1", "1|1")) {
        out_vars[i, basecall2 := cell_lines[1]]
      } else if (genotype_1 %in% c("1/1", "1|1") & genotype_2 %in% c("1/0", "0/1", "0|1", "1|0")) {
        out_vars[i, basecall2 := cell_lines[2]]
      }
    } else if (current_basecall == "ALT") {
      if (genotype_1 %in% c("./.", ".|.", "0/0", "0|0") & genotype_2 %in% c("1/1", "1|1")) {
        out_vars[i, basecall2 := cell_lines[2]]
      } else if (genotype_1 %in% c("1/1", "1|1") & genotype_2 %in% c("./.", ".|.", "0/0", "0|0")) {
        out_vars[i, basecall2 := cell_lines[1]]
      } else if (genotype_1 %in% c("./.", ".|.", "0/0", "0|0") & genotype_2 %in% c("1/0", "0/1", "0|1", "1|0")) {
        out_vars[i, basecall2 := cell_lines[2]]
      } else if (genotype_1 %in% c("1/0", "0/1", "0|1", "1|0") & genotype_2 %in% c("./.", ".|.", "0/0", "0|0")) {
        out_vars[i, basecall2 := cell_lines[1]]
      }
    }
  }

  
  out_reads <- out_vars[, .(readcall = read_decision(basecall2), nVars =length(unique(POS))), by = c("RG","UB","GE","GEin","readID")]
  
  return(out_reads)
}



makeWide <- function(allele_dat, metric = c("reads","UMIs"), level = c("exon","intron","exonintron")){
  dat <- allele_dat
  if(metric == "UMIs"){
    dat <- allele_dat[!UB==""]
  }
  if(level == "exon"){
    dat[, GeneID := GE]
  }else if(level == "intron"){
    dat[, GeneID := GEin]
  }else{
    dat[, GeneID := GE][
        is.na(GeneID), GeneID := GEin]
  }
  dat <- dat[!is.na(GeneID)]
  dat[, c("GE","GEin") := NULL]
  
  cell_lines <- unique(dat$readcall)
  cell_lines <- cell_lines[cell_lines != "other"]
  
  dat_summarized <- dat[, .(reads = length(unique(readID)), UMIs = length(unique(UB))), by =c("RG","GeneID","readcall")]
  
  out_tabs <- list(
    dcast(dat_summarized[readcall == cell_lines[1]], formula = GeneID ~ RG, fill=0, value.var = metric),
    dcast(dat_summarized[readcall == cell_lines[2]], formula = GeneID ~ RG, fill=0, value.var = metric)
  )
  names(out_tabs) <- cell_lines
  return(out_tabs)
}


# 4.Input variables -------------------------------------------------------
option_list <- list(
  make_option(c("-y", "--yaml"), type="character",
              help="zUMIs config file. Mandatory"),
  make_option(c("-v", "--vcf"), type="character",
              help="SNP list (VCF file) with variant annotation. Mandatory"),
  make_option(c("-m","--minCount"), type="integer",
              help="Cutoff for minimum coverage in a Cell/Gene pair. Default: 0",
              default=0),
  make_option(c("-b","--bam"), type="character",
              help="BAM file with reads. Mandatory"),
  make_option(c("-c","--cellBCs"), type="character",
              help="Cell types to be analyzed. Mandatory"),
  make_option(c("-f","--subfolder"), type="character",
              help="subfolder of output directory", default = "separation_results"),
  make_option(c("-a","--annotation"), type="character",
              help="Path of barcode annotation file")
)
option_list <- parse_args(OptionParser(option_list=option_list))

if (any(is.null(option_list$yaml),is.null(option_list$vcf),is.null(option_list$bam))) {
  stop("All mandatory parameters must be provided. See script usage (--help)")
}


# parse args: YAML and subfolder
opt   <- read_yaml(option_list$yaml)
ncores <- opt$num_threads
subfolder_path <- option_list$subfolder

outpath <- paste0(opt$out_dir,"/",subfolder_path,"/")
log_dir <- outpath


if(!dir.exists(outpath)){
  try(system(paste("mkdir",outpath)))
}

outpath <- paste0(outpath,opt$project,".")

# parse args: annotation
barcode_annotation_path <- option_list$annotation


# parse args: VCF
path_snps <- option_list$vcf

# parse args: BAM
path_bam <- option_list$bam

# parse args: cellBC
cellBCs_path <- option_list$cellBCs
cellBCs <- read.table(cellBCs_path, header = TRUE)
cellBCs <- cellBCs$XC_DNBPE

# Setting up
setwd(opt$out_dir)
setDTthreads(ncores)

# Initialize log file
logfile <- paste0(log_dir, "/processing_log.txt")
logwrite <- function(msg) {
  write(msg, file = logfile, append = TRUE)
}

# Record start time
logwrite(paste("Process started at:", Sys.time(), "\n"))

# Record input parameters
logwrite("--------------------Input Parameters--------------------:\n")
logwrite(paste("yaml file:", option_list$yaml, "\n"))
logwrite(paste("VCF file:", option_list$vcf, "\n"))
logwrite(paste("BAM file:", option_list$bam, "\n"))
logwrite(paste("CellBCs file:", option_list$cellBCs, "\n"))
logwrite(paste("Number of threads:", ncores, "\n"))

# 5.VCF loading --------------------------------------------------------------


print("Reading Variants...")
if(grepl(path_snps, pattern = ".gz$")){
  vcf <- fread(cmd = paste("zcat",path_snps," | grep -v '^##'"), header = T)
}else{
  vcf <- fread(cmd = paste("grep -v '^##'",path_snps), header = T)
}

setnames(vcf,c("CHROM",colnames(vcf)[2:length(colnames(vcf))]))

vcf <- vcf[nchar(REF)==1][nchar(ALT)==1]

print("Done!")
Sys.time()

chroms_todo <- unique(vcf$CHROM)
chroms_todo <- chroms_todo[! grepl(pattern = "KI", chroms_todo)]
chroms_todo <- chroms_todo[! grepl(pattern = "GL", chroms_todo)]
chroms_todo <- chroms_todo[! grepl(pattern = "MT", chroms_todo)]

logwrite("VCF file loaded successfully.\n")
logwrite("--------------------First few lines of VCF--------------------:\n")
logwrite(paste(head(vcf, n=3), collapse = "\n"))


# 6.Bam loading and read separation ------------------------------------------------------------
vcf_list <- split(vcf[CHROM %in% chroms_todo], by = "CHROM")
out_list <- lapply(chroms_todo, function(chr){
  print(paste("Working on ",chr,"..."))
  logwrite(paste("Working on chromosome", chr, "..."))
  #load reads
  print("Loading reads from bam file...")
  this_reads <- load_reads_from_bam(path_bam = path_bam, cellBCs = cellBCs, ncores = ncores, chroms_todo = chr)
  this_reads <- rbindlist(this_reads)
  logwrite("Reads loaded and processed.\n")
  logwrite("First few lines of loaded reads:\n")
  logwrite(paste(head(this_reads, n=3), collapse = "\n"))
  print("Done!")
  print("Genotyping reads based on SNPs...")

  
  out_reads <- SNPdeconv(vcf_chunk = vcf_list[[chr]], reads_chunk = this_reads, cellBCs = cellBCs, read_layout = opt$read_layout )
  logwrite("Cell line assignment.\n")
  logwrite("First few lines of Cell line assignment:\n")
  logwrite(paste(head(out_reads, n=3), collapse = "\n"))
  outread_filename <- paste0(log_dir, "/", chr, "_out_read.txt")
  write.table(out_reads, file = outread_filename, sep = "\t", row.names = FALSE, quote = FALSE)
  # 提取 out_reads 中的 readID 列
  readID_readcall <- out_reads[, .(readID, readcall)]
  # 保存到文本文件
  outread_readID_filename <- paste0(log_dir, "/", chr, "_out_read_readID.txt")
  fwrite(readID_readcall, file = outread_readID_filename, sep = "\t", quote = FALSE)

  print("Done!")
  return(out_reads)

})



# 8.Results -------------------------------------------------------------------------------------

tmp_rds_path <- paste0(outpath,"tmp.rds")
saveRDS(out_list, file = tmp_rds_path, compress = FALSE)
print("Combining output")
out <- rbindlist(out_list)
print("Done!")

gc(verbose = FALSE)

print("nVars_distribution_plot")

nVars_distribution_plot <- ggplot(out, aes(x = nVars)) +
  geom_bar(stat = "count", fill = "skyblue", color = "black") +
  labs(title = "Distribution of nVars", x = "nVars", y = "Frequency") +
  theme_minimal() +
  geom_text(stat = "count", aes(label = ..count..), vjust = -0.5)  # 添加数字标签


# 指定保存路径
nVars_distribution_plot_path <- paste0(log_dir, "/nVars_distribution_plot.png")

# 保存柱状图到指定路径
ggsave(nVars_distribution_plot_path, plot = nVars_distribution_plot, width = 8, height = 6)



print("Making count tables in different flavors...")
final <- list(
  "counts" = list(
    "exon" = test <- makeWide(out, metric = "reads", level = "exon"),
    "intron" = test <- makeWide(out, metric = "reads", level = "intron"),
    "exonintron" = test <- makeWide(out, metric = "reads", level = "exonintron")
  ),
  
  "umicounts" = list(
    "exon" = test <- makeWide(out, metric = "UMIs", level = "exon"),
    "intron" = test <- makeWide(out, metric = "UMIs", level = "intron"),
    "exonintron" = test <- makeWide(out, metric = "UMIs", level = "exonintron")
  )
)
print("Done!")

print("Saving output...")
saveRDS(final, file = paste0(outpath,"separated_counts.rds"), compress = FALSE)
print("Done!")
print(timestamp())


# 8.Statistical analysis -------------------------------------------------------------------------------------


# load RDA file
dgecounts <- readRDS(tmp_rds_path)
# load barcode annotation
barcode_annotation <- read.table(barcode_annotation_path, header = TRUE, sep = "\t")

logwrite("Barcode annotation file loaded successfully.\n")
logwrite("--------------------First few lines of Barcode annotation file--------------------:\n")
logwrite(paste(head(barcode_annotation, n=3), collapse = "\n"))

# Initialize an empty read_assignment_result data.table
read_assignment_result <- data.table(RG = character(), TdT_count = integer(), Nalm6_count = integer(), other_count = integer())

# Iterate through all data.tables
for (dt in dgecounts) {
  # Use data.table's .N to count TdT and Nalm6 occurrences
  counts <- dt[, .(TdT_count = sum(readcall == "TdT"), 
                   Nalm6_count = sum(readcall == "Nalm6"), 
                   other_count = sum(readcall == "other")),
                by = RG]
  
  # Bind the current data.table's results to the overall result
  read_assignment_result <- rbind(read_assignment_result, counts, use.names = TRUE, fill = TRUE)
}

# Group by RG and sum the counts
final_read_assignment_result <- read_assignment_result[, .(TdT_count = sum(TdT_count, na.rm = TRUE), 
                            Nalm6_count = sum(Nalm6_count, na.rm = TRUE),
                            other_count = sum(other_count, na.rm = TRUE)), 
                        by = RG]
barcode_annotation_subset <- barcode_annotation[, c("condition", "RG")]
final_read_assignment_result_with_barcodeinfo <- merge(final_read_assignment_result, barcode_annotation_subset, by = "RG")
# Sort by 'condition' first, and then by 'TdT_count'
final_read_assignment_result_sorted <- final_read_assignment_result_with_barcodeinfo[order(final_read_assignment_result_with_barcodeinfo$condition, 
                                                           final_read_assignment_result_with_barcodeinfo$TdT_count), ]


# Save the final result as a text file
fwrite(final_read_assignment_result_sorted, paste0(outpath,"final_read_assignment_result_sorted.txt"), sep = "\t")

#colnames(final_read_assignment_result_sorted)

# Initialize counters
total_obs <- 0
total_TdT <- 0
total_Nalm6 <- 0
total_other <- 0

# Iterate through all data.tables
for (dt in dgecounts) {
  # Add the number of observations in the current data.table to total_obs
  total_obs <- total_obs + nrow(dt)
  
  # Count TdT and Nalm6 in the current data.table
  total_TdT <- total_TdT + sum(dt$readcall == "TdT")
  total_Nalm6 <- total_Nalm6 + sum(dt$readcall == "Nalm6")
  total_other <- total_other + sum(dt$readcall == "other")
}

# Print the results
logwrite("--------------------Statistics of separated reads--------------------:\n")
logwrite(paste("Total number of observations:", total_obs, "\n"))
logwrite(paste("Total number of TdT:", total_TdT, "\n"))
logwrite(paste("Total number of Nalm6:", total_Nalm6, "\n"))
logwrite(paste("Total number of other:", total_other, "\n"))


# --------------------------------------------------------------------------------------------
# Scatter plot

#print(head(barcode_annotation))
colnames(barcode_annotation)

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
    TdT_count = as.numeric(TdT_count)
  )

# 获取坐标轴范围
x_min <- min(final_read_assignment_result_sorted$Nalm6_count, na.rm = TRUE)
x_max <- max(final_read_assignment_result_sorted$Nalm6_count, na.rm = TRUE)
y_min <- min(final_read_assignment_result_sorted$TdT_count, na.rm = TRUE)
y_max <- max(final_read_assignment_result_sorted$TdT_count, na.rm = TRUE)

# 获取所有 celltype
celltypes <- unique(final_read_assignment_result_sorted$celltype)

# 初始化 HTML 列表
html_list <- list()

# 遍历每个 celltype 并生成图片
for (celltype in celltypes) {
  subset <- final_read_assignment_result_sorted %>% filter(celltype == !!celltype)
  
  # 创建图
  plot <- ggplot(subset, aes(x = Nalm6_count, y = TdT_count, color = time)) +
    geom_point(size = 2, alpha = 0.6) +
    labs(title = paste("Celltype:", celltype), x = "Nalm6 Count", y = "TdT Count") +
    xlim(x_min, x_max) +
    ylim(y_min, y_max) +
    theme_minimal()
  
  # 保存图片
  filename <- paste0(outpath, "scatterplot_", celltype, ".png")
  ggsave(filename, plot, width = 5, height = 5)
  
  # 添加图片到 HTML 列表
  html_list <- append(html_list, tags$div(
    tags$h3(paste("Celltype:", celltype)),
    tags$img(src = paste0("scatterplot_", celltype, ".png"), width = "500px", height = "500px")
  ))
}

# 生成 HTML 文件
html_page <- tags$html(
  tags$head(tags$title("Scatterplots by Celltype")),
  tags$body(
    tags$h1("Scatterplots by Celltype"),
    html_list
  )
)

# 保存 HTML 文件
html_file <- paste0(outpath, "scatterplots_summary.html")
write(as.character(html_page), html_file)

cat("Plots and HTML file saved to:", outpath, "\n")

# ---------------------------------------------------------------------

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
    Read_with_SNP_counts = TdT_count + Nalm6_count + other_count,
    Read_with_SNP_proportion = Read_with_SNP_counts / Read_Count,
    Read_TdT_Nalm6 = TdT_count + Nalm6_count,
    Read_TdT_Nalm6_proportion = Read_TdT_Nalm6 / Read_Count
  )

# Filter data for each cell type
data_TdT_sc <- filter(final_read_assignment_result_sorted, celltype == 'TdT_sc')
data_Nalm6_mono <- filter(final_read_assignment_result_sorted, celltype == 'Nalm6_mono')
data_TdT_doublets <- filter(final_read_assignment_result_sorted, celltype == 'TdT_doublets')
data_TdT_mono <- filter(final_read_assignment_result_sorted, celltype == 'TdT_mono')

# Initialize an HTML list to store plots
html_list <- list()

# Define a function to create and save violin plots
plot_violin <- function(input_data, celltype) {
  # Define plot colors
  colors <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd")
  
  # Define columns to be plotted
  cols_to_plot <- c("TdT_count", "Nalm6_count", "Read_Count", "Read_with_SNP_proportion","Read_TdT_Nalm6_proportion")
  
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
    filename <- paste0(outpath, "violin_plot_", celltype, "_", col, ".png")
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
  combined_filename <- paste0(outpath, "combined_violin_plot_", celltype, ".png")
  ggsave(combined_filename, combined_plot, width = 12, height = 10)
  
  # Append the combined plot to the HTML list
  html_list <<- append(html_list, tags$div(
    tags$h3(paste("Combined Violin Plot for Celltype:", celltype)),
    tags$img(src = paste0("combined_violin_plot_", celltype, ".png"), width = "1000px", height = "800px")
  ))
}

# Generate plots for each subset of data
plot_violin(data_Nalm6_mono, "Nalm6_mono")
plot_violin(data_TdT_mono, "TdT_mono")
plot_violin(data_TdT_doublets, "TdT_doublets")
plot_violin(data_TdT_sc, "TdT_sc")
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
html_file <- paste0(outpath, "violin_plots_summary.html")
write(as.character(html_page), html_file)

cat("Plots and HTML report saved to:", outpath, "\n")