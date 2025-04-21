# 1.Description
#   This R script is designed for separating the reads from the bam file based on the vcf file
#   (including different heterozygous and homozygous SNPs between two cell lines) 
#   so that the mixed reads in the doublets can be accurately assigned to the corresponding cells

# 2.New functionalities of read-separation algorithm(v0.6 version)
#   1) ReadID record of separated reads
#     - The logic is that, if one of read1 and read2 is assigned to a specific type of cell line, 
#       the other will also be included. Therefore, in this algorithm, we record all the ReadID of 
#       separated reads and their corresponding celltype, so we can extract sub_bam files in the 
#       following modules to get both read1 and read2 based on the readID.

#   2) Separation of Reads Based on both heterozygous and homozygous SNPs
#     - The script also includes functionality to separate reads based on the presence of both 
#       heterozygous and homozygous SNPs rather than just homozygous SNPs

#   3) New parameter("-f","--subfolder") to specify the path of the subfolder

#   4) New parameter("-a","--annotation") to specify the path of barcode annotation file

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
#      - annotation: the path of barcode annotation file
#      - subfolder: the path of the subfolder

#   2)Outputs
#     - processing_log.txt: A log file containing information about the processing steps.
#     - "chr"_out_reads.txt: A file containing the results of read separation with readcall(Tsc or Nalm6) and	nVars as annotations.
#     - "chr"_out_vars.txt: A file containing the results of read separation
#     - tmp.rds: A file containing the records of separated read(Sum of all "chr"_out_read_pair.txt files)
#     - separated_counts.rds: A file containing the count tables for the reads.
#     - Tsc_readIDs.txt & Nalm6_readIDs.txt: ReadIDs of known cell lines
#      - nVars_distribution_plot.png: Statistics on the distribution of SNPs on each read

#   3)Example
#   Rscript DoubletSeparation_SNPbased_v0.4_processing_read_pair.R \
#   --yaml a.yaml \
#   --vcf b.vcf.gz \
#   --bam c.bam \
#   -c d.txt 
#   -f subfolder
#   -a barcode_annotation.txt

# 4.Rscript structure
# 1) Packages loading
# 2) log info
# 3) Function
# 4) Input variables
# 5) VCF loading
# 6) Bam loading and read separation
# 7) Sort ReadID of seperated reads
# 8）Results

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
logfile <- paste0(log_dir, "/doublets_separation.processing_log.txt")
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
  # Extract the readID column from out_reads
  readID_readcall <- out_reads[, .(readID, readcall)]
  outread_readID_filename <- paste0(log_dir, "/", chr, "_out_read_readID.txt")
  fwrite(readID_readcall, file = outread_readID_filename, sep = "\t", quote = FALSE)

  print("Done!")
  return(out_reads)

})



# 7) Sort ReadID of seperated reads ---------------------------------------------------------

# Merge all classified readIDs to generate Tsc_readIDs.txt and Nalm6_readIDs.txt
file_pattern <- "_out_read_readID.txt"  # Match file name pattern

readID_files <- list.files(log_dir, pattern = file_pattern, full.names = TRUE)
# Read and merge all files
readID_files_combined <- readID_files %>%
  lapply(read.delim, header = TRUE, sep = "\t") %>%
  bind_rows()

# Filter out rows where readcall is not "other"
filtered_readID_files_combined <- readID_files_combined %>%
  filter(readcall != "other")

# Group by readcall and save readID to different files
tdt_readIDs <- filtered_readID_files_combined %>%
  filter(readcall == "Tsc") %>%
  pull(readID)  # Extract the readID column

nalm6_readIDs <- filtered_readID_files_combined %>%
  filter(readcall == "Nalm6") %>%
  pull(readID)  

writeLines(tdt_readIDs, paste0(outpath, "Tsc_readIDs.txt"))
writeLines(nalm6_readIDs, paste0(outpath, "Nalm6_readIDs.txt"))

cat("Processing completed! The results have been saved as 'Tsc_readIDs.txt' and 'Nalm6_readIDs.txt'!\n")




# 8.Results -------------------------------------------------------------------------------------

tmp_rds_path <- paste0(outpath,"tmp.rds")
saveRDS(out_list, file = tmp_rds_path, compress = FALSE)
print("Combining output")
out <- rbindlist(out_list)
print("Done!")

gc(verbose = FALSE)

#print("nVars_distribution_plot")
# Statistics on the distribution of SNPs on each read
nVars_distribution_plot <- ggplot(out, aes(x = nVars)) +
  geom_bar(stat = "count", fill = "skyblue", color = "black") +
  labs(title = "Distribution of nVars", x = "nVars", y = "Frequency") +
  theme_minimal() +
  geom_text(stat = "count", aes(label = ..count..), vjust = -0.5)  # Add numeric labels

nVars_distribution_plot_path <- paste0(log_dir, "/nVars_distribution_plot.png")
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
