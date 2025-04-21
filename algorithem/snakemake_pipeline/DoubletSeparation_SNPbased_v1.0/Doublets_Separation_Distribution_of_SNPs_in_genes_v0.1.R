suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(rtracklayer)) # 用于读取 GTF 文件
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(ggplot2))


vcf_file <- "/mnt/nvme/home/zhouhui/project_workspace/Tcell_NALM6_project/Freebayes_multi_samples_calling/Freebayes_multi_samples_calling_snponly_simple_Nalm6_TdT_only01_different_homo_hete.vcf.gz"   # 你的 VCF 文件路径

print("Reading Variants...")
if(grepl(vcf_file, pattern = ".gz$")){
  vcf_data <- fread(cmd = paste("zcat",vcf_file," | grep -v '^##'"), header = T)
}else{
  vcf_data <- fread(cmd = paste("grep -v '^##'",vcf_file), header = T)
}

print("Reading GTF")
gtf_file <- "/mnt/nvme/home/chris/resources/genomes/Human/gencode.v42.primary_assembly.annotation.nochr.rRNA.gtf"  # 你的 GTF 文件路径
gtf <- import(gtf_file)

# ----------------------------------------------------------
print("Processing gtf...")
# 提取基因信息 (假设是"gene_id"列)
genes <- gtf[gtf$type == "gene"]  # 只保留基因信息

print("gene")
print(head(genes))  # 检查基因信息



# ----------------------------------------------------------
print("Processing Variants...")
# 提取 VCF 文件中的染色体和位置
vcf_positions <- vcf_data[, c("#CHROM", "POS")]
# 使用 setnames 修改列名
setnames(vcf_positions, old = "#CHROM", new = "CHROM")

# 使用 grepl 匹配含有 "KI" 或 "GL" 的染色体
vcf_positions_filtered <- vcf_positions[!grepl("^(KI|GL)", vcf_positions$CHROM), ]


print("vcf_positions")
print(head(vcf_positions_filtered))  # 检查 SNP 信息

print("Structure of vcf_positions")
str(vcf_positions_filtered)





# 创建 SNP 范围
snp_gr <- GRanges(seqnames = vcf_positions_filtered$CHROM,
                  ranges = IRanges(start = vcf_positions_filtered$POS, end = vcf_positions_filtered$POS))

print("snp_gr")
print(snp_gr)


# ----------------------------------------------------------


# 查找 SNP 所在的基因
snp_in_genes <- findOverlaps(snp_gr, genes)
print("snp_in_genes")
str(snp_in_genes)

# 获取重叠关系的 SNP 和基因索引
snp_idx <- queryHits(snp_in_genes)  # SNP 索引
gene_idx <- subjectHits(snp_in_genes)  # 基因索引

# 从元数据提取基因名
gene_names <- mcols(genes)$gene_name

# 统计每个基因的 SNP 数量
snp_per_gene <- table(gene_idx)

# 创建包含基因名的数据框
snp_counts_df <- data.frame(
  Gene_Index = as.numeric(names(snp_per_gene)),  # 基因索引
  Gene_Name = gene_names[as.numeric(names(snp_per_gene))],  # 基因名
  SNP_Count = as.vector(snp_per_gene)           # 对应 SNP 数量
)

# 按照 SNP_Count 降序排列，并选择前 20 个基因
top_20_genes <- snp_counts_df[order(-snp_counts_df$SNP_Count), ][1:20, ]

# 查看前 20 个基因的信息
print(top_20_genes)


# 计算 SNP_Count 列的总和
total_snp_count <- sum(snp_counts_df$SNP_Count)

# 输出总和
print("total_snp_count")
print(total_snp_count)

# 保存为 CSV 文件
write.csv(snp_counts_df, "/mnt/nvme/home/zhouhui/project_workspace/Tcell_NALM6_project/Deconvolute_algorithm/DoubletSeparation_SNPbased_v0.6.0/snp_counts_per_gene.csv", row.names = FALSE)

# 对 SNP_Count 列进行频率统计
freq_table <- table(snp_counts_df$SNP_Count)
freq_percentage <- prop.table(freq_table) * 100

# 将统计结果转换为数据框
freq_df <- data.frame(
  SNP_Count = as.numeric(names(freq_table)),  # 不同的 SNP 数量
  Frequency = as.vector(freq_table),         # 出现的次数
  Percentage = as.vector(freq_percentage)    # 出现的百分比
)

# 绘制频率分布柱状图
ggplot(freq_df, aes(x = SNP_Count, y = Frequency)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  labs(
    title = "Distribution of SNP Counts per Gene",
    x = "SNP Count",
    y = "Number of Genes"
  ) +
  theme_minimal()
# 保存为 PDF 文件
ggsave("/mnt/nvme/home/zhouhui/project_workspace/Tcell_NALM6_project/Deconvolute_algorithm/DoubletSeparation_SNPbased_v0.6.0/snp_count_frequency_distribution.pdf", width = 10, height = 6)
