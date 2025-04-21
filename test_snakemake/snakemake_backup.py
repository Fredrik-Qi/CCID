# 用于存储此前snakemake的所有参数，现在从新开始添加规则
# 1.Read the config file ------------------------------------------------------

import yaml
import os

Configfile = "/mnt/data/user/zhouhui/T_cells_dymamics_pipeline/config/Tcell_dynamics.run.yaml"  
with open(Configfile, 'r') as f:
    config = yaml.safe_load(f)



# 1.1 Input files from config 基础信息 ----------------------------------------------
# 输入文件和数据
PROJECT_NAME = config['project']

GTF_FILE = config['reference']['GTF_file']

THREADS = config['num_threads']


OUT_DIR = config['out_dir']
ZUMIS_PATH = config['zUMIs_directory']

READ_SEPARATION_RESULTS_DIR = os.path.join(OUT_DIR, "read_separation_results")
os.makedirs(READ_SEPARATION_RESULTS_DIR, exist_ok=True)

SPLIT_BAM_OUTPUT_DIR = os.path.join(OUT_DIR, "split_bam_output")
os.makedirs(SPLIT_BAM_OUTPUT_DIR, exist_ok=True)

CELLBCs_FOR_READSEPARATION = config['cellBCs_FOR_READSEPARATION']
ANNOTATION = config['annotation']


# 代码
SNP_distribution_per_gene_script = config['Script']['SNP_distribution_per_gene_script']
split_bam_script = config['Script']['split_bam_script']
Doublets_separation_script = config['Script']['Doublets_separation_script']
Plotting_script = config['Script']['Plotting_script']
Sub_bam_script = config['Script']['Sub_bam_script']
#SPLIT_BAM_SCRIPT = os.path.join(ZUMIS_PATH, 'split_bam_by_cellBCs.sh') #还需要写一下


# 1.2 Input files from config 新生成的信息 ----------------------------------------------

# 1）zUMIs 产生的bam文件的路径
ZUMIS_BAM = os.path.join(OUT_DIR, f"{PROJECT_NAME}.filtered.Aligned.GeneTagged.UBcorrected.sorted.bam")

# 2）bam_split 产生的文件
SPLIT_BAM_OUTPUT_Tcell = os.path.join(SPLIT_BAM_OUTPUT_DIR, f"{PROJECT_NAME}.Tsc.bam")
SPLIT_BAM_OUTPUT_NALM6 = os.path.join(SPLIT_BAM_OUTPUT_DIR, "Nalm6sc_subset_stats_all_0.txt")

# 3）SNP calling
#VCF_FILE = ... 待写

# 4）Parameter processing and combination ---------------------------------------


## 4）Output files for Rule SNP_distribution_per_gene ----------------------
SNP_COUNTS_PER_GENE_PATH = os.path.join(READ_SEPARATION_RESULTS_DIR, "snp_counts_per_gene.csv")
COUNTS_DISTRIBUTION_PATH = os.path.join(READ_SEPARATION_RESULTS_DIR, "snp_count_frequency_distribution.pdf")
LOG_PATH_1 = os.path.join(READ_SEPARATION_RESULTS_DIR, "SNP_distribution_per_gene.processing_log.txt")

# 5）Output files for Rule doublets_separation -----------------------------
TEM_RDS_PATH = os.path.join(READ_SEPARATION_RESULTS_DIR, "Tcell_dynamics.tmp.rds")
COUNTS_RDS_PATH = os.path.join(READ_SEPARATION_RESULTS_DIR, "Tcell_dynamics.separated_counts.rds")
NALM6_READID_PATH = os.path.join(READ_SEPARATION_RESULTS_DIR, "Tcell_dynamics.Nalm6_readIDs.txt")
TDT_READID_PATH = os.path.join(READ_SEPARATION_RESULTS_DIR, "Tcell_dynamics.TdT_readIDs.txt")
LOG_PATH_2 = os.path.join(READ_SEPARATION_RESULTS_DIR, "doublets_separation.processing_log.txt")

# 6）Output files for Rule doublets_separation_plotting -----------------------------
READ_ASSIGNMENT_RESULT = os.path.join(READ_SEPARATION_RESULTS_DIR, "final_read_assignment_result_sorted.txt")
SCATTER_PLOT_HTML = os.path.join(READ_SEPARATION_RESULTS_DIR, "scatterplots_summary.html")
VIOLLIN_PLOT_HTML = os.path.join(READ_SEPARATION_RESULTS_DIR, "violin_plots_summary.html")
LOG_PATH_3 = os.path.join(READ_SEPARATION_RESULTS_DIR, "doublets_separation_plotting.processing_log.txt")

# 7) Output files for Rule read_separation_for_bam
NALM6_SEPARATED_BAM_PATH = os.path.join(READ_SEPARATION_RESULTS_DIR, "Tcell_dynamics.Nalm6_readIDs_separated_bam.bam")
TDT_SEPARATED_BAM_PATH = os.path.join(READ_SEPARATION_RESULTS_DIR, "Tcell_dynamics.TdT_readIDs_separated_bam.bam")



# 2.Rule to process the VCF and BAM files and associate SNPs with genes
rule all:
    input:
        run_zUMIs,
        extract_barcodes_and_split_bam,
        SNP_calling,
        SNP_COUNTS_PER_GENE_PATH,
        COUNTS_DISTRIBUTION_PATH,
        TEM_RDS_PATH,
        COUNTS_RDS_PATH,
        NALM6_READID_PATH,
        TDT_READID_PATH,
        NALM6_SEPARATED_BAM_PATH,
        TDT_SEPARATED_BAM_PATH,
        READ_ASSIGNMENT_RESULT,
        SCATTER_PLOT_HTML,
        VIOLLIN_PLOT_HTML

rule run_zUMIs:
    input:
        config_file = Configfile
    params:
        script = ZUMIS_SCRIPT
    output:
        "zUMIs_output.done",
        ZUMIS_BAM
    log:
        "logs/zUMIs.log"
    shell:
        "{params.script} -c -y {input.config_file} > {log} 2>&1 && touch {output}"
# 在这一步中我们得到的是zUMIs_output的结果，所以下一步就是根据这一结果的路径，以及bam的结果，对于bam文件进行
# 拆分，然后才能送到后面的文件中进行SNP calling
       


rule extract_barcodes_and_split_bam:
    input:
        barcode = ANNOTATION  
    output:
        SPLIT_BAM_OUTPUT_NALM6
    params:
        script = split_bam_script,  
        output_path = SPLIT_BAM_OUTPUT_DIR  
    log:
        "logs/extract_barcodelist.log"
    shell:
        """
        chmod +x {params.script}
        bash {params.script} {input.barcode} {params.output_path}
        """



# 这一个规则正在修改
#rule SNP_calling:
#    input:
#        VCF_FILE
#    output:
#        SNP_calling_results  # 这个是结果，也是其它的开始，后续的vcf这个可以取消掉
#    log:
#        "logs/SNP_calling.log"
#    shell:
#        这部"分需要重新写一个东西"

rule SNP_distribution_per_gene:
    input:
        vcf = VCF_FILE,
        gtf = GTF_FILE
    params:
        output_dir = READ_SEPARATION_RESULTS_DIR,
        script = SNP_distribution_per_gene_script
    output:
        SNP_COUNTS_PER_GENE_PATH,
        COUNTS_DISTRIBUTION_PATH
    log:
        LOG_PATH_1
    shell:
        """
        Rscript {params.script} --vcf {input.vcf} --gtf {input.gtf} \
        --outputDir {params.output_dir} 
        """        

rule doublets_separation:
    input:
        yaml = Configfile,
        vcf = VCF_FILE,
        bam = ZUMIS_BAM,
        barcodelist = CELLBCs_FOR_READSEPARATION,
        annotation = ANNOTATION
    params:
        subfolder = READ_SEPARATION_RESULTS_DIR, 
        script = Doublets_separation_script
    output:
        TEM_RDS_PATH,
        COUNTS_RDS_PATH,
        NALM6_READID_PATH,
        TDT_READID_PATH
    log:
        LOG_PATH_2
    shell:
        """
        Rscript {params.script} --yaml {input.yaml} --vcf {input.vcf} --bam {input.bam} \
        -c {input.barcodelist} -f {params.subfolder} -a {input.annotation}
        """

rule doublets_separation_plotting:
    input:
        rds = TEM_RDS_PATH,
        annotation = ANNOTATION
    params:
        output_dir = READ_SEPARATION_RESULTS_DIR,
        script = Plotting_script
    output:
        READ_ASSIGNMENT_RESULT,
        SCATTER_PLOT_HTML,
        VIOLLIN_PLOT_HTML
    log:
        LOG_PATH_3
    shell:
        """
        Rscript {params.script} --rds {input.rds} --annotation {input.annotation} \
        --outdir {params.output_dir} 
        """      

rule read_separation_for_bam:
    input:
        nalm6_readID = NALM6_READID_PATH,
        tdt_readID = TDT_READID_PATH,
        bam = ZUMIS_BAM
    params:
        output_dir = READ_SEPARATION_RESULTS_DIR,
        script = Sub_bam_script
    output:
        NALM6_SEPARATED_BAM_PATH,
        TDT_SEPARATED_BAM_PATH
    shell:
        """
        chmod +x {params.script}
        {params.script} {input.nalm6_readID} {input.bam} {params.output_dir} 50
        {params.script} {input.tdt_readID} {input.bam} {params.output_dir} 50
        """      
