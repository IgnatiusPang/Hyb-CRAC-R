import os.path
import pandas

# module load python/3.8.3


## Example run commands 
# cd /srv/scratch/z3371724/StaphCLASH2020/Source/Demultiplex
# qsub -I -N DMMiseq -l nodes=1:ppn=12,walltime=12:00:00,mem=60gb
# snakemake --dry-run --snakefile dm_miseq.smk
# snakemake --printshellcmds -j 12 --snakefile merge_demultiplexv3.smk  get_max_reads
# snakemake -j 12 --printshellcmds --snakefile dm_miseq.smk --cluster-config cluster.json \
# --cluster "qsub -N {cluster.N} -l nodes={cluster.nodes}:ppn={cluster.ppn},walltime={cluster.walltime},mem={cluster.mem} -M {cluster.email} -m ae -j oe"

## Cluster parameters 
BASE_DIR="/home/ignatius/PostDoc/2020/Hyb-CRAC-R"
BBDUK_NCPUS=12 # Number of CPUs for BBduk trimming

## GitHub base dir

## Input Samples Parameter
FASTQ_DIR=os.path.join( BASE_DIR, "Data/CLASH_FASTQ")
OUTPUT_DIR=os.path.join( BASE_DIR, "Results/DM_MiSeq") 

FILES_R1=[ "DGM-RNase-CLASH_S1_L001_R1_001.fastq.gz"]

FILES_R2=[ "DGM-RNase-CLASH_S1_L001_R2_001.fastq.gz"]

ORGANISM="Sa_JKD6009"
REF_DIR= os.path.join(BASE_DIR, "Data", "Ref_Seq") 
LINKER_FILE = os.path.join(REF_DIR, "linker.fa")
BARCODE_FILE = os.path.join(REF_DIR, "DM_MiSeq_CLASH_barcodes.tab")
BARCODE_NAMES = [ "NNNAGAGC_L5Ba_control_1" , "NNNGTGAGC_L5Bb_control_2", 
                  "NNNCACTAGC_L5Bc_rnc-HTF_1",  "NNNTCTCTAGC_L5Bd_rnc-HTF_2" ]  

# BBmerge percentage threshold above which BBmerge is used
PERC_MERGED_THRESHOLD=10

## Genomic parameters 
GENOME_DIR=os.path.join( BASE_DIR, "Data/Genomic/JKD6009")
GTF_FILE_DIR =  os.path.join(GENOME_DIR, "Annotation/Hyb_pyCRAC")
pyCRAC_GTF_FILE = os.path.join(GTF_FILE_DIR,  "GCF_900607245.1_JKD6009_genomic_merge_features_simple_annot_checked_edited.gtf")
GENOME_SIZE =  os.path.join(GENOME_DIR,  "Sa_JKD6009.genome")
NOVO_INDEX =   os.path.join(GENOME_DIR,  "Sa_JKD6009.nix")  
HYB_DB_NAME = "Sa_JKD6009"
RNA_TYPE_PRIORITIES =  os.path.join(GTF_FILE_DIR,  "rna_type_priority.tab")

## Source Directories 
SOURCE_DIR= os.path.join( BASE_DIR, "Source" ) 
R_SCRIPTS_DIR=os.path.join(SOURCE_DIR, "Demultiplex", "scripts" ) 
ANNOTATION_SCRIPT_DIR = os.path.join(SOURCE_DIR, "Jai_Script", "Annotation" ) 
STATS_SCRIPT_DIR = os.path.join(SOURCE_DIR, "Jai_Script", "Hyb_stats" ) 
PYCRAC_SRC_DIR = os.path.join(SOURCE_DIR, "Jai_Script", "pyCRAC" ) 

## Location of CRAC installation environment
CRAC_ENVIRONMENT_ACTIVATE="/srv/scratch/z3371724/my_python_dir/my_python_env/bin/activate"

rule all:
  input:
    sample_r1 = expand(os.path.join( FASTQ_DIR, "{fastq}"), fastq=FILES_R1),
    sample_r2 = expand(os.path.join( FASTQ_DIR, "{fastq}"), fastq=FILES_R2),
    linker = LINKER_FILE,
    all_fa = os.path.join(OUTPUT_DIR,  "Merged", "all.fq"),
    ihist = os.path.join(OUTPUT_DIR, "Merged", "ihist.txt"),
    barcodes_file = BARCODE_FILE,
    fastq_file= expand(  os.path.join(OUTPUT_DIR, "Demultiplexed", "all_{barcode}.fastq"), barcode=BARCODE_NAMES),
    trimmed_files = expand( os.path.join(OUTPUT_DIR,  "Trimmed", "all_{barcode}.trimmed.pyCRAC.fq"), barcode=BARCODE_NAMES),
    fastqc_files = expand( os.path.join(OUTPUT_DIR, "FastQC", "all_{barcode}.trimmed.pyCRAC" + "_fastqc.html"), barcode=BARCODE_NAMES),
    zip_fastqc_files = expand( os.path.join(OUTPUT_DIR, "FastQC", "all_{barcode}.trimmed.pyCRAC" + "_fastqc.zip"), barcode=BARCODE_NAMES),
    multiqc_report = os.path.join(OUTPUT_DIR, "MultiQC", "multiqc_report.html"),
    compressed_files = expand( os.path.join(OUTPUT_DIR, "Compressed", "all_{barcode}.compressed.pyCRAC.fasta"), barcode=BARCODE_NAMES) ,
    hyb_detect_output = expand( os.path.join(OUTPUT_DIR, "Hyb_Detect", "all_{barcode}.compressed.pyCRAC_comp_" + ORGANISM + "_hybrids.hyb"), barcode=BARCODE_NAMES),
    hyb_merge_output = expand( os.path.join(OUTPUT_DIR, "Hyb_Merge", "all_{barcode}.compressed.pyCRAC_comp_" + ORGANISM + "_hybrids.merged.hyb"), barcode=BARCODE_NAMES),
    linked_file= expand( os.path.join(OUTPUT_DIR,  "pyCRAC", "all_{barcode}.trimmed.pyCRAC.fastq"), barcode=BARCODE_NAMES),
    num_mapped_reads_files = expand(os.path.join( OUTPUT_DIR,
                  "pyCRAC",
                  "all_{barcode}",
                  "pyReadCounters_blocks_nomuts_analyses",
                  "all_{barcode}.trimmed.pyCRAC_trimmed_blocks_nomuts_file_statistics_cDNAs.txt"),
                  barcode=BARCODE_NAMES ),
    # pyCRAC_log = expand( os.path.join(OUTPUT_DIR, "pyCRAC", "all_{barcode}.crac_log_file.txt"), barcode=BARCODE_NAMES),
    sgr_positive = expand( os.path.join(OUTPUT_DIR,  "pyCRAC",
                           "all_{barcode}", "sgr_files",
                           "all_{barcode}.trimmed.pyCRAC_trimmed_blocks_nomuts_plus_strand_reads.sgr"),
                           barcode=BARCODE_NAMES  ),
    sgr_negative = expand( os.path.join(OUTPUT_DIR,  "pyCRAC",
                           "all_{barcode}", "sgr_files",
                           "all_{barcode}.trimmed.pyCRAC_trimmed_blocks_nomuts_minus_strand_reads.sgr"),
                           barcode=BARCODE_NAMES  ),
    pws_files = expand( os.path.join(OUTPUT_DIR,  "pyCRAC",
                        "all_{barcode}", "pws_files",
                        "all_{barcode}.pws"),
                        barcode=BARCODE_NAMES  ),
    hyb_annot_output =  expand( os.path.join(OUTPUT_DIR, "Hyb_Annot",
                                "all_{barcode}.compressed.pyCRAC_comp_" + ORGANISM + "_hybrids_annot.hyb" ),
                                barcode=BARCODE_NAMES  ),
    max_reads_files = expand( os.path.join( OUTPUT_DIR, "Hyb_Annot",
                                            "all_{barcode}.hyb.max_reads.txt" ),
                                            barcode=BARCODE_NAMES  ),
    p_values_file= expand( os.path.join( OUTPUT_DIR, "Hyb_Annot", "all_{barcode}.hyb.max_reads.txt.p_values.txt"),
                               barcode=BARCODE_NAMES  ),
    unique_hybrids_file = expand( os.path.join(OUTPUT_DIR, "Hyb_Annot", "all_{barcode}.unique.hyb.p_values.tab" ),
                               barcode=BARCODE_NAMES  ),
    multimapped_hybrids = expand( os.path.join(OUTPUT_DIR, "Hyb_Annot", "all_{barcode}.multimapped.hyb.p_values.tab" ),
                               barcode=BARCODE_NAMES  )

  
 
## Include the shared rules 
include: "common_rules.smk"
