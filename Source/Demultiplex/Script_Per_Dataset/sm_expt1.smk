import os.path
import pandas

# module load python/3.8.3


## Example run commands 
# cd /srv/scratch/z3371724/StaphCLASH2020/Source/Demultiplex
# qsub -I -N DMMiseq -l nodes=1:ppn=12,walltime=12:00:00,mem=60gb
# snakemake --dry-run --snakefile merge_demultiplexv4.smk
# snakemake --printshellcmds -j 12 --snakefile merge_demultiplexv3.smk  get_max_reads
# snakemake -j 12 --printshellcmds --snakefile merge_demultiplexv3.smk --cluster-config cluster.json \
# --cluster "qsub -N {cluster.N} -l nodes={cluster.nodes}:ppn={cluster.ppn},walltime={cluster.walltime},mem={cluster.mem} -M {cluster.email} -m ae -j oe"


# flexbar -r forward.fastq -p reverse.fastq -qf i1.8 -n 10 -ao 7 --adapters three_adapters.fasta --adapters2 five_adapters_reverse_complemented.fa --output-reads forward_trimmed.fastq --output-reads2 reverse_trimmed.fastq -qt 30 -ap ON
 

## Cluster parameters 
BASE_DIR="/srv/scratch/treelab/igy/Projects/StaphCLASH2020"
BBDUK_NCPUS=12

## Input Samples Parameter
FASTQ_DIR="/srv/scratch/treelab/JJT_fastq/stuart_mackellar_staph_clash_2017/JKD6009/r1"
OUTPUT_DIR="/srv/scratch/treelab/igy/Data/CLASH_JKD6009/SM1/Processed_Data" # "/srv/scratch/treelab/igy/Data/CLASH_JKD6009/DM_MiSeq/test_snakemake"

FILES_R1=[ "170221_K00166_0184_BHHMMVBBXX_8_SG-PE-013_1.fastq"]
FILES_R2=[ "170221_K00166_0184_BHHMMVBBXX_8_SG-PE-013_2.fastq"]

# /srv/scratch/treelab/JJT_fastq/stuart_mackellar_staph_clash_2017/JKD6009/r1/170221_K00166_0184_BHHMMVBBXX_8_SG-PE-013_1.fastq
# /srv/scratch/treelab/JJT_fastq/stuart_mackellar_staph_clash_2017/JKD6009/r1/170221_K00166_0184_BHHMMVBBXX_8_SG-PE-013_2.fastq 
# /srv/scratch/treelab/JJT_fastq/stuart_mackellar_staph_clash_2017/JKD6009/r2/170221_K00166_0184_BHHMMVBBXX_8_SG-PE-019_1.fastq
# /srv/scratch/treelab/JJT_fastq/stuart_mackellar_staph_clash_2017/JKD6009/r2/170221_K00166_0184_BHHMMVBBXX_8_SG-PE-019_2.fastq
# /srv/scratch/treelab/JJT_fastq/stuart_mackellar_staph_clash_2017/JKD6009/r3/CLASH5_1.fastq
# /srv/scratch/treelab/JJT_fastq/stuart_mackellar_staph_clash_2017/JKD6009/r3/CLASH5_2.fastq

# 
# # R1
# NNNTGTCAC	RNaseIII_TSB
# NNNACAGTG	RNaseIII_RPMI
# # R2
# NNNCGTGAT	RNase_III_TSB
# NNNGCACTA	RNase_III_RPMI_0_ligation
# NNNTAGTGC	RNase_III_RPMI_60_ligation
# NNNATCACG	RNase_III_RPMI_120_ligation
# # R3
# NNNTAGTGC	RNaseIII_TSB
# NNNATCACG	RNaseIII_RPMI

ORGANISM="Sa_JKD6009"
REF_DIR= os.path.join(BASE_DIR, "Data", "Ref_Seq") 
LINKER_FILE = os.path.join(REF_DIR, "linker.fa")
BARCODE_FILE = os.path.join(REF_DIR, "SM1_CLASH_barcodes.tab")
BARCODE_NAMES = [ "NNNTGTCAC_RNaseIII_TSB", "NNNACAGTG_RNaseIII_RPMI" ]

## Genomic parameters 
GENOME_DIR="/srv/scratch/treelab/genome/saus/JKD6009"
GTF_FILE_DIR =  os.path.join(GENOME_DIR, "Hyb_pyCRAC")
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
    compressed_files = expand( os.path.join(OUTPUT_DIR, "Compressed", "all_{barcode}.compressed.pyCRAC.fasta"), barcode=BARCODE_NAMES),
    hyb_detect_output = expand( os.path.join(OUTPUT_DIR, "Hyb_Detect", "all_{barcode}.compressed.pyCRAC_comp_Sa_JKD6009_hybrids.hyb"), barcode=BARCODE_NAMES),
    hyb_merge_output = expand( os.path.join(OUTPUT_DIR, "Hyb_Merge", "all_{barcode}.compressed.pyCRAC_comp_Sa_JKD6009_hybrids.merged.hyb"), barcode=BARCODE_NAMES),
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
                                "all_{barcode}.compressed.pyCRAC_comp_Sa_JKD6009_hybrids_annot.hyb" ),
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
