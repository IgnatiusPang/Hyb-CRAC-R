
# Possible Improvements:                              
# See here on how to best load modules within each rule, rather than embedding 'module load' into the shell section
# envmodules:
#     "bio/bwa/0.7.9",
#     "bio/samtools/1.9"
# https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html#using-environment-modules                              

# * merge_p_values.Rmd to read in a parameter file, the first column is the data type label, the second column is the file path
# * Annotation code rewritten in R and compare output with old code to make sure they are the same
# * Fix the pyCRAC ruffus pipeline. 

## To Fix:
# Genome name in file names need to be a parameter

rule run_bbmerge:
  input:
    sample_r1 = expand(os.path.join( FASTQ_DIR, "{fastq}"), fastq=FILES_R1),
    sample_r2 = expand(os.path.join( FASTQ_DIR, "{fastq}"), fastq=FILES_R2),
    linker = LINKER_FILE
  output:
    all_fa = os.path.join(OUTPUT_DIR, "Merged", "all.fq"),
    ihist = os.path.join(OUTPUT_DIR, "Merged", "ihist.txt"),
    perc_merged = os.path.join(OUTPUT_DIR, "Merged", "perc_merged.log")
  shell:
    """
    module load java/8u231-jre
    module load bbmap/38.51
    
    cd $TMPDIR/
    
    if [[ ! -d {OUTPUT_DIR}/Merged ]]; then
      mkdir -p {OUTPUT_DIR}/Merged
    fi
    
    ## Check zipping of r1 
    if file --mime-type "{input.sample_r1}" | grep -q gzip$; then
      echo "{input.sample_r1} is gzipped"
      cat "{input.sample_r1}" > $TMPDIR/sample_r1.gz
    else
      echo "{input.sample_r1} is not gzipped"
      cp "{input.sample_r1}" "$TMPDIR/sample_r1"
      gzip "$TMPDIR/sample_r1"
    fi
    
    
    ## Check zipping of r2 
    if file --mime-type "{input.sample_r2}" | grep -q gzip$; then
      echo "{input.sample_r2} is gzipped"
      cat "{input.sample_r2}" > $TMPDIR/sample_r2.gz
    else
      echo "{input.sample_r2} is not gzipped"
      cp "{input.sample_r2}" "$TMPDIR/sample_r2"
      gzip "$TMPDIR/sample_r2"
    fi

    ## Merge the reads
    bbmerge.sh \
    --tossbrokenreads \
    t={BBDUK_NCPUS} \
    in1=$TMPDIR/sample_r1.gz \
    in2=$TMPDIR/sample_r2.gz \
    adapters={input.linker} \
    out=merged.fq \
    outu=unmerged.fq ihist={output.ihist} > {output.perc_merged} 2>&1
    
    
    # This is the example of the string we are trying to parse: Joined:                 447526          2.656%
    perc_merged=$(cat {output.perc_merged}  |  grep 'Joined' | sed -E 's/.*\s+([0-9\.]+)%/\\1/')
    
    PERC_MERGED_THRESHOLD={PERC_MERGED_THRESHOLD}

    ## If percentage merge threshold is not met, just use R1 reads as the input data for demultiplexing
    if (( $(echo "$perc_merged > $PERC_MERGED_THRESHOLD" | bc -l) )); then
      echo "Use bbmerge"
      # Inspirations from this gist code:
      # https://gist.github.com/nathanhaigh/3521724
      cat unmerged.fq | paste - - - - | awk '$2 !~ /^2:/ {{print $0}}' | tr "\\t" "\\n" > unmerged_r1.fq 
  
      cat merged.fq unmerged_r1.fq > {output.all_fa}
  
      rm merged.fq
      rm unmerged.fq
      rm unmerged_r1.fq
      
      rm $TMPDIR/sample_r1.gz
      rm $TMPDIR/sample_r2.gz
    else   
      gunzip  $TMPDIR/sample_r1.gz
      cp sample_r1 {output.all_fa}
      
      rm $TMPDIR/sample_r1
      rm $TMPDIR/sample_r2.gz
    fi

    """


rule run_dmplx:
  input:
    all_fa = os.path.join(OUTPUT_DIR, "Merged", "all.fq"),
    barcodes_file = BARCODE_FILE,
  output:
    os.path.join(OUTPUT_DIR, "Demultiplexed", "all_others.fastq"),
    expand( os.path.join(OUTPUT_DIR, "Demultiplexed", "all_{barcode}.fastq"), barcode=BARCODE_NAMES),
    os.path.join(OUTPUT_DIR, "Demultiplexed", "all_barcode_statistics.txt"),
    os.path.join(OUTPUT_DIR, "Demultiplexed", "all_random_nucleotide_statistics.txt"),
  shell:
    """
    module load python/2.7.15
    
    ## Activate python environment for pyCRAC scripts
    . {CRAC_ENVIRONMENT_ACTIVATE}
    
    pip -V
    
    if [[ ! -d {OUTPUT_DIR}/Demultiplexed ]]; then
      mkdir -p {OUTPUT_DIR}/Demultiplexed
    fi
    
    ### Activate the python library
    cd {OUTPUT_DIR}/Demultiplexed
    pyBarcodeFilter.py -f {input.all_fa} -b {input.barcodes_file}
   """



rule run_bbduk:
  input:
     fastq_file=  os.path.join(OUTPUT_DIR, "Demultiplexed", "all_{barcode}.fastq")
  output:
    out_file = os.path.join(OUTPUT_DIR,  "Trimmed", "all_{barcode}.trimmed.pyCRAC.fq"),
    log_file = os.path.join(OUTPUT_DIR,  "Trimmed", "all_{barcode}.trimmed.pyCRAC.log"),
    linked_file=os.path.join(OUTPUT_DIR,  "pyCRAC", "all_{barcode}.trimmed.pyCRAC.fastq")
  shell:
    """
    module load java/8u231-jre
    module load bbmap/38.51
    
    if [[ ! -d {OUTPUT_DIR}/Trimmed ]]; then
      mkdir -p {OUTPUT_DIR}/Trimmed
    fi
    
    bbduk.sh in={input.fastq_file} \
    out={output.out_file} \
    ref={REF_DIR}/phix174_ill.ref.fa,\
{REF_DIR}/linker.fa,\
{REF_DIR}/phix_adapters.fa,\
{REF_DIR}/illumina_adapter.fa,\
{REF_DIR}/poly.fa \
    ktrim=r \
    mink=11 \
    hdist=1 hdist2=1 \
    k=23 \
    mlf=0.13 minlen=20 ordered zl=4 \
    qtrim=rl trimq=28 \
    entropy=0.2 entropywindow=20 entropyk=5 \
    t={BBDUK_NCPUS} > {output.log_file} 2>&1
    
      if [[ ! -d {OUTPUT_DIR}/pyCRAC ]]; then
        mkdir -p {OUTPUT_DIR}/pyCRAC
      fi   
      
      
      if [[ ! -L {output.linked_file} ]]; then
        echo "Symbolic link do not exists, creating it now."
        cp {output.out_file} {output.linked_file}
      fi      

    """


rule run_fastqc:  
  input:
    trimmed_file =os.path.join(OUTPUT_DIR, "Trimmed", "all_{barcode}.trimmed.pyCRAC.fq"),
  output:
    fastqc_file = os.path.join(OUTPUT_DIR, "FastQC", "all_{barcode}.trimmed.pyCRAC" + "_fastqc.html"),
    zip_fastqc_file = os.path.join(OUTPUT_DIR, "FastQC", "all_{barcode}.trimmed.pyCRAC" + "_fastqc.zip"),
  shell:
    """
        module load fastqc

        if [[ ! -d {OUTPUT_DIR}/FastQC ]]; then
          mkdir -p {OUTPUT_DIR}/FastQC
        fi
        
        cd  {OUTPUT_DIR}/FastQC

        fastqc {input.trimmed_file} --outdir={OUTPUT_DIR}/FastQC
    """
  

rule run_multiqc:  
  input:
    fastqc_files=expand( os.path.join(OUTPUT_DIR, "FastQC", "all_{barcode}.trimmed.pyCRAC" + "_fastqc.html"), barcode=BARCODE_NAMES),
  output:
    multiqc_report=  os.path.join(OUTPUT_DIR, "MultiQC", "multiqc_report.html"),
    multiqc_data_dir= directory(os.path.join(OUTPUT_DIR, "MultiQC", "multiqc_data")),
  shell:
    """
      module load multiqc

      if [[ ! -d {OUTPUT_DIR}/MultiQC ]]; then
        mkdir -p {OUTPUT_DIR}/MultiQC
      fi    
    
      cd {OUTPUT_DIR}/MultiQC
      multiqc {OUTPUT_DIR}/FastQC
    """

  
rule run_pcr_de_dup:
  input:
    trimmed_file=os.path.join( OUTPUT_DIR, "Trimmed", "all_{barcode}.trimmed.pyCRAC.fq" )
  output:
    compressed_fasta= os.path.join( OUTPUT_DIR, 
                                   "Compressed", 
                                   "all_{barcode}.compressed.pyCRAC.fasta" )
  shell:
    """
        module load python/2.7.15
    
        ## Activate python environment for pyCRAC scripts
        . {CRAC_ENVIRONMENT_ACTIVATE}
    
        if [[ ! -d {OUTPUT_DIR}/Compressed ]]; then
          mkdir -p {OUTPUT_DIR}/Compressed
        fi
    
      pyFastqDuplicateRemover.py -f {input.trimmed_file} -o {output.compressed_fasta}
    """
  
## hyb_detect can generate empty *_ua.hyb output file, especially for control files where intput are small.
rule run_hyb_detect:
  input: 
    compressed_fasta = os.path.join( OUTPUT_DIR, 
                                   "Compressed", 
                                   "all_{barcode}.compressed.pyCRAC.fasta" )
  output:
    hyb_detect_output =  os.path.join( OUTPUT_DIR, 
                                       "Hyb_Detect", 
                                       "all_{barcode}.compressed.pyCRAC_comp_" + ORGANISM + "_hybrids.hyb")                           
  shell:
    """
      module unload python
      module load python/2.7.15 ##?? Not sure if this is right
      module add blat/35
      module add blast+/2.8.1
      module add bowtie/2.3.4.2
      module add fastx/0.0.14
      module add fastqc/0.11.7
      module add pblat/2.0
      module add viennarna/2.4.9
      module add unafold/3.9-1
      module unload intel/2020.0.088
      module add intel/19.0.0.117 ## Igy added, required for flexbar
      module add flexbar/3.5.0
      module add perl/5.28.0
      module add hyb/20160304 ## Igy added, required for hyb    
    
      if [[ ! -d {OUTPUT_DIR}/Hyb_Detect ]]; then
        mkdir -p {OUTPUT_DIR}/Hyb_Detect
      fi       
    
      echo "## Make blastdb index"
      if [[ ! -f {GENOME_DIR}/Sa_JKD6009.nhr ]]; then
        make_hyb_db {GENOME_DIR}/Sa_JKD6009.fasta
      fi 
      
      ## Export the blastdb index location
      export HYB_DB={GENOME_DIR}
      
      ## Run Hyb detect on all the files 
      echo "## Run Hyb Detect"
      
      cd {OUTPUT_DIR}/Hyb_Detect

      hyb detect db={HYB_DB_NAME} anti=1 format=fasta in={input.compressed_fasta} > {OUTPUT_DIR}/Hyb_Detect/hyb_detect.{wildcards.barcode}.log 2>&1
      

    """  


rule run_hyb_merge:
  input:
    hyb_input_file =  os.path.join(OUTPUT_DIR, "Hyb_Detect", 
                      "all_{barcode}.compressed.pyCRAC_comp_" + ORGANISM + "_hybrids.hyb" )
  output:
    hyb_merge_output = os.path.join(OUTPUT_DIR, "Hyb_Merge", 
                       "all_{barcode}.compressed.pyCRAC_comp_" + ORGANISM + "_hybrids.merged.hyb")
  shell:
    """
          
      module unload python
      module load python/2.7.15 ##?? Not sure if this is right
      module add blat/35
      module add blast+/2.8.1
      module add bowtie/2.3.4.2
      module add fastx/0.0.14
      module add fastqc/0.11.7
      module add pblat/2.0
      module add viennarna/2.4.9
      module add unafold/3.9-1
      module unload intel/2020.0.088
      module add intel/19.0.0.117 ## Igy added, required for flexbar
      module add flexbar/3.5.0
      module add perl/5.28.0
      module add hyb/20160304 ## Igy added, required for hyb    
    
      if [[ ! -d {OUTPUT_DIR}/Hyb_Merge ]]; then
          mkdir -p {OUTPUT_DIR}/Hyb_Merge
      fi   

    /srv/scratch/z3371724/StaphCLASH2020/Source/Jai_Script/Annotation/combine_hyb_merge_1.3.pl TWO_WAY_MERGE=1 \
      {input.hyb_input_file}  > {output.hyb_merge_output}
      
    """

rule run_pyCRAC:
  input:
    linked_file=os.path.join(OUTPUT_DIR, "pyCRAC", "all_{barcode}.trimmed.pyCRAC.fastq"),
    pyCRAC_GTF_FILE=pyCRAC_GTF_FILE,
    GENOME_SIZE=GENOME_SIZE,
    NOVO_INDEX=NOVO_INDEX
  output:
    sgr_positive =  os.path.join(OUTPUT_DIR,   "pyCRAC",
                   "all_{barcode}" , "sgr_files",
                    "all_{barcode}.trimmed.pyCRAC_trimmed_blocks_nomuts_plus_strand_reads.sgr"),
    sgr_negative = os.path.join(OUTPUT_DIR,   "pyCRAC",
                   "all_{barcode}", "sgr_files",
                    "all_{barcode}.trimmed.pyCRAC_trimmed_blocks_nomuts_minus_strand_reads.sgr"),
    pyCRAC_log = os.path.join(OUTPUT_DIR,  "pyCRAC" , "all_{barcode}.crac_log_file.txt"), ## This will keep the log file 
    num_mapped_reads_file = os.path.join( OUTPUT_DIR,
                  "pyCRAC",
                  "all_{barcode}",
                  "pyReadCounters_blocks_nomuts_analyses",
                  "all_{barcode}.trimmed.pyCRAC_trimmed_blocks_nomuts_file_statistics_cDNAs.txt"),
  shell:
    """
    
     # module load python/2.7.15
      module unload intel/2020.0.088
      module load intel/19.0.0.117
      module load novocraft/3.04.06
      module load flexbar/3.5.0

    ## Activate python environment for pyCRAC scripts
    . {CRAC_ENVIRONMENT_ACTIVATE}
    

      if [[ ! -d {OUTPUT_DIR}/pyCRAC ]]; then
        mkdir -p {OUTPUT_DIR}/pyCRAC
      fi

      cd {OUTPUT_DIR}/pyCRAC
      
      INPUT_FILE={input.linked_file}
      FILE_NAME=${{INPUT_FILE##*/}}


      ## Input file must be in current directory and have the .fastq suffix
      ## The CRAC_pipeline_SE script itself re-run PCR de-duplication step with 'pyFastqDuplicateRemover.py'
      {PYCRAC_SRC_DIR}/CRAC_pipeline_SE.py \
      -f $FILE_NAME \
      -g {pyCRAC_GTF_FILE} \
      -c {GENOME_SIZE} \
      --novoindex={NOVO_INDEX} \
      --name="all_{wildcards.barcode}" \
      -a AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT \
      -p 1 > {output.pyCRAC_log} 2>&1

    """


rule  run_sgr2pws:
  input:
      plus_strand_file=os.path.join(OUTPUT_DIR, "pyCRAC", "all_{barcode}", "sgr_files",
              "all_{barcode}.trimmed.pyCRAC_trimmed_blocks_nomuts_plus_strand_reads.sgr"),
      minus_strand_file=os.path.join(OUTPUT_DIR, "pyCRAC", "all_{barcode}", "sgr_files",
              "all_{barcode}.trimmed.pyCRAC_trimmed_blocks_nomuts_minus_strand_reads.sgr")
  output:
      pws_files = os.path.join(OUTPUT_DIR, "pyCRAC", "all_{barcode}", "pws_files", "all_{barcode}.pws"),
  shell:
    """
      # module load python/3.8.2

     ## Activate python environment for pyCRAC scripts
     . {CRAC_ENVIRONMENT_ACTIVATE}
    
      MY_OUTPUT_DIR="{OUTPUT_DIR}/pyCRAC/all_{wildcards.barcode}/pws_files"

      if [[ ! -d $MY_OUTPUT_DIR ]]; then
        mkdir -p $MY_OUTPUT_DIR
      fi

      /srv/scratch/z3371724/StaphCLASH2020/Source/Jai_Script/Hyb_stats/sgr2pws.py \
       -p  {input.plus_strand_file} \
       -n  {input.minus_strand_file} > \
       {output.pws_files}
    """
    
 
## Annotation  
rule run_annotation:
  input:
     HYB_FILE = os.path.join(OUTPUT_DIR, "Hyb_Merge", 
                  "all_{barcode}.compressed.pyCRAC_comp_" + ORGANISM + "_hybrids.merged.hyb"),
     GTF_FILE = HYB_GTF_FILE,
     RNA_TYPE_FILE = RNA_TYPE_PRIORITIES
     CHROMOSOME_FEATURE_FILE = CHROMOSOME_FEATURE_FILE
  output: 
    hyb_annot_output = os.path.join(OUTPUT_DIR, "Hyb_Annot", "all_{barcode}.compressed.pyCRAC_comp_" + ORGANISM + "_hybrids_annot.hyb"),
  envmodules:
      "R/3.6.3"
  script:
      "/srv/scratch/treelab/igy/Projects/StaphCLASH2020/Source/Demultiplex/scripts/annotation.R"          
  #   
  # shell:
  #   """
  #     module load python/3.8.2
  # 
  #     if [[ ! -d {OUTPUT_DIR}/Hyb_Annot ]]; then
  #       mkdir -p {OUTPUT_DIR}/Hyb_Annot
  #     fi  
  #     
  #     HYB={input.HYB_FILE}
  # 
  #     name=${{HYB%%.hyb}}
  #     prefix_merged=${{name##*/}}
  #     prefix_annot=${{prefix_merged%%.merged}}
  #     
  #     ## {GTF_FILE_DIR} and {ANNOTATION_SCRIPT_DIR} are hard coded variables 
  #     
  #     # ncRNA > CDS > 5UTR = 3UTR > rRNA > tRNA > others > operon > transcript?? 
  #     #  
  #     # CRISPR
  #     # pseudogene
  #     # riboswitch
  #     # RNase_P_RNA
  #     # RNA_thermometer
  #     # SRP_RNA
  #     # tmRNA
  #     
  #     {ANNOTATION_SCRIPT_DIR}/annotate_hyb_2_edited.py $HYB -g {GTF_FILE_DIR}/other_features_ncRNA.gtf > $TMPDIR/$prefix_annot\_file1
  #     {ANNOTATION_SCRIPT_DIR}/annotate_hyb_2_edited.py $TMPDIR/$prefix_annot\_file1 -g {GTF_FILE_DIR}/main_features_CDS.gtf > $TMPDIR/$prefix_annot\_file2
  #     {ANNOTATION_SCRIPT_DIR}/annotate_hyb_2_edited.py $TMPDIR/$prefix_annot\_file2 -g {GTF_FILE_DIR}/main_features_5UTR.gtf > $TMPDIR/$prefix_annot\_file3
  #     {ANNOTATION_SCRIPT_DIR}/annotate_hyb_2_edited.py $TMPDIR/$prefix_annot\_file3 -g {GTF_FILE_DIR}/main_features_3UTR.gtf > $TMPDIR/$prefix_annot\_file4
  #     {ANNOTATION_SCRIPT_DIR}/annotate_hyb_2_edited.py $TMPDIR/$prefix_annot\_file4 -g {GTF_FILE_DIR}/main_features_3UTR.gtf > $TMPDIR/$prefix_annot\_file5
  #     {ANNOTATION_SCRIPT_DIR}/annotate_hyb_2_edited.py $TMPDIR/$prefix_annot\_file5 -g {GTF_FILE_DIR}/other_features_rRNA.gtf > $TMPDIR/$prefix_annot\_file6
  #     {ANNOTATION_SCRIPT_DIR}/annotate_hyb_2_edited.py $TMPDIR/$prefix_annot\_file6 -g {GTF_FILE_DIR}/other_features_tRNA.gtf > $TMPDIR/$prefix_annot\_file7
  #     {ANNOTATION_SCRIPT_DIR}/annotate_hyb_2_edited.py $TMPDIR/$prefix_annot\_file7 -g {GTF_FILE_DIR}/other_features_CRISPR.gtf > $TMPDIR/$prefix_annot\_file8
  #     {ANNOTATION_SCRIPT_DIR}/annotate_hyb_2_edited.py $TMPDIR/$prefix_annot\_file8 -g {GTF_FILE_DIR}/other_features_pseudogene.gtf > $TMPDIR/$prefix_annot\_file9
  #     {ANNOTATION_SCRIPT_DIR}/annotate_hyb_2_edited.py $TMPDIR/$prefix_annot\_file9 -g {GTF_FILE_DIR}/other_features_riboswitch.gtf > $TMPDIR/$prefix_annot\_file10
  #     {ANNOTATION_SCRIPT_DIR}/annotate_hyb_2_edited.py $TMPDIR/$prefix_annot\_file10 -g {GTF_FILE_DIR}/other_features_RNase_P_RNA.gtf > $TMPDIR/$prefix_annot\_file11
  #     {ANNOTATION_SCRIPT_DIR}/annotate_hyb_2_edited.py $TMPDIR/$prefix_annot\_file11 -g {GTF_FILE_DIR}/other_features_RNA_thermometer.gtf > $TMPDIR/$prefix_annot\_file12
  #     {ANNOTATION_SCRIPT_DIR}/annotate_hyb_2_edited.py $TMPDIR/$prefix_annot\_file12 -g {GTF_FILE_DIR}/other_features_SRP_RNA.gtf > $TMPDIR/$prefix_annot\_file13
  #     {ANNOTATION_SCRIPT_DIR}/annotate_hyb_2_edited.py $TMPDIR/$prefix_annot\_file13 -g {GTF_FILE_DIR}/other_features_tmRNA.gtf > $TMPDIR/$prefix_annot\_file14
  #     {ANNOTATION_SCRIPT_DIR}/annotate_hyb_2_edited.py $TMPDIR/$prefix_annot\_file14 -g {GTF_FILE_DIR}/other_features_operon.gtf > {OUTPUT_DIR}/Hyb_Annot/$prefix_annot\_annot.hyb
  # 
  #     # The last pass of annotate_hyb_2_edited.py of the annotation shell activates a function that
  #     # outputs individuals files for each class of interaction eg: scRNA-mRNA
  #     # it also outputs a .counter file that tabulates how many of each interaction you have.
  # 
  #     # {ANNOTATION_SCRIPT_DIR}/annotate_hyb_2_edited.py {OUTPUT_DIR}/Hyb_Annot/$prefix_annot\_annot.hyb -g {GTF_FILE_DIR}/Sa_JKD6008_1.5_intergenic_region.gtf -o False > \
  #     #  {OUTPUT_DIR}/Hyb_Annot/$prefix_annot\.miscel
  # 
  #     echo annotation finished for $prefix_annot
  # 
  #   """



## Create datasets.txt
rule get_dataset_stats:
  input:
    num_mapped_reads_files = expand(os.path.join( OUTPUT_DIR,
                  "pyCRAC",
                  "all_{barcode}",
                  "pyReadCounters_blocks_nomuts_analyses",
                  "all_{barcode}.trimmed.pyCRAC_trimmed_blocks_nomuts_file_statistics_cDNAs.txt"),
                  barcode=BARCODE_NAMES ),
    hyb_files =  expand( os.path.join( OUTPUT_DIR, "Hyb_Annot", 
                                       "all_{barcode}.compressed.pyCRAC_comp_" + ORGANISM + "_hybrids_annot.hyb"),
                  barcode=BARCODE_NAMES ),
    pws_files = expand( os.path.join( OUTPUT_DIR, "pyCRAC", 
                                      "all_{barcode}", "pws_files", "all_{barcode}.pws"),
                  barcode=BARCODE_NAMES )
  output: 
    dataset_stats_file=os.path.join(OUTPUT_DIR, "Stats" , "datasets.txt" )
  envmodules:
      "R/3.6.3"
  script:
    "/srv/scratch/treelab/igy/Projects/StaphCLASH2020/Source/Demultiplex/scripts/format_stats_dataset.R"
    
    
## Get the max reads     
rule get_max_reads:
  input:
      ## For improvement: HYB, PWS and MAX_READS are defined below, how can I make it only declared once?
      HYB = os.path.join( OUTPUT_DIR, "Hyb_Annot", 
                          "all_{barcode}.compressed.pyCRAC_comp_" + ORGANISM + "_hybrids_annot.hyb"), 
      PWS = os.path.join( OUTPUT_DIR, "pyCRAC", 
                          "all_{barcode}", "pws_files", "all_{barcode}.pws")
  output: 
      MAX_READS = os.path.join( OUTPUT_DIR, "Hyb_Annot",
                                "all_{barcode}.hyb.max_reads.txt" )
  envmodules:
      "R/3.6.3"
  script:
    "/srv/scratch/treelab/igy/Projects/StaphCLASH2020/Source/Demultiplex/scripts/max_reads.R"                                



## Calculate the generic hybrid p-values
rule calc_hyb_p_val:
    input:
          MAX_READS= os.path.join( OUTPUT_DIR, "Hyb_Annot", "all_{barcode}.hyb.max_reads.txt"),
          ## For improvement: Best to go through the table and check all files listed in there are present, not yet implemented 
          dataset_stats_file=os.path.join(OUTPUT_DIR, "Stats" , "datasets.txt" )
    output: 
          p_values_file= os.path.join(OUTPUT_DIR, "Hyb_Annot", "all_{barcode}.hyb.max_reads.txt.p_values.txt"),
    run:
          #dataset_stats_file=os.path.join(OUTPUT_DIR, "Stats" , "datasets.txt" )
          dataset_table = pandas.read_table(input.dataset_stats_file, header=None, sep="\t")
          matching_line = (dataset_table[1]==  ("all_" + wildcards.barcode + ".pws") )
          
          print( "all_" + wildcards.barcode + ".pws")
          print( matching_line )
          
          num_mapped_reads = dataset_table[matching_line][2].values[0]
          
          print( num_mapped_reads)
          #MAX_READS=os.path.join( OUTPUT_DIR, "Hyb_Annot", "".join(["all_", wildcards.barcode, ".hyb.max_reads.txt"]) )
          rscipt_cmd = ['module load R/3.6.3; Rscript ' ,  STATS_SCRIPT_DIR ,
                              '/calculate_generic_hybrid_p_values.R ' ,
                              input.MAX_READS , ' ' ,
                              str(num_mapped_reads) , ';'] 
          calc_p_value_cmd = "".join( rscipt_cmd )
          print( calc_p_value_cmd )                              
          os.system( calc_p_value_cmd )
          
          
## Separate unique hybrid and multimapped hybrid. 
## A multimapepd hybrid is when at least one half of the hybrid is mapped to multiple location in the genome. 
## For example, duplicated gene or regions where the sequence is exactly the same.           
rule run_sep_unique_multimapped:
    input:
      p_values_file= os.path.join(OUTPUT_DIR, "Hyb_Annot", "all_{barcode}.hyb.max_reads.txt.p_values.txt")
    output:
      unique_hybrids_file  =  os.path.join(OUTPUT_DIR, "Hyb_Annot", "all_{barcode}.unique.hyb.p_values.tab" ),
      multimapped_hybrids_file = os.path.join(OUTPUT_DIR, "Hyb_Annot", "all_{barcode}.multimapped.hyb.p_values.tab" )  
    envmodules:
      "R/3.6.3"
    script:
      "/srv/scratch/treelab/igy/Projects/StaphCLASH2020/Source/Demultiplex/scripts/separate_unique_multimapped_hybrids.R"                                


