In silico analysis of RNase III-RNA binding sites and regulatory RNA-RNA interactions.
1. Analysis of CLASH dataset. Analysis of the CLASH data was like those described in Waters et al. (2017), Tree et al. (2018), and Mediati et al. (2022). The pipeline for CLASH data analysis, called Hyb-CRAC-R, has been implemented in Snakemake workflow management language. Installation instructions and codes are available in the following GitHub repository (https://github.com/IgnatiusPang/Hyb-CRAC-R). An example Snakemake script ‘dm_miseq.smk’ could be found in the folder ‘Source/Demultiplex/Script_Per_Dataset’ in the GitHub repository. A command to run the pipeline is as follows, where the ‘-j 2’ options assign 2 two CPUs to perform the computational analysis.

snakemake -j 1 --printshellcmds --snakefile dm_miseq.smk 

The parameters in the ‘dm_miseq.smk’ file that need to be adjusted are listed a follows:

Cluster parameters 
BASE_DIR Is the main directory where all results and scripts are kept
BBDUK_NCPUS is the number of compute CPUs that will be used by the bbduck tool
Input Samples Parameter
FASTQ_DIR is the directory where the FASTQ sequence are stored.
OUTPUT_DIR is the directory where all the output files and directories are saved. 
FILES_R1 is a Python list with the name of the FASTQ R1 files. E.g. [ "DGM-RNase-CLASH_S1_L001_R1_001.fastq.gz"]
FILES_R2 is a Python list with the name of the FASTQ R2 files. E.g. [ "DGM-RNase-CLASH_S1_L001_R2_001.fastq.gz"]
ORGANISM is the name of the organism in the FASTA genome file e.g. Sa_JKD6009
REF_DIR is the directory containing the reference sequences such as the linker sequences and the barcode sequences
LINKER_FILE is the FASTA file containing the linker sequences e.g. os.path.join(REF_DIR, "linker.fa")
BARCODE_FILE is the FASTA file containing the barcode sequences e.g.  os.path.join(REF_DIR, "DM_MiSeq_CLASH_barcodes.tab")
BARCODE_NAMES is a Python list listing the name of the barcodes e.g. [ "NNNAGAGC_L5Ba_control_1" , "NNNGTGAGC_L5Bb_control_2", 
                  "NNNCACTAGC_L5Bc_rnc-HTF_1",  "NNNTCTCTAGC_L5Bd_rnc-HTF_2" ]  
PERC_MERGED_THRESHOLD is the BBmerge percentage threshold above which BBmerge is used (e.g. 10 )

Genomic parameters 
GENOME_DIR is the directory where the Genomic data and GTF files are stored. 
GTF_FILE_DIR is the directory where the GTF file is stored (e.g.  os.path.join(GENOME_DIR, "Hyb_pyCRAC") )
pyCRAC_GTF_FILE is the path to the GTF file used for pyCRAC analysis (e.g. os.path.join(GTF_FILE_DIR,  "GCF_900607245.1_JKD6009_genomic_merge_features_simple_annot_checked_edited.gtf") )
HYB_GTF_FILE is the path to the GTF file used for ‘hyb’ analysis (e.g. os.path.join(GTF_FILE_DIR,  "compiled_edited_features.gtf"))
GENOME_SIZE is a file containing the header of the FASTA file containing the genome (without the ‘>’ sign) followed by a tab separator and an integer representing the size of the organism’s genome  
NOVO_INDEX is the path to the Novoalign index file (e.g. os.path.join(GENOME_DIR,  "Sa_JKD6009.nix")  ) 

HYB_DB_NAME is a string representing the name of the database used in ‘hyb’ software (e.g. "Sa_JKD6009")
RNA_TYPE_PRIORITIES is a file listing the different types of RNA in order of priority for annotation (e.g. os.path.join(GTF_FILE_DIR,  "rna_type_priority.tab"))
CHROMOSOME_FEATURE_FILE A tab separated file with two columns, the column on the left list the chromosome feature names in the GTF file, the column on the right are the chromosome feature names referenced in the Hyb and genomic FASTA file 
Source Directories 
SOURCE_DIR is the path of the directory in which omost of the source codes for this project is being kept (e.g. os.path.join( BASE_DIR, "Source" ) )
R_SCRIPTS_DIR is the directory containing the R scripts (e.g. os.path.join(SOURCE_DIR, "Demultiplex", "scripts" ) )
ANNOTATION_SCRIPT_DIR is the pathway to a file containing the script for performing the annotation step 
STATS_SCRIPT_DIR is a directory containing scripts for analysis of the output files from ‘hyb’ (e.g. os.path.join(SOURCE_DIR, "Jai_Script", "Hyb_stats" ) )
PYCRAC_SRC_DIR is the local directory in which the pyCRAC source codes are kept (e. os.path.join(SOURCE_DIR, "Jai_Script", "pyCRAC" ) 

Location of CRAC installation environment
CRAC_ENVIRONMENT_ACTIVATE is the path of the Python virtual environment activation script 
 

2. Merging of pair-end reads was performed using ‘bbmerge’. The two halves of the paired-end reads were merged together using ‘bbmerge’ (Bushnell, 2014) with the following linker sequence (NAGATCGGAAGAGCACACGTCTG). Both the merged and unmerged reads were collated for demultiplexing. The following is the command line for running bbmerge.

bbmerge.sh \
    --tossbrokenreads \
    t=1 \
    in1=sample_r1.gz \
    in2=sample_r2.gz \
    adapters=linker.fa \
    out=merged.fq \
outu=unmerged.fq \
ihist=ihist.txt

The parameter ‘—tossbrokenreads’ remove broken reads, the parameters ‘in1’ and ‘in2’ indicates the two FASTQ files for the paired-end sequencing analysis in gzipped compressed format. The ‘adapter’ parameter represents the ‘linker.fa’ file in FASTA format listing all the adaptor sequences, the ‘out’ parameter indicates the output file the merged reads and the ‘outu’ parameter indicates the output file for the unmerged FASTQ reads. The ‘ihist’ parameter indicates a file in which the insert size histogram will be written.

In the following commands, the R1 reads from the ‘unmerged.fa’ FASTQ file was extracted into a new file ‘unmerged_r1.fq’ .                                           
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      
cat unmerged.fq | paste - - - - | awk '$2 !~ /^2:/ {{print $0}}' | tr "\\t" "\\n" > unmerged_r1.fq 

The next command then concatenates the file containing the merged ‘merged.fq’ and unmerged ‘unmerged_r1.fq’ reads together into one file named ‘all.fq’.

 cat merged.fq unmerged_r1.fq > all.fq

3. Demultiplexing of samples was performed using ‘pyCRAC’. Demultiplexing was performed using the ‘pyBarcodeFilter.py’ script from ‘pyCRAC’ (Webb et al., 2014), this resulted in one FASTQ format file per sample barcode. The following analysis steps were performed on each demultiplexed FASTQ file, with one FASTQ file per sample. In the following command line call to ‘pyBarcodeFilter.py’, the ‘-f’ parameter indicates the input file ‘all.fq’ and the ‘-b’ parameter indicates a tab-separated file listing the barcodes ‘barcodes.tab’. 

pyBarcodeFilter.py -f all.fq -b barcodes.tab

Each of the ouput file will have the format ‘all_{barcode}.fastq’ in which {barcode} represents the barcode sequence for all sequences tagged with the same barcode.

The barcode file is in the following tap-separated format, in which the first column represents the barcode in which the letter ‘N’ represents any one of the nucleotides C, G, T or A, and the second column represents the name of the output FASTA file containing all sequences with the barcode. 
NNNTAAGCN	BC1_L5Aa_untag_1
NNNATTAGCN	BC1_L5Ab_untag_2
NNNGCGCAGCN	BC1_L5Ac_untag_3
NNNAGAGCN	BC1_L5Ba_Hfq_1
NNNGTGAGCN	BC1_L5Bb_Hfq_2
NNNCACTAGCN	BC1_L5Bc_Hfq_3

4. Quality and adapter trimming was performed using ‘bbduk’. Adapter and quality trimming. The following ‘bbduk’ command line code (Bushnell, 2014) was used to trim the sequences. 
bbduk.sh in= all_{barcode}.fastq \
    out= all_{barcode}.trimmed.pyCRAC.fq \
    ref=phix174_ill.ref.fa,\
linker.fa,\
phix_adapters.fa,\
illumina_adapter.fa,\
poly.fa \
    ktrim=r \
    mink=11 \
    hdist=1 hdist2=1 \
    k=23 mlf=0.13 minlen=20 ordered zl=4 \
    qtrim=rl trimq=28 \
    entropy=0.2 entropywindow=20 entropyk=5 \
    t={BBDUK_NCPUS}
The parameters used with the ‘bbduck’ the command line is briefly explained here. The ‘in’ parameter indicates the input file, which was the output from the previous ‘pyBarcodeFilter.py’ step, with the {barcode} representing the nucleotide sequence of the barcode in which all the sequences in the file was previously tagged with. The ‘out’ parameter represents the output file. The 3’-adapter sequences were right-trimmed using the option ‘ktrim=r’, where all sequences on the right-hand side of the sequences are removed from subsequent analyses. To specify the size of the kmers, the parameter ‘k=23’ specifies the kmer size to use for the majority of the length of the sequence and the ‘mink=11’ parameter specifies a shorter kmer size of 11 at the last 11 bases of the reads. The ‘ref’ parameter specifies a list of FASTA file, separated by commas, that contains the adapter sequences to trim. At most one mismatch was allowed per sequence, this is specified with the ‘hdist=1’ and ‘hdist2=1’ parameters. For length trimming, the ‘minlen=20’ parameter would discard reads less than 20 nucleotides long after trimming. For quality trimming, the ‘qtrim=rl’ parameter specified trimming on both the right and left side of the sequence and the ‘trimq=28’ parameter will perform quality trimming to Q28 using the Phred algorithm. The “ordered” parameter keeps the reads in the same order as the input file when multiple CPU cores were used. To remove sequences that are of low complexity, the parameters ‘entropy=0.2 entropywindow=20 entropyk=5’ will filter out reads that have an average entropy of under 0.1 with a sliding window size of 20 and a kmer size of 5. The parameter ‘zl=4’ indicates the use of level 4 compression algorithm to the output file.  A FASTQC (Andrews et al., 2012) analysis performed to enable visual inspection to ensure the trimmed output sequences are of high quality.

These 3’-adapter sequences include the following: adapter sequence (NAGATCGGAAGAGCACACGTCTG), Illumina universal adapter (AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCG ATCT), indexed adapter (AGATCGGAAGAGCACACGTCTGAACTCCAGTCACATC), PhiX adapters (AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTGAAA and AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATTAAAAAA), the PhiX genome (RefsSq: NC_001422.1), and any poly-A or poly-T sequences.

5. De-duplication of PCR amplified reads was performed using ‘pyCRAC’. Reads as a result of PCR duplicates were removed using the following pyFastqDuplicateRemover.py command line tool from the pyCRAC library (Webb et al., 2014).

      pyFastqDuplicateRemover.py -f all_{barcode}.trimmed.pyCRAC.fq -o all_{barcode}.compressed.pyCRAC.fasta

In the above command, the ‘-f’ parameter indicates the input file all_{barcode}.trimmed.pyCRAC.fq, and the ‘-o’ parameter indicates the output all_{barcode}.compressed.pyCRAC.fasta. The ‘{barcode}’ string represents the barcode associated with the file. 

6. Identification of the two RNAs that form the sequence hybrid representing the
RNA–RNA interaction was performed using ‘Hyb’. The ‘hyb detect’ function of the hyb pipeline (Travis et al., 2014) was used to detect the genomic co-ordinates of the two RNA sequences from the linked sequence hybrid of the RNA-RNA interaction. 

hyb detect db={HYB_DB_NAME} anti=1 format=fasta in={input.compressed_fasta}

Multiple reads for the same RNA-RNA interaction were counted and the sequences merged using the combine_hyb_merge_1.3.pl custom script from the Hyb-CRAC-R (https://github.com/IgnatiusPang/Hyb-CRAC-R ) GitHub directory (Hyb-CRAC-R/Source/Jai_Script/Annotation/combine_hyb_merge_1.3.pl) , with the parameter ‘TWO_WAY_MERGE=1’, which enable merging of hybrids in both orientation, and the input file name ‘input_file.hyb’, which is the output file from the previous ‘hyb detect’ step.

combine_hyb_merge_1.3.pl TWO_WAY_MERGE=1  input_file.hyb > merged_file.hyb

This resulted in a table of RNA-RNA interactions with a count for the number of reads supporting each RNA-RNA interaction. Only hybrid sequences in which the two RNA sequences were unambiguously mapped to two corresponding genomic locations were used for further analysis. Hybrid sequences in which two RNA sequences were less than 100 or equal to nucleotides apart from each other were excluded from further analyses.

7. Each RNA transcript from the RNA-RNA interaction were annotated by comparing the genomic co-ordinates with the transcriptome annotations. Annotations were assigned based on priority, ranked from highest to lowest priority. These priorities are ncRNA, CDS, 5' UTR, 3' UTR, rRNA, tRNA, CRISPR, pseudogene, riboswitch, Ribonuclease P RNA, RNA thermometer, Signal recognition particle RNA, Transfer-messenger RNA, and operon. The annotation of RNA-RNA interaction was performed using the following R script from the Hyb-CRAC-R GitHub repository (Source/Demultiplex/scripts/annotation.R). The user will need to edit the following parameters in the R script. The ‘hyb_file’ parameter is the input file ‘merged_file.hyb’ generated in the previous step from running the ‘combine_hyb_merge_1.3.pl’ script, the ‘rna_type_file’ represents the RNA types and their priority of being reported, and the ‘output_file’ indicates the output file. The ‘rna_type_file’ has the following format, in which the first column is the ‘rank’ which states which priority the ’rna_type’ in the second column will be reported.
rank	rna_type
1	ncRNA
2	CDS
3	5UTR
4	3UTR
5	rRNA
6	tRNA
7	CRISPR
8	pseudogene
9	riboswitch
10	Rnase_P_RNA
11	RNA_termometer
12	SRP_RNA
13	tmRNA
14	operon

8. Obtaining read counts of RNA transcripts. The read counts corresponding to each RNA species from the hybrid sequences was identified using ‘pyCRAC’. The read counts of individual RNA transcripts corresponding to those identified in the hybrid sequences was obtained from the merged and trimmed paired-end reads data using the ‘CRAC_pipeline_SE.py’ script from the ‘pyCRAC’ pipeline. This script also runs its own PCR reads de-duplication step within the script itself so the input file does not have to be de-duplicated for PCR repeats. For each of the two RNA sequences from a hybrid sequence, the read count corresponds to the highest read depth within the region spanned by the RNA sequence. 
The sequence alignment within the ‘pyCRAC’ pipeline was performed using ‘NovoAlign’ (http://www.novocraft.com/products/novoalign/). The command to run the genome index function of novoindex software includes the output index file ‘Sa_JKD6009.nix’ and the input file ‘Sa_JKD6009.fasta’.
novoindex Sa_JKD6009.nix Sa_JKD6009.fasta
The ‘CRAC_pipeline_SE.py’ command  uses the following options, ‘-f’ the input file which was produced by the ‘bbduck’ sequence trimming step, the ‘-g’ parameter indicates the path location of the GTF genome annotation file (from the GitHub directory Data/Genomic/JKD6009/Annotation/Hyb_pyCRAC/). The ‘-c’ parameter indicates a tab-separated file with two columns, the name of the genome as used in the GTF file ‘Saa_6009_C1’ and the size of the genome 2930441. The ‘—novoindex’ parameter represents the indexed genome file ‘Sa_JKD6009.nix’ generated above. The ‘-a’ parameter specifies the adaptor nucleotide sequence. The ‘-name’ parameter provides a unique name for the analysis, with the {barcode} tag representing one of the multiplexing barcode from this experiment. The ‘-p’ parameter indicates the number of CPU processes that would be used to perform the read mapping. The script outputs two ‘.sgr’ file, one for the plus strand and one for the negative strand. 

CRAC_pipeline_SE.py \
      -f all_{barcode}.trimmed.pyCRAC.fastq \
      -g merge_features_simple_annot_checked_edited.gtf \
      -c Sa_JKD6009.genome \
      --novoindex= Sa_JKD6009.nix \
      --name="all_{barcode}" \
      -a AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT \
      -p 1 > {output.pyCRAC_log}

The output from the ‘pyCRAC’ pipeline is in SGR format which was converted into PWS format using a custom python script (sgr2pws.py) from the GitHub repository (Source/Jai_Script/Hyb_stats/sgr2pws.py). The ‘sgr2pws.py’ script takes as input the .sgr files for the plus strand using the ‘-p’ parameter and negative strand using the ‘-n’ parameter. The output file is the .pws file.  
sgr2pws.py \
       -p  input.plus_strand_file.sgr \
       -n  input.minus_strand_file.sgr > \
       output.pws

9. Statistical analysis of the significance of the RNA–RNA interactions was calculated using custom R scripts as described in Waters et al. (2017). Probabilistic analysis of RNA-RNA interactions. In order to identify RNA-RNA interactions that are of biological origin, we compared each observed interaction with the background probability that the same RNA-RNA interaction event could occur randomly in vitro. The analysis was based on the method described in Sharma et al. (2016) and Waters et al. (2017). Briefly, the probability of drawing a gene at random, P(gx), was estimated by the maximum read depth for the gene divided by the total number of reads in the replicate (N). The background probability for observing a RNA-RNA interaction, pdf (gx, gy), was estimated by multiplying the probabilities P(gx) and P(gy). Ligation events that were not observed were set to a probability of zero, as per Sharma et al. (2016). The probability distribution pdf(gx, gy) was re-normalized to the sum of 1. The background probability of observing k occurrences of interaction between gx and gy was modeled using the binomial distribution, k ~ Binomial(p=pdf(gx, gy), N). These calculations were used to assign a p-value (ie: P(X < k)) to each experimentally observed pair of RNA-RNA interaction. To perform the calculations, the following four R scripts from the GitHub respository would be ran sequentially in these order:
Source/Demultiplex/scripts/format_stats_dataset.R
Source/Demultiplex/scripts/max_reads.R
Source/Demultiplex/scripts/calculate_generic_hybrid_p_values.R
Demultiplex/scripts/separate_unique_multimapped_hybrids.R

First, we run ‘format_stats_dataset.R’, changing the parameters inside the R script. The input parameter ‘num_mapped_reads_files` is a file generated in the previous step ‘CRAC_pipeline_SE.py’. The input parameter ‘hyb files‘ is a list of all the ‘.hyb’ files associated with each barcode analyzed after then have been analyzed for known RNA targets by the ‘annotation.R’ script. The pws files were the files generated by the ‘sgr2pws.py’ script. The ‘dataset_stats_file’ is the name of the output file ‘datasets.txt’, which looks like below, with the first column denoting the file name and the second column represents the number of reads in the file. 

replicate_1.hyb replicate_1.pws 1234567 
replicate_2.hyb replicate_2.pws 1443584 
replicate_3.hyb replicate_3.pws 135463

The max_reads.R script identifies the maximum read depth for all RNA-RNA in one file. The following parameters were sued with ‘max_reads.R’. This includes the name of the output file in the “MAX_READS” parameter, the hyb output file from after the annotation of the RNA step, and the .pws files from running the ‘CRAC_pipeline_SE.py’ script. 

To run the ‘calculate_generic_hybrid_p_values.R’ script, the script takes in command line output files from running ‘max_reads.R’ lines and the ‘dataset.txt’ file generated by ‘sgr2pws.py’. 

The separate_unique_multimapped_hybrids.R script identifies uniquely mapped RNA-RNA hybrids or multiply mapped RNA-RNA hybrids. It takes into input the output from the ‘calculate_generic_hybrid_p_values.R’ for the analyses.

10. Merging of probabilities from multiple biological replicate experiments. For RNA-RNA interactions that were observed in multiple biological replicate experiments, the corrected pvalue was combined using the Fisher’s method. If the RNA-RNA interaction was only identified in one biological replicate, the original p-value was used and Fisher's method was not applied. After processing with Fisher’s method, the p-values were corrected for false discovery rates using the Benjamini-Hochberg multiple testing correction. Statistically significant RNA-RNA interactions were defined as those with a false discovery rate threshold of p-value < 0.05 (Tree et al., 2018). This step is performed using the ‘merge_hybrid_p_values.R’ script from the following GitHub directory (https://github.com/IgnatiusPang/Hyb-CRAC-R/tree/main/Source/Jai_Script/Hyb_stats). The input parameters were adjusted inside the file itself. It takes a list of files generated by running ‘calculate_generic_hybrid_p_values.R’ and also a list of names representing each sample to be merged. The output is a file which contains the merged p-values from combining data from a number of biological replicates analysis.
