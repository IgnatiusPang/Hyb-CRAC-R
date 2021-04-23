# README #

## This repository contains the complete data processing pipeline for analyzing paired-end CRAC data (such as kinetic-CRAC data).##

## The CRAC_pipeline_PE.py  script runs the following steps over multiple processors: ##

1) Demultiplexing of reads using 5' adapter barcode sequences (pyBarcodeFilter.py)

2) Adapter removal and trimming of reads (Flexbar and TrimNucs.py)

3) collapsing the data (pyFastqDuplicateRemover.py)

4) aligning the data to the reference genome (yeast in our case) (novoalign)

5) counting the number of reads that map to each annotated feature (pyReadCounters.py)

6) made bedgraph files for visualization in genome browsers (pyGTF2Bedgraph.py)

7) making sgr files for each file (pyGTF2sgr.py)


## Important! ##
Use the -h flag to get a detailed help menu. 

The CRAC_pipeline_PE.py script needs to be run from the folder that contains the fastq files

The barcode list file should contain two tab-separated columns in which the first column is the barcode sequence and the second column is the name of the experiment 

The chromosome_lengths file should contain two tab-separated columns in which the first column has the chromosome name and the second the chromosome length.

## To be able to run this pipeline you need to install the following dependencies: ##

[Python 2.7](https://www.python.org/download/releases/2.7/)

ruffus python package (installable using easy_install or pip)

[novoalign](http://www.novocraft.com)


[pyCRAC 1.2.2.7](https://bitbucket.org/sgrann/pycrac)


[Flexbar](https://github.com/seqan/flexbar)

## Additionally: ##

I also included an R script (written by Alina Selega) that performs the Gaussian process analysis that we used for selecting transcripts with differential cross-linking in the time series

Also included is the TrimNucs.py script that the pipeline uses for trimming the last random nucleotide from the 5' end of the transcript after demultiplexing.

This needs to be added to your shell path.

If you have any problems with the script (it is still in testing phase), please do not hesitate to contact me:
sgrannem@staffmail.ed.ac.uk