#!/bin/bash
# hyb_p_values.sh
# shell to iteratively process CLASH hyb datasets to calculate cumulative p-values
#
# copyright 2017 Jai J. Tree 
# GNU Lesser General Public License
#
# input files:
# HYB file for each dataset
# PWS file for each dataset (matched CRAC data)
#
# workflow:
# 	1. max_reads.sh - report CRAC data at hyb positions
# 	2. calculate_generic_hybrid_p_values.R - report p-value and append to table
# iterate over pairs of p-value table to generate cumulative counts of p
# 	3. merge_hybrid_p_values.R 
# 	4. reformat table merged table for another iteration
#
# generate PWS files using SGR outputs from CRAC_pipeline and sgr2pws.py - not part of the pipeline
#
# usage: hyb_p_values.sh INPUT_DATASETS.txt
#
# input is a space delimited file of HYB PWS SIZE
# where size is the total number of mapped reads for a dataset

#############################################

module rm R
module add R/3.3.1 # the R dplyr library requires R 3.2.1 or higher

# software locations
merge_hybrids="/srv/scratch/z3371724/StaphCLASH2020/Source/Jai_Script/Hyb_stats/merge_hybrid_p_values.R"						# update with local location 
calc_generic_p_value="/srv/scratch/z3371724/StaphCLASH2020/Source/Jai_Script/Hyb_stats/calculate_generic_hybrid_p_values.R"	# update with local location

# generate p_values and connection scores for each hybrid.

while read -r hyb pws size remainder; 
	do
#	echo skipping...
	echo running max_reads.sh on $hyb $pws $size 
	max_reads.sh $hyb $pws 
	Rscript $calc_generic_p_value $hyb\.max_reads.txt $size 
	done < $1

# output files are in the format $hyb\.max_reads.txt.p_values.txt
# pairwise merge of hyb files

echo finsihed making p-val files
echo making cumulative table

while read -r hyb pws remainder;
	do
	if [ $last_file  ] 
	then
		echo entering while loop with $last_file and $hyb													
		name=$(basename $hyb _scRNA.hyb) 
		Rscript $merge_hybrids $last_file $hyb\.max_reads.txt.p_values.txt							# use sqldf to merge data table together | union chiastic table 
		echo outfile name: $hyb\.max_reads.txt.p_values.txt_merged_output_table.txt
		head -n 1 $hyb\.max_reads.txt.p_values.txt_merged_output_table.txt > file1												# strip off the header
		sed 's/"//g' $hyb\.max_reads.txt.p_values.txt_merged_output_table.txt | tail -n +2 | cut -f 2- | sort | uniq >> file1		# remove first field and quotes | sort | remove duplicates 
		format_merged_hybrids_from_sqldf.py file1 > cumulative_hybrid_p_values.txt												# reformat merged table to add together p-values | return single table
		last_file=cumulative_hybrid_p_values.txt																				# use new merged and formatted table to merge with a new dataset
	else
		last_file=$hyb\.max_reads.txt.p_values.txt
	fi
	done < $1
	
