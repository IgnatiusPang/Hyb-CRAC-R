#!/usr/bin/bash
# staph progressive annotation shell
# annotates hybrids with CDS > scRNA > rRNA > tRNA > misc_features > 5UTR > 3 UTR > intergenic regions
# 5UTR and 3UTR are set ranges

HYB=$1
name=${HYB%%.hyb}  #$(basename $1 .hyb)
curr_path=$(dirname $1)
loc=/srv/scratch/treelab/scripts/jaitree/gtf   

SCRIPT_LOC=/srv/scratch/z3371724/StaphCLASH2020/Source/Jai_Script/Annotation

$SCRIPT_LOC/annotate_hyb_2.py $1 -g $loc/Sa_JKD6008_1.5_protein_coding.gtf > $curr_path/file1
$SCRIPT_LOC/annotate_hyb_2.py $curr_path/file1 -g $loc/Sa_JKD6008_1.5_scRNA.gtf > $curr_path/file2
$SCRIPT_LOC/annotate_hyb_2.py $curr_path/file2 -g $loc/Sa_JKD6008_1.5_rRNA.gtf > $curr_path/file3
$SCRIPT_LOC/annotate_hyb_2.py $curr_path/file3 -g $loc/Sa_JKD6008_1.5_tRNA.gtf > $curr_path/file4
$SCRIPT_LOC/annotate_hyb_2.py $curr_path/file4 -g $loc/Sa_JKD6008_1.5_misc_feature.gtf > $curr_path/file5
$SCRIPT_LOC/annotate_hyb_2.py $curr_path/file5 -g $loc/Sa_JKD6008_1.5_5UTR.gtf > $curr_path/file6
$SCRIPT_LOC/annotate_hyb_2.py $curr_path/file6 -g $loc/Sa_JKD6008_1.5_3UTR.gtf > $curr_path/file7
$SCRIPT_LOC/annotate_hyb_2.py $curr_path/file7 -g $loc/Sa_JKD6008_1.5_intergenic_region.gtf > $name\_annot.hyb
$SCRIPT_LOC/annotate_hyb_2.py $name\_annot.hyb -g $loc/Sa_JKD6008_1.5_intergenic_region.gtf -o False

rm $curr_path/file*
echo annotation finished for $name 

