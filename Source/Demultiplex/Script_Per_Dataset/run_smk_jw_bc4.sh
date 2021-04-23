#!/bin/bash

#PBS -N JW_BC4
#PBS -l nodes=1:ppn=1
#PBS -l mem=60gb
#PBS -l walltime=11:55:00
#PBS -j oe
#PBS -M i.pang@unsw.edu.au
#PBS -m ae

cd ${PBS_O_WORKDIR}


## Start python
module load python/3.8.3

# snakemake --dry-run --printshellcmds -j 12 --snakefile jw_bc4.smk

snakemake -j 12 --printshellcmds --snakefile jw_bc4.smk --cluster-config  cluster.json \
 --latency-wait 20 --use-envmodules \
 --cluster "qsub -N {cluster.N} -l nodes={cluster.nodes}:ppn={cluster.ppn},walltime={cluster.walltime},mem={cluster.mem} -M {cluster.email} -m ae -j oe"  

# 
# snakemake \
# --snakefile jw_bc4.smk \
#  --dag \
#  | dot -Tpdf \
#  > dag.pdf
# 
# 
# 
# snakemake \
#     --snakefile jw_bc4.smk \
#     --rulegraph \
#     | dot -Tpdf \
#     > rulegraph.pdf
#     


# qsub -I -N Smk_Staging -l nodes=1:ppn=12,walltime=12:00:00,mem=60gb
# qsub -I -N Smk_Staging -l nodes=1:ppn=1,walltime=12:00:00,mem=6gb

#  cd /srv/scratch/treelab/igy/Data/CLASH_JKD6009/DM_MiSeq/Processed_Data/
# cd /srv/scratch/treelab/igy/Projects/StaphCLASH2020/Source/Demultiplex/

