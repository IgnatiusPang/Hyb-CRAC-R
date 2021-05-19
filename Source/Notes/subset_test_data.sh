
cd /Users/ipang/PostDoc/2020/Hyb-CRAC-R/Data/CLASH_FASTQ

/Users/ipang/Programs/2021/seqtk/./seqtk sample -s100 DGM-RNase-CLASH_S1_L001_R1_001.fastq.gz 500000 > sub_DM_R1.fastq
/Users/ipang/Programs/2021/seqtk/./seqtk sample -s100 DGM-RNase-CLASH_S1_L001_R2_001.fastq.gz 500000 > sub_DM_R2.fastq

gzip sub_DM_R1.fastq
gzip sub_DM_R2.fastq


	
