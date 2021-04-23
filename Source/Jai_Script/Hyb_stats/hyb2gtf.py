#!/usr/bin/python
# hyb2gtf.py. Script to convert hyb file to gtf format
#
# copyright 2017 Jai J. Tree 
# GNU Lesser General Public License

import os
import sys
os.chdir(os.getcwd())

hyb = open(sys.argv[1])

for line in hyb:
	x = line.strip().split("\t")
	name = x[3].split("|")
	chromo = name[0]
	gene_name = name[1]
	RNA_type = name[3]
	name2 = x[9].split("|")
	chromo2 = name2[0]
	gene_name2 = name2[1]
	RNA_type2 = name2[3]
	cnt = x[15].split(";")			#have added hybrid counts to the GTF - column 6 (before strand)
	count = cnt[0].split("=")		#
    	if x[6] < x[7]:
        	strand_1 = "+"
    	else:
        	strand_1 = "-"
    	if x[12] < x[13]:
        	strand_2 = "+"
    	else:
        	strand_2 = "-"
    	out = chromo + '\t' + RNA_type + '\t' + "hybrid" + '\t' + str(min(x[6:8])) + '\t'  + str(max(x[6:8])) + '\t' + count[1] + '\t' +  strand_1 + '\t' + x[0] + '\t' +  'gene_id "' + gene_name + '"; other_bit "' + gene_name2 + '";' + '\n' + chromo2 + '\t' + RNA_type2 + '\t' + "hybrid" + '\t' + str(min(x[12:14])) + '\t' + str(max(x[12:14])) + '\t' + count[1] + '\t' +  strand_2 + '\t' + x[0] + '\t' +  'gene_id "' + gene_name2 + '"; other_bit "' + gene_name + '";'
    	print out
hyb.close()

