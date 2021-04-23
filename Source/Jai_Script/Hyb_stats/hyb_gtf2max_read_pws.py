#!/usr/bin/python
#
# copyright 2017 Jai J. Tree 
# GNU Lesser General Public License
#
# script to retrieve the maximal read count at a hybrid half
# Script returns hybrid ID, #hybrids, maximal read count for both halves
# Inputs need to be sorted by chromosome 
# remove header from PWS then: sort -k1,1 -k4,4 -k2n,2
# Divide HYB file into each feature using hyb2gtf.py (make sure counts are included)
# sort HYG.GTF sort -k1,1 -k7,7 -k4n,4 

import sys
import os

HYB_GTF = open(sys.argv[1])
PWS = open(sys.argv[2])
OUT = open(sys.argv[1] + "_correlate.txt", "w")

OUT.write( "#" + os.getcwd() + "\n#" + sys.argv[0] + "\n#HYB_GTF=" + sys.argv[1] + "\n#PWS=" + sys.argv[2] + "\n" )

def get_max(start, end, chrom, strand): 
	MAX = float(0)
	for a in PWS:
		if a[0] != "#":
			b = a.split("\t")
			if b[0] == chrom:
				if b[3] == strand:
					if int(b[1]) >= int(start):
						if int(b[1]) <= int(end):
							if float(b[4]) > MAX:
								MAX = float(b[4])
					if int(b[1]) > int(end):
						break
	PWS.seek(0)
	return MAX

OUT.write( "#ID\tmax_reads\thyb_count\trna_class\tchromo\tname\tstart\tend\n" )

last_end = int(0)

for i in HYB_GTF: 
	j = i.split("\t")
	count = j[5]
#	merged = cnt[2].split("=")
	
	nme = j[8].split('"')
	name = nme[1]
	chromo = j[0]
	rna_class_1 = j[1]
	strand = j[6]
	
	ID = j[7]
	start_1 = j[3]
	end_1 = j[4]
	
	hyb1_max = get_max(min(start_1, end_1), max(start_1, end_1), chromo, strand )
	
	out = ID + "\t" + str(hyb1_max) + "\t" + str(count) + "\t" + rna_class_1 + "\t" + chromo + "\t" + name + "\t" + start_1 + "\t" + end_1 + "\t" + strand
	
#	print out
	OUT.write( out + "\n" )
	
