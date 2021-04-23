#!/usr/bin/python
# annotate_hyb_2.py
#
# copyright 2017 Jai J. Tree 
# GNU Lesser General Public License
#
# Have added some code that allows the script to detect whether thereeis an existing annotation. 
# Will look for annototation only if it was unknown or hasn't been annotated, 
# this should allow rerunning annotation with new feature GTFs without earasing old annotations.
#
# Use -o True to write each interaction class to a seperate file (eg: sRNA-mRNA) and generate the 
# file containing a count of each interactions type (.counter). 

import fileinput
import sys
import os
from GTF_stripper import *
from optparse import OptionParser

os.chdir(os.getcwd())

parser = OptionParser()
parser.add_option("-g", "--gtf", dest="annot_gtf", help="gtf to use for hyb annotation")
parser.add_option("-o", "--all_outputs", dest="outputs", help="writes individual interaction class files and a coutner", default = False)
(options, args) = parser.parse_args()

annot_gtf = options.annot_gtf
outputs = options.outputs

def feature_overlap(starts, ends, chromo):
	"""Method to search the sakai dictionary of features to find an overlap. NB: stops at the first overlapping feature - those that map across boundaries will only return the first found NBB: dictionaries are searched in random order. Returns a list containing feature annoation data: chromo, gene_name, gene_name, source"""
	if starts < ends:
		strand = "+"
	else:
		strand = "-"
	if isinstance(chromo, list):
		chromo = chromo[0]
	start = min(int(starts), int(ends))				 		
	end = max(int(starts), int(ends))
	out = []
		
	for i in sakai:
		if chromo == sakai[i][5]:
			if start >= sakai[i][0] and start <= sakai[i][1]:
				if strand == sakai[i][2]:
					out = [sakai[i][5], sakai[i][3], sakai[i][3], sakai[i][4]]
					break
			elif end >= sakai[i][0]:
				if end <= sakai[i][1]:
					if sakai[i][2] == strand:
						out = [sakai[i][5], sakai[i][3], sakai[i][3], sakai[i][4]]
						break
	if len(out) < 1:
		out = [chromo, "none", "none", "unknown"]
	return out

def check_annot(chromosome):
	"""Checks to see if hyb file already annotated, returns chromosomes. Hand the tab split line - list of annotations or tuple of chromosomes"""
	if "|" in chromosome[3]:
		c1 = chromosome[3].split("|")
		c2 = chromosome[9].split("|")
		return [c1, c2]
	else:
		return chromosome[3], chromosome[9]

#make a dictionary of features from the Sakai GTF

sakai = dict()								
gtf_foo = open(annot_gtf)
for i in gtf_foo:
	if not i.startswith("#"):
		gtf = GTF_stripper(i)
		if gtf.GTF["gene_name"] in sakai:					
			gene_name = gtf.GTF["gene_name"]
			while gene_name in sakai:
				iteration = gene_name.split("_")
				if len(iteration) == 1:
					gene_name = gene_name + "_1"
				else:
					gene_name = iteration[0] + "_" + str(int(iteration[1]) + 1)
		else:
			gene_name = gtf.GTF["gene_name"]

		sakai[gene_name] = [int(gtf.GTF["start"]), int(gtf.GTF["end"]), gtf.GTF["strand"], gene_name.replace("_", "-"), gtf.GTF["source"], gtf.GTF["chromosome"]] #.replace("_", "-")]

gtf_foo.close()
	
count = int(1)
converter = dict()

def make_key(RNA, converter):
	if RNA not in converter.keys():
		global count						
		count = count*10
		converter[RNA] = count
		return count
	else:
		return converter[RNA]

counter = dict()
decrypter = dict()
type_writer = dict()

for line in open(args[0]):
	l = line.split("\t")
	rna1, rna2 = None, None
	a = check_annot(l)							#looks for existing annotation and checks whether it's unknown - runs feature_overlap if none or no annotation
	if not isinstance(a[0], list) or a[0][3] == "unknown":		
		rna1 = feature_overlap( int(l[6]), int(l[7]), a[0] )
	else:
		rna1 = a[0]
	if not isinstance(a[1], list) or a[1][3] == "unknown":
		rna2 = feature_overlap( int(l[12]), int(l[13]), a[1] )
	else:
		rna2 = a[1]		

	l[3] = "|".join(rna1)
	l[9] = "|".join(rna2)
	out = "\t".join(l).strip()
	print out
	key = make_key(rna1[3], converter) + make_key(rna2[3], converter)				#make a key that is the same for AB and BA interactions
	if key not in decrypter:						
		decrypter[key] = rna1[3] + "-" + rna2[3]				#make a human readable key
	
	if outputs:
		if key in counter:									#count number of interaction types
			counter[key] += 1		
			type_writer[key] += out + "\n"					#write hyb line to a dictionary using type as a key
		else:
			counter[key] = 1
			type_writer[key] = out + "\n"

if outputs:
	fo = open(args[0] + ".counter", "w")						#write the numbers of interactions to a .counter file
	fo.write("interacting_RNAs" +"\t" + "number_found" + "\n")
	for k, v in sorted(counter.items()):
		fo.write(decrypter[k] + "\t" + str(v) + "\n")
	fo.close()

	for l, m in type_writer.items():
		type_file = open(args[0] + "_" + decrypter[l] + ".hyb", "w")	#write each interaction type to a seperate file
		type_file.write(m)
		type_file.close()


