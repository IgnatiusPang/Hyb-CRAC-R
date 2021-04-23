#!/usr/bin/python
# not compatible with python 3
__author__ = "Sander Granneman"
__copyright__ = "Copyright 2014"

import sys
import time
from optparse import *

def headersplitter(string,filetype="fastq"):
	""" splits the header on the @ character"""
	if filetype is "fastq":
		Fld = string.split("@")[1:]
		assert len(Fld) == 2, "Can not split the %s header correctly! %s @ characters in the header\n" % (string, len(Fld))
		header_f = "@%s" % Fld[0]
		header_r = "@%s" % Fld[1]
		return header_f,header_r
	elif filetype is "fasta":
		Fld = string.split(">")[1:]
		if len(Fld) is 2:
			header_f = ">%s" % Fld[0]
			header_r = ">%s" % Fld[1]
		elif len(Fld) is 1:
			header_f = ">%s/1" % Fld[0]
			header_r = ">%s/2" % Fld[0]
		return header_f,header_r
	else:
		raise TypeError("I do not recognize the file type %s. Please use 'fasta' or 'fastq'\n")

def trimNucs(data,numberofnucs=1,end="left",file_type="fastq",addtoheader=False,pairedend=False):
	""" trims remaining ramdon nucleotides from the 5' ends of the reads. 
	These need to be removed prior to removing duplicate reads """
	if file_type == "fasta":
		for line in data:
			if line[0] == ">":
				sys.stdout.write("%s" % line)
			else:
				if end == "left":
					sys.stdout.write("%s" % line[numberofnucs:])	   # remove the nucleotides from the left
				else:
					sys.stdout.write("%s\n" % line[:-numberofnucs])
	elif file_type == "fastq":
		while True:
			try:
				name	 = data.next().strip()
				sequence = data.next().strip()
				plus	 = data.next().strip()
				quality	 = data.next().strip()
			except StopIteration:
				break
			if end == "left":
				trimmedseq = sequence[:numberofnucs]
				sequence = sequence[numberofnucs:]
				quality	 = quality[numberofnucs:]
			else:
				trimmedseq = sequence[-numberofnucs:]
				sequence = sequence[:-numberofnucs]
				quality = quality[:-numberofnucs]
			if addtoheader:
				pairedend = False
				Fld = name.split("@")[1:]
				if len(Fld) == 2:
					pairedend = True
					header_f = "@%s" % Fld[0]
					header_r = "@%s" % Fld[1]
					sys.stdout.write("%s%s%s\n%s\n%s\n%s\n" % (header_f,trimmedseq,header_r,sequence,plus,quality))
				else:
					sys.stdout.write("%s%s\n%s\n%s\n%s\n" % (name,trimmedseq,sequence,plus,quality))
			else:
				sys.stdout.write("%s\n%s\n%s\n%s\n" % (name,sequence,plus,quality))
									
def main():
	parser = OptionParser(usage="usage: %prog [options] -f filename -b barcode_list", version="%prog 0.0.1")
	parser.add_option("-f", "--input_file", dest="filename",help="name of the FASTQ or FASTA input file. Default is standard input", metavar="FILE")
	parser.add_option("--file_type", dest="filetype",help="type of file, uncompressed (fasta or fastq). Default is fastq", metavar="FASTQ",default="fastq")
	parser.add_option("-n", dest="nucs", type="int",help="indicate the number of nucleotides you want to have trimmed from the 5' end. Default = 1 ", default=1)
	parser.add_option("--right", dest="right",action="store_true",help="to trim from the right end of the sequence, default is left.",default=False)
	parser.add_option("--addtoheader",dest="addtoheader",action="store_true",help="to add the sequence of the trimmed nucleotides to the header. Default is False",default=False)
	(options, args) = parser.parse_args()
	data = sys.stdin
	trimfrom = "left"
	if options.filename:
		data = open(options.filename,"r")
	if options.right:
		trimfrom = "right"
	trimNucs(data,numberofnucs=options.nucs,end=trimfrom,file_type=options.filetype,addtoheader=options.addtoheader)
	
if __name__ == "__main__":
	main()
