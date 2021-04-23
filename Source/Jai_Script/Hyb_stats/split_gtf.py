#!/usr/bin/python
#file splitter to reduce computation time 

import sys
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-s", "--split", dest="split", help="genome integer to divide file")
parser.add_option("-f", "--file", dest="file", help="input file")
parser.add_option("-t", "--type", dest="type", help="pws or gtf file", default="gtf")
parser.add_option("-l", "--len", dest="pws_len", help="the length of the pws file", default="3000000")
(options, args) = parser.parse_args()

split = int(options.split)
file = options.file
type = options.type
pws_len = int(options.pws_len)

if options.type == "gtf":
	field = int(3)
if options.type == "pws":
	field = int(1)
	
GTF = open(file)

split_int = split

for i in range(1,(pws_len/split_int)+1):
	open(file + str(split_int*i) + "." + type , "w")
#	print "opened file: " + str(file + str(split_int*i) + "." + type)

OUT1 = open(file + str(split_int) + "." + type , "w")

for a in GTF:
	b = a.strip().split("\t")
	if int(b[field]) < split_int:
		OUT1.write(a)
	else:
		OUT1.close()
		split_int += split
		OUT1 = open(file + str(split_int) + "." + type , "w")
		OUT1.write(a)

OUT1.close()
