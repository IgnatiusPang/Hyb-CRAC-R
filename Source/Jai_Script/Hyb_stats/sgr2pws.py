#!/usr/bin/python
#sgr2pws.py
#converts sgr into pws format
#
# copyright 2017 Jai J. Tree 
# GNU Lesser General Public License

from optparse import OptionParser

parser = OptionParser()
parser.add_option("-p", "--positive_strand", dest="plus", help="Postive strand SGR file")
parser.add_option("-n", "--minus_strand", dest="minus", help="Minus strand SGR file")
(options, args) = parser.parse_args()

plus = open(options.plus)
minus = open(options.minus)

def parse_sgr(foo, strand):
	"""Parses file into PWS"""
	for i in foo:
		l = i.split("\t")
		print l[0] + "\t" + l[1] + "\t.\t" + strand + "\t" + l[2].strip() + "\t0\t0\t0\t0"

parse_sgr(plus, "+")
parse_sgr(minus, "-")

plus.close()
minus.close()

