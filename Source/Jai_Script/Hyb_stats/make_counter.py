#!/usr/bin/python

import sys

counter = int(sys.argv[1])
pws_len = int(sys.argv[2])

for i in range(1, (pws_len/counter)+2 ):
	print i*counter

