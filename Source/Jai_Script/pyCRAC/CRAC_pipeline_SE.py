#!/usr/bin/env python
__author__		= "Sander Granneman"
__copyright__	= "Copyright 2016"
__version__		= "0.1.6"
__credits__		= ["Sander Granneman"]
__email__		= "sgrannem@staffmail.ed.ac.uk"
__status__		= "beta"

from ruffus import *
from ruffus.cmdline import MESSAGE
from collections import defaultdict
import ruffus.cmdline as cmdline
import subprocess
import glob
import re
import os
import argparse
import sys
import time

parser = cmdline.get_argparse(description="CRAC pipeline for processing single-end CRAC data")
parser.add_argument("-f",dest="forwardreads",help="the path to your fastq file with reads.",metavar="data_1.fastq",default=None)
parser.add_argument("-g","--gtf",dest="gtf",help="the path to your gtf annotation file",metavar="rRNA.gtf",default=None)
parser.add_argument("-c","--chromosome",dest="chromosome",help="the path to your chromosome length file",metavar="chromosome_lengths.txt",default=None)
parser.add_argument("-i","--ignoremutations",dest="ignoremuts",action="store_true",help="to tell the novoalign parser to ignore mutations. Useful when analyzing ChemModSeq data. Default is OFF",default=False)
parser.add_argument("--novoindex",dest="novoindex",help="the path to your novoindex file",metavar="yeast.novoindex",default=None)
parser.add_argument("--name",dest="name",help="provide a single word describing the run. Default is 'analysis' with a time stamp",default="analysis_%s" % time.strftime("%d%m%Y"))
parser.add_argument("-b","--barcodes",dest="barcodes",help="the path to the file containing the list of barcodes",metavar="barcodes.txt",default=None)
parser.add_argument("-a","--adapter",dest="adapter",help="provide your 3' adapter sequence. Default is TGGAATTCTCGGGTGCCAAGGC",default="TGGAATTCTCGGGTGCCAAGGC")
parser.add_argument("-m","--mismatches",dest="mismatches",type=int,help="indicate how many mismatches you allow for demultiplexing. Default is 1",default=1)
parser.add_argument("-p","--processors",dest="processors",type=int,help="indicate how many processors you want to use for analyses. Default is 8",default=8)
args = parser.parse_args()

def getFileBaseName(filename):
	""" removes path and extension from file name """
	return os.path.splitext(os.path.basename(filename))[0]

def getFullFilePath(filename):
	""" grabs the full path name for the file """
	return os.path.abspath(filename)

def getBarcodeInfo(barcodefile):
	return ["%s.fastq" % "_".join(line.strip().split()) for line in open(barcodefile,"r").readlines()]
	
def runFlexBar(inputfile,outputfile):
	""" runs Flexbar on the data to remove the adapter sequence from the forward reads """
	outputfilename = "%s/%s_trimmed" % (os.path.join(root_dir,"flexbar_trimmed"),re.search("(.*).(san)?fastq",inputfile).group(1))
	cmd = "flexbar -r '%s' -qf i1.8 -t '%s' -n 10 -ao 7 -as %s -qt 30" % (inputfile,outputfilename,args.adapter)
	logger.info(cmd)
	os.system(cmd)
	#os.system("touch %s" % outputfile)

def demultiplexSamples(inputfile,outputfile):
	""" demultiplexes all the samples """
	os.chdir(os.path.join(root_dir,"demultiplexed"))
	cmd = "pyBarcodeFilter.py -f '%s' -b '%s' -m '%s'" % (inputfile,os.path.join(home_dir,args.barcodes),args.mismatches)
	logger.info(cmd)
	os.system(cmd)
	#os.system("touch %s" % outputfile)
	os.chdir(root_dir)

def trimAndCollapse(inputfile,outputfile):
	""" Removes the last random nucleotide from the forward read, collapses the data and then splits it again into fasta files """
	cmd = "/srv/scratch/z3371724/StaphCLASH2020/Source/Jai_Script/pyCRAC/TrimNucs.py -n 1 --addtoheader -f '%s' | pyFastqDuplicateRemover.py -o '%s'" % (inputfile,outputfile)
	logger.info(cmd)
	os.system(cmd)
	#os.system("touch %s" % outputfile)

	
def alignReads(inputfile,outputfile):
	""" Runs novoalign on all the collapsed files"""
	cmd = "novoalign -d '%s' -f '%s' -r Random > '%s'" % (args.novoindex,inputfile,outputfile)
	logger.info(cmd)
	os.system(cmd)
	#os.system("touch %s" % outputfile)

def runPyReadCounters(inputfile,outputfile):
	""" runs pyReadCounters on all the files """
	outputfile = re.search("(.*)_count_output_reads.gtf",outputfile).group(1)
	if args.ignoremuts:
		string = "--mutations=nomuts"
	else:
		string = ""
	cmd = "pyReadCounters.py -f '%s' --gtf '%s' -v --rpkm -o '%s' '%s'" % (inputfile,args.gtf,outputfile,string)
	logger.info(cmd)
	os.system(cmd)
	#os.system("touch %s_count_output_reads.gtf" % outputfile)
	
def runSecondPyReadCounters(inputfile,outputfile):
	""" runs pyReadCounters on all the files """
	outputfile = re.search("(.*)_count_output_cDNAs.gtf",outputfile).group(1)
	cmd = "pyReadCounters.py -f '%s' --gtf '%s' -v --rpkm -o '%s' --mutations  nomuts --blocks" % (inputfile,args.gtf,outputfile)
	logger.info(cmd)
	os.system(cmd)
	#os.system("touch %s_count_output_reads.gtf" % outputfile)
	
def makeCoverageSgrFiles(inputfile,outputfiles):
	""" takes the pyReadCounters gtf output files and generates coverage sgr files """
	outputfile = re.search("(.*)_plus_strand_reads.sgr",outputfiles[0]).group(1)
	cmd = "pyGTF2sgr.py --gtf '%s' --zeros --count -v -o '%s' -c '%s'" % (inputfile,outputfile,args.chromosome)
	logger.info(cmd)
	os.system(cmd)
	#for i in outputfiles:
	#	os.system("touch %s" % i)
		
def makeCoverageBedgraphFiles(inputfile,outputfiles):
	""" takes the pyReadCounters gtf output files and generates bedgraph files for viewing of the data in genome browsers """
	outputfile = re.search("(.*)_plus_strand_reads.bedgraph",outputfiles[0]).group(1)
	cmd = "pyGTF2bedGraph.py --gtf '%s' --count -v --permillion -o '%s' -c '%s'" % (inputfile,outputfile,args.chromosome)
	logger.info(cmd)
	os.system(cmd)
	#for i in outputfiles:
	#	os.system("touch %s" % i)

### start of pipeline

# setting the main working directory
home_dir = os.getcwd()
root_dir = "%s/%s" % (home_dir,args.name)
try:
	os.mkdir(root_dir)
except OSError:
	pass

### setting up logging
if not args.log_file:
	logfile = "%s/%s" % (root_dir,"log_file.txt")
else:
	logfile = "%s/%s" % (root_dir,args.log_file)
logger, logger_mutex = cmdline.setup_logging ("ChemModSeqPipeline",logfile,10)

### starting the pipeline commands

pipeline = Pipeline(name = "CRAC_Pipeline_Single_end")

startingfiles = args.forwardreads


pipeline.transform( task_func = runFlexBar,
					input  = startingfiles,
					filter = formatter("^.+/([^/]+).(san)?fastq"),						
					output = "%s/flexbar_trimmed/{1[0]}_trimmed.fastq" % root_dir,
				).follows(pipeline.mkdir(os.path.join(root_dir,"flexbar_trimmed")))
				
if args.barcodes:
	pipeline.subdivide(task_func = demultiplexSamples,
						input  = runFlexBar,
						filter = formatter("^.+/([^/]+).fastq"),
						output = ["%s/demultiplexed/{1[0]}_%s" % (root_dir,barcodestring) for barcodestring in getBarcodeInfo(args.barcodes)]
					).follows(pipeline.mkdir(os.path.join(root_dir,"demultiplexed")))

	pipeline.transform(task_func = trimAndCollapse,
						input  = demultiplexSamples,
						filter = formatter("^.+/([^/]+)_trimmed_(?P<BARCODE>.*).fastq"),						
						output = "%s/collapsed/{1[0]}_trimmed_{BARCODE[0]}.fasta" % root_dir
					).follows(pipeline.mkdir(os.path.join(root_dir,"collapsed"))).follows(demultiplexSamples)

else:
	pipeline.transform(task_func = trimAndCollapse,
						input  = runFlexBar,
						filter = formatter("^.+/([^/]+).fastq"),						
						output = "%s/collapsed/{1[0]}.fasta" % root_dir,
					).follows(pipeline.mkdir(os.path.join(root_dir,"collapsed")))

pipeline.transform( task_func = alignReads,
					input  = trimAndCollapse,
					filter = formatter("^.+/([^/]+).fasta"),
					output = "%s/novo_files/{1[0]}.novo" % root_dir,
				).follows(pipeline.mkdir(os.path.join(root_dir,"novo_files")))

pipeline.transform( task_func = runPyReadCounters,
					input  = alignReads,
					filter = formatter("^.+/([^/]+).novo"),
					output = "%s/pyReadCounters_analyses/{1[0]}_count_output_reads.gtf" % root_dir,
				).follows(pipeline.mkdir(os.path.join(root_dir,"pyReadCounters_analyses")))

pipeline.transform(task_func = runSecondPyReadCounters,
					input  = alignReads,
					filter = formatter("^.+/([^/]+).novo"),
					output = "%s/pyReadCounters_blocks_nomuts_analyses/{1[0]}_blocks_nomuts_count_output_cDNAs.gtf" % root_dir,
				).follows(pipeline.mkdir(os.path.join(root_dir,"pyReadCounters_blocks_nomuts_analyses"))).follows(alignReads)

pipeline.subdivide(task_func = makeCoverageSgrFiles,
					input  = runSecondPyReadCounters,
					filter = formatter("^.+/([^/]+)_count_output_cDNAs.gtf"),
					output = ["%s/sgr_files/{1[0]}_plus_strand_reads.sgr" % root_dir,"%s/sgr_files/{1[0]}_minus_strand_reads.sgr" % root_dir]
				).follows(pipeline.mkdir(os.path.join(root_dir,"sgr_files"))).follows(runSecondPyReadCounters)
				
pipeline.subdivide(task_func = makeCoverageBedgraphFiles,
					input  = runSecondPyReadCounters,
					filter = formatter("^.+/([^/]+)_count_output_cDNAs.gtf"),
					output = ["%s/bedgraph_files/{1[0]}_plus_strand_reads.bedgraph" % root_dir,"%s/bedgraph_files/{1[0]}_minus_strand_reads.bedgraph" % root_dir]
				).follows(pipeline.mkdir(os.path.join(root_dir,"bedgraph_files"))).follows(runSecondPyReadCounters)

pipeline_run(multiprocess=args.processors,verbose=5,checksum_level=3)
