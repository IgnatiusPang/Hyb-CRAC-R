#!/usr/bin/python
#
# copyright 2017 Jai J. Tree 
# GNU Lesser General Public License
#
# script to reformat merged hybrid table from R sqldf join 
# into a format that can be recycled through sqldf 

import math
import fileinput 

def which_overlap(hybIN, hybA, hybB):
	""" method to find overlapping hybrid half supply tuples of (start, end). Returns 0 (hybA overlap), 1 (hybB overlap)"""
	in_start, in_end = hybIN
	a_start, a_end = hybA
	b_start, b_end = hybB
	if in_start <= a_start and in_end >= a_start:
		return False
	elif in_start <= a_end and in_end >= a_end:
		return False
	elif in_start <= b_start and in_end >= b_start:
		return True 
	elif in_start <= b_end and in_end >= b_end:
		return True 

def connection_score(hyb_count, gx, gy):
	""" method to return a cumulative connection score from sum of hybrids and max_reads"""		
	return hyb_count / (math.sqrt(gx * gy))
	
## reformat input table

l = ["ID"]*32											# hold last line in l

for line in fileinput.input():
	x = line.split("\t")								# new line = x
	if line.startswith("ID"):
		print "ID\tg_x\thyb_count\trna_class\tchromo\tname\tstart\tend\tstrand\tID.1\tg_y\thyb_count.1\trna_class.1\tchromo.1\tname.1\tstart.1\tend.1\tstrand.1\tmin\tmax\tp_x_y\tpdf_x_y\tp_value\tbh_adj_p_value\tconnection_score\ttotal_hybrids\tcumulative_p_value\tcumulative_connection_score\tnumber_expts"
	elif l[1] == "ID":
		l=x
	elif x[0:29] == l[0:29] and x[30] != "NA":			# sqldf duplicate entry found - assign x to l 
		l[49] = str( min(float(x[49]), float(l[49]) ) )	# update p_x_y
		l[50] = str( min(float(x[50]), float(l[50]) ) )	# update pdf_x_y
		l[51] = str( min(float(x[51]), float(l[51]) ) )	# update p_value
		l[52] = str( min(float(x[52]), float(l[52]) ) )	# updaate bh_adj_p_valu
		l[53] = str( max(float(x[53]), float(l[53]) ) )	# update con_score_max
		l[31] = str(int(x[31]) + int(l[31]))			# update total_hybrids

		l[56] = str( float(x[56]) * float(l[56]) ) 		# update cum_p_val
		l[57] = str( l[57] )			# don't update num expts - +1 is added once new line doesn't match
			
	elif not l[0].startswith("ID") and x[0:29] != l[0:29] and l[30] != "NA":
		chiastic = which_overlap( (int(l[6]),int(l[7])) , (int(l[35]),int(l[36])) , (int(l[44]),int(l[45])) )
		if chiastic:
			g_x = str( max(int(l[1]) , int(l[39]) ) )
			g_y = str( max(int(l[10]) , int(l[30]) ) )
			start = str( min(int(l[6]), int(l[44]) ) )
			end = str( max(int(l[7]), int(l[45]) ) )
			start_1 = str( min(int(l[15]), int(l[35]) ) )
			end_1 = str( max(int(l[16]), int(l[36]) ) ) 
		else:
			g_x = str( int(l[1]) + int(l[30]) )
			g_y = str( int(l[10]) + int(l[39]) )
			start = str( min(int(l[6]), int(l[35]) ) )
			end = str( max(int(l[7]), int(l[36]) ) )
			start_1 = str( min(int(l[15]), int(l[44]) ) )
			end_1 = str( max(int(l[16]), int(l[45]) ) )

		id = l[0] + "," + l[29]
		hyb_count = str(max(int(l[2]), int(l[31]) ))
		rna_class = l[3]
		chromo = l[4]
		name = l[5]
		strand = l[8]
		rna_class_1 = l[12]
		chromo_1 = l[13]
		name_1 = l[14]
		strand_1 = l[17]
		min_reads = str( min(int(g_x), int(g_y) ) )
		max_reads = str( max(int(g_x), int(g_y) ) )
	
		p_x_y = str( min(float(l[20]), float(l[49]) ) )
		pdf_x_y = str( min(float(l[21]), float(l[50]) ) )
		p_value = str( min(float(l[22]), float(l[51]) ) )
		bh_adj_p_value = str( min(float(l[23]), float(l[52]) ) )
		con_score_max = str( max(float(l[24]), float(l[53]) ) )
		total_hybrids = str(int(l[2]) + int(l[31]))
		if int(total_hybrids) == 0 or int(g_x) == 0 or int(g_y) == 0:
			cuml_con_score = con_score_max
		else:
			cuml_con_score = str(connection_score(float(total_hybrids), int(g_x), int(g_y) ))
		cuml_p_val = str( float(l[27]) * float(l[56]) )
		num_expt = str( int(l[28]) + 1 )
	
		out = [id, g_x, hyb_count, rna_class, chromo, name, start, end, strand, id, g_y, hyb_count, rna_class_1, chromo_1, name_1, start_1, end_1, strand_1, min_reads, max_reads, p_x_y, pdf_x_y, p_value, bh_adj_p_value, con_score_max, total_hybrids, cuml_p_val, cuml_con_score, num_expt] 
		print "\t".join(out)
		l=x

	elif x[0:29] != l[0:29] and l[30] == "NA":
		print "\t".join(l[0:29])
		l=x

	
