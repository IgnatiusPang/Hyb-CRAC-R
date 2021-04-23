#!/usr/bin/python

import sys

gff = open(sys.argv[1])
GTF = open(sys.argv[1] + ".gtf", "w")

for line in gff:
	if not line.startswith("#"):
		feat_info = {}
		i = line.split('\t')
		feat = i[8].split(';')
		for x in feat:
			y = x.split('=')
			feat_info[y[0]] = y[1]
		id = feat_info["ID"].strip()
		if "Name" in feat_info.keys():
			name = feat_info["Name"].strip()
		else:
			name = feat_info["ID"].strip()
		rna_type = i[2].strip()
		allowedfeatures = ["exon","CDS","5UTR","3UTR","start_codon","stop_codon"]
		if i[2] not in allowedfeatures:
			i[2] = "exon"
		info = []
		for k in feat_info.keys():
			info.append(k + ' "' + feat_info[k].strip() + '";')
		GTF.write("\t".join(i[0:8]) + '\tgene_id "' + id + '"; transcript_id "' + name + '"; rna_type "' + rna_type + '";' + "".join(info) + '\n')
