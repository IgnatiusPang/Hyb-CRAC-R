#!/usr/bin/python
# Class that can be imported to strip out the features of a GTF
# copyright 2017 Jai J. Tree 
# GNU Lesser General Public License

class GTF_stripper():
	"""Class that allows the features of a GTF line to be returned as a dictionary. Feature names are keys: GTF[chromosome], GTF[source], GTF[feature], GTF[start], GTF[end], GTF[strand], GTF[transcript_id]"""
	def __init__(self, GTF_line):
		"""init method to put features into dicitonary"""
		GTF_features = GTF_line.split("\t")
		
		self.GTF = dict()
		self.GTF["chromosome"] = GTF_features[0]
		self.GTF["source"] = GTF_features[1]
		self.GTF["feature"] = GTF_features[2]
		self.GTF["start"] = int(GTF_features[3])
		self.GTF["end"] = int(GTF_features[4])
		self.GTF["score"] = GTF_features[5]
		self.GTF["strand"] = GTF_features[6]
		self.GTF["frame"] = GTF_features[7]
		
		description = GTF_features[8].split('"')
		self.GTF["gene_name"] = description[1]

