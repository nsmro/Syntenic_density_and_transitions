#! usr/local/env python3
import sys
import re

"""
Usage: gff_preparation my_file.gff 
for SkowalevskiiJGIv3.0.longestTrs.gff3
SACKO names are gene names (i.e. not on CDS lines) We want them to be here.
"""


d_ID = {}
search_string_pacid = 'pacid=(\d+);*'
search_string_protID = 'Name=(\w+);*'

input_gff = sys.argv[1]

with open(input_gff, "r") as f:
		for line in f:
			if line.startswith("#"):
				pass
			else:
				line = line.strip()
				line = line.split('\t')
				if 'mRNA' in line[2]:
					try:
						d_value = re.search(search_string_protID, line[8]).group(1) #prot ID is value
						d_key = re.search(search_string_pacid, line[8]).group(1) #pacid is key
						if d_key not in d_ID.keys():
							d_ID[d_key] = d_value
					except:
						pass #skips line if any re.search is "None"

with open("output.txt", "w") as g:
	with open(input_gff, "r") as f:
		for line in f:
			line = line.strip()
			parsed_line = line.split('\t')
			if len(parsed_line) == 9:
				if parsed_line[2] == 'CDS':
					line_pacid = re.search(search_string_pacid, parsed_line[8]).group(1)
					output_line = line+';protein_id='+d_ID[line_pacid]
					print(output_line)
				else:
					output_line = line
					print(output_line)
			else:
				print(line)
