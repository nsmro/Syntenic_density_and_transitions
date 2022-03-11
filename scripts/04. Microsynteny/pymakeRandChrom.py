#! /usr/bin/env python3

import argparse
import sys
import os
import numpy as np
import pandas as pd

parser = argparse.ArgumentParser(description = "Randomizes the gene order of a chrom file. Outputs randomized chrom with extension .rand.n in the working directory")
parser.add_argument('input', type = str, nargs = '?', default = None)
parser.add_argument('-f', '--filetype', choices = ['chrom','list'], help = 'type of input: either chrom file directly or a newline-separated list of chrom files. If -f is not used, will try to get filenames from stdin.')
parser.add_argument('-n', '--number', help = 'number of randomized chroms to generate', default = 1, type = int)
parser.add_argument('-o', '--output_directory', help = 'working directory to use. If none is provided, randomized chrom files will be written in current directory.', default = os.getcwd())


args = parser.parse_args()

number_randomizations = args.number + 1
accession_column = 1

filelist = []
if args.filetype == 'list':
	with open(args.input, 'r') as f:
		for filename in f:
			filename = filename.rstrip()
			filelist.append(filename)
elif args.filetype == 'chrom':
	filelist.append(args.input)
elif sys.stdin.isatty() is False:
	for line in sys.stdin:
		line = line.rstrip()
		filelist.append(line)
else:
	print('WARNING! Chrom files to randomize should be provided through a filelist, directly in the command, or through stdin.')

for filename in filelist:
	for n in range(1, number_randomizations):
		output_filename = '{}/{}.rand.{}'.format(args.output_directory, os.path.basename(filename), str(n))
		current_chrom = pd.read_csv(filename, header = None, sep = '\t')
		current_chrom[accession_column] = np.random.permutation(current_chrom[accession_column])
		with open(output_filename, 'w') as f:
			f.write(current_chrom.to_csv(sep='\t', index=False, header=False))
