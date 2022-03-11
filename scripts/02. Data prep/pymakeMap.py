#! /usr/bin/env python3

import argparse
import sys
import os
import re
from operator import itemgetter
import csv
import itertools


def parsegff(inputfile):
	"""parsegff ouputs nested list architecture. List of lists of CDS sharing same protein of lists sharing the same chromosome.
	e.g. Protein1 has two locus assigned to it:
	[[[Chr1, Protein1, CDS1][Chr1, Protein1, CDS2][Chr2, Protein1, CDS1][Chr2, Protein1, CDS2]]]
	We'll keep only the first group of each sublevel (i.e. one locus per protein).
	In the case of there is only one locus per protein, we just end up with a protein of length 1, and sublist_tmp[0] extracts it.
	In NCBI, NC_ taken in priority (locus on a chromosome), then NT_ then NW_ ))  
	"""
	list_raw_chrom = []
	with open(inputfile, 'r', encoding="ISO-8859-1") as f:
		for line in f:
			if (line.startswith('#') or len(line.split('\t')) < 9):
				pass
			else:
				line = line.strip()
				line = line.split('\t')
				if line[2] == feature:
					try:
						seq_name = re.search(seq_name_search_string, line[8]).group(1)
						seq_chrom = line[0]
						seq_start = line[3]
						seq_stop =  line[4]
						seq_orientation = line[6]
						seq_name_w_prefix = seq_prefix+'_'+seq_name
						list_temp = [seq_prefix,seq_name_w_prefix, seq_chrom, seq_orientation, seq_start, seq_stop]
						list_raw_chrom.append(list_temp)
					#in certain cases, CDS do not have protein product (e.g. in NCBI gffs "exception=rearrangement required for product")
					except:
						pass
				else:
					pass
	list_groupby_protein_acc = []
	for _, g_same_protein in itertools.groupby(sorted(list_raw_chrom, key = itemgetter(1)), key = itemgetter(1)):
		list_groupby_protein_acc.append(list(g_same_protein))
	for sublist in list_groupby_protein_acc:
		sublist_tmp = []
		for _, g_same_chrom in itertools.groupby(sorted(sublist, key = itemgetter(2)), key = itemgetter(2)):
			sublist_tmp.append(list(g_same_chrom))
		yield sublist_tmp[0]#Just the first locus. The sublist is sorted alphabetically so in the case of NCBI proteins, locus located on DNA with accessions starting with NC_ (longest, if there is conflict with alternative loci) are used.
	print('gff file loaded!')

def storechrom(parsed_gff_file):
	y = 0
	for sublist in parsed_gff_file:
		y += 1
		ls_coordinates = []
		for x in range(0,len(sublist)):
			ls_coordinates.append([int(sublist[x][4]),int(sublist[x][5])])
		yield[sublist[0][0], sublist[0][1]+suffix, sublist[0][2], sublist[0][3], ls_coordinates]	#keep exon information.
	print('{} sequences found in the gff.'.format(str(y)))

def fastafilteredchrom(unfilteredchrom, fastafile):
	ls_accessions = []
	with open(fastafile, 'r') as f:
		for line in f:
			if line.startswith('>'):
				line = line.lstrip('>')
				line = line.rstrip()
				line = line.split(' ')
				ls_accessions.append(line[0])
			else:
				pass
	print('{} sequences found in the fasta file.'.format(len(ls_accessions)))
	y = 0
	for seq_record in unfilteredchrom:
		if seq_record[1] in ls_accessions:
			y += 1
			yield seq_record
			del ls_accessions[ls_accessions.index(seq_record[1])] #reduce the list size, so that when we iterate through again, the list is smaller
	print('{} sequences from the gff will be retained.'.format(y))


def clusfilteredchrom(unfilteredchrom, clusfile):
	ls_accessions = []
	with open(clusfile, 'r') as f:
		for line in f:
			line = line.rstrip()
			line = line.split('\t')
			species_count = 0
			for x in line [2:]:
				species_count += seq_prefix in x
			if species_count < len(line[2:]): # can only be lower or equal, not more.
				ls_accessions.extend(line[2:])			
			else:
				pass
	ls_accessions_prefix = [accession for accession in ls_accessions if accession.startswith(seq_prefix)]  
	y = 0
	for seq_record in unfilteredchrom:
		if seq_record[1] in ls_accessions_prefix:
			y += 1
			yield seq_record
			del ls_accessions_prefix[ls_accessions_prefix.index(seq_record[1])]
	print('{} sequences in the clus file start with {} and are not species-specific.'.format(str(y), seq_prefix))

def test_transcript_overlap(ls_exons_transcript1,ls_exons_transcript2):
	"""
	testing the overlap of two transcripts.
	We also check that the overlap is of at least 20 bp as minimum overlap to avoid exons such as [100,3000] and [3000,400] to be considered overlapping.
	"""
	start1 = min(list(itertools.chain.from_iterable(ls_exons_transcript1)))
	end1 = max(list(itertools.chain.from_iterable(ls_exons_transcript1)))
	start2 = min(list(itertools.chain.from_iterable(ls_exons_transcript2)))
	end2 = max(list(itertools.chain.from_iterable(ls_exons_transcript2)))
	length_overlap = min(abs(end1 - start1), abs(end1 - start2), abs(end2 - start1), abs(end2 - start2))
	overlap = max(start1, start2) > min(end1, end2)
	if overlap and length_overlap > 20:
		return True
	else:
		return False

def test_shared_exon(ls_exons_transcript1,ls_exons_transcript2):
	ls_exons_transcript1_sets = [set(x) for x in ls_exons_transcript1] # transform exons into sets. this way, even if exons are not in the same order in two transcripts, they'll be found to be identical (as {a,b} = {b,a})
	ls_exons_transcript2_sets = [set(x) for x in ls_exons_transcript2]
	return any([x in ls_exons_transcript1_sets for x in ls_exons_transcript2_sets]) #any([]) returns False. The list comprehension looks if one set of coordinates(exon). If no subset is shared, it returns False.

def test_shared_coordinates(ls_exons_transcript1,ls_exons_transcript2):
	"""
	make a list of maxipmum and minimums of exons. If a maximum or a minimum of the exons is shared, there is an overlap.
	this distinguishes between a simple shared coordinate, where if the start and stop of an exon are shared it'd consider thgere is an overlap.
	e.g. exons[1,3] and [2,3] are overlapping, but exons [1,3] and [3,4] are not
	"""
	transcript1_minls = [min(x) for x in ls_exons_transcript1]
	transcript1_maxls = [max(x) for x in ls_exons_transcript1]
	transcript2_minls = [min(x) for x in ls_exons_transcript2]
	transcript2_maxls = [max(x) for x in ls_exons_transcript2]
	max_shared = any([x in transcript1_maxls for x in transcript2_maxls])
	min_shared = any([x in transcript1_minls for x in transcript2_minls])
	return any([min_shared,max_shared])

def assign_ids_to_pair(AN1, AN2, Dict_AN):
	if AN1 in Dict_AN.keys() and AN2 in Dict_AN.keys():
		if Dict_AN[AN1] == Dict_AN[AN2]: # Check if the two gene IDs are the same
			pass
		elif Dict_AN[AN1] != Dict_AN[AN2]: #if gene IDs are distinct, we'll fuse the two groups
			for AN, ID in Dict_AN.items():
				if ID == Dict_AN[AN2]:
					Dict_AN[AN] = Dict_AN[AN1]
					Dict_AN[AN2] = Dict_AN[AN1]
	elif AN1 not in Dict_AN.keys() and AN2 not in Dict_AN.keys():
		new_id = max(Dict_AN.values()) + 1
		Dict_AN[AN1] = new_id
		Dict_AN[AN2] = new_id
	elif AN1 not in Dict_AN.keys():
		Dict_AN[AN1] = Dict_AN[AN2]
	elif AN2 not in Dict_AN.keys():
		Dict_AN[AN2] = Dict_AN[AN1]

#almostchrom is the output of storechrom, fastafilteredchrom, clusfilteredchrom functions.
#Nested list architecture, racapitulates chrom for the first 4 colums. Only a fifth column is a nested list of exons [[exon1_start, exon1_stop],[exon2_start, exon2_stop]].
def filterchrombysharedexons(almostchrom):
	tmp_chrom = []
	filteredchrom = []
	d_AN_gene_ids = {'':0}#for later, this'll come in handy
	for _,scaffold in itertools.groupby(sorted(almostchrom, key = itemgetter(2)), key = itemgetter(2)): #sublist of this groupy are still chrom, but grouped by scaffold
		sorted_scaffold = sorted(scaffold, key = lambda x: len(x[4]), reverse = True)
		for transcript1 in sorted_scaffold: #sort the genes within the scaffold according to their number of exons. This
			accession_transcript1 = transcript1[1]
			transcript1_exon_ls = transcript1[4]
			for transcript2 in sorted_scaffold[(sorted_scaffold.index(transcript1)+1):]: # loop only starting from the othergene located after the gene already compared genes before to everything else. Also, if indexes are out of bounds, this won't raise an error, as the list will be empty, it won't iterate through it and this do nothing.
				transcript2_exon_ls = transcript2[4]
				if test_transcript_overlap(transcript1_exon_ls, transcript2_exon_ls) is True:
					nb_exons_transcript1 = len(transcript1_exon_ls) #total number of exons = length of list
					nb_exons_transcript2 = len(transcript2_exon_ls)
					length_transcript1 = max(itertools.chain.from_iterable(transcript1_exon_ls)) - min(itertools.chain.from_iterable(transcript1_exon_ls))   # chain to isolate coordinates. Max - min = CDS length
					length_transcript2 = max(itertools.chain.from_iterable(transcript2_exon_ls)) - min(itertools.chain.from_iterable(transcript2_exon_ls)) 
					accession_transcript2 = transcript2[1]
					if nb_exons_transcript1 >= 3:
						if nb_exons_transcript2 >= 3:
							transcript1_internal_exons_ls = transcript1_exon_ls[1:-1]
							transcript2_internal_exons_ls = transcript2_exon_ls[1:-1]
							if test_shared_exon(transcript1_internal_exons_ls, transcript2_internal_exons_ls) is True:
								assign_ids_to_pair(accession_transcript1, accession_transcript2, d_AN_gene_ids)
						elif nb_exons_transcript2 == 2:
							if test_shared_coordinates(transcript1_exon_ls, transcript2_exon_ls) is True:
								assign_ids_to_pair(accession_transcript1, accession_transcript2, d_AN_gene_ids)
						elif nb_exons_transcript2 == 1:
							ls_test = []
							for exon in transcript1_exon_ls:
								ls_test.append(test_transcript_overlap([exon],[transcript2[4][0]]))
							if any(ls_test) is True:
								assign_ids_to_pair(accession_transcript1, accession_transcript2, d_AN_gene_ids)
					elif nb_exons_transcript1 == 2:
						if nb_exons_transcript2 == 2:
							if test_shared_coordinates(transcript1_exon_ls, transcript2_exon_ls) is True:
								assign_ids_to_pair(accession_transcript1, accession_transcript2, d_AN_gene_ids)
						elif nb_exons_transcript2 == 1:
							ls_test = []
							for exon in transcript1_exon_ls:
								ls_test.append(test_transcript_overlap([exon],[transcript2[4][0]]))
							if any(ls_test) is True:
								assign_ids_to_pair(accession_transcript1, accession_transcript2, d_AN_gene_ids)
					elif nb_exons_transcript1 == 1:
						if test_transcript_overlap([transcript1[4][0]],[transcript2[4][0]]) is True:
							assign_ids_to_pair(accession_transcript1, accession_transcript2, d_AN_gene_ids)
			if accession_transcript1 not in d_AN_gene_ids.keys():
				new_id = max(d_AN_gene_ids.values()) + 1
				d_AN_gene_ids[accession_transcript1] = new_id
	for record in almostchrom:
		exon_ls = record[4]
		sum_length_of_exons = sum([max(stop - start, start - stop) for start,stop in exon_ls]) #max is not necessary, but it is yet another check, if exons happen to be reversed
		modified_sublist = record[0:4]#Prefix, accession, scaffold and strand
		start = min(itertools.chain.from_iterable(exon_ls))
		stop = max(itertools.chain.from_iterable(exon_ls))
		modified_sublist.append(start)
		modified_sublist.append(stop)
		modified_sublist.append(d_AN_gene_ids[record[1]])
		modified_sublist.append(sum_length_of_exons)
		tmp_chrom.append(modified_sublist)
	y = 0
	for _,genes in (itertools.groupby(sorted(tmp_chrom, key = itemgetter(6)), key = itemgetter(6))): #groupby gene id
		genes = sorted(genes, key = itemgetter(7), reverse = True) # sort eachy gene the length of exon sum
		y += 1
		yield [str(x) for x in genes[0][0:6]]
	print('{} sequences have been identified as the longest isoform of their respective genes.'.format(y))




if __name__ == "__main__"::
	parser = argparse.ArgumentParser(description="Converts gff into chrom, only one start/stop ouputted per protein .\
		For NCBI, in the case of several loci for one protein script prioritizes chromosome coordinates before scaffolds.\
		Tested on NCBI, ENSEMBL, JGI, and B. lanceolatum gtf files. \
		For any other GFF file, check your outputted chrom file.")
	parser.add_argument("-gff", "--gff_input", help = "Gff input file.", required = True)
	parser.add_argument("-p", "--prefix", help = "PREFIX is the PREFIX used in your fasta file.", required = True)
	parser.add_argument("-f", "--feature", help = "FEATURE is the feature column value in the GFF file (e.g. exon, mRNA, CDS).", required = True)
	parser.add_argument("-k", "--key", help = "KEY is in the last column of the GFF and appears as `KEY=SEQUENCENAME` in the input gff.", required = True)
	parser.add_argument("-o", "--output_chrom", help = "Name of output chrom file.",required = True)
	parser.add_argument("-d", "--delete_redundancy", help = " If two transcripts share one exon (same start and stop, same scaffold), only the longest is kept.", action = 'store_true', default = False)
	parser.add_argument("-r", "--ref_filter", help = "OPTIONAL: reference format for filtering. files are either 'fasta' file or 'clus'. use only with -F", choices = ['fasta','clus'], default = None)
	parser.add_argument("-F", "--filter_file", help = "OPTIONAL: File to be used for filtering. Use only with -r", default = None)
	parser.add_argument("-S", "--suffix", help = "suffix found in the accessions of the fasta but absent of the gff, e.g. _1. Use is your proteins are obtained from transdecoder.", default = '')
	args = parser.parse_args()
	input_gff = args.gff_input
	seq_prefix = args.prefix
	feature = args.feature
	seq_key = args.key
	output_chrom = args.output_chrom
	suffix = args.suffix
	filter_by_exons = args.delete_redundancy
	#if the chrom file does not recover the search string, this is the one to modify
	#It is also possible that the KEY you specified is not located on the FEATURE-positive lines
	seq_name_search_string = (seq_key+'\s*[:="\']*([^,|;"\']+)[;|"\s\n]*')

	unfiltered_chrom = list(storechrom(parsegff(input_gff)))
	
	none_conditions_filtering = [args.ref_filter == None, args.filter_file == None]
	some_conditions_filtering = [args.ref_filter != None, args.filter_file != None]
	
	#I no filtering by exon is required, we test which filters are applied.
	if filter_by_exons is False:
		with open(output_chrom, 'w') as output:
			if all(none_conditions_filtering) is True:
					file_to_use = unfiltered_chrom
			elif all(some_conditions_filtering) is True:
				if args.ref_filter == 'fasta':
					file_to_use = list(fastafilteredchrom(unfiltered_chrom, args.filter_file))
				elif args.ref_filter == 'clus':
					file_to_use = list(clusfilteredchrom(unfiltered_chrom, args.filter_file))	
			else:
				print('Unspecified filtering method or file for accession filtering using fasta/clus file')
			for line in file_to_use:
				output_line = line[0:4]
				output_line.append(str(min(itertools.chain.from_iterable(line[4]))))
				output_line.append(str(max(itertools.chain.from_iterable(line[4]))))
				output.write('\t'.join(output_line)+'\n')
	elif filter_by_exons is True:
		with open(output_chrom, 'w') as output:
			if all(none_conditions_filtering) is True:
				file_to_use = filterchrombysharedexons(unfiltered_chrom)
			elif all(some_conditions_filtering) is True:
				if args.ref_filter == 'fasta':
					file_to_use = filterchrombysharedexons(list(fastafilteredchrom(unfiltered_chrom, args.filter_file)))
				elif args.ref_filter == 'clus':
					file_to_use = filterchrombysharedexons(list(clusfilteredchrom(unfiltered_chrom, args.filter_file)))
			else:
				print('Unspecified filtering method or file for accession filtering using fasta/clus file')
			for line in file_to_use:
				output.write('\t'.join(line)+'\n')







