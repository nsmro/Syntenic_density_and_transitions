#! /usr/bin/env python3

import argparse
import csv
import glob
import itertools
import os
import subprocess


parser = argparse.ArgumentParser(description = 'takes chom and genome folders, lists the files, returns a df with whole genome length gene counts and species')
parser.add_argument('-c', '--chrom', help = 'chrom folder', required = True)
parser.add_argument('-g', '--genome', help = 'genome folder', required = True)
parser.add_argument('-o', '--output', help = 'output name', required = True)
args = parser.parse_args()


def parse_fasta(handle):
    """
    parses fasta file
    :param handle: filehandle of a fasta file
    :return: a generator, pairs of header/sequence items. Sequence devoid of newline characters.
    """
    fasta_iter = (list(g) for _,g in itertools.groupby(handle, lambda l: l.startswith('>')))
    for header_group, sequence_group in zip(*[fasta_iter]*2):
        header_string = ''.join(header_group).lstrip('>').rstrip()
        header = ''.join(header_group).lstrip('>').rstrip().split()[0]
        seq = ''.join([line.rstrip() for line in sequence_group])
        yield header, seq


def count_file_lines(file_path):
    """
    Counts the number of lines in a file using wc utility.
    :param file_path: path to file
    :return: int, no of lines
    """
    num = subprocess.check_output(['wc', '-l', file_path])
    num = num.split()
    return int(num[0])



chrom_ls = sorted(glob.glob(f'{args.chrom}/*chrom'))
genome_ls = sorted([x for x in glob.glob(f'{args.genome}/*genome') if os.path.isdir(x) is False])

print(chrom_ls)

print(genome_ls)

results = [['species', 'genome_length', 'gene_count']]
for chrom, genome in zip(chrom_ls, genome_ls):
    prefix = os.path.basename(chrom).split('.')[0]
    print(f'prefix is: {prefix}')
    genome_length = 0
    print(f'parsing genome file {genome}...')
    with open(genome) as f:
        for _, seq in parse_fasta(f):
            genome_length += len(seq)
    print(f'Done!')
    print(f'parsing chrom file {genome}...')
    print(f'Done!')
    gene_count = count_file_lines(chrom)
    results.append([prefix, genome_length, gene_count])

with open(f'{args.output}.csv', 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerows(results)
