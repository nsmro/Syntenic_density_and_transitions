#! /usr/bin/env python3

import sys
import os
import re
import itertools
import argparse
import collections


parser = argparse.ArgumentParser(description = 'Take as input a chrom file, a gff and a fasta of transcripts')
parser.add_argument('-c',
                    '--chrom',
                    help = 'chrom file',
                    required = True)

parser.add_argument('-a',
                    '--annot',
                    help = 'gff file. script works best if mitochondrial genomes are out',
                    required = True)
parser.add_argument('-f',
                    '--fasta',
                    help = 'transcript fasta file',
                    required = True)
parser.add_argument('-ft',
                    '--feature_transcript',
                    help ='feature line where transcripts are mapped to parents',
                    type = str,
                    required = True)
parser.add_argument('-fp',
                    '--feature_protein',
                    help ='feature line where proteibns are mapped to parents',
                    type = str,
                    required = True)
parser.add_argument('-kpp',
                    '--key_parent_protein',
                    help ='gff key, where "key_parent=parent on lines yhere ye look at protein"',
                    type = str,
                    required = True)
parser.add_argument('-kpt',
                    '--key_parent_transcript',
                    help ='gff key, where "key_parent=parent on lines yhere ye look at protein"',
                    type = str,
                    required = True)
parser.add_argument('-kp',
                    '--key_protein_id',
                    help ='gff key, where "key_protein_id=accession"',
                    type = str,
                    required = True)
parser.add_argument('-kt',
                    '--key_transcript_id',
                    help ='gff key, where "key_transcript=accession"',
                    type = str,
                    required = True)
parser.add_argument('-o',
                    '--output',
                    help = 'prefix of the output file, default is prefix parsed form the chrom. results saved in PREFIX.binsize.density',
                    default = False)
args = parser.parse_args()

def parse_chrom(handle):
    """
    :param handle: handle of the chrom file
    :returns: accession list, no prefix
    """
    output_list = []
    for line in handle:
        prefix, acc, *_ = line.rstrip().split('\t')
        finalacc = acc.replace(f'{prefix}_', '')
        output_list.append(finalacc)
    return output_list

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

def parse_gff(handle):
    """
    parse gff file in to dict
    :param handle:
    :returns: dict transcript_id: protein_id
    """
    parent_transcript_dict = {}
    protein_parent_dict = {}
    re_key = '\s*[:="\']*([^,|;"\']+)[;|"\s\n]*'
    p_parent_re = args.key_parent_protein + re_key
    t_parent_re = args.key_parent_transcript + re_key
    transcript_re = args.key_transcript_id + re_key
    protein_re = args.key_protein_id + re_key
    for line in handle:
        if (line.startswith('#') or len(line.split('\t')) < 9):
            pass
        else:
            fields = line.rstrip().split('\t')
            _, _, feature, *_,comments = fields
            if feature == args.feature_transcript:
                try:
                    tran_id = re.search(transcript_re, comments).group(1)
                    parent_id = re.search(t_parent_re, comments).group(1)
                    parent_transcript_dict[parent_id] = tran_id
                except AttributeError:
                    pass
            elif feature == args.feature_protein:
                try:
                    prot_id = re.search(protein_re, comments).group(1)
                    parent_id = re.search(p_parent_re, comments).group(1)
                    protein_parent_dict[prot_id] = parent_id
                except AttributeError:
                    pass
    print(list(protein_parent_dict.items()))
    output_dict = {parent_transcript_dict[parent]:prot for prot, parent in protein_parent_dict.items()}
    return output_dict

with open(args.annot, 'r') as annot:
    acc_dict = parse_gff(annot)

with open(args.chrom, 'r') as chrom:
    non_redundant_prots = parse_chrom(chrom)

with open(args.fasta, 'r') as fasta, open(args.output, 'w') as out:
    for header, seq in parse_fasta(fasta):
        try:
            header_prot = acc_dict[header]
            if header_prot in non_redundant_prots:
                out.write(f'>{header_prot}\n{seq}\n')
        except KeyError:
            pass







