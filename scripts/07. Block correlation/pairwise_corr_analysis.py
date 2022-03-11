#! /usr/bin/env python3


import pandas as pd
import numpy as np
import scipy.stats
import itertools
import collections
import argparse

parser = argparse.ArgumentParser(description = 'takes a tidy dataframe as input, made with make_tidy_density_df.py, outputs block corr and bp distance between every 2 genes.')
parser.add_argument('-t', '--results_df', help = 'tidy dataframe, one observation by line', required = True)
parser.add_argument('-s', '--full_synt_file', help = 'original (unfiltered synt file)', required = True)
parser.add_argument('-og', '--orthology', help = 'orthology file, clus format', required = True)
parser.add_argument('-e', '--expression_data', help = 'full tpm table')
parser.add_argument('-c', '--chrom', help = 'prefix to add to accessions in expression data')
parser.add_argument('-o', '--output', help = 'output prefix', required = True)

args = parser.parse_args()

"""
tmp_args = collections.namedtuple('args', ['results_df', 'full_synt_file', 'expression_data', 'chrom','orthology'])

args = tmp_args('../../02_REDUX_gene_density_analysis/density_whole_genome/key_nodes.tidydf.csv',
                 '../../01_microsynteny/chrom_of/5.blocks.3.syn_corrected.synt',
                 '/scratch/robert/2019_12_Neuropeptides/TPM_normalizations/CRAGI/CRAGI_transcript_tpms_all_samples.tsv',
                 '../../01_microsynteny/chrom/CRAGI.chrom',
                 '../../01_microsynteny/Orthofinder.clus')
"""

def load_expr(expr, prefix):
    """
    load expression data
    unexpressed genes (0 tpms in all stages) are deleted
    :param expr: expression table, first line is header, firts column is accession
    second column length, rest of fields are expression by stage
    """
    output_dict = {}
    with open(expr, 'r') as f:
        header = f.readline()
        for line in f:
            transcript_id, _, *tpms = line.rstrip().split() #second column is length
            tpms = np.array([np.float64(x) for x in tpms])
            if prefix != '':
                transcript_id = f'{prefix}_{transcript_id}'
            if tpms.max() == 0:
                pass
            else:
                output_dict[transcript_id] = tpms
    return output_dict


def cleanup_acc_ls(acc_str, exp_dict):
    """
    deletes acc for acc ls if they're not in the filtered exp dict
    :param acc_ls: a list of accessions, as strings, comma separated
    :param exp_dict: expresison dictionary, output of exp_dict
    :returns: filtered acc_ls (as list object)
    """
    acc_ls = acc_str.split(',')
    to_remove = []
    for gene in acc_ls:
        if exp_dict.get(gene) is None:
            to_remove.append(gene)
    acc_ls = list(set(acc_ls) - set(to_remove))
    return acc_ls


def load_ortho(ortho):
    """
    load orthology
    :param ortho: Orthofinder output, clus format
    :returns: dict, acc as keys, OG as values
    """
    output_dict = {}
    with open(ortho, 'r') as f:
        for line in f:
            og, _, *acc_ls = line.rstrip().split('\t')
            for acc in acc_ls:
                output_dict[acc] = og
    return output_dict


order = ['Vertebrate', 'Tunicate', 'Cephalochordate', 'Ambulacrarian', 'Lophotrochozoan', 'Ecdysozoan', 'Acoel', 'Cnidarian', 'Placozoan', 'Ctenophore', 'Poriferan', 'Metazoa_outgroup']
species_dict = {'Vertebrate' : ['HOMSA','MUSMU','CHEMY','GALGA','XENTR','LATCH','MAYZE','HIPCO','DANRE','LEPOC','CALMI'],
                'Tunicate': ['CIOIN'],
                'Cephalochordate' : ['BRALA'],
                'Ambulacrarian' : ['SACKO','PTYFL','STRPU','ACAPL'],
                'Lophotrochozoan' : ['CAPTE', 'EUPSC','LOTGI','MIZYE','CRAGI','HELRO','ADIVA','LINAN','SCHME'],
                'Ecdysozoan' : ['DROME','ANOGA','TRICA','DAPPU','STRMA','IXOSC','PARTE','CAEEL'],
                'Acoel' : ['HOFMI'],
                'Cnidarian' : ['NEMVE','EXAPA','ACRMI','HYDVU','CLYHE','AURAU'],
                'Placozoan' : ['HOIHO','TRIAD'],
                'Ctenophore' : ['PLEBA','MNELE'],
                'Poriferan' : ['SYCCI','AMPQU'],
                'Metazoa_outgroup' : ['SALRO','CAPOW'],
}


order_sp = [species for taxon in order for species in species_dict[taxon]]



results_df = pd.read_csv(args.results_df)

chrom_df = pd.read_csv(args.chrom, sep = '\t', names = ['prefix', 'accession', 'chromosome', 'strand', 'start', 'end'])
prefix = chrom_df.prefix.values[0]
ortho_dict = load_ortho(args.orthology)
exp_dict = load_expr(args.expression_data, prefix)

full_df = pd.read_csv(args.results_df)
species_df = full_df.query('species == @prefix')

with open(args.output, 'w') as g:
    g.write(f'genea\tgeneb\tblock_id\tnode\tpara\tcorr\tdist\trandom\n')
    for row in species_df.values.tolist():
        node, _, species, random, block_id, _, _, acc_str, *_ = row
        acc_ls = cleanup_acc_ls(acc_str, exp_dict)
        for genea, geneb in itertools.combinations(acc_ls,2):
            exp_genea = exp_dict[genea]
            exp_geneb = exp_dict[geneb]
            coords_df = chrom_df.query('accession in [@genea, @geneb]')[['start', 'end']]
            coords_ls = sorted([int(y) for x in coords_df.values for y in x])
            corr = scipy.stats.spearmanr(exp_genea, exp_geneb).correlation
            dist = coords_ls[2] - coords_ls[1]
            if ortho_dict[genea] == ortho_dict[geneb]:
                para = 'para'
            else:
                para = 'not_para'
            g.write(f'{genea}\t{geneb}\t{block_id}\t{node}\t{para}\t{corr}\t{dist}\t{random}\n')


