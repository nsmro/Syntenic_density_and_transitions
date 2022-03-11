#! /usr/bin/env python3

import pandas as pd
import numpy as np
import scipy.stats
import itertools
import collections
import argparse


parser = argparse.ArgumentParser(description = 'takes a tidy dataframe as input, made with make_tidy_density_df.py')
parser.add_argument('-t', '--results_df', help = 'tidy dataframe, one observation by line', required = True)
parser.add_argument('-s', '--full_synt_file', help = 'original (unfiltered synt file)', required = True)
parser.add_argument('-e', '--expression_data', help = 'full tpm table')
parser.add_argument('-p', '--prefix', help = 'prefix to add to accessions in expression data')
parser.add_argument('-o', '--output', help = 'output prefix', required = True)

args = parser.parse_args()


def load_expr(expr):
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
            if args.prefix != '':
                transcript_id = f'{args.prefix}_{transcript_id}'
            if tpms.max() == 0:
                pass
            else:
                output_dict[transcript_id] = tpms
    return output_dict


def block_correlation(acc_str):
    """
    computes block correlation from block_id
    uses exp_dict(names: numpy 1D arrays with expression per stage)
    and block_dict(block_id: list of expressed genes of the block)
    :param myid: block_id
    :returns: block correlation
    """
    acc_ls = acc_str.split(',')
    corrs = []
    to_remove = []
    for gene in acc_ls:
        if exp_dict.get(gene) is None:
            to_remove.append(gene)
    acc_ls = list(set(acc_ls) - set(to_remove))
    if len(acc_ls) < 3:
        return np.nan
    else:
        for genea, geneb in itertools.combinations(acc_ls, 2):
            exp_genea = exp_dict[genea]
            exp_geneb = exp_dict[geneb]
            corr = scipy.stats.spearmanr(exp_genea, exp_geneb).correlation
            if corr == 1:
                corr = np.float64(0.9999999999999999) #arctanh(1) = inf, mean of any array with one infinity value is infinity, and tanh(inf) = 1. This avoids block correlation of one when one pair has a value of 1.
            if corr == -1:
                corr = np.float64(-0.9999999999999999)            
            corrs.append(corr)
        corrs = np.array(corrs)
        block_corr = np.tanh(np.mean(np.arctanh(corrs)))
        return block_corr


"""

tmp_args = collections.namedtuple('tmp_args', ['results_df', 'full_synt_file', 'expression_data', 'prefix', 'output'])

args = tmp_args('../../02_REDUX_gene_density_analysis/density_whole_genome/key_nodes.tidydf.csv',
                 '/scratch/robert/2019_microsynteny_size_constraints/01_microsynteny/chrom_of/5.blocks.3.syn_corrected.synt',
                 '../CALMI/CALMI_transcript_tpms_all_samples.tsv',
                 'CALMI',
                 'density_blockcorrelation')


"""

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


exp_dict = load_expr(args.expression_data)

species_df_obs = results_df.query("random == 'observed' & species == @args.prefix")
species_df_obs['block_corr'] = species_df_obs['acc_ls'].map(lambda x: block_correlation(x))
species_df_obs = species_df_obs.dropna()

blocks_to_keep = species_df_obs.block_id.tolist()

species_df_rand = results_df.query("random == 'random' & species == @args.prefix & block_id in @blocks_to_keep" )
species_df_rand['block_corr'] = species_df_rand['acc_ls'].map(lambda x: block_correlation(x))
species_df_rand = species_df_rand.dropna()

species_df = pd.concat([species_df_obs, species_df_rand])

species_df.to_csv(args.output, index = False)


