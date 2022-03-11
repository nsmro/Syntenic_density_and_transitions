#! /usr/bin/env python3

import argparse
import collections
import glob
import itertools
import numpy as np
import os
import pandas as pd
import sys

parser = argparse.ArgumentParser(description = """
makes a single tidy df with gene densities by block (one line per iteration)
finds para/not para blocks. random blocks para attributes are inferrend from the observed counterpart
provides two lists off accs. acc_ls is syntenic blocks, all_acc_ls also includes intervening genes
""")
parser.add_argument('-s', '--synt_folder', help = """folder where synt files, filename should start by nodename, dot-separated, extension should be *synt. random should appear in the name.
                                                    e.g. Bilateria.novel.synt, and the random blocks Bilateria.novel.random.synt""", default = os.getcwd())
parser.add_argument('-g', '--genome_folder', help = 'folder where genome files are located, prefix in the name dot separated. e.g. EUPSC.lachesis201904.genome.fa, CRAGI.NCBIgenome.fasta', required = True)
parser.add_argument('-c', '--chrom_folder', help = 'folder where chrom files are located, files should have the *chrom extension', required = True)
parser.add_argument('-m', '--multi_sp', help = 'Multi species block (total), used to get multi_sp ID', required = True)
parser.add_argument('-og', '--ortho', help = 'orthology file, clus format', required = True)
parser.add_argument('-o', '--output', help = 'output prefix', default = 'density')
args = parser.parse_args()


#tmp_args = collections.namedtuple('tmp_args', ['synt_folder', 'genome_folder', 'chrom_folder', 'multi_sp', 'ortho', 'output'])
#args = tmp_args(os.getcwd(), '../../genomes/', '../../../01_microsynteny/chrom/', '../../../01_microsynteny/chrom_of/5.blocks.3.syn.clusters', '../../../01_microsynteny/Orthofinder.clus', 'key_nodes_w_consecutive_pairs')

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


def parse_blocks(syntfile):
    """
    Get coordinates of all the blocks in a file
    :param syntfile : synt file
    :return: dictionary, keys are species prefixes, values are a list of coordinates
    """
    block_dict = {}
    block = collections.namedtuple('block',
                                    ['id','chromosome', 'start', 'end', 'acc_ls'])
    with open(syntfile, 'r') as f:
        for line in f:
            line = line.rstrip().split('\t')
            block_id, species, _, _, _, _, _, coords, _, *rest = line #works with output of pick_random_blocks
            acc_ls = rest[0] #acc_ls is the first item there
            chromosome, start, stop = coords.replace('..', ':').split(':')
            if block_dict.get(species, False) is False:
                block_dict[species] = [block(block_id, chromosome, int(start), int(stop), acc_ls)]
            else:
                block_dict[species].append(block(block_id, chromosome, int(start), int(stop), acc_ls))
    return block_dict


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


def block_para(acc_ls):
    """
    Determines whether a block is para or not para based on OG composition
    if more than 50 % of genes belong to same OGm it's para
    :param myid: acc_ls, string with comma-separated accessions
    :returns: para/not_para
    """
    acc_ls = acc_ls.split(',')
    og_ls = [ortho_dict[acc] for acc in acc_ls]
    counts = []
    for og in set(og_ls):
        counts.append(og_ls.count(og))
    para_conds = [count > (0.4 * len(og_ls)) for count in counts]
    if any(para_conds) is True:
        return 'para'
    else:
        return 'not_para'


def median_consecutivedist(complete_acc_ls):
    """
    when provided a list of accessions, calculates the mean intergenic distance between genes of the list 
    """
    dist_ls = collections.deque()
    acc_ls = complete_acc_ls.split(",")
    prefix = acc_ls[0].split('_')[0]
    df_coords = chrom_dict[prefix].query('accession in @acc_ls').sort_values(by = "start")
    acc_ls = df_coords.accession.tolist()
    start_ls = df_coords.start.tolist()
    end_ls = df_coords.end.tolist()
    for i in range(len(acc_ls) -1):
        j = i+1
        acc1,end_1 = acc_ls[i], end_ls[i]
        acc2, start2 = acc_ls[j], start_ls[j]
        intergenic_dist = start2 - end_1
        if intergenic_dist < 0: intergenic_dist = 0
        dist_ls.append(intergenic_dist)
    return statistics.median(dist_ls)



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

species_dict_r = {species: taxon for taxon,species_ls in species_dict.items() for species in species_ls}


synt_obs_ls = [file for file in glob.glob(f'{args.synt_folder}/*synt') if 'random' not in file]
synt_rand_ls = [file for file in glob.glob(f'{args.synt_folder}/*synt') if 'random' in file]

nodes_obs_dict = {os.path.basename(syntfile).split('.')[0]: parse_blocks(syntfile) for syntfile in synt_obs_ls}
nodes_rand_dict = {os.path.basename(syntfile).split('.')[0]: parse_blocks(syntfile) for syntfile in synt_rand_ls}


print('Loading genome info...', file = sys.stderr)

genomes_ls = [file for file in glob.glob(f'{args.genome_folder}/*') if 'genome' in file]
genomes_lengths_dict = {}
for genome_file in genomes_ls:
    prefix = os.path.basename(genome_file).split('.')[0]
    print(prefix)
    with open(genome_file, 'r') as f:
        lengths_d = {header: len(sequence) for header, sequence in parse_fasta(f)}
        genomes_lengths_dict[prefix] = lengths_d

print('Done!\n', file = sys.stderr)


print('Processing chrom...', file = sys.stderr)

chrom_ls = [file for file in glob.glob(f'{args.chrom_folder}/*chrom')]
chrom_dict = {}
total_density_dict = {}
total_lengths_dict = {}
for chrom in chrom_ls:
    prefix = os.path.basename(chrom).split('.')[0]
    chrom_dict[prefix] = pd.read_csv(chrom, sep = '\t', names = ['prefix', 'accession', 'chromosome', 'strand', 'start', 'end'])
    total_nb_genes = len(chrom_dict[prefix])
    total_length = sum(genomes_lengths_dict[prefix].values())
    total_density_dict[prefix] = total_nb_genes/total_length
    total_lengths_dict[prefix] = total_length

print(f'Done!\n', file = sys.stderr)

print('Processing multi_sp...', file = sys.stderr)
multi_sp_dict = {}
with open(args.multi_sp, 'r') as f:
    for line in f:
        multi_sp, *block_id_ls = line.rstrip().split('\t')
        multi_sp_dict.update({block_id:multi_sp for block_id in block_id_ls})

print(f'Done!\n', file = sys.stderr)

print('Processing orthology...', file = sys.stderr)
ortho_dict = load_ortho(args.ortho)
print(f'Done!\n', file = sys.stderr)


results_ls = collections.deque()
for node in nodes_obs_dict.keys():
    for species in nodes_obs_dict[node].keys():
        taxon = species_dict_r[species]
        chrom_df = chrom_dict[species]
        print('obs', node, species)
        for block in nodes_obs_dict[node][species]:
            tmp_df = chrom_df.query(('chromosome == @block.chromosome'))
            overlapping_genes_df = tmp_df.query('@block.start <= start < @block.end |@block.start < end <= @block.end| start <= @block.start < end|start < @block.end <= end')
            density = len(overlapping_genes_df) / (block.end - block.start)
            all_acc_ls = overlapping_genes_df.accession.tolist()
            all_acc_ls = ",".join(all_acc_ls)
            iteration = 1
            results_ls.append([node, taxon, species, 'observed', block.id, iteration, density, block.acc_ls, all_acc_ls])


for node in nodes_rand_dict.keys():
    for species in nodes_rand_dict[node].keys():
        taxon = species_dict_r[species]
        chrom_df = chrom_dict[species]
        print('rand', node, species)
        for block in nodes_rand_dict[node][species]:
            tmp_df = chrom_df.query(('chromosome == @block.chromosome'))
            overlapping_genes_df = tmp_df.query('@block.start <= start < @block.end |@block.start < end <= @block.end| start <= @block.start < end|start < @block.end <= end')
            density = len(overlapping_genes_df) / (block.end - block.start)
            all_acc_ls = overlapping_genes_df.accession.tolist()
            all_acc_ls = ",".join(all_acc_ls)
            block_id, iteration = block.id.split('.')
            results_ls.append([node, taxon, species, 'random', block_id, iteration, density, block.acc_ls, all_acc_ls])


results_df = pd.DataFrame(results_ls, columns =  ['node', 'taxon', 'species', 'random', 'block_id', 'iteration', 'density', 'acc_ls', 'all_acc_ls'])

#df for having rand and obs on the same plots
results_df['total_density'] = results_df['species'].map(total_density_dict)
results_df['total_genome_length'] = results_df['species'].map(total_lengths_dict)
results_df['density_ratio'] = results_df['density'] / results_df['total_density']
results_df['multi_sp'] = results_df['block_id'].map(multi_sp_dict)

#Assign para based on which obs the random block was sampled from
obs_df = results_df.query("random == 'observed'")
obs_df['para'] = obs_df['acc_ls'].map(lambda x: block_para(x))
para_dict = dict(zip(obs_df['block_id'], obs_df['para']))

results_df['para'] = results_df['block_id'].map(para_dict)
results_df['median_dist_pair'] = results_df['all_acc_ls'].map(lambda x: median_consecutivedist(x))
results_df['median_dist_pair_norm'] = results_df['median_dist_pair'] / results_df['total_genome_length']


results_df.to_csv(f'{args.output}.tidydf.csv', index = False)

