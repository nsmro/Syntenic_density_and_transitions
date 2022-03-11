import argparse
import collections
import glob
import itertools
import numpy as np
import os
import pandas as pd
import sys

parser = argparse.ArgumentParser(description = """
Each block is decomposed into N bins (bin size = 1/N block width)
the upstream and downstream 2N bins are also included
Gene density is calculated for 5N bins
Outputs tidy df where each line is a bin
block_id, multi_sp, species obs and bin location are outputted
""")
parser.add_argument('-s', '--synt_file', help = "synt file where to pick the blocks from", type = str, required = True)
parser.add_argument('-r', '--rand_synt_file', help = "randomly sampled blocks", type = str, required = True)
parser.add_argument('-m', '--multi_sp', help = 'Multi species block (total), used to get multi_sp ID', type = str, required = True)
parser.add_argument('-g', '--genome_folder', help = 'folder where genome files are located, prefix in the name dot separated. e.g. EUPSC.lachesis201904.genome.fa, CRAGI.NCBIgenome.fasta',type = str, required = True)
parser.add_argument('-c', '--chrom_folder', help = 'folder where chrom files are located, files should have the *chrom extension', required = True)
parser.add_argument('-n', '--n_bins', help = 'number of bins by block', default = 5, type = int)
parser.add_argument('-o', '--output', help = 'output name')
args = parser.parse_args()

#tmp_args = collections.namedtuple('tmp_args', 'synt_file rand_synt_file multi_sp genome_folder chrom_folder n_bins output')
#args = tmp_args("/scratch/robert/2019_microsynteny_size_constraints/03_REDUX_density_HOX_WNT/wnt/wnt.syn.synt",
#                "/scratch/robert/2019_microsynteny_size_constraints/03_REDUX_density_HOX_WNT/wnt/wnt.syn.random.synt",
#                "/scratch/robert/2019_microsynteny_size_constraints/03_REDUX_density_HOX_WNT/wnt/wnt.syn.clusters",
#                "/scratch/robert/2019_microsynteny_size_constraints/02_REDUX_gene_density_analysis/genomes/",
#                "/scratch/robert/2019_microsynteny_size_constraints/01_microsynteny/chrom/",
#                5,
#                "test")

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


def get_ranges(coords_string, nb_bins, species, lengths_dict):
    """
    Return pairs of values start end of bins, given a block, and a nb of bins
    :param coords_string: coords as written in synt file (scaffold:start...end)
    :param nb_bins: nb of bins the block should be split into
    :param species: prefix of species where the block is found
    :param lengths_dict: nested length dict where k,k,v are prefix, scaffold, length
    :returns: a list of lists of length 4,
    start end of the bins that can be sampled (won't sample outside of block boundaries),
    and the normalized coords so that 0 = block start, 1 = block end 
    """
    output = []
    bin_factors = [round(x, 1) for x in np.linspace(-2, 3, num = nb_bins * 5 + 1)]
    bin_factors = zip(bin_factors, bin_factors[1:])
    chromosome, start, end = coords_string.replace('..', ':').split(':')
    start = int(start)
    end = int(end)
    len_block = end - start
    bin_size = round(len_block / nb_bins)
    bin_lower = start - 2 * len_block
    bin_upper = end + 2 * len_block
    len_chrom = lengths_dict[species][chromosome]
    bin_breaks = [i for i in range(bin_lower, bin_upper, bin_size)] + [bin_upper]
    bin_breaks = zip(bin_breaks, bin_breaks[1:])
    bin_bounds = [itertools.chain.from_iterable((n, bp)) for n, bp in zip(bin_factors, bin_breaks) if (any([x < 0 for x in bp]) is False and any([x > len_chrom for x in bp]) is False)]
    return bin_bounds


print('Loading genome scaffold lengths...', file = sys.stderr)
genomes_ls = [file for file in glob.glob(f'{args.genome_folder}/*') if 'genome' in file]
genomes_lengths_dict = {}
for genome_file in genomes_ls:
    prefix = os.path.basename(genome_file).split('.')[0]
    print(prefix)
    with open(genome_file, 'r') as f:
        lengths_d = {header: len(sequence) for header, sequence in parse_fasta(f)}
        genomes_lengths_dict[prefix] = lengths_d

print('Done!\n', file = sys.stderr)

print('Loading chromfiles...', file = sys.stderr)

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

print('Loading multi_sp...', file = sys.stderr)
multi_sp_dict = {}
with open(args.multi_sp, 'r') as f:
    for line in f:
        multi_sp, *block_id_ls = line.rstrip().split('\t')
        multi_sp_dict.update({block_id:multi_sp for block_id in block_id_ls})

print(f'Done!\n', file = sys.stderr)

#bin_loc is a factor, N bins of normalized bp size 1/n
header = ['block_id', 'iteration', 'multi_sp', 'species', 'obs', 'bin_location', 'breaks_location', 'density']
results = []
for file, observed in [[args.synt_file, 'observed'], [args.rand_synt_file, 'random']]:
    with open(file, 'r') as f:
        for line in f:
            block_id, species, _, _, _, _, _, coords, *_ = line = line.rstrip().split('\t')
            print(block_id)
            block_chromosome = coords.split(':')[0]
            iteration = 0
            multi_sp = multi_sp_dict.get(block_id)
            if multi_sp == None:
                block_id, iteration = block_id.split('.')
                multi_sp = multi_sp_dict.get(block_id)
            chrom_df = chrom_dict[species]
            tmp_df = chrom_df.query(('chromosome == @block_chromosome'))
            for start_factor, end_factor, start_breaks, end_breaks in \
                get_ranges(coords, args.n_bins, species, genomes_lengths_dict):
                factor_location = round((start_factor + end_factor) / 2, 1)
                break_location = f'[{start_breaks}:{end_breaks})'
                end_breaks = end_breaks -1
                overlapping_genes_df = tmp_df.query('@start_breaks <= start < @end_breaks |@start_breaks < end <= @end_breaks| start <= @start_breaks < end|start < @end_breaks <= end')            
                density = len(overlapping_genes_df) / (end_breaks - start_breaks)
                results.append([block_id,
                                iteration,
                                multi_sp,
                                species,
                                observed,
                                factor_location,
                                break_location,
                                density])


df = pd.DataFrame(results, columns = header)

df.to_csv(args.output, sep = '\t', index = False)

