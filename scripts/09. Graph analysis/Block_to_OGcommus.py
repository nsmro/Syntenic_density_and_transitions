#! /usr/bin/env python3

import sys
import collections
import itertools
import glob
import argparse
import statistics

parser = argparse.ArgumentParser(description = 'Provides reports of block content of specified noppdes')
parser.add_argument("-c", "--clusters_id", help = "Tsv file, the multi-species clusters of microsyntenic blocks, output of the microsynteny pipeline (makeClusters3.pl).", required = True)
parser.add_argument("-b", "--block_list", help = "Tsv file, microsyntenic blocks details, output of the microsynteny pipeline (makeClusters3.pl).", required = True)
parser.add_argument("-r", "--random_blocks", help = "Tsv file, microsyntenic blocks details, output of Bob's script, corrected to fit the same format as the observed blocks", required = True)
parser.add_argument("-g", "--chromfiles", help = "List of the chromfiles to use", type = str, required = True, nargs = "+")
parser.add_argument("-og", "--orthology", help = "Orthology file, clus format", type = str, required = True)
parser.add_argument("-o", "--output", help = "Prefix of the output files.", type = str, default = "output")
parser.add_argument("--custom_orthology", help = "use flag if custom_orthology.", action = "store_true")
args = parser.parse_args()


def nested_dict():
    """ Create a nested dict. Function is picklable """
    return collections.defaultdict(nested_dict)


def read_chromfiles(chromfiles_list, orthology_dict):
    """
    Reads a list of files in the .chrom format
    Returns a nested dict as follows:
    [OG] [Species] [chromosome/scaffold] [accession] [coordinates]
    Coordinates is a tuple of len 3:
        (position on the chrom (1-based) , start, stop)
    Second dict will just contain accessions
    if args.custim_orthology is true, all the values will be stored in the chrom
    if it's false, orthology will be used to filter the coordinates
    """
    output_dict = nested_dict()
    accession_pos_dict = {}
    species_ls = [name.split('/')[-1].split('.')[0] for name in chromfiles_list]
    for file in chromfiles_list:
        species = file.split('/')[-1].split('.')[0]
        sys.stderr.write(f'Reading file {file}\n')
        tmp_coords = nested_dict()
        with open(file, 'r') as chrom_file:
            for line in chrom_file:
                _, accession, scaffold, _, start, stop = line.rstrip(None).split('\t')
                if args.custom_orthology is True:
                    tmp_coords[scaffold][accession] = [start,stop]
                elif accession in orthology_dict.keys():
                    tmp_coords[scaffold][accession] = [start,stop]
        for chromosome in tmp_coords.keys():
            chrsort = sorted((tmp_coords[chromosome]).keys(),key = lambda x: int((tmp_coords[chromosome])[x][0]))
            pos = 0
            for gene in chrsort:
                pos = pos + 1
                if args.custom_orthology is False:
                    output_dict[orthology_dict[gene]][species][chromosome][gene] = [pos] + tmp_coords[chromosome][gene]
                else:
                    try:
                        output_dict[orthology_dict[gene]][species][chromosome][gene] = [pos] + tmp_coords[chromosome][gene]
                    except KeyError:
                        output_dict['other'][species][chromosome][gene] = [pos] + tmp_coords[chromosome][gene]
                accession_pos_dict[gene] = [pos] + tmp_coords[chromosome][gene]
    sys.stderr.write('Chroms loaded!\n')
    return output_dict, species_ls, accession_pos_dict


def block_id_dict(handle, og_dictionary):
    """
    returns a dict with block id as keys,
    named tuples as values: species, chromosome, acc_ls, og_pairs
    og_pairs are the OG_pairs found in a given block.
    another dict is the reverse dict (acc as keys, block as values)
    aother dict is the index of the accession within the block (used for determinining random dists)
    """
    output_dict = {}
    output_dict_acc = {}
    output_dict_acc_index = {}
    block = collections.namedtuple('block', 'species accessions chromosome og_pairs')
    for line in handle:
        line = line.rstrip().split()
        block_id = line[0]
        species = line[1]
        chromosome = line[7].split(':')[0]
        acc_ls = line[9].split(',')
        og_pairs = set([frozenset([og_dictionary[x],og_dictionary[y]]) for x,y in itertools.combinations(acc_ls, 2)])
        output_dict[block_id] = block(species, acc_ls, chromosome, og_pairs)
        for x in  range(0,len(acc_ls)):
            acc = acc_ls[x]
            output_dict_acc[acc] = block_id
            output_dict_acc_index[acc] = x
    return output_dict, output_dict_acc, output_dict_acc_index

def rand_blocks_positions(handle, accession_pos_dict):
    output_dict = {}
    for line in handle:
        line = line.rstrip().split()
        block_id, iteration_number = line[0].split('.')
        acc_ls = line[9].split(',')  
        if iteration_number == "1":
            for x in range(0, len(acc_ls)):
                pos, start, end = [int(x) for x in accession_pos_dict[acc_ls[x]]]
                output_dict[frozenset({block_id, x})] = [[pos],[start],[end]]
        else:
            for x in range(0, len(acc_ls)):
                pos, start, end = [int(x) for x in accession_pos_dict[acc_ls[x]]]
                output_dict[frozenset({block_id, x})][0].append(pos)
                output_dict[frozenset({block_id, x})][1].append(start)
                output_dict[frozenset({block_id, x})][2].append(end)
    return output_dict

#nested lists frozenset{block_string, index} [[pos,pos][starts, start][end,end]]

def og_dict(handle):
    """
    dict with OGs or accessions as keys
    accession return OGs,
    OGs returns acc_ls
    """
    output_dict = {}
    for line in handle:
        line = line.rstrip().split()
        OG = line[0]
        for acc in line[2:]:
            output_dict[acc] = OG
        output_dict[OG] = line[2:]
    return output_dict

def clus_id_dict(handle):
    """
    dict
    """
    clus_id_dict = {}
    for line in handle:
        line = line.rstrip().split()
        clus_id_dict[line[0]] = line[1:]
    return clus_id_dict




with open(args.orthology, 'r') as filehandle:
    OG_d = og_dict(filehandle)

with open(args.clusters_id, 'r') as filehandle:
    clus_d = clus_id_dict(filehandle)

with open(args.block_list) as filehandle:
    block_d, acc_d, acc_index_d = block_id_dict(filehandle, OG_d)

chrom_dict, total_sp_ls, accession_pos_coords_d = read_chromfiles(args.chromfiles,OG_d)

class SpeciesDist(object):
    """
    this is just to have a fixed length list of items created in each dictionary entry in results_dict
    We'll populate it with NA at the start
    And if the OG pair is found in the cluster, we'll change this value to the actual distance
    """
    __slots__ = total_sp_ls
    def __init__(self):# we'll declare values as NAs in all attributes after creating the instance
        pass

with open(args.random_blocks) as filehandle:
    rand_blocks_positions_dict = rand_blocks_positions(filehandle, accession_pos_coords_d)

dist_dict, bp_dict, dist_dict_rand, bp_dict_rand, OG_commu_dict = {}, {}, {}, {}, {}

for multi_sp_block in clus_d.keys():
    species_ls = [block_d[block].species for block in clus_d[multi_sp_block]]
    acc_ls_total = itertools.chain.from_iterable([block_d[block].accessions for block in clus_d[multi_sp_block]])
    OG_commu = list(set([OG_d[acc] for acc in acc_ls_total]))
    OG_commu_dict[multi_sp_block] = OG_commu
    for species in species_ls:
        og_pair_set_sp = set(itertools.chain.from_iterable([block_d[block].og_pairs for block in clus_d[multi_sp_block] if block_d[block].species == species]))
        acc_ls_set_sp = list(itertools.chain.from_iterable([block_d[block].accessions for block in clus_d[multi_sp_block] if block_d[block].species == species]))
        for pair_fzset in og_pair_set_sp:
            if len(pair_fzset) == 1: #OG to self
                og1,og2 = list(pair_fzset)*2
            else:
                og1,og2 = pair_fzset
            key_string_ogpair = f'{multi_sp_block}_{og1}_{og2}'
            og1_filt = [acc for acc in OG_d[og1] if acc in acc_ls_set_sp]
            og2_filt = [acc for acc in OG_d[og2] if acc in acc_ls_set_sp]
            acc_pairs = [[x,y] for x in og1_filt for y in og2_filt if x !=y and acc_d[x] == acc_d[y]]
            for mydict in dist_dict, bp_dict, dist_dict_rand, bp_dict_rand:
                try:
                    mydict[key_string_ogpair]
                except KeyError: #if key not in the dict, we initialize it
                    mydict[key_string_ogpair] = SpeciesDist()
                    for item in total_sp_ls: #initialize Species_dist object with NAs, as we iterate only through the existing OG pairs in species
                        setattr(mydict[key_string_ogpair], item, 'NA')
            for gene1,gene2 in acc_pairs:
                block_id_pair = acc_d[gene1]
                chromosome = block_d[block_id_pair].chromosome
                pos1, start1, end1 = [int(x) for x in chrom_dict[og1][species][chromosome][gene1]]
                pos2, start2, end2 = [int(x) for x in chrom_dict[og2][species][chromosome][gene2]]
                overlap_genes = min(start2,end2) <= max(start1, end1) and max(start2,end2) >= min(start1, end1)
                if overlap_genes is True:
                    dist_bp = 0
                else:
                    dist_bp = min(abs(end1 - start2), abs(end1 - end2), abs(start1 - start2), abs(start1 - end2))
                dist_genes = abs(pos2 - pos1)
                gene1_index, gene2_index = acc_index_d[gene1],acc_index_d[gene2]
                pos1_r_ls, start1_r_ls, end1_r_ls = rand_blocks_positions_dict[frozenset({block_id_pair,gene1_index})]
                pos2_r_ls, start2_r_ls, end2_r_ls = rand_blocks_positions_dict[frozenset({block_id_pair,gene2_index})]
                dist_gene_rand = statistics.median([abs(pos2_r - pos1_r) for pos1_r, pos2_r in zip(pos1_r_ls, pos2_r_ls)])
                dist_bp_rand_ls = []
                for start1, end1, start2, end2 in zip(start1_r_ls, end1_r_ls, start2_r_ls, end2_r_ls):
                    overlap_genes = min(start2,end2) <= max(start1, end1) and max(start2,end2) >= min(start1, end1)
                    if overlap_genes is True:
                        dist_bp_rand_ls.append(0)
                    else:
                        dist_bp_rand_ls.append(min(abs(end1 - start2), abs(end1 - end2), abs(start1 - start2), abs(start1 - end2)))
                dist_bp_rand = statistics.median(dist_bp_rand_ls)
                if getattr(dist_dict[key_string_ogpair], species) == 'NA' or getattr(dist_dict[key_string_ogpair], species) > dist_genes:
                    setattr(dist_dict[key_string_ogpair], species, dist_genes)
                    setattr(dist_dict_rand[key_string_ogpair], species, dist_gene_rand)
                if getattr(bp_dict[key_string_ogpair], species) == 'NA' or getattr(bp_dict[key_string_ogpair], species) > dist_bp:
                    setattr(bp_dict[key_string_ogpair], species, dist_bp)
                    setattr(bp_dict_rand[key_string_ogpair], species, dist_bp_rand)   


with open(f'{args.output}.OG_commus', 'w') as f:
    for cluster_id, og_commu in OG_commu_dict.items():
        f.write(cluster_id + '\t' + '\t'.join(og_commu) + '\n')

with open(f'{args.output}.dist', 'w') as f:
    f.write('OG1\tOG2\t' + '\t'.join(total_sp_ls) + '\n')
    for key, value in dist_dict.items():
        multi_sp_block, OG1, OG2 = key.split('_')
        dist_ls_string = '\t'.join([str(getattr(dist_dict[key], species)) for species in total_sp_ls])
        f.write(f'{multi_sp_block}\t{OG1}\t{OG2}\t{dist_ls_string}\n')

with open(f'{args.output}.bp', 'w') as f:
    f.write('OG1\tOG2\t' + '\t'.join(total_sp_ls) + '\n')
    for key, value in bp_dict.items():
        multi_sp_block, OG1, OG2 = key.split('_')
        dist_ls_string = '\t'.join([str(getattr(bp_dict[key], species)) for species in total_sp_ls])
        f.write(f'{multi_sp_block}\t{OG1}\t{OG2}\t{dist_ls_string}\n')

with open(f'{args.output}.rand.dist', 'w') as f:
    f.write('OG1\tOG2\t' + '\t'.join(total_sp_ls) + '\n')
    for key, value in dist_dict_rand.items():
        multi_sp_block, OG1, OG2 = key.split('_')
        dist_ls_string = '\t'.join([str(getattr(dist_dict_rand[key], species)) for species in total_sp_ls])
        f.write(f'{multi_sp_block}\t{OG1}\t{OG2}\t{dist_ls_string}\n')

with open(f'{args.output}.rand.bp', 'w') as f:
    f.write('OG1\tOG2\t' + '\t'.join(total_sp_ls) + '\n')
    for key, value in bp_dict_rand.items():
        multi_sp_block, OG1, OG2 = key.split('_')
        dist_ls_string = '\t'.join([str(getattr(bp_dict_rand[key], species)) for species in total_sp_ls])
        f.write(f'{multi_sp_block}\t{OG1}\t{OG2}\t{dist_ls_string}\n')
