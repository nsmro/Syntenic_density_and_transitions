#!/usr/bin/env python3

import argparse
import pandas as pd
from goatools.obo_parser import GODag
from goatools.anno.idtogos_reader import IdToGosReader
from goatools.goea.go_enrichment_ns import GOEnrichmentStudyNS
from goatools.godag_plot import plot_gos, plot_results, plot_goid2goobj


parser = argparse.ArgumentParser(description = """Plots goterms enriched GO terms. Pvalue of 0,01, BH procedure correction for multiple testing.
    Will test all 3 namespaces (Biological Process, Molecular Function, Cellular Compartment""")
parser.add_argument('-i', '--ids2go', help = 'GO annotation, ids2go. Two fields, first one is accessions, second one is semicolon-delimited lis of GO annotations. ', required = True)
parser.add_argument('-go', '--obo', help = 'obo file to use in the GO enrichment analysis', required = True)
parser.add_argument('-s', '--study', help = 'list of sequences from the study, newline-delimited', required = True)
parser.add_argument('-b', '--background', help = 'List of the background sequences, i.e. all the sequences from the animal, newline-delimited', required = True)
parser.add_argument('-p', '--propagate_counts', help = 'propagate GO_counts, default: False', action = 'store_true')
parser.add_argument('-o', '--output', help = 'Name of the output filename prefix')
args = parser.parse_args()


def file2list(sequences_list):
    acc_ls = []
    with open(sequences_list, 'r') as f:
        for line in f:
            acc = line.rstrip()
            acc_ls.append(acc)
    return acc_ls

if args.output is None:
    output_fname = f'{args.study.split(".")[0]}_GOAE'
else:
    output_fname = args.output

#Import the GO hierarchy, downloaded from http://geneontology.org/ontology/go-basic.obo
obodag = GODag(args.obo)

# Read ids2go format. Store annotations in a list of named tuples. Specify obo graph to use
objanno = IdToGosReader(args.ids2go, godag = obodag)

"""
Get namespace2association. Basically, a nested dict as follows:
    ns2assoc[namespace][association]: GO
    namespace is: BP (biological_process), MF (molecular_function), CC (cellular_component)
    assocation is: a protein_id
    GO is: the set of GO IDs associated with that protein_ID
"""
ns2assoc = objanno.get_ns2assc()

background = file2list(args.background)

goeaobj = GOEnrichmentStudyNS(
        background, # List of background proteins
        ns2assoc, # geneid/GO associations
        obodag, # Ontologies
        propagate_counts = args.propagate_counts,
        alpha = 0.05, # we'll change the pvalue cut-off to 0.01
        methods = ['fdr_bh']) # defult multipletest correction method

study_list = file2list(args.study)

# 'p_' means "pvalue". 'fdr_bh' is the multiple test correction, Benjamini-Hochberg.
goea_results_all = goeaobj.run_study(study_list)
goea_results_sig = [r for r in goea_results_all if r.p_fdr_bh < 0.05]

goeaobj.wr_xlsx(f'{output_fname}.xlsx', goea_results_sig)
goeaobj.wr_txt(f'{output_fname}.txt', goea_results_sig)

#Isolate only the significant GOs we use the attribute GO of the GOEnrichmentRecord object located in the significant goea_results
GO_list = [GOEnrichmentRecord.GO for GOEnrichmentRecord in goea_results_sig]

#This'll make one single plot with the 3 GO namespaces.
plot_gos(f'{output_fname}.pdf',
    GO_list, # Source GO ids
    obodag,
    goea_results = goea_results_sig) # Use pvals for coloring
