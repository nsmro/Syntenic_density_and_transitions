#! /usr/bin/env python3

import re
import sys
import pandas as pd
import glob

"""
Usage:
MUSMU_reformat.py protlist expression_data MUSMU_gff
"""

if len(4 > sys.argv > 1):
    _, protein_file, gff_file = sys.argv
else:
    raise ArgumentError("""
        Usage:
        MUSMU_reformat.py protlist expression_data MUSMU_gff
        protlist: list of neyline separated protein accessions
        expression data: mouse encode file to use
        gff file:a gff file
        """)


#load proteins
protein_ls = {}
with open(protein_file, 'r') as f:
    for line in f:
        line = line.rstrip()
        prot = line.split('_', 1)[1]
        protein_ls[prot] = ''

#load gff
accessions_dict = {} 
re_transcript = 'Parent'+'\s*[:="\']*([^,|;"\']+)[;|"\s\n]*'
re_protein = 'protein_id'+'\s*[:="\']*([^,|;"\']+)[;|"\s\n]*'
with open(gff_file) as f:
    for line in f:
        if (line.startswith('#') or len(line.split('\t')) < 9):
            pass
        else:
            fields = line.rstrip().split('\t')
            _, _, feature, *_,comments = fields
            if feature == 'CDS':
                if 'exception' not in comments:
                    protein_id = re.search(re_protein, comments).group(1)
                    if protein_ls.get(protein_id) != None:
                        transcript_id = re.search(re_transcript, comments).group(1)
                        transcript_id = transcript_id.split('-')[1].split('.')[0]
                        accessions_dict[transcript_id] = protein_id
                        del protein_ls[protein_id]

results =  {}
for file in glob.glob('19-tissues-expr/*.expr'):
    basename = os.path.basename(file)
    with open(file, 'r') as f:
        _ = f.readline()
        for line in f:
            line = line.rstrip().split()
            gene_id , *_, FPKM, _, _ = line
            protein_id = accessions_dict.get(gene_id)
            if protein_id != None:
                if results.get(basename) == None:
                    results[basename] = {protein_id:FPKM}
                else:
                    results[basename][protein_id] = FPKM

df = pd.DataFrame(results)
df = df.fillna(0.0)


stages = {'bone_marrow' : ['boneMarrow1-zy24.gene.expr', 'boneMarrow2-zy26.gene.expr'],
'E14.5_brain': ['brain-E14.5-1.expr', 'brain-E14.5-2.expr'],
'cerebellum': ['cerebellum1-zy21.gene.expr', 'cerebellum2-zy22.gene.expr'],
'cortex': ['cortex1-zy13.gene.expr', 'cortex2-zy14.gene.expr'],
'E14.5_heart': ['heart-E14.5-1.expr', 'heart-E14.5-2.expr'],
'heart':['heart1-zy6.gene.expr', 'heart2-zy7.gene.expr'],
'intestine':['intestine-2.expr', 'intestine-3.expr'],
'kidney':['kidney1-zy15.gene.expr', 'kidney2-zy16.gene.expr'],
'E14.5_limb': ['limb-E14.5-1.expr', 'limb-E14.5-2.expr'],
'E14.5-liver':['liver-E14.5-1.expr', 'liver-E14.5-2.expr'],
'liver':['liver1-zy4.gene.expr', 'liver2-zy5.gene.expr'],
'lung':['lung1-zy10.gene.expr', 'lung2-zy11.gene.expr'],
'mESC':['mESC-zy27.gene.expr', 'mESC-zy28.gene.expr'],
'MEF':['mef-male1-zy17.gene.expr', 'mef-male2-zy18.gene.expr'],
'olfactory_bulb':['olfactory-1.expr', 'olfactory-2.expr'],
'placenta':['placenta-1.expr', 'placenta-2.expr'],
'spleen':['spleen1-zy8.gene.expr', 'spleen2-zy9.gene.expr'],
'testes':['testes-1.expr', 'testes-2.expr'],
'thymus':['thymus-1.expr', 'thymus-2.expr']}

df_med = pd.DataFrame()

for stage in stages.keys():
    df_med[stage] = df[stages[stage]].median(axis = 1)