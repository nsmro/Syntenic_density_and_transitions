#! /usr/bin/env python3

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.stats
import itertools



#plot preferences, move onto other file
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


results_df = pd.read_csv('key_nodes.tidydf.csv')

results_df_fig2 = results_df.groupby(['multi_sp', 'random', 'taxon', 'node']).median() # median of all the values
results_df_fig2.reset_index(inplace = True)
del results_df_fig2['block_id']
del results_df_fig2['iteration']
del results_df_fig2['density']
del results_df_fig2['total_density']


df = results_df_fig2.copy()
df = df.pivot_table(index=['multi_sp', 'taxon', 'node'],
                    columns='random', 
                    values=['density_ratio'])
df.columns = ['obs_density_ratio', 'rand_density_ratio']
df.reset_index(inplace = True)
df['deviation_to_random'] = (df['obs_density_ratio'] - df['rand_density_ratio']) / df['obs_density_ratio']


df2 = df.copy()
df2['order'] = df2['taxon'].map(lambda x: order.index(x))
df2 = df2.pivot_table(index=['multi_sp', 'node'],
                      columns='taxon', 
                      values=['deviation_to_random'])
df2.columns = [x[1] for x in df2.columns]
df2.reset_index(inplace = True)

order_taxons = [tax for tax in order if tax in df.taxon.to_list()]
order_columns = [x for x in df2.columns if x not in order_taxons] + order_taxons

df2 = df2[order_columns]

df2.to_csv('raw_data_scatter.csv', index = False)