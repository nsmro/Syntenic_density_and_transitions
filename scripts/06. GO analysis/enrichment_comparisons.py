#!/usr/bin/env python3


import os
import pandas as pd
import glob

taxo_dict = {'ACRMI':'Cnidarian','ADIVA':'Lophotrochozoan','AMPQU':'Poriferan','ACAPL':'InvDeut','ANOGA':'Ecdysozoan','AURAU':'Cnidarian','BRALA':'InvDeut','CAEEL':'Ecdysozoan','CALMI':'Vertebrate','CAPTE':'Lophotrochozoan','CAPOW':'Metazoa_outgroup','CHEMY':'Vertebrate','CIOIN':'InvDeut','CLYHE':'Cnidarian','CRAGI':'Lophotrochozoan','DANRE':'Vertebrate','DAPPU':'Ecdysozoan','DROME':'Ecdysozoan','EUPSC':'Lophotrochozoan','EXAPA':'Cnidarian','GALGA':'Vertebrate','HELRO':'Lophotrochozoan','HIPCO':'Vertebrate','HOFMI':'Acoel','HOIHO':'Placozoan','HOMSA':'Vertebrate','HYDVU':'Cnidarian','IXOSC':'Ecdysozoan','LATCH':'Vertebrate','LEPOC':'Vertebrate','LINAN':'Lophotrochozoan','LOTGI':'Lophotrochozoan','MAYZE':'Vertebrate','MIZYE':'Lophotrochozoan','MNELE':'Ctenophore','MUSMU':'Vertebrate','NEMVE':'Cnidarian','PARTE':'Ecdysozoan','PLEBA':'Ctenophore','PTYFL':'InvDeut','SACKO':'InvDeut','SALRO':'Metazoa_outgroup','SCHME':'Lophotrochozoan','STRMA':'Ecdysozoan','STRPU':'InvDeut','SYCCI':'Poriferan','TRICA':'Ecdysozoan','TRIAD':'Placozoan','XENTR':'Vertebrate'}
filelist_meta = glob.glob('Metazoa/*.xlsx')
filelist_planu = glob.glob('Planulozoa/*.xlsx')
filelist_bila = glob.glob('Bilateria/*.xlsx')
taxonomy_ls = ['Metazoa_outgroup', 'Ctenophore', 'Poriferan', 'Placozoan', 'Cnidarian', 'Acoel', 'Ecdysozoan', 'Lophotrochozoan', 'InvDeut', 'Vertebrate']

def make_dicts(filelist):
    GO_dict = {}
    GO_strings = {}
    for file in filelist:
        df = pd.read_excel(file)
        species = os.path.basename(file).split('_')[0]
        df_enriched = df.query("enrichment =='e'")
        GO_ls = list(df_enriched.get('GO'))
        for GO in GO_ls:
            if GO in GO_dict.keys():
                GO_dict[GO].append(species)
            else:
                GO_dict[GO] = [species]
                GO_strings[GO] = df_enriched.query('GO == @GO').iloc[0]['name']
    return GO_dict, GO_strings

#Metazoa terms
GO_dict_meta, GO_strings_meta = make_dicts(filelist_meta)

basal_metazoan = ['Ctenophore','Poriferan']
other_metazoans = ['Placozoan', 'Cnidarian', 'Acoel', 'Ecdysozoan', 'Lophotrochozoan', 'InvDeut', 'Vertebrate']

GO_ls_bytaxo = []
for GO,species_ls in GO_dict_meta.items():
	taxo_ls = [taxo_dict[species] for species in species_ls]
	nb_basal_in_taxon = len([species for species in taxo_ls if species in basal_metazoan])
	nb_other_in_taxon = len([species for species in taxo_ls if species in other_metazoans])
	total_nb_species_with_GO = len(species_ls)
	if nb_basal_in_taxon > 0 and nb_other_in_taxon > 0 and total_nb_species_with_GO >= 8:
		GO_ls_bytaxo.append([GO, GO_strings_meta[GO], ','.join(species_ls)])

with open('Metazoa_report_GO_enriched.tsv', 'w') as f:
	f.write(pd.DataFrame(GO_ls_bytaxo, columns = ['GO', 'name', 'species list']).to_csv(sep = '\t', index = False))

#Planulozoa novel terms
GO_dict_planu, GO_strings_planu = make_dicts(filelist_planu)

Cnidaria = ['Cnidarian']
Bilateria = ['Acoel', 'Ecdysozoan', 'Lophotrochozoan', 'InvDeut', 'Vertebrate']
GO_ls_bytaxo = []
for GO,species_ls in GO_dict_planu.items():
        taxo_ls = [taxo_dict[species] for species in species_ls]
        nb_Cnidaria_in_taxon = len([species for species in taxo_ls if species in Cnidaria])
        nb_Bilateria_in_taxon = len([species for species in taxo_ls if species in Bilateria])
        total_nb_species_with_GO = len(species_ls)
        if nb_Cnidaria_in_taxon >= 3 and nb_Bilateria_in_taxon >= 8:
                GO_ls_bytaxo.append([GO, GO_strings_planu[GO], ','.join(species_ls)])

GO_ls_bytaxo = list(sorted(GO_ls_bytaxo, key = lambda x:x[2]))

with open('Planu_report_GO_enriched.tsv', 'w') as f:
        f.write(pd.DataFrame(GO_ls_bytaxo, columns = ['GO', 'name', 'species list']).to_csv(sep = '\t', index = False))

#Bilateria novel terms
GO_dict_planu, GO_strings_planu = make_dicts(filelist_bila)

protostomia = ['Ecdysozoan', 'Lophotrochozoan']
deuterostomia = ['InvDeut', 'Vertebrate']
GO_ls_bytaxo = []
for GO,species_ls in GO_dict_planu.items():
	taxo_ls = [taxo_dict[species] for species in species_ls]
	nb_protostomia_in_taxon = len([species for species in taxo_ls if species in protostomia])
	nb_deuterostomia_in_taxon = len([species for species in taxo_ls if species in deuterostomia])
	total_nb_species_with_GO = len(species_ls)
	if nb_protostomia_in_taxon >= 4 and nb_deuterostomia_in_taxon >= 4:
		GO_ls_bytaxo.append([GO, GO_strings_planu[GO], ','.join(species_ls)])

GO_ls_bytaxo = list(sorted(GO_ls_bytaxo, key = lambda x:x[2]))

with open('Bila_report_GO_enriched.tsv', 'w') as f:
	f.write(pd.DataFrame(GO_ls_bytaxo, columns = ['GO', 'name', 'species list']).to_csv(sep = '\t', index = False))
