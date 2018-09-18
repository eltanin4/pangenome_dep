import requests
import pandas as pd
from bs4 import BeautifulSoup
from tqdm import tqdm
import urllib.request
import numpy as np
from uniqify import uniqify
from unlistify import unlistify
import os
import json
from load_kegg import *

# Getting the set of bacterial abbreviations for each organism in KEGG.
pangenome_df = pd.read_csv('fuller_pangenome_df.csv')
old_pangenome_df = pd.read_csv('pangenome_df.csv')
species = uniqify(list(pangenome_df['species'].values))
old_species = uniqify(list(old_pangenome_df['species'].values))
all_strains = uniqify(list(pangenome_df['strain'].values))
strain_abbrs = uniqify(list(pangenome_df['kegg_abbr'].values))
index_array = np.array(list(range(len(all_strains))))

# Getting all SEED organisms.
seed_orgs = pd.read_excel('seed_orgs.xlsx')
list_seed_names = list(seed_orgs['seed_name'].values)
seed_onames = [s[s.find("(")+1:s.find(")")] for s in list_seed_names]
smap = dict(zip(seed_onames, list_seed_names))

exact_in_seed_dict = {}
exact_species = {}
for s in list(set(seed_onames) & set(all_strains)):
    exact_in_seed_dict[smap[s]] = pangenome_df.loc[
                                  pangenome_df['strain'] == s, :]['kegg_abbr'].values[0]
    exact_species[smap[s]] = pangenome_df.loc[
                                  pangenome_df['kegg_abbr'] == exact_in_seed_dict[smap[s]], :]['species'].values[0]

# Getting all SEED autocompletion reactions.
seed_rxns = pd.read_excel('seed_rxns.xlsx')
list_seed_rxns = list(seed_rxns['seed_rxns'].values)

holder = 'https://www.patricbrc.org/api/model_reaction/?http_accept=application/json&eq(id,'
path = 'seed_rxn_files/'
# Getting all EC numbers for each species.
for rxn in tqdm(list_seed_rxns):
    if not os.path.isfile(path + rxn + '.txt'):
        # K01977 is the KO identifier for 16s rRNA.
        urllib.request.urlretrieve( holder + rxn + ')', path + rxn + '.txt' )



# Getting all abbrs for each reaction.
seed_to_kegg_rxn_map = {}
kegg_to_seed_rxn_map = {}
for rxn in tqdm(list_seed_rxns):
    s = open(path + rxn + '.txt', 'r').readlines()[0]
    if s != '[]':
        seed_to_kegg_rxn_map[rxn] = json.loads(s[1:-1])['abbreviation']
        kegg_to_seed_rxn_map[json.loads(s[1:-1])['abbreviation']] = rxn

# Finding what's in KEGG.
common_rxns = set(list(rxn_map.keys())) & set(list(kegg_to_seed_rxn_map.keys()))
filt_kegg_to_seed_rxn_map = {k:kegg_to_seed_rxn_map[k] for k in common_rxns}
filt_seed_to_kegg_rxn_map = {value: key for key, value in filt_kegg_to_seed_rxn_map.items()}

# Getting gap-filling information.
full_gap_df = pd.read_excel('seed_gaps_filled.xlsx').fillna('present')
study_gap_fill = []
# r = list(filt_seed_to_kegg_rxn_map.keys())[0]
for c in tqdm(list(exact_in_seed_dict.keys())):
    for r in tqdm(list(filt_seed_to_kegg_rxn_map.keys())):
        thisval = full_gap_df.loc[full_gap_df['seed_rxns'] == r, c].values[0]
        if thisval == 'missing':
            study_gap_fill.append((exact_species[c], exact_in_seed_dict[c], filt_seed_to_kegg_rxn_map[r], thisval))

study_gap_fill_df = pd.DataFrame(study_gap_fill,
                                columns=['species', 'kegg_abbr', 'kegg_rxn', 'is_missing'])
study_gap_fill_df.to_csv('all_gap_fills.csv', index=None)

study_gap_fill_df = pd.read_csv('all_gap_fills.csv',
                                columns=['species', 'kegg_abbr', 'kegg_rxn', 'is_missing'])
strains_affected = set(study_gap_fill_df['kegg_abbr'].tolist())
species_affected = []
for this_species in tqdm(species):
    these_strains = list(pangenome_df.loc[pangenome_df['species'] == this_species]['kegg_abbr'])
    for this_strain in these_strains:
        if this_strain in strains_affected:
            species_affected.append(this_species)
