import numpy as np
import pickle
import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm
from uniqify import uniqify
from unlistify import unlistify
from load_kegg import *
from give_scope import giveScope
import pandas as pd
import itertools
from sampler import collect_samples
from copy import deepcopy
import random

def propagate_single_for_medium(org, medium):
    # Defining the seed set to be the medium and the currency,
    seedVec = np.zeros(len(rxnMat.T))
    seedVec[[kegg_to_id[e] for e in Currency + medium]] = 1

    # Getting all the reactions performable by this organism.
    can = ''.join(open('strain_reactions/' + org + '.txt', 'r').readlines()).split()
    can = list((set(can)) & set(rxns))
    avrxns = [rxn_kegg_to_id[e] for e in can]

    # Calculating the metabolites within the scope of 
    # this organism's reaction network. 
    scopeMets = giveScope(rxnMat[avrxns], prodMat[avrxns], seedVec, sumRxnVec[avrxns])[0]
    
    # Finding how much of the core is within the network's scope.
    return uniqify([id_to_kegg[e] for e in np.where(scopeMets)[0]])


# Getting the set of bacterial abbreviations for each organism in KEGG.
pangenome_df = pd.read_csv('fuller_pangenome_df.csv')
old_pangenome_df = pd.read_csv('pangenome_df.csv')
species = uniqify(list(pangenome_df['species'].values))
old_species = uniqify(list(old_pangenome_df['species'].values))
all_strains = uniqify(list(pangenome_df['kegg_abbr'].values))
index_array = np.array(list(range(len(all_strains))))


# Generating a vector of all metabolites initially provided, i.e. seeds.
seeds_df = pd.read_csv('../black_queen_critique/seeds_from_vitkup.csv')
media = list(seeds_df['kegg_id'].values)
media_sets = list(itertools.combinations(media, 1))

# species = [sp for sp in species if sp not in old_species]

def propagate_rxns_for_medium(orgrxns, medium):
    # Defining the seed set to be the medium and the currency,
    seedVec = np.zeros(len(rxnMat.T))
    seedVec[[kegg_to_id[e] for e in Currency + medium]] = 1

    # Getting all the reactions performable by this organism.
    avrxns = orgrxns[:]

    # Calculating the metabolites within the scope of 
    # this organism's reaction network. 
    scopeMets = giveScope(rxnMat[avrxns], prodMat[avrxns], seedVec, sumRxnVec[avrxns])[0]
    
    # Finding how much of the core is within the network's scope.
    return uniqify([id_to_kegg[e] for e in np.where(scopeMets)[0]])


# Getting sets as synthetic communities.
species_set = species[40:50]
print(species_set)
for this_species in tqdm(species_set):
    # print(this_species)
    these_strains = list(pangenome_df.loc[pangenome_df['species'] == this_species]['kegg_abbr'])
    other_strains = list(pangenome_df.loc[pangenome_df['species'] != this_species]['kegg_abbr'])
    if len(these_strains) > 25:
        these_strains = random.sample(these_strains, 25)
    other_strains = random.sample(other_strains, 25)
    pairs = list(itertools.product(these_strains, other_strains))
    
    core_monos, core_cos = [], []
    for this_pair in tqdm(pairs):
        # abiotic_sigma = {e: {} for e in these_strains}
        for medium in tqdm(media_sets):
            # Simulating mono-culture growth.
            scope_mono = [propagate_single_for_medium(org, list(medium)) 
                          for org in this_pair]
            core_mono = [list(set(sc) & set(Core)) for sc in scope_mono]
            
            secs = set(unlistify(scope_mono))
            
            # Simulating iterative secretions.
            while True:
                scope_co = [propagate_single_for_medium(org, list(medium) + list(secs)) 
                            for org in this_pair]
                new_secs = set(unlistify(scope_co))
                if len(new_secs) == len(secs):
                    break
                secs = deepcopy(new_secs)
            
            # Simulating co-culture biosynthetic capacity.
            core_co = [list(set(sc) & set(Core)) for sc in scope_co]

            core_monos.append(core_mono)
            core_cos.append(core_co)

    pickle.dump(core_monos, open('inter_species_growth/' + this_species + '_core_monos.dat', 'wb'))
    pickle.dump(core_cos, open('inter_species_growth/' + this_species + '_core_cos.dat', 'wb'))

