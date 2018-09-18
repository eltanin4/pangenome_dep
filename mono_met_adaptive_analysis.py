import numpy as np
import pickle
import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm
from uniqify import uniqify
from unlistify import unlistify
from load_kegg import *
from give_scope import giveScope
import itertools
from sampler import collect_samples
from copy import deepcopy
import random
import operator

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

# species = [sp for sp in species if sp in old_species]
# pangenome_df = old_pangenome_df[:]
simDir = 'growth_tests/'

known_metrics_df = pd.read_csv('core_acc/fuller_known_metrics_kegg.csv')
known_metrics_df = known_metrics_df.loc[:, ['species', 'ko_phi']]

def core_rxns(strains, CUTOFF=0.95):
    all_sets = []
    for org in strains:
        # Getting all the reactions performable by this organism.
        can = ''.join(open('strain_reactions/' + org + '.txt', 'r').readlines()).split()
        can = list((set(can)) & set(rxns))
        all_sets.append(set([rxn_kegg_to_id[e] for e in can]))

    pangenes = uniqify(unlistify(all_sets))
    core_genes = []
    for gene in pangenes:
        g_frac = 0.0
        for oi, org in enumerate(strains):
            if gene in all_sets[oi]:
                g_frac += 1
        if g_frac / len(strains) >= CUTOFF:
            core_genes.append(gene)

    return core_genes

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

def propagate_core_for_medium(corerxns, medium):
    # Defining the seed set to be the medium and the currency,
    seedVec = np.zeros(len(rxnMat.T))
    seedVec[[kegg_to_id[e] for e in Currency + medium]] = 1

    # Getting all the reactions performable by these corerxns.
    can = ['R' + (5 - len(str(int(e)))) * '0' + str(int(e)) for e in corerxns]
    can = list((set(can)) & set(rxns))
    avrxns = [rxn_kegg_to_id[e] for e in can]

    # Calculating the metabolites within the scope of 
    # this organism's reaction network. 
    scopeMets = giveScope(rxnMat[avrxns], prodMat[avrxns], seedVec, sumRxnVec[avrxns])[0]
    
    # Finding how much of the core is within the network's scope.
    return uniqify([id_to_kegg[e] for e in np.where(scopeMets)[0]])

known_fluidity = []
mono_met_adapt = []
core_prod = {}
for this_species in tqdm(species):
    if this_species == 'Mannheimia haemolytica':
        continue
    these_strains = list(pangenome_df.loc[pangenome_df['species'] == this_species]['kegg_abbr'])

    core_core = {}
    this_core = core_rxns(these_strains, 1.0)
    for i, medium in enumerate(media_sets):
        core_core[i] = propagate_core_for_medium(this_core, list(medium))
    core_set = set(unlistify(core_core.values())) & set(Core)

    core_prod[this_species] = list(core_set)
    extra_acc = {}
    for this_strain in tqdm(these_strains):
        monos = []
        for medium in tqdm(media_sets):
            monos.append(propagate_single_for_medium(this_strain, list(medium)))
        extra_acc[this_strain] = [(set(e) & set(Core)).difference(core_core[i]) for e, i in zip(monos, range(len(media_sets)))]
        # extra_acc[this_strain] = [bool((set(e) & set(Core).difference(core_set)) for e in monos]

    known_fluidity.append(known_metrics_df.loc[
                         known_metrics_df['species'] == this_species, :]['ko_phi'].values[0])
    mono_met_adapt.append(extra_acc.items())

known_fluidity = []
mono_met_adapt = []
core_prod = {}
for this_species in tqdm(species):
    if this_species == 'Mannheimia haemolytica':
        continue
    these_strains = list(pangenome_df.loc[pangenome_df['species'] == this_species]['kegg_abbr'])

    core_core = {}
    this_core = core_rxns(these_strains, 1.0)
    for i, medium in enumerate(media_sets):
        core_core[i] = propagate_core_for_medium(this_core, list(medium))
    core_set = set(unlistify(core_core.values()))

    core_prod[this_species] = list(core_set)
    extra_acc = {}
    for this_strain in tqdm(these_strains):
        monos = []
        for medium in tqdm(media_sets):
            monos.append(propagate_single_for_medium(this_strain, list(medium)))
        extra_acc[this_strain] = [(set(e)).difference(core_core[i]) for e, i in zip(monos, range(len(media_sets)))]
        # extra_acc[this_strain] = [bool((set(e) & set(Core).difference(core_set)) for e in monos]

    known_fluidity.append(known_metrics_df.loc[
                         known_metrics_df['species'] == this_species, :]['ko_phi'].values[0])
    mono_met_adapt.append(extra_acc.items())

# Differences in the biosynthetic abilities of different strains becomes more pronounced as their overall accessory genome content increases, i.e. with increasing genome fluidity.
# mean: Spearman rho=0.62, p=1e-11
known_fluidity = []
met_fluidity = []
for this_species in tqdm(species):
    these_strains = list(pangenome_df.loc[pangenome_df['species'] == this_species]['kegg_abbr'])
    pairs = list(itertools.combinations(these_strains, 2))
    this_extra_acc = dict(mono_met_adapt[species.index(this_species)])
    
    mean_pair_overlap = []

    for (o1, o2) in pairs:
        # this_pair_overlap = []
        o1_can = set(unlistify(this_extra_acc[o1]))
        o2_can = set(unlistify(this_extra_acc[o2]))

        # for i, medium in enumerate(tqdm(media_sets)):
        if o1_can and o2_can:
            this_pair_overlap = len(o1_can ^ o2_can) / len(o1_can | o2_can)
            mean_pair_overlap.append(this_pair_overlap)
    
    known_fluidity.append(known_metrics_df.loc[
                         known_metrics_df['species'] == this_species, :]['ko_phi'].values[0])
    met_fluidity.append(np.mean(mean_pair_overlap))

mergedf = pd.DataFrame(list(zip(known_fluidity, met_fluidity)), columns=['phi', 'eta'])
import seaborn as sns
sns.jointplot(x="eta", y="phi", data=mergedf, kind='reg')
plt.xlim(-0.02, 0.3)
plt.ylim(-0.015, 0.25)
plt.savefig('plots/mono_met_overlap_corr_phi_eta_core.svg')
plt.show()

# Tentative: # of potentially advantageous strain-specific (not in to all strains)
# biosynthetic (core) molecules producible per strain (note: only by the accessory content).
# To test if "ostensibly dispensible genes can be commonly advantageous" and if that advantage 
# scales with the extent of "dispensible" gene content.
known_fluidity = []
niche_spec_adapts = []
spec_names = []
for this_species in tqdm(species):
    these_strains = list(pangenome_df.loc[pangenome_df['species'] == this_species]['kegg_abbr'])
    this_extra_acc = dict(mono_met_adapt[species.index(this_species)])
    
    ts = []
    ns = []
    for i, medium in enumerate(media_sets):
        s = [this_extra_acc[this_strain][i] for this_strain in these_strains]
        ts.append(set.union(*s) - set.intersection(*s))
        ns.append(np.median([len(e - set.intersection(*s)) for e in s]))

    # niche_spec_adapts.append(len(set(unlistify(ts))) / len(these_strains))

    spec_names.append(this_species)
    known_fluidity.append(known_metrics_df.loc[
                         known_metrics_df['species'] == this_species, :]['ko_phi'].values[0])
    niche_spec_adapts.append(np.median(ns))

# mean_mean: Spearman rho=0.44, p=7e-6
mergedf = pd.DataFrame(list(zip(spec_names, known_fluidity, niche_spec_adapts)), columns=['species', 'phi', 'mono_niche'])
import seaborn as sns
sns.jointplot(x="mono_niche", y="phi", data=mergedf, kind='reg')
plt.xlim(-1, max(niche_spec_adapts) + 2)
plt.ylim(-0.015, 0.25)
plt.savefig('plots/medians_mono_niche_specific_corr_phi_core.svg')
plt.show()
