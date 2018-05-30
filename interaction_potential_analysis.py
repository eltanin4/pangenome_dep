import numpy as np
import pickle
import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm
from uniqify import uniqify
from unlistify import unlistify
# from load_kegg import *
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

try:
    species = [sp for sp in species if sp != 'Mannheimia haemolytica']
except:
    pass
simDir = 'growth_tests/'

known_metrics_df = pd.read_csv('core_acc/fuller_known_metrics_kegg.csv')
known_metrics_df = known_metrics_df.loc[:, ['species', 'ko_phi']]

kdeps = {}
known_fluidity = []
mean_int_pot = []
for this_species in tqdm(species):
    kdeps[this_species] = []
    core_monos = pickle.load(open(simDir + this_species + '_core_monos.dat', 'rb'))
    core_cos = pickle.load(open(simDir + this_species + '_core_cos.dat', 'rb'))
    
    for i in tqdm(range(0, len(core_monos), len(media_sets))):
        monos = unlistify(core_monos[i : i + len(media_sets)])
        cos = unlistify(core_cos[i : i + len(media_sets)])
        num_monos, num_cos = list(map(len, monos)), list(map(len, cos))
        mean_int_of_pair = np.mean(list(map(operator.sub, num_cos, num_monos)))

        kdeps[this_species].append(mean_int_of_pair)

    known_fluidity.append(known_metrics_df.loc[
                         known_metrics_df['species'] == this_species, :]['ko_phi'].values[0])
    mean_int_pot.append(np.mean(kdeps[this_species]))

# mean_mean: Spearman rho=0.56, p=4e-9
mergedf = pd.DataFrame(list(zip(known_fluidity, mean_int_pot)), columns=['phi', 'kdep'])
import seaborn as sns
sns.jointplot(x="kdep", y="phi", data=mergedf, kind='reg')
plt.xlim(-0.3, 4)
plt.ylim(-0.02, 0.22)
# plt.savefig('plots/fuller_max_pair_mean_media_corr_phi_kdep.svg')
plt.show()

# new = pd.DataFrame(list(zip(species, mean_int_pot)), columns=['species', 'kdep'])
# known_metrics_df = pd.read_csv('core_acc/fuller_known_metrics_kegg.csv')
# new_metrics_df = known_metrics_df.merge(new, on='species', how='inner')
# new_metrics_df.to_csv('core_acc/fuller_measured_metrics_df.csv', index=None)

# measured_metrics = pd.read_csv('core_acc/chosen_fuller_measured_metrics.csv')
# import seaborn as sns
# sns.jointplot(x="kdep", y="ko_phi", data=measured_metrics, kind='reg')
# plt.savefig('plots/chosen_fuller_mean_pair_mean_media_corr_phi_kdep.svg')
# # plt.xlim(0.0, 4.0)
# plt.show()

# mergedf = pd.DataFrame(list(zip(known_fluidity, mean_int_pot)), columns=['phi', 'kdep'])
# from sklearn import datasets, linear_model
# from sklearn.metrics import mean_squared_error, r2_score
# # Create linear regression object
# regr = linear_model.LinearRegression()
# regr.fit(mean_int_pot, known_fluidity)
# y_fit = regr.predict(mean_int_pot)