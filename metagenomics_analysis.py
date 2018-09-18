import requests
import pandas as pd
from bs4 import BeautifulSoup
from tqdm import tqdm
import urllib.request
import numpy as np
from unlistify import unlistify
from uniqify import uniqify
import itertools

###############################################################################
# Sample retrieval from Chaffron et al data
###############################################################################
samples = pd.read_csv('gg_sample_details_otus_filtered_file.d.01.tsv', sep='\t', header=None).fillna('')
events = samples[1].tolist()
otu_col = samples[6].tolist()
otus_semiraw = [e.split()[-1] for e in otu_col]
extracted_otu_cluster_dict = {}
for e, otu in zip(events, otus_semiraw):
    try:
        extracted_otu_cluster_dict[e].append(otu)
    except:
        extracted_otu_cluster_dict[e] = [otu]

extracted_otu_clusters = list(extracted_otu_cluster_dict.values())

###############################################################################
# Processing BLAST query file to match otu names to 16s seq_ids
###############################################################################
blast_infile = open('gg_otus_allseqs.d.01.fas.txt', 'r')
seqs_raw = list(filter(None, [s.strip() for s in blast_infile.readlines() 
                              if not s.isupper()]))[2:]
blast_infile.close()

seq_clusters = {}
seq_to_otu_dict = {}
for entry in seqs_raw:
    if entry[0] == '#':
        this_otu = entry[2:].split()[0]
        seq_clusters[this_otu] = []
        continue
    seq_clusters[this_otu].append(int(entry.split()[0][1:]))
    seq_to_otu_dict[int(entry.split()[0][1:])] = this_otu

###############################################################################
# Reading BLAST output and mapping strains to OTUs
###############################################################################
strain_to_seq_map = pd.read_csv('test_out.txt', sep='\t', 
                    usecols=['strain_abbr', 'seq_id']).set_index(
                    'strain_abbr').to_dict()['seq_id']
strain_to_otu_map = {}
for s in strain_to_seq_map:
    try:
        strain_to_otu_map[s] = seq_to_otu_dict[strain_to_seq_map[s]]
    except:
        pass
# strain_to_otu_map = {s : seq_to_otu_dict[strain_to_seq_map[s]] 
#                      for s in strain_to_seq_map}
otu_to_strain_map = {value: key for key, value in strain_to_otu_map.items()}

unfound_strains, found_strains = [], []
for s in list(strain_to_otu_map.values()):
    if s not in unlistify(extracted_otu_clusters):
        unfound_strains.append(otu_to_strain_map[s])
    else:
        found_strains.append(otu_to_strain_map[s])

###############################################################################
# Measuring strain-to-strain co-occurrence
###############################################################################
pangenome_df = pd.read_csv('fuller_pangenome_df.csv')
old_pangenome_df = pd.read_csv('pangenome_df.csv')
species = uniqify(list(pangenome_df['species'].values))
old_species = uniqify(list(old_pangenome_df['species'].values))
all_strains = uniqify(list(pangenome_df['strain'].values))
strain_abbrs = uniqify(list(pangenome_df['kegg_abbr'].values))

occur_propensity = {}
cooccuring_strains = set()
for this_species in tqdm(species):
    these_strains = list(pangenome_df.loc[pangenome_df['species'] == this_species]['kegg_abbr'])
    available_strains = [s for s in these_strains if s in found_strains]
    pairs = list(itertools.combinations(available_strains, 2))

    if not pairs:
        continue

    cooccurrence = 0
    for this_pair in tqdm(pairs):
        for this_clust in extracted_otu_clusters:
            if (strain_to_otu_map[this_pair[0]] in this_clust and
                strain_to_otu_map[this_pair[1]] in this_clust):
                cooccurrence += 1
                cooccuring_strains.add(this_pair[0])
                cooccuring_strains.add(this_pair[1])
    occur_propensity[this_species] = cooccurrence / len(pairs)

cooccuring_species = list(occur_propensity.keys())

###############################################################################
# Plotting the MDP pattern for co-occurring species
###############################################################################
# known_metrics_df = pd.read_csv('core_acc/fuller_measured_metrics.csv')
# cooccuring_df = known_metrics_df[known_metrics_df['species'].isin(cooccuring_species)]
# import seaborn as sns
# sns.jointplot(x="MDP", y="phi", data=cooccuring_df, kind='reg')
# plt.xlim(-0.3, 4)
# plt.ylim(-0.01, 0.11)
# # plt.savefig('plots/cooccuring_species_metagenomics.svg')
# plt.show()

strain_spec_list = [(st, pangenome_df.loc[pangenome_df['kegg_abbr'] == st,:]['species'].values[0])
                    for st in cooccuring_strains]
strain_spec_occur_df = pd.DataFrame(strain_spec_list,
                                columns=['kegg_abbr', 'species'])
strain_spec_occur_df.to_csv('all_cooccurs.csv', index=None)
