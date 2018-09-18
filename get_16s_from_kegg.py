import requests
import pandas as pd
from bs4 import BeautifulSoup
from tqdm import tqdm
import urllib.request
import numpy as np
from uniqify import uniqify
from unlistify import unlistify
import os

# Getting the set of bacterial abbreviations for each organism in KEGG.
pangenome_df = pd.read_csv('fuller_pangenome_df.csv')
old_pangenome_df = pd.read_csv('pangenome_df.csv')
species = uniqify(list(pangenome_df['species'].values))
old_species = uniqify(list(old_pangenome_df['species'].values))
all_strains = uniqify(list(pangenome_df['kegg_abbr'].values))
index_array = np.array(list(range(len(all_strains))))


holder = 'http://rest.kegg.jp/link/'
path = 'strain_16s_files/'
# Getting all EC numbers for each species.
for strain in tqdm(all_strains):
    if not os.path.isfile(path + strain + '.txt'):
        # K01977 is the KO identifier for 16s rRNA.
        urllib.request.urlretrieve( holder + strain + '/K01977', path + strain + '.txt' )

def give_16s(this_gene):
    nucleotide_holder = 'https://www.genome.jp/dbget-bin/www_bget?-f+-n+n+'    
    url = nucleotide_holder + this_gene
    response = requests.get(url)
    soup = BeautifulSoup(response.text, 'lxml')
    nuc_tag = soup.find_all('pre')[0]
    return find_between_r(str(nuc_tag), '(N)', '</pre>').replace('\n', '')

# Getting all KEGG reaction IDs.
rrna_dict = {}
for strain in tqdm(all_strains):
    gn_df = pd.read_csv( path + strain + '.txt', sep='\t', header=None, names=['ko', 'gene'] )
    genes = list(gn_df['gene'])
    all_16s = [give_16s(tgene) for tgene in genes]
    rrna_dict[strain] = all_16s[:]

# Mapping the longest 16s sequence to each strain as its unique identity.
rrna_dict = {strain:max(rrna_dict[strain], key=len) for strain in all_strains if rrna_dict[strain]}

# Checking if all sequences correspond uniquely to strains.

# Saving the rRNA sequences as a pandas dataframe
rrna_df = pd.DataFrame([(strain, rrna_dict[strain]) 
                        for strain in rrna_dict], 
                        columns=['strain', 'rrna_16s'])
rrna_df.to_csv('16s_rrna_df.csv', index=None)

# Writing .faa file in FASTA format to make a database.
rrna_df = pd.read_csv('16s_rrna_df.csv')
db_file = open('16s_rrna_db.faa', 'w')
for strain in tqdm(rrna_dict):
    db_file.writelines('>' + strain + '\n')
    current_rrnas = rrna_df.loc[rrna_df['strain'] == strain]['rrna_16s'].values
    longest_assigned = max(current_rrnas, key=len)
    db_file.writelines(longest_assigned + '\n')

db_file.close()
