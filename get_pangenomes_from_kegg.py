import requests
import pandas as pd
from bs4 import BeautifulSoup
from tqdm import tqdm
import urllib.request
import numpy as np

url = 'http://www.genome.jp/kegg/catalog/org_list4.html'
response = requests.get(url)
soup = BeautifulSoup(response.text, 'lxml')
bacTable = soup.find_all('table')[2]

td_tags = bacTable.find_all('td')
strain_dict = {}
i = 0
while i < len(td_tags):
    tag = td_tags[i]
    if len(tag.attrs) != 1:
        i += 1
    elif tag.attrs['align'] == 'left':
        strain_dict[tag.findAll(text=True)[0]] = int(td_tags[i+2].findAll(text=True)[0])
        i += 3

names, nums = zip(*strain_dict.items())
names, nums = np.array(names), np.array(nums)

REP_CUTOFF = 5
# Getting a list of species that pass the representation cutoff.
allowed_species = names[np.where(nums > REP_CUTOFF)]

# Now getting links to each organism.
holder = 'http://www.genome.jp'
links = {}
i = 0
while i < len(td_tags):
    tag = td_tags[i]
    if len(tag.attrs) != 1: 
       i += 1
    elif tag.attrs['align'] == 'left':
        links[tag.findAll(text=True)[0]] = tag.find('a').attrs['href']
        i += 3

# Downloading the species pangenome pages for the allowed species from KEGG.
pangenome_dict = {}
pangenome_full_dict = {}
for sp in tqdm(allowed_species):
    response = requests.get(holder + links[sp])
    soup = BeautifulSoup(response.text, 'lxml')
    sp_table = soup.find_all('table')[3]

    tags = sp_table.find_all('tr')[0].find_all('tr')
    pangenome_dict[sp] = [t.find('a').findAll(text=True)[0] for t in tags]
    pangenome_full_dict[sp] = [t.find_all('td')[-1].findAll(text=True)[0] for t in tags]

holder = 'http://rest.kegg.jp/link/'
path = 'strain_ecfiles/'
# Getting all EC numbers for each species.
for ts in tqdm(pangenome_dict):
    for strain in pangenome_dict[ts]:
        if not os.path.isfile(path + strain + '.txt'):
            urllib.request.urlretrieve( holder + strain + '/ec', path + strain + '.txt' )

# Getting all KEGG reaction IDs.
rxn_dict = {}
refDF = pd.read_csv('ECREF.txt', delimiter='\t', sep='delimiter', header=None, names=['ec', 'reaction'])
for ts in tqdm(pangenome_dict):
    for strain in pangenome_dict[ts]:
        ecs = pd.read_csv( path + strain + '.txt', sep='\t', header=None, names=['ec', 'gene'] )
        reactions = list(set(pd.merge( ecs, refDF, how='right', on='ec' ).dropna( )[ 'reaction' ].tolist( )))
        keggRxns = [ thisRxn[ 3: ] for thisRxn in reactions ]
        rxn_dict[strain] = keggRxns
        with open('strain_reactions/' + strain + '.txt', 'w') as f:
            f.write('\n'.join(rxn_dict[strain]) + '\n')
        f.close()

# Saving the pangenome as a pandas dataframe
pangenome_df = pd.DataFrame([(e,pangenome_full_dict[e][i], pangenome_dict[e][i]) 
                             for e in pangenome_dict 
                             for i in range(len(pangenome_dict[e]))], 
                             columns=['species', 'strain', 'kegg_abbr'])
# pangenome_df.to_csv('fuller_pangenome_df.csv', index=None)
