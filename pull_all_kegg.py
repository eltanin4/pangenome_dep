import os
import numpy as np
from bs4 import BeautifulSoup
import urllib3
import urllib.request
import shutil
import pickle
import re
import collections
from print_progress_bar import *
from copy import deepcopy
from tqdm import tqdm

#-------------------------------------------------------------------------
# Getting all compound and reaction maps from KEGG.
#-------------------------------------------------------------------------
import pandas as pd
# refDF = pd.read_csv('compound_to_reaction.txt', delimiter='\t', sep='delimiter', header=None, names=['cpd', 'rxn'])
# newDF = pd.read_csv('glycan_to_reaction.txt', delimiter='\t', sep='delimiter', header=None, names=['gly', 'rxn'])
# rxnids = sorted(set([int(e[4:]) for e in set(list(refDF['rxn'].values))]))
# rxnids = [ 'R' + (5 - len(str(e))) * '0' + str(e) for e in rxnids ]
# new_rxnids = sorted(set([int(e[4:]) for e in set(list(newDF['rxn'].values))]))
# new_rxnids = [ 'R' + (5 - len(str(e))) * '0' + str(e) for e in new_rxnids ]
# rxnids = sorted(uniqify(rxnids + new_rxnids))
# with open('akshit_rxns.txt', 'w') as f:
#     for r in rxnids:
#         f.write(r + '\n')

# cpdids = sorted(set([int(e[5:]) for e in set(list(refDF['cpd'].values))]))
# cpdids = [ 'C' + (5 - len(str(e))) * '0' + str(e) for e in cpdids ]
# with open('akshit_mets.txt', 'w') as f:
#     for c in cpdids:
#         f.write(c + '\n')

# glyDF = pd.read_csv('glycan_to_reaction.txt', delimiter='\t', sep='delimiter', header=None, names=['gly', 'rxn'])
# glyids = sorted(set([int(e[4:]) for e in set(list(glyDF['gly'].values))]))
# glyids = [ 'G' + (5 - len(str(e))) * '0' + str(e) for e in glyids ]
# with open('akshit_mets.txt', 'a') as f:
#     for c in glyids:
#         f.write(c + '\n')
# f.close()

rxnids = ''.join(open('akshit_rxns.txt', 'r').readlines()).split()
cpdids = ''.join(open('akshit_mets.txt', 'r').readlines()).split()

#-------------------------------------------------------------------------
# This code extracts all KEGG webpages for reactions.
#-------------------------------------------------------------------------
rxn_url_holder = 'http://www.genome.jp/dbget-bin/www_bget?rn:'

handler = urllib.request.ProxyHandler({'http': 'proxy.ncbs.res.in:3128'})
opener = urllib.request.build_opener(handler)
urllib.request.install_opener(opener)
urllib.request.ProxyHandler({'ftp': 'proxy.ncbs.res.in:3128'})

for i in range(len(rxnids)):
    print_progress_bar(i, len(rxnids), 'Pulling all reaction files')
    r_id = rxnids[i]
    current_url = rxn_url_holder + r_id
    path = 'reaction_files/' + r_id + '.html'
    if not os.path.isfile(path):
        urllib.request.urlretrieve(current_url, path)

def modify_stoich_for_reaction(side, st, stoich_rid, stoich_matrix):
    for clause in side.split(' + '):
        # Checking if non-unit stoichiometry:
        if len(clause.split()) == 2:
            c_id = cpd_kegg_to_id[clause.split()[1][:6]]
            if clause.split()[0].isdigit():
                stoich_matrix[stoich_rid, c_id] = st * int(clause.split()[0])
            else:
                stoich_matrix[stoich_rid, c_id] = st
        elif len(clause.split()) == 1:
            stoich_matrix[stoich_rid, cpd_kegg_to_id[clause.split()[0][:6]]] = st
        else:
            raise ValueError('Unidentified clause structure noticed.')

rxn_kegg_to_id = {rxnids[i]: i for i in range(len(rxnids))}
cpd_kegg_to_id = {cpdids[i]: i for i in range(len(cpdids))}
stoich_matrix = np.zeros((len(rxnids) * 2, len(cpdids)))
#-------------------------------------------------------------------------
# This converts the extracted html files into reactions and then into the
# stoichiometric matrix.
#-------------------------------------------------------------------------
irreversibles = []
for i in tqdm(range(len(rxnids))):
    # Setting up things for this iteration.
    r_id = rxnids[i]
    self_rid = rxn_kegg_to_id[r_id]
    path = 'reaction_files/' + r_id + '.html'

    # Extracting the relevant reaction string from the webpage.
    soup = BeautifulSoup(open(path, 'rb'), 'lxml')
    letters = soup.find_all('td', class_='td21')
    letters_2 = soup.find_all('td', class_='td20')
    for el in letters + letters_2:
        if 'cpd:' in str(el) or 'gl:' in str(el):
            eq = el
            break

    eq_string = ''.join(el.findAll(text=True))[:-1]


    # Checking if the reaction is reversible or irreversible.
    is_reversible = True
    try:
        reactant_side, product_side = re.split(r'\s*<=>\s*', eq_string)
    except:
        is_reversible = False
        irreversibles.append(r_id)
        reactant_side, product_side = re.split(r'\s*=>\s*', eq_string)

    # Modifying the stoichiometric matrix.
    self_rid = rxn_kegg_to_id[r_id]
    modify_stoich_for_reaction(reactant_side, -1, self_rid, stoich_matrix)
    modify_stoich_for_reaction(product_side, +1, self_rid, stoich_matrix)

    # Adding the reverse reaction to the matrix if needed.
    if is_reversible:
        self_rid = rxn_kegg_to_id[r_id] + len(rxnids)
        modify_stoich_for_reaction(reactant_side, +1, self_rid, stoich_matrix)
        modify_stoich_for_reaction(product_side, -1, self_rid, stoich_matrix)

#-------------------------------------------------------------------------
# This code extracts all KEGG webpages for compounds.
#-------------------------------------------------------------------------
cpd_url_holder = 'http://www.genome.jp/dbget-bin/www_bget?cpd:'
gly_url_holder = 'http://www.genome.jp/dbget-bin/www_bget?gl:'

for i in tqdm(range(len(cpdids))):
    c_id = cpdids[i]
    if c_id[0] == 'C':
        current_url = cpd_url_holder + c_id
    else:
        current_url = gly_url_holder + c_id
    path = 'compound_files/' + c_id + '.html'
    if not os.path.isfile(path):
        urllib.request.urlretrieve(current_url, path)

#-------------------------------------------------------------------------
# This converts the extracted html files into compound names
#-------------------------------------------------------------------------
names, exceptions = [], []
for i in tqdm(range(len(cpdids))):
    # Setting up things for this iteration.
    c_id = cpdids[i]
    self_cid = cpd_kegg_to_id[c_id]
    path = 'compound_files/' + c_id + '.html'

    # Extracting the relevant reaction string from the webpage.
    soup = BeautifulSoup(open(path, 'rb'), 'lxml')
    try:
        containter = soup.find_all('td', class_='td21')[0]
        names.append(''.join(containter.findAll(text=True)).strip().replace('\n', ' '))
    except:
        exceptions.append(c_id)
        names.append('Unknown name')

# # Saving as a file.
# with open('akshit_names.txt', 'w') as f:
#     for n in names:
#         f.write(n + '\n')
# f.close()
