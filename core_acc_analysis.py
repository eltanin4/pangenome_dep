###############################################################################
# Fluidity analysis.
###############################################################################
def genome_fluidity(strains):
    pairs = list(itertools.combinations(strains, 2))
    fluidity = []
    for (o1, o2) in pairs:
        # Getting all the reactions performable by this organism.
        can1 = ''.join(open('strain_reactions/' + o1 + '.txt', 'r').readlines()).split()
        can2 = ''.join(open('strain_reactions/' + o2 + '.txt', 'r').readlines()).split()
        can1 = (set(can1)) & set(rxns)
        can2 = (set(can2)) & set(rxns)
        fluidity.append(len(can1 ^ can2) / len(can1 | can2))

    return np.mean(fluidity)

genome_fluidity_dict = {}
for this_species in tqdm(species):
    these_strains = list(pangenome_df.loc[pangenome_df['species'] == this_species]['kegg_abbr'])
    genome_fluidity_dict[this_species] = genome_fluidity(these_strains)

fluidity_df = pd.DataFrame(list(genome_fluidity_dict.items()), columns=['species', 'fluidity'])
fluidity_df = fluidity_df.sort_values(['fluidity'])

###############################################################################
# Core analysis
###############################################################################
def return_core_fraction(strains, CUTOFF=0.95):
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
                g_frac += 1 / len(strains)
        if g_frac >= CUTOFF:
            core_genes.append(gene)

    return len(core_genes) / len(pangenes)

core_acc_dict = {}
total_pairs = 0
for this_species in species:
    these_strains = list(pangenome_df.loc[pangenome_df['species'] == this_species]['kegg_abbr'])
    # total_pairs += len(list(itertools.combinations(these_strains, 2)))
    core_acc_dict[this_species] = return_core_fraction(these_strains)

core_df = pd.DataFrame(list(core_acc_dict.items()), columns=['species', 'core_frac'])
core_df = core_df.sort_values(['core_frac'])
core_df.to_csv('core_acc/core_frac_metabolism_kegg.csv', index=None)
core_fluid_df = core_df.merge(fluidity_df, on='species').sort_values(['fluidity'])
core_fluid_df.to_csv('core_acc/core_fluid_metabolism_kegg.csv', index=None)

mergedf = met_core_df.merge(gen_core_df, on='species', how='inner')
import seaborn as sns
sns.jointplot(x="phi", y="fluidity", data=mergedf)
plt.show()
