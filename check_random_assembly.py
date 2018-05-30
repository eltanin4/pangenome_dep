###############################################################################
# Core analysis
###############################################################################
def give_pangenes(strains):
    all_sets = []
    for org in strains:
        # Getting all the reactions performable by this organism.
        kos = pd.read_csv( 'strain_kofiles/' + org + '.txt', sep='\t', header=None, names=['ko', 'gene'] )['ko']
        all_sets.append(set(list(kos.values)))

    return uniqify(unlistify(all_sets))

def return_core_fraction(strains, CUTOFF=0.95):
    all_sets = []
    for org in strains:
        # Getting all the reactions performable by this organism.
        kos = pd.read_csv( 'strain_kofiles/' + org + '.txt', sep='\t', header=None, names=['ko', 'gene'] )['ko']
        all_sets.append(set(list(kos.values)))

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

fig, ax = plt.subplots(1)
keeper = []
for this_species in tqdm(species):
    acc_frac = {}
    these_strains = list(pangenome_df.loc[pangenome_df['species'] == this_species]['kegg_abbr'])
    for this_cutoff in tqdm(np.linspace(0.0, 1.0, 100)):
        acc_frac[this_cutoff] = 1.0 - return_core_fraction(these_strains, this_cutoff)

    x, y = zip(*sorted(acc_frac.items()))
    ax.plot(x, y)
    keeper.append((x, y))
plt.show()

def give_new_core(pangenes, strains, CUTOFF=0.95):
    core_genes = []
    for gene in pangenes:
        g_frac = 0.0
        for oi, org in enumerate(strains):
            if gene in org:
                g_frac += 1 / len(strains)
        if g_frac >= CUTOFF:
            core_genes.append(gene)

    return len(core_genes) / len(pangenes)

fig, ax = plt.subplots(1)
acc_frac_rand, keeper = {}, []
for this_species in tqdm(species):
    these_strains = list(pangenome_df.loc[pangenome_df['species'] == this_species]['kegg_abbr'])
    pangenes = give_pangenes(these_strains)
    new_genomes = []
    for org in these_strains:
        kos = pd.read_csv( 'strain_kofiles/' + org + '.txt', sep='\t', header=None, names=['ko', 'gene'] )
        new_genomes.append(random.sample(pangenes, len(list(kos['ko'].drop_duplicates().values))))
    for this_cutoff in tqdm(np.linspace(0.0, 1.0, 100)):
        acc_frac_rand[this_cutoff] = 1.0 - give_new_core(pangenes, new_genomes, this_cutoff)

    x, y = zip(*sorted(acc_frac_rand.items()))
    ax.plot(x, y)
    keeper.append((x, y))
plt.show()

fig, ax = plt.subplots(1)
for i in range(len(species)):
    ax.plot(rand_keeper[i][0], [1 - e for e in rand_keeper[i][1]], c='black')
    ax.plot(keeper[i][0], [1-e for e in keeper[i][1]], c='r')
plt.show()
