holder = 'http://rest.kegg.jp/link/'
path = 'strain_kofiles/'
# Getting all EC numbers for each species.
for ts in tqdm(pangenome_dict):
    for strain in pangenome_dict[ts]:
        if not os.path.isfile(path + strain + '.txt'):
            urllib.request.urlretrieve( holder + strain + '/ko', path + strain + '.txt' )

ko_dict = {}
for ts in tqdm(pangenome_dict):
    for strain in pangenome_dict[ts]:
        kos = pd.read_csv( path + strain + '.txt', sep='\t', header=None, names=['ko', 'gene'] )
        ko_dict[strain] = list(kos['ko'].drop_duplicates().values)

###############################################################################
# Fluidity analysis.
###############################################################################
def ko_genome_fluidity(strains):
    pairs = list(itertools.combinations(strains, 2))
    fluidity = []
    N = len(strains)
    for (o1, o2) in pairs:
        # Getting all the reactions performable by this organism.
        can1 = set(ko_dict[o1])
        can2 = set(ko_dict[o2])
        fluidity.append(len(can1 ^ can2) / len(can1 | can2))


    return np.mean(fluidity)

pangenome_df = pd.read_csv('fuller_pangenome_df.csv')
species = uniqify(list(pangenome_df['species'].values))
ko_genome_fluidity_dict = {}
for this_species in tqdm(species):
    these_strains = list(pangenome_df.loc[pangenome_df['species'] == this_species]['kegg_abbr'])
    ko_genome_fluidity_dict[this_species] = ko_genome_fluidity(these_strains)

ko_fluidity_df = pd.DataFrame(list(ko_genome_fluidity_dict.items()), columns=['species', 'ko_phi'])
ko_fluidity_df = ko_fluidity_df.sort_values(['ko_phi'])
gen_core_df = pd.read_csv('core_acc/andreani_pangenome_df.csv')
mergedf = ko_fluidity_df.merge(gen_core_df, on='species', how='inner')
import seaborn as sns
sns.jointplot(x="phi", y="ko_phi", data=mergedf)
plt.show()
