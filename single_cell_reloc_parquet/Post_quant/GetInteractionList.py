#%%
import pandas as pd
import os

if __name__ == '__main__':
	Global_Variables = {
		"microfluidics_results": "F:/Microfluidics/RES_N_ULTS",
		"information_path": "C:/Users/pcnba/Grant Brown's Lab Dropbox/Peter Bartlett/Peter Bartlett Data/Code/Data_copies/Information_files",
		"post_path": "D:/ALL_FINAL"}
	os.chdir(Global_Variables['information_path'])

	string_interactions = pd.read_csv("C:/Users/pcnba/Grant Brown's Lab Dropbox/Peter Bartlett/Peter Bartlett Data/Code/Data_copies/Information_files/string_interactions.tsv", sep = '\t')
	string_interactions.rename(columns={'#node1': 'node1'}, inplace = True)

	string_nodes = pd.read_csv("C:/Users/pcnba/Grant Brown's Lab Dropbox/Peter Bartlett/Peter Bartlett Data/Code/Data_copies/Information_files/STRING_nodes.csv")
	string_nodes['Protein'] = string_nodes['@id'].apply(lambda x: x[14:])
	string_nodes['@id'] = string_nodes['@id'].apply(lambda x: x[9:])

	cluster_interactions = pd.read_csv("C:/Users/pcnba/Grant Brown's Lab Dropbox/Peter Bartlett/Peter Bartlett Data/Code/Data_copies/Information_files/Clustered_stringUniverse.csv")

	cluster_interactions[['Prot1', 'inter', 'Prot2']] = cluster_interactions['shared name'].str.split(pat = '(ppp)', expand = True)
	cluster_interactions.drop(columns = ['shared name', 'inter'], inplace = True)
	cluster_interactions['Prot1'] = cluster_interactions['Prot1'].apply(lambda x: x[5:-2])
	cluster_interactions['Prot2'] = cluster_interactions['Prot2'].apply(lambda x: x[7:])

	Prots1 = string_interactions['node1_string_id'].apply(lambda x: x[5:]).to_list()
	Prot2 = string_interactions['node2_string_id'].apply(lambda x: x[5:]).to_list()

	# Prot1 = cluster_interactions['Prot1'].to_list()
	# Prot2 = cluster_interactions['Prot2'].to_list()

	n_interactions = string_interactions.groupby('node1').node2.agg('count')
	n_interactions = pd.DataFrame(n_interactions).reset_index(drop = False)
	n_interactions.rename(columns = {'node1': 'Protein' ,'node2': 'n_interactions'}, inplace = True)
	n_interactions.sort_values(by = 'n_interactions', inplace= True)

	Prots_all = Prot1 + Prot2
	Prots_all = list(set(Prots_all))

	linked = string_nodes[string_nodes['Protein'].isin(Prots_all)]['Protein']
	unlinked = string_nodes[~string_nodes['Protein'].isin(Prots_all)]['Protein']



#%%
mcode = pd.read_csv("C:/Users/pcnba/Grant Brown's Lab Dropbox/Peter Bartlett/Peter Bartlett Data/Code/Data_copies/Information_files/Mcode_test.txt", sep = '\t').set_index('Cluster')


string_nodes.loc[string_nodes['@id'].isin(mcode.loc[2,'Node IDs'].split(', '))]['updated_yet_perc'].mean()
string_nodes.loc[string_nodes['@id'].isin(mcode.loc[2,'Node IDs'].split(', '))]['updated_yet_perc'].std()

#%%
string_nodes['C1']
string_nodes['C2']
string_nodes['C3']
string_nodes['C4']
string_nodes['C5']
string_nodes['C6']
string_nodes['C7']
string_nodes['C8']
string_nodes['C9']
string_nodes['C10']
string_nodes['C11']
string_nodes['C12']



#%%
