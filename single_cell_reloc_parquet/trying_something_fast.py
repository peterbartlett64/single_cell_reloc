#%%
import pandas as pd
import os
#%%
df = pd.read_csv("C:/Users/pcnba/OneDrive - University of Toronto/Desktop/SGD_features.tab", sep = "\t")
df.rename(columns= {"1":   "Primary SGDID", "2":   "Feature type", "3":   "Feature qualifier", "4":   "Feature name", "5":   "Standard gene name", "6":   "Alias (multiples separated by |)", "7":   "Parent feature name", "8":   "Secondary SGDID (optional, multiples separated by |)", "9":   "Chromosome", "10":  "Start_coordinate", "11":  "Stop_coordinate", "12":  "Strand", "13":  "Genetic position", "14":  "Coordinate version", "15":  "Sequence version", "16":  "Description"}, inplace=True)
dfORF_CDS = df.loc[(df['Feature type'] == 'ORF') |(df['Feature type'] == 'ORF')]# | (df['Feature type'] == 'CDS')]
#%%
def f_closest_same(chrom, start, stop, strand, df):
	check_strand = strand
	search_res = df.loc[(df["Chromosome"] == chrom) & (df["Strand"] == check_strand)].reset_index(drop = True)
	i = search_res.loc[(search_res['Start_coordinate'] == start) & (search_res['Stop_coordinate'] == stop)].index
	upstream = search_res.iloc[i-1,:]
	downstream = search_res.iloc[i+1:]

	proximal = (upstream['Primary SGDID'], downstream['Primary SGDID'])
	distance = (start - upstream['Stop_cooridnate'], downstream['Start_coordinate'] - stop)
	return(proximal, distance)

def f_closest_opposite(chrom, start, stop, strand,df):
	if strand == "W":
		check_strand = "C"
	elif start == "C":
		check_strand = "W"
	search_res = df.loc[(df["Chromosome"] == chrom) & (df["Strand"] == check_strand)].reset_index(drop = True)
	i = search_res.loc[(search_res['Start_coordinate'] == start) & (search_res['Stop_coordinate'] == stop)].index[0]
	try:
		upstream_r = search_res.iloc[i-1,:]
		upstream = upstream_r['Primary SGDID']
		d_u = start - upstream['Start_coordinate']
	except IndexError:
		upstream = None
		d_u = None

	try:
		downstream_r= search_res.iloc[i+1:]
		downstream = downstream_r['Primary SGDID']
		d_d = downstream['Stop_cooridinate'] - stop
	except IndexError:
		downstream = None
		d_d = None
	proximal = (upstream, downstream)
	distance = (d_u, d_d)
	return(proximal, distance)
#%%
dfORF_CDS['Closest_same'] = dfORF_CDS.apply(lambda x: f_closest_same(x['Chromosome'], x['Start_coordinate'], x['Stop_coordinate'], x['Strand'], dfORF_CDS), axis = 1)
#%%
dfORF_CDS['Closest_opposite'] = dfORF_CDS.apply(lambda x: f_closest_opposite(x['Chromosome'], x['Start_coordinate'], x['Stop_coordinate'], x['Strand'], dfORF_CDS))
# %%
