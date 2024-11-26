import os
import seaborn as sns
import plotly.express as px
from joblib import Parallel, delayed
from glob import glob

sns.set(rc = {'figure.figsize':(15,8)}) #* Set the figure size to be larger

#, Single Prearson Outputs
# prots = full_data["Protein"].unique()
def graph_pearson_single(i, full_data):
	prot = max_agg_Loc_score.iloc[i,:]["Protein"]
	f = max_agg_Loc_score.iloc[i,:]["Frame_x"]
	try: 
		sub_data = full_data.loc[(full_data["Protein"] == prot) & (full_data["Frame_x"] == f)].dropna(subset=['z_score_Loc', 'z_score_Abund'])
	except ValueError:
		sub_data["z_score_Abund"] = stats.zscore(sub_data.loc[:,"Abundance"])
		sub_data["z_score_Loc"] = stats.zscore(sub_data.loc[:,"Loc_score"])
	z_Loc_Abun_r = stats.pearsonr(sub_data["z_score_Loc"], sub_data["z_score_Abund"])	
	scatter = sns.scatterplot(data= sub_data, x = "z_score_Abund", y = "z_score_Loc")
	fig = scatter.get_figure()
	fig.savefig(f"Pearson{z_Loc_Abun_r}--{prot}_{f}.svg")
	print(f"{prot} is done pearson'n")
	return(prot, z_Loc_Abun_r)

def graph_pearson_plus_minus(i, full_data):
	try:
		prot = max_agg_Loc_score.iloc[i,:]["Protein"]
		f = max_agg_Loc_score.iloc[i,:]["Frame_x"]
		try: 
			sub_data = full_data.loc[(full_data["Protein"] == prot) & ((full_data["Frame_x"] == f)| (full_data["Frame_x"] == f+1)| (full_data["Frame_x"] == f-1))].dropna(subset=['z_score_Loc', 'z_score_Abund'])
		except ValueError:
			# sub_data["z_score_Abund"] = stats.zscore(sub_data.loc[:,"Abundance"])
			# sub_data["z_score_Loc"] = stats.zscore(sub_data.loc[:,"Loc_score"])
			sub_data["z_score_Abund"] = sub_data.groupby("Frame_x")["Abundance"].transform(stats.zscore)
			sub_data["z_score_Loc"] = sub_data.groupby("Frame_x")["Loc_score"].transform(stats.zscore)
		z_Loc_Abun_r = stats.pearsonr(sub_data["z_score_Loc"], sub_data["z_score_Abund"])	
		scatter = sns.scatterplot(data= sub_data, x = "z_score_Abund", y = "z_score_Loc", hue = "Frame_x")
		fig = scatter.get_figure()
		fig.savefig(f"Pearson{z_Loc_Abun_r}--{prot}_{f}.svg")
		print(f"{prot} is done pearson'n")
		return(prot, z_Loc_Abun_r)
	except Exception as e:
		return(f"{e} on {prot}")

pr = os.cpu_count()
# full_data = full_data.loc[full_data["Protein"] == "ECO1"]

if __name__ == "__main__":
	full_data_name = input("FileName")
	full_data = pd.read_csv(full_data_name)
	x = Parallel(n_jobs= pr, verbose= 100)(delayed(graph_pearson_plus_minus)(i = i, full_data = full_data) for i in range(len(max_agg_Loc_score)))
	print(x)
	

