#Heatmap_plot
#%% Import the plotnine functions
import pandas as pd
import os
from plotnine.data import economics
from plotnine import ggplot, aes, geom_pointdensity, aes, stat_bin, facet_wrap, geom_p
import polars as pl
#%%
dataframe = pd.read_parquet("E:/ALL_FINAL/Combined_by_perc/SAE2_selected.parquet")
# dataframe_pl = pl.read_parquet("E:/ALL_FINAL/Combined_by_perc/XRS2d0210_selected.parquet")
#%%

datatypes_dict = {'Date': 'category',
				  'Frame': 'int16',
					"Unique_Frame": 'category'}
# 					"x80thPercentile_norm_OBJ_Median_GFP",
# 					"x90thPercentile_norm_OBJ_Median_GFP",
# 					"x95thPercentile_norm_OBJ_Median_GFP",
# 					"x99thPercentile_norm_OBJ_Median_GFP",
# 					"x90thPercentile_norm_OBJ_Mean_Intensity_GFP",
# 					"x95thPercentile_norm_OBJ_Mean_Intensity_GFP",
# 					"x99thPercentile_norm_OBJ_Mean_Intensity_GFP",
# 					"x90thPercentile_norm_OBJ_Total_Intensity_GFP",
# 					"x95thPercentile_norm_OBJ_Total_Intensity_GFP",
# 					"x99thPercentile_norm_OBJ_Total_Intensity_GFP",
# 					"x90thPercentile_norm_BKGRND_Median_GFP",
# 					"x95thPercentile_norm_BKGRND_Median_GFP",
# 					"x99thPercentile_norm_BKGRND_Median_GFP",
# 					"x90thPercentile_norm_BKGRND_Mean_GFP",
# 					"x95thPercentile_norm_BKGRND_Mean_GFP",
# 					"x99thPercentile_norm_BKGRND_Mean_GFP",
# 					"x90thPercentile_norm_BKGRND_Total_intensity_GFP",
# 					"x95thPercentile_norm_BKGRND_Total_intensity_GFP",
# 					"x99thPercentile_norm_BKGRND_Total_intensity_GFP",
# 					"averageIntensity_GFP_Frame",
# 					"averageIntesntiy_GFP_Background",
# 					"averageIntensity_GFP_Object",
# 					"GFP_spread",
# 					"Progen_bud",
# 					"x80thPercentile_Diff_background_mKate",
# 					"x90thPercentile_Diff_background_mKate",
# 					"x99thPercentile_Diff_background_mKate",
# 					"x80thPercentile_Diff_background_mKO",
# 					"x90thPercentile_Diff_background_mKO",
# 					"x99thPercentile_Diff_background_mKO",
# 					"averageIntensity_mKO_Frame",
# 					"averageIntesntiy_mKO_Background",
# 					"averageIntensity_mKO_Object",
# 					"mKO_spread",
# 					"averageIntensity_mKate_Frame",
# 					"averageIntesntiy_mKate_Background",
# 					"averageIntensity_mKate_Object",
# 					"mKate_spread",
# 					"factor_median_OBJ_GFP",
# 					"factor_mean_OBJ_GFP",
# 					"factor_total_OBJ_GFP",
# 					"factor_GFP_background_Med",
# 					"factor_GFP_background_Avg",
# 					"factor_GFP_background_Tot",
# 					"factor_median _OBJ_KO",
# 					"factor_mean_OBJ_KO",
# 					"factor_total_OBJ_KO",
# 					"factor_mKO_background_Med",
# 					"factor_mKO_background_Avg",
# 					"factor_mKO_background_Tot",
# 					"factor_median_OBJ_mKate",
# 					"factor_mean_OBJ_mKate",
# 					"factor_total_OBJ_mKate",
# 					"factor_mKate_background_Med",
# 					"factor_mKate_background_Avg",
# 					"factor_mKate_background_Tot",
# 					"x80thPercentile_GFP_RAW",
# 					"x60thPercentile_GFP_RAW",
# 					"x90thPercentile_GFP_RAW",
# 					"x95thPercentile_GFP_RAW",
# 					"x99thPercentile_GFP_RAW",
# 					"max_GFP_RAW",
# 					"x60thPercentile_mKa_RAW",
# 					"x80thPercentile_mKa_RAW",
# 					"x90thPercentile_mKa_RAW",
# 					"x95thPercentile_mKa_RAW",
# 					"x99thPercentile_mKa_RAW",
# 					"max_mKa_RAW",
# 					"x60thPercentile_mKO_RAW",
# 					"x80thPercentile_mKO_RAW",
# 					"x90thPercentile_mKO_RAW",
# 					"x95thPercentile_mKO_RAW",
# 					"x99thPercentile_mKO_RAW",
# 					"max_mKO_RAW",
# 					"TrackID_valid",
# 					"Myo1Identity",
# 					"mKO_foldChange",
# 					"mKO_direction",
# 					"mKA_foldChange",
# 					"mKa_direction",
# 					"byProgen_bud",
# 					"byRange",
# 					"Col_info",
# 					"Protein",
# 					"Is_treated",
# 					"Frames_post_treatment",
# 					"Unique_pos",
# 					"Loc_score",
# 					"Upper",
# 					"Lower",
# 					"Relocalized",
# 					"mKa_factor_upper",
# 					"mKa_factor_lower",
# 					"mKa_RAWfactor_upper",
# 					"mKa_RAWfactor_lower",
# 					"<lambda_0>_x",
# 					"cell_count",
# 					"mKO_factor_upper",
# 					"mKO_factor_lower",
# 					"mKO_RAWfactor_upper",
# 					"mKO_RAWfactor_lower",
# 					"<lambda_0>_y",
# 					"count",
# 					"log_Abundance",
# 					"log_Loc_score",
# 					"z_score_logLoc",
# 					"z_score_logAbund":
# }
#%%
dataframe.astype(datatypes_dict)
#%%
(ggplot(dataframe, aes(x = "Loc_score", y = "Frame"))
	+ geom_density_ridges(fill = 'gray90'))
    # + facet_wrap("~Frame"))

#%% Crate the plot( as a function)
def heatmap_points(dataframe:pd.DataFrame, location_save:str):
	plot = (ggplot(dataframe, aes())
		+geom
		+facet_wrap("~frame"))

	return(plot)

lsize = 0.65
fill_alpha = 0.7

(ggplot(dataframe, aes('Frame', 'Loc_score', fill='Frame'))
 + geom_violin(m1, alpha=fill_alpha, size=lsize, show_legend=False)
 + geom_point(m2, color='none', alpha=fill_alpha, size=2, show_legend=False)
 + geom_line(m2, color='gray', size=lsize, alpha=0.6)
 + geom_boxplot(width=shift, alpha=fill_alpha, size=lsize, show_legend=False)
 + scale_fill_manual(values=['dodgerblue', 'darkorange'])
 + theme_classic()
 + theme(figure_size=(8, 6))
)
