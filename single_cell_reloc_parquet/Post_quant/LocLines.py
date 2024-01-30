#%%
#, This is an attempt at using the lifelines package to do survival analysis on the data
import pandas as pd
from lifelines import KaplanMeierFitter

full_data = pd.read_parquet("D:\ALL_FINAL\Combined_by_perc\merged_data_final.parquet")
#%%
data_after = full_data.loc[full_data['Frames_post_treatment'] >= 0]
data_after.set_index('Cell_Barcode', inplace = True)
data_cats = data_after[['Protein']]
#%%

#%%
#, This is the simplest impolementation of checking the first frame which is True. It will not be portable to more high end covariates
moved_frame = data_after.groupby('Cell_Barcode').apply(lambda x: x[x['Yet'] == 1]['Frames_post_treatment'].iloc[0] if x['Yet'].any() else None)

#Check the postions where there was a true duration, and then set the Observed for those cells to 1
moved_frame = moved_frame.rename(columns = {0: 'Max_neg_time'})
moved_frame.loc[~(moved_frame['Duration'].isna()), 'Observed'] = 1 #* The barcodes which are not NA are the ones where the cell moved
moved_frame['Observed'] = moved_frame['Observed'].fillna('0')
max_post = data_after.groupby('Cell_Barcode').Frames_post_treatment.agg('max').to_frame()

merg = moved_frame.merge(max_post, left_index = True, right_index = True)
merg['max_neg_time'] = merg['max_neg_time'].fillna(merg['Frames_post_treatment'])


#* Using the progen bud binary column to check first time cell exists, check what the corresponding Frame_post_treatment is
first_frame = data_after.groupby('Cell_Barcode').apply(lambda x: x[x['Progen_bud_binary'] == 1]['Frames_post_treatment'].iloc[0] if x['Progen_bud_binary'].any() else None).rename(columns = {0: 'Max_neg_time'})
merg = merg.merge(first_frame, left_index = True, right_index = True)
merg['Duration'] = merg['Max_neg_time'] - merg['Frames_post_treatment']

#%%
#* Merge in categorical covariates from data_cats
merg_g = merg.merge(data_cats, right_index = True, left_index = True, how = 'left')

#%%
kmf = KaplanMeierFitter()
T = merg['Duration']
C = merg['Observed']
kmf.fit(T,C)
#%%
#* Plotting the survival curve for each protein
for r in merg_g['Protein'].unique():
	ix = merg_g['Protein'] == r
	kmf.fit(T.loc[ix], C.loc[ix], label = r)
	kmf.plot(ax = ix)
