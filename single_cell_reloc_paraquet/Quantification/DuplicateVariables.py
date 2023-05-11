#,,,,,,,,,,,,,,,,, kdalkfd aklfakf a
DennormGFPbyMedian_object = DenGFP_object/factor_median_OBJ_GFP
		DennormGFPbyMean_object = DenGFP_object/factor_mean_OBJ_GFP
		DennormGFPbyTotal_intensity_object = DenGFP_object/Denfactor_total_OBJ_GFP

		## Normalizing by the median mean and sum of intensity of the non-cell space
		DennormGFPbyBackground_median = DenGFP_object/Denfactor_GFP_background_Med
		DennormGFPbyBackground_mean = DenGFP_object/Denfactor_GFP_background_Avg
		DennormGFPbyBackground_total = DenGFP_object/Denfactor_GFP_background_Tot


		## Decile bins for normalized against OBJECT (cell) MEDIAN
		# x10thPercentile_norm_OBJ_Median_GFP = np.percentile(DennormGFPbyMedian_object, 10, axis = 0)
		# x20thPercentile_norm_OBJ_Median_GFP = np.percentile(DennormGFPbyMedian_object, 20, axis = 0)
		# x30thPercentile_norm_OBJ_Median_GFP = np.percentile(DennormGFPbyMedian_object, 30, axis = 0)
		# x40thPercentile_norm_OBJ_Median_GFP = np.percentile(DennormGFPbyMedian_object, 40, axis = 0)
		# x50thPercentile_norm_OBJ_Median_GFP = np.percentile(DennormGFPbyMedian_object, 50, axis = 0)
		# x60thPercentile_norm_OBJ_Median_GFP = np.percentile(DennormGFPbyMedian_object, 60, axis = 0)
		# x70thPercentile_norm_OBJ_Median_GFP = np.percentile(DennormGFPbyMedian_object, 70, axis = 0)
		x80thPercentile_norm_OBJ_Median_denGFP = np.percentile(DennormGFPbyMedian_object, 80, axis = 0)
		x90thPercentile_norm_OBJ_Median_denGFP = np.percentile(DennormGFPbyMedian_object, 90, axis = 0)
		x95thPercentile_norm_OBJ_Median_denGFP = np.percentile(DennormGFPbyMedian_object, 95, axis = 0)
		x99thPercentile_norm_OBJ_Median_denGFP = np.percentile(DennormGFPbyMedian_object, 99, axis = 0) # This is for finding foci

		## Decile bins for normalized against OBJECT (cell) MEAN
		# x10thPercentile_norm_OBJ_Mean_Intensity_GFP = np.percentile(DennormGFPbyMean_object, 10, axis = 0)
		# x20thPercentile_norm_OBJ_Mean_Intensity_GFP = np.percentile(DennormGFPbyMean_object, 20, axis = 0)
		# x30thPercentile_norm_OBJ_Mean_Intensity_GFP = np.percentile(DennormGFPbyMean_object, 30, axis = 0)
		# x40thPercentile_norm_OBJ_Mean_Intensity_GFP = np.percentile(DennormGFPbyMean_object, 40, axis = 0)
		# x50thPercentile_norm_OBJ_Mean_Intensity_GFP = np.percentile(DennormGFPbyMean_object, 50, axis = 0)
		# x60thPercentile_norm_OBJ_Mean_Intensity_GFP = np.percentile(DennormGFPbyMean_object, 60, axis = 0)
		# x70thPercentile_norm_OBJ_Mean_Intensity_GFP = np.percentile(DennormGFPbyMean_object, 70, axis = 0)
		# x80thPercentile_norm_OBJ_Mean_Intensity_GFP = np.percentile(DennormGFPbyMean_object, 80, axis = 0)
		x90thPercentile_norm_OBJ_Mean_Intensity_denGFP = np.percentile(DennormGFPbyMean_object, 90, axis = 0)
		x95thPercentile_norm_OBJ_Mean_Intensity_denGFP = np.percentile(DennormGFPbyMean_object, 95, axis = 0)
		x99thPercentile_norm_OBJ_Mean_Intensity_denGFP = np.percentile(DennormGFPbyMean_object, 99, axis = 0)  # This is for finding foci

		## Decile bins for normalized against OBJECT (cell) TOTAL
		# x20thPercentile_norm_OBJ_Total_Intensity_GFP = np.percentile(DennormGFPbyTotal_intensity_object, 20, axis = 0)
		# x10thPercentile_norm_OBJ_Total_Intensity_GFP = np.percentile(DennormGFPbyTotal_intensity_object, 10, axis = 0)
		# x30thPercentile_norm_OBJ_Total_Intensity_GFP = np.percentile(DennormGFPbyTotal_intensity_object, 30, axis = 0)
		# x40thPercentile_norm_OBJ_Total_Intensity_GFP = np.percentile(DennormGFPbyTotal_intensity_object, 40, axis = 0)
		# x50thPercentile_norm_OBJ_Total_Intensity_GFP = np.percentile(DennormGFPbyTotal_intensity_object, 50, axis = 0)
		# x60thPercentile_norm_OBJ_Total_Intensity_GFP = np.percentile(DennormGFPbyTotal_intensity_object, 60, axis = 0)
		# x70thPercentile_norm_OBJ_Total_Intensity_GFP = np.percentile(DennormGFPbyTotal_intensity_object, 70, axis = 0)
		# x80thPercentile_norm_OBJ_Total_Intensity_GFP = np.percentile(DennormGFPbyTotal_intensity_object, 80, axis = 0)
		x90thPercentile_norm_OBJ_Total_Intensity_denGFP = np.percentile(DennormGFPbyTotal_intensity_object, 90, axis = 0)
		x95thPercentile_norm_OBJ_Total_Intensity_denGFP = np.percentile(DennormGFPbyTotal_intensity_object, 95, axis = 0)
		x99thPercentile_norm_OBJ_Total_Intensity_denGFP = np.percentile(DennormGFPbyTotal_intensity_object, 99, axis = 0)  # This is for finding foci

		## Decile bins normalized against BACKGROUND intesntity MEDIAN
		# x10thPercentile_norm_BKGRND_Median_GFP = np.percentile(DennormGFPbyBackground_median, 10, axis = 0)
		# x20thPercentile_norm_BKGRND_Median_GFP = np.percentile(DennormGFPbyBackground_median, 20, axis = 0)
		# x30thPercentile_norm_BKGRND_Median_GFP = np.percentile(DennormGFPbyBackground_median, 30, axis = 0)
		# x40thPercentile_norm_BKGRND_Median_GFP = np.percentile(DennormGFPbyBackground_median, 40, axis = 0)
		# x50thPercentile_norm_BKGRND_Median_GFP = np.percentile(DennormGFPbyBackground_median, 50, axis = 0)
		# x60thPercentile_norm_BKGRND_Median_GFP = np.percentile(DennormGFPbyBackground_median, 60, axis = 0)
		# x70thPercentile_norm_BKGRND_Median_GFP = np.percentile(DennormGFPbyBackground_median, 70, axis = 0)
		# x80thPercentile_norm_BKGRND_Median_GFP = np.percentile(DennormGFPbyBackground_median, 80, axis = 0)
		x90thPercentile_norm_BKGRND_Median_denGFP = np.percentile(DennormGFPbyBackground_median, 90, axis = 0)
		x95thPercentile_norm_BKGRND_Median_denGFP = np.percentile(DennormGFPbyBackground_median, 95, axis = 0)
		x99thPercentile_norm_BKGRND_Median_denGFP = np.percentile(DennormGFPbyBackground_median, 99, axis = 0)  # This is for finding foci

		##  Decile bins normalized against BACKGROUND intenstiy MEAN
		# x10thPercentile_norm_BKGRND_Mean_GFP = np.percentile(DennormGFPbyBackground_mean, 10, axis = 0)
		# x20thPercentile_norm_BKGRND_Mean_GFP = np.percentile(DennormGFPbyBackground_mean, 20, axis = 0)
		# x30thPercentile_norm_BKGRND_Mean_GFP = np.percentile(DennormGFPbyBackground_mean, 30, axis = 0)
		# x40thPercentile_norm_BKGRND_Mean_GFP = np.percentile(DennormGFPbyBackground_mean, 40, axis = 0)
		# x50thPercentile_norm_BKGRND_Mean_GFP = np.percentile(DennormGFPbyBackground_mean, 50, axis = 0)
		# x60thPercentile_norm_BKGRND_Mean_GFP = np.percentile(DennormGFPbyBackground_mean, 60, axis = 0)
		# x70thPercentile_norm_BKGRND_Mean_GFP = np.percentile(DennormGFPbyBackground_mean, 70, axis = 0)
		# x80thPercentile_norm_BKGRND_Mean_GFP = np.percentile(DennormGFPbyBackground_mean, 80, axis = 0)
		x90thPercentile_norm_BKGRND_Mean_denGFP = np.percentile(DennormGFPbyBackground_mean, 90, axis = 0)
		x95thPercentile_norm_BKGRND_Mean_denGFP = np.percentile(DennormGFPbyBackground_mean, 95, axis = 0)
		x99thPercentile_norm_BKGRND_Mean_denGFP = np.percentile(DennormGFPbyBackground_mean, 99, axis = 0)  # This is for finding foci

		##  Decile bins normalized against BACKGROUND intensity TOTAL
		# x10thPercentile_norm_BKGRND_Total_intensity_GFP = np.percentile(DennormGFPbyBackground_total, 10, axis = 0)
		# x20thPercentile_norm_BKGRND_Total_intensity_GFP = np.percentile(DennormGFPbyBackground_total, 20, axis = 0)
		# x30thPercentile_norm_BKGRND_Total_intensity_GFP = np.percentile(DennormGFPbyBackground_total, 30, axis = 0)
		# x40thPercentile_norm_BKGRND_Total_intensity_GFP = np.percentile(DennormGFPbyBackground_total, 40, axis = 0)
		# x50thPercentile_norm_BKGRND_Total_intensity_GFP = np.percentile(DennormGFPbyBackground_total, 50, axis = 0)
		# x60thPercentile_norm_BKGRND_Total_intensity_GFP = np.percentile(DennormGFPbyBackground_total, 60, axis = 0)
		# x70thPercentile_norm_BKGRND_Total_intensity_GFP = np.percentile(DennormGFPbyBackground_total, 70, axis = 0)
		# x80thPercentile_norm_BKGRND_Total_intensity_GFP = np.percentile(DennormGFPbyBackground_total, 80, axis = 0)
		x90thPercentile_norm_BKGRND_Total_intensity_denGFP = np.percentile(DennormGFPbyBackground_total, 90, axis = 0)
		x95thPercentile_norm_BKGRND_Total_intensity_denGFP = np.percentile(DennormGFPbyBackground_total, 95, axis = 0)
		x99thPercentile_norm_BKGRND_Total_intensity_denGFP = np.percentile(DennormGFPbyBackground_total, 99, axis = 0)  # This is for finding foci

		denGFP_spread = len(DenGFP_object)

		##### Metrics from KO intensity
		averageIntensity_denmKO_Frame = np.mean(Denraw_mKO) #average intesnity of frame including the CELLS and NON-CEll SPACES
		averageIntesntiy_mKO_Background = np.mean(DenmKO_background_Raw) #average intensity of NON-CELL space
		averageIntensity_mKO_Object= np.mean(DenmKO_object) #average intensity of GFP within cell

		segCellmKOvalues = DenmKO_object - backgroundKO_val # subtract background
		finalKOvalues = segCellmKOvalues.copy()#*segCellmKOvalues*segCellmKOvalues # improve signal-noise ratio

		normKObyMedian_object = DenmKO_object/factor_median_OBJ_mKate
		normKObyMean_object = DenmKO_object/factor_mean_OBJ_mKate
		normKObyTotal_intensity_object = DenmKO_object/factor_total_OBJ_mKate

		# normKObyBackground_median = DenmKO_object/factor_mKO_background_Med
		# normKObyBackground_mean = DenmKO_object/factor_mKO_background_Avg
		# normKObyBackground_total = DenmKO_object/factor_mKO_background_Tot

		x80thPercentile_Diff_background_mKO = np.percentile(finalKOvalues, 80)
		x90thPercentile_Diff_background_mKO = np.percentile(finalKOvalues, 90)
		x99thPercentile_Diff_background_mKO = np.percentile(finalKOvalues, 99)


		x60thPercentile_denGFP_RAW = np.percentile(DenGFP_object,60)
		x80thPercentile_denGFP_RAW = np.percentile(DenGFP_object,80)
		x90thPercentile_denGFP_RAW = np.percentile(DenGFP_object,90)
		x95thPercentile_denGFP_RAW = np.percentile(DenGFP_object,95)
		x99thPercentile_denGFP_RAW = np.percentile(DenGFP_object,99)
		max_denGFP_RAW = max(DenGFP_object)

		x60thPercentile_denmKa_RAW = np.percentile(DenmKate_object,60)
		x80thPercentile_denmKa_RAW = np.percentile(DenmKate_object,80)
		x90thPercentile_denmKa_RAW = np.percentile(DenmKate_object,90)
		x95thPercentile_denmKa_RAW = np.percentile(DenmKate_object,95)
		x99thPercentile_denmKa_RAW = np.percentile(DenmKate_object,99)
		max_denmKa_RAW = max(DenmKate_object)

		x60thPercentile_denmKO_RAW = np.percentile(DenmKO_object,60)
		x80thPercentile_denmKO_RAW = np.percentile(DenmKO_object,80)
		x90thPercentile_denmKO_RAW = np.percentile(DenmKO_object,90)
		x95thPercentile_denmKO_RAW = np.percentile(DenmKO_object,95)
		x99thPercentile_denmKO_RAW = np.percentile(DenmKO_object,99)
		max_denmKO_RAW = max(DenmKO_object)

		# x90thPercentile_norm_OBJ_Median_mKO = np.percentile(normKObyMedian_object, 90)
		# x99thPercentile_norm_OBJ_Median_mKO = np.percentile(normKObyMedian_object, 99)
		# x90thPercentile_norm_OBJ_Mean_mKO = np.percentile(normKObyMean_object, 90)
		# x99thPercentile_norm_OBJ_Mean_mKO = np.percentile(normKObyMean_object, 99)
		# x90thPercentile_norm_OBJ_Total_intensity_mKO = np.percentile(normKObyTotal_intensity_object, 90)
		# x99thPercentile_norm_OBJ_Total_intensity_mKO = np.percentile(normKObyTotal_intensity_object, 99)

		# x90thPercentile_norm_BKGRND_Median_mKO = np.percentile(normKObyBackground_median, 90)
		# x99thPercentile_norm_BKGRND_Median_mKO = np.percentile(normKObyBackground_median, 99)
		# x90thPercentile_norm_BKGRND_Mean_mKO = np.percentile(normKObyBackground_mean, 90)
		# x99thPercentile_norm_BKGRND_Mean_mKO = np.percentile(normKObyBackground_mean, 99)
		# x90thPercentile_norm_BKGRND_Total_intensity_mKO = np.percentile(normKObyBackground_total, 90)
		# x99thPercentile_norm_BKGRND_Total_intensity_mKO = np.percentile(normKObyBackground_total, 99)


		mKO_spread = len(DenmKO_object)

		#### Metrics from mKate intensity
		averageIntensity_mKate_Frame = np.mean(raw_mKate) #average intesnity of frame including the CELLS and NON-CEll SPACES
		averageIntesntiy_mKate_Background = np.mean(mKate_background_Raw) #average intensity of NON-CELL space
		averageIntensity_mKate_Object= np.mean(DenmKate_object) #average intensity of GFP within cell

		# normKatebyMedian_object = DenmKO_object/factor_median_OBJ_mKate
		# normKatebyMean_object = DenmKO_object/factor_mean_OBJ_mKate
		# normKatebyTotal_intensity_object = DenmKO_object/factor_total_OBJ_mKate

		segCellmKAvalues = DenmKate_object - backgroundKA_val
		finalKAvalues = segCellmKAvalues.copy()#*segCellmKAvalues*segCellmKAvalues

		normKatebyBackground_median = DenmKate_object/factor_mKate_background_Med
		normKatebyBackground_mean = DenmKate_object/factor_mKate_background_Avg
		normKatebyBackground_total = np.sum(DenmKate_object)/factor_mKate_background_Tot

		#Percentiles of Cubic
		x80thPercentile_Diff_background_mKate = np.percentile(finalKAvalues, 80)
		x90thPercentile_Diff_background_mKate = np.percentile(finalKAvalues, 90)
		x99thPercentile_Diff_background_mKate = np.percentile(finalKAvalues, 99)

		# x90thPercentile_norm_OBJ_Median_mKate = np.percentile(normKatebyMedian_object, 90)
		# x99thPercentile_norm_OBJ_Median_mKate = np.percentile(normKatebyMedian_object, 99)
		# x90thPercentile_norm_OBJ_Mean_mKate = np.percentile(normKatebyMean_object, 90)
		# x99thPercentile_norm_OBJ_Mean_mKate = np.percentile(normKatebyMean_object, 99)
		# x90thPercentile_norm_OBJ_Total_intensity_mKate = np.percentile(normKatebyTotal_intensity_object, 90)
		# x99thPercentile_norm_OBJ_Total_intensity_mKate = np.percentile(normKatebyTotal_intensity_object, 99)

		# x90thPercentile_norm_BKGRND_Median_mKate = np.percentile(normKatebyBackground_median, 90)
		# x99thPercentile_norm_BKGRND_Median_mKate = np.percentile(normKatebyBackground_median, 99)
		# x90thPercentile_norm_BKGRND_Mean_mKate = np.percentile(normKatebyBackground_mean, 90)
		# x99thPercentile_norm_BKGRND_Mean_mKate = np.percentile(normKatebyBackground_mean, 99)
		# x90thPercentile_norm_BKGRND_Total_intesnity_mKate = np.percentile(normKatebyBackground_total, 90)
		# x99thPercentile_norm_BKGRND_Total_intensity_mKate = np.percentile(normKatebyBackground_total, 99)

		DenmKate_spread = len(DenmKate_object)




		added on Novemebr 30th because it is required for merging with info
								# 'x10thPercentile_norm_OBJ_Median_GFP' : [x10thPercentile_norm_OBJ_Median_GFP],
								# 'x20thPercentile_norm_OBJ_Median_GFP' : [x20thPercentile_norm_OBJ_Median_GFP],
								# 'x30thPercentile_norm_OBJ_Median_GFP' : [x30thPercentile_norm_OBJ_Median_GFP],
								# 'x40thPercentile_norm_OBJ_Median_GFP' : [x40thPercentile_norm_OBJ_Median_GFP],
								# 'x50thPercentile_norm_OBJ_Median_GFP' : [x50thPercentile_norm_OBJ_Median_GFP],
								# 'x60thPercentile_norm_OBJ_Median_GFP' : [x60thPercentile_norm_OBJ_Median_GFP],
								# 'x70thPercentile_norm_OBJ_Median_GFP' : [x70thPercentile_norm_OBJ_Median_GFP],
								'x80thPercentile_norm_OBJ_Median_GFP' : [x80thPercentile_norm_OBJ_Median_GFP],
								'x90thPercentile_norm_OBJ_Median_GFP' : [x90thPercentile_norm_OBJ_Median_GFP],
								'x95thPercentile_norm_OBJ_Median_GFP' : [x95thPercentile_norm_OBJ_Median_GFP],
								'x99thPercentile_norm_OBJ_Median_GFP' : [x99thPercentile_norm_OBJ_Median_GFP],

								# 'x10thPercentile_norm_OBJ_Mean_Intensity_GFP' : [x10thPercentile_norm_OBJ_Mean_Intensity_GFP],
								# 'x20thPercentile_norm_OBJ_Mean_Intensity_GFP' : [x20thPercentile_norm_OBJ_Mean_Intensity_GFP],
								# 'x30thPercentile_norm_OBJ_Mean_Intensity_GFP' : [x30thPercentile_norm_OBJ_Mean_Intensity_GFP],
								# 'x40thPercentile_norm_OBJ_Mean_Intensity_GFP' : [x40thPercentile_norm_OBJ_Mean_Intensity_GFP],
								# 'x50thPercentile_norm_OBJ_Mean_Intensity_GFP' : [x50thPercentile_norm_OBJ_Mean_Intensity_GFP],
								# 'x60thPercentile_norm_OBJ_Mean_Intensity_GFP' : [x60thPercentile_norm_OBJ_Mean_Intensity_GFP],
								# 'x70thPercentile_norm_OBJ_Mean_Intensity_GFP' : [x70thPercentile_norm_OBJ_Mean_Intensity_GFP],
								# 'x80thPercentile_norm_OBJ_Mean_Intensity_GFP' : [x80thPercentile_norm_OBJ_Mean_Intensity_GFP],
								'x90thPercentile_norm_OBJ_Mean_Intensity_GFP' : [x90thPercentile_norm_OBJ_Mean_Intensity_GFP],
								'x95thPercentile_norm_OBJ_Mean_Intensity_GFP' : [x95thPercentile_norm_OBJ_Mean_Intensity_GFP],
								'x99thPercentile_norm_OBJ_Mean_Intensity_denGFP' : [x99thPercentile_norm_OBJ_Mean_Intensity_denGFP],

								# 'x10thPercentile_norm_OBJ_Total_Intensity_GFP' : [x10thPercentile_norm_OBJ_Total_Intensity_GFP],
								# 'x20thPercentile_norm_OBJ_Total_Intensity_GFP' : [x20thPercentile_norm_OBJ_Total_Intensity_GFP],
								# 'x30thPercentile_norm_OBJ_Total_Intensity_GFP' : [x30thPercentile_norm_OBJ_Total_Intensity_GFP],
								# 'x40thPercentile_norm_OBJ_Total_Intensity_GFP' : [x40thPercentile_norm_OBJ_Total_Intensity_GFP],
								# 'x50thPercentile_norm_OBJ_Total_Intensity_GFP' : [x50thPercentile_norm_OBJ_Total_Intensity_GFP],
								# 'x60thPercentile_norm_OBJ_Total_Intensity_GFP' : [x60thPercentile_norm_OBJ_Total_Intensity_GFP],
								# 'x70thPercentile_norm_OBJ_Total_Intensity_GFP' : [x70thPercentile_norm_OBJ_Total_Intensity_GFP],
								# 'x80thPercentile_norm_OBJ_Total_Intensity_GFP' : [x80thPercentile_norm_OBJ_Total_Intensity_GFP],
								'x90thPercentile_norm_OBJ_Total_Intensity_denGFP' : [x90thPercentile_norm_OBJ_Total_Intensity_denGFP],
								'x95thPercentile_norm_OBJ_Total_Intensity_denGFP' : [x95thPercentile_norm_OBJ_Total_Intensity_denGFP],
								'x99thPercentile_norm_OBJ_Total_Intensity_denGFP' : [x99thPercentile_norm_OBJ_Total_Intensity_denGFP],

								# 'x10thPercentile_norm_BKGRND_Median_GFP' : [x10thPercentile_norm_BKGRND_Median_GFP],
								# 'x20thPercentile_norm_BKGRND_Median_GFP' : [x20thPercentile_norm_BKGRND_Median_GFP],
								# 'x30thPercentile_norm_BKGRND_Median_GFP' : [x30thPercentile_norm_BKGRND_Median_GFP],
								# 'x40thPercentile_norm_BKGRND_Median_GFP' : [x40thPercentile_norm_BKGRND_Median_GFP],
								# 'x50thPercentile_norm_BKGRND_Median_GFP' : [x50thPercentile_norm_BKGRND_Median_GFP],
								# 'x60thPercentile_norm_BKGRND_Median_GFP' : [x60thPercentile_norm_BKGRND_Median_GFP],
								# 'x70thPercentile_norm_BKGRND_Median_GFP' : [x70thPercentile_norm_BKGRND_Median_GFP],
								# 'x80thPercentile_norm_BKGRND_Median_GFP' : [x80thPercentile_norm_BKGRND_Median_GFP],
								'x90thPercentile_norm_BKGRND_Median_denGFP' : [x90thPercentile_norm_BKGRND_Median_denGFP],
								'x95thPercentile_norm_BKGRND_Median_denGFP' : [x95thPercentile_norm_BKGRND_Median_denGFP],
								'x99thPercentile_norm_BKGRND_Median_denGFP' : [x99thPercentile_norm_BKGRND_Median_denGFP],

								# 'x10thPercentile_norm_BKGRND_Mean_GFP' : [x10thPercentile_norm_BKGRND_Mean_GFP],
								# 'x20thPercentile_norm_BKGRND_Mean_GFP' : [x20thPercentile_norm_BKGRND_Mean_GFP],
								# 'x30thPercentile_norm_BKGRND_Mean_GFP' : [x30thPercentile_norm_BKGRND_Mean_GFP],
								# 'x40thPercentile_norm_BKGRND_Mean_GFP' : [x40thPercentile_norm_BKGRND_Mean_GFP],
								# 'x50thPercentile_norm_BKGRND_Mean_GFP' : [x50thPercentile_norm_BKGRND_Mean_GFP],
								# 'x60thPercentile_norm_BKGRND_Mean_GFP' : [x60thPercentile_norm_BKGRND_Mean_GFP],
								# 'x70thPercentile_norm_BKGRND_Mean_GFP' : [x70thPercentile_norm_BKGRND_Mean_GFP],
								# 'x80thPercentile_norm_BKGRND_Mean_GFP' : [x80thPercentile_norm_BKGRND_Mean_GFP],
								'x90thPercentile_norm_BKGRND_Mean_denGFP' : [x90thPercentile_norm_BKGRND_Mean_denGFP],
								'x95thPercentile_norm_BKGRND_Mean_denGFP' : [x95thPercentile_norm_BKGRND_Mean_denGFP],
								'x99thPercentile_norm_BKGRND_Mean_denGFP' : [x99thPercentile_norm_BKGRND_Mean_denGFP],

								# 'x10thPercentile_norm_BKGRND_Total_intensity_GFP' : [x10thPercentile_norm_BKGRND_Total_intensity_GFP],
								# 'x20thPercentile_norm_BKGRND_Total_intensity_GFP' : [x20thPercentile_norm_BKGRND_Total_intensity_GFP],
								# 'x30thPercentile_norm_BKGRND_Total_intensity_GFP' : [x30thPercentile_norm_BKGRND_Total_intensity_GFP],
								# 'x40thPercentile_norm_BKGRND_Total_intensity_GFP' : [x40thPercentile_norm_BKGRND_Total_intensity_GFP],
								# 'x50thPercentile_norm_BKGRND_Total_intensity_GFP' : [x50thPercentile_norm_BKGRND_Total_intensity_GFP],
								# 'x60thPercentile_norm_BKGRND_Total_intensity_GFP' : [x60thPercentile_norm_BKGRND_Total_intensity_GFP],
								# 'x70thPercentile_norm_BKGRND_Total_intensity_GFP' : [x70thPercentile_norm_BKGRND_Total_intensity_GFP],
								# 'x80thPercentile_norm_BKGRND_Total_intensity_GFP' : [x80thPercentile_norm_BKGRND_Total_intensity_GFP],
								'x90thPercentile_norm_BKGRND_Total_intensity_denGFP' : [x90thPercentile_norm_BKGRND_Total_intensity_denGFP],
								'x95thPercentile_norm_BKGRND_Total_intensity_denGFP' : [x95thPercentile_norm_BKGRND_Total_intensity_denGFP],
								'x99thPercentile_norm_BKGRND_Total_intensity_denGFP' : [x99thPercentile_norm_BKGRND_Total_intensity_denGFP],

								'averageIntensity_GFP_Frame' : [averageIntensity_GFP_Frame],
								'averageIntesntiy_GFP_Background' : [averageIntesntiy_GFP_Background],
								'averageIntensity_GFP_Object' : [averageIntensity_GFP_Object],
								'denGFP_spread': [denGFP_spread],
								# 'GFP_background': [GFP_background_Raw],

								'Progen_bud': [progen_bud],

								'x80thPercentile_Diff_background_mKate' : [x80thPercentile_Diff_background_mKate],
								'x90thPercentile_Diff_background_mKate' : [x90thPercentile_Diff_background_mKate],
								'x99thPercentile_Diff_background_mKate' : [x99thPercentile_Diff_background_mKate],

								# 'x90thPercentile_norm_OBJ_Median_mKO' : [x90thPercentile_norm_OBJ_Median_mKO],
								# 'x99thPercentile_norm_OBJ_Median_mKO' : [x99thPercentile_norm_OBJ_Median_mKO],
								# 'x90thPercentile_norm_OBJ_Mean_mKO' : [x90thPercentile_norm_OBJ_Mean_mKO],
								# 'x99thPercentile_norm_OBJ_Mean_mKO' : [x99thPercentile_norm_OBJ_Mean_mKO],
								# 'x90thPercentile_norm_OBJ_Total_intensity_mKO' : [x90thPercentile_norm_OBJ_Total_intensity_mKO],
								# 'x99thPercentile_norm_OBJ_Total_intensity_mKO' : [x99thPercentile_norm_OBJ_Total_intensity_mKO],
								# 'x90thPercentile_norm_BKGRND_Median_mKO' : [x90thPercentile_norm_BKGRND_Median_mKO],
								# 'x99thPercentile_norm_BKGRND_Median_mKO' : [x99thPercentile_norm_BKGRND_Median_mKO],
								# 'x90thPercentile_norm_BKGRND_Mean_mKO' : [x90thPercentile_norm_BKGRND_Mean_mKO],
								# 'x99thPercentile_norm_BKGRND_Mean_mKO' : [x99thPercentile_norm_BKGRND_Mean_mKO],
								# 'x90thPercentile_norm_BKGRND_Total_intensity_mKO' : [x90thPercentile_norm_BKGRND_Total_intensity_mKO],
								# 'x99thPercentile_norm_BKGRND_Total_intensity_mKO' : [x99thPercentile_norm_BKGRND_Total_intensity_mKO],

								'x80thPercentile_Diff_background_mKO' : [x80thPercentile_Diff_background_mKO],
								'x90thPercentile_Diff_background_mKO' : [x90thPercentile_Diff_background_mKO],
								'x99thPercentile_Diff_background_mKO' : [x99thPercentile_Diff_background_mKO],

								'averageIntensity_denmKO_Frame' : [averageIntensity_denmKO_Frame],
								'averageIntesntiy_mKO_Background' : [averageIntesntiy_mKO_Background],
								'averageIntensity_mKO_Object' : [averageIntensity_mKO_Object],
								'mKO_spread': [mKO_spread],
								# 'mKO_background': [DenmKO_background_Raw],

								# 'x90thPercentile_norm_OBJ_Median_mKate' : [x90thPercentile_norm_OBJ_Median_mKate],
								# 'x99thPercentile_norm_OBJ_Median_mKate' : [x99thPercentile_norm_OBJ_Median_mKate],
								# 'x90thPercentile_norm_OBJ_Mean_mKate' : [x90thPercentile_norm_OBJ_Mean_mKate],
								# 'x99thPercentile_norm_OBJ_Mean_mKate' : [x99thPercentile_norm_OBJ_Mean_mKate],
								# 'x90thPercentile_norm_OBJ_Total_intensity_mKate' : [x90thPercentile_norm_OBJ_Total_intensity_mKate],
								# 'x99thPercentile_norm_OBJ_Total_intensity_mKate' : [x99thPercentile_norm_OBJ_Total_intensity_mKate],
								# 'x90thPercentile_norm_BKGRND_Median_mKate' : [x90thPercentile_norm_BKGRND_Median_mKate],
								# 'x99thPercentile_norm_BKGRND_Median_mKate' : [x99thPercentile_norm_BKGRND_Median_mKate],
								# 'x90thPercentile_norm_BKGRND_Mean_mKate' : [x90thPercentile_norm_BKGRND_Mean_mKate],
								# 'x99thPercentile_norm_BKGRND_Mean_mKate' : [x99thPercentile_norm_BKGRND_Mean_mKate],
								# 'x90thPercentile_norm_BKGRND_Total_Intesnity_mKate' : [x90thPercentile_norm_BKGRND_Total_intesnity_mKate],
								# 'x99thPercentile_norm_BKGRND_Total_intensity_mKate' : [x99thPercentile_norm_BKGRND_Total_intensity_mKate],

								# 'x80thPercentile_Diff_background_mKa' : [x80thPercentile_Diff_background_mKate],
								# 'x90thPercentile_Diff_background_mKa' : [x90thPercentile_Diff_background_mKate],
								# 'x99thPercentile_Diff_background_mKa' : [x99thPercentile_Diff_background_mKate],

								'averageIntensity_mKate_Frame' : [averageIntensity_mKate_Frame],
								'averageIntesntiy_mKate_Background' : [averageIntesntiy_mKate_Background],
								'averageIntensity_mKate_Object' : [averageIntensity_mKate_Object],
								'DenmKate_spread': [DenmKate_spread],
								# 'mKate_background': [mKate_background_Raw],

								# Show calcultion factors specific to the cell and frame number (timepoint)
								'factor_median_OBJ_GFP' : [factor_median_OBJ_GFP],
								'factor_mean_OBJ_GFP' : [factor_mean_OBJ_GFP],
								'factor_total_OBJ_GFP' : [factor_total_OBJ_GFP],
								'factor_GFP_background_Med' : [factor_GFP_background_Med],
								'factor_GFP_background_Avg' : [factor_GFP_background_Avg],
								'factor_GFP_background_Tot' : [factor_GFP_background_Tot],

								'factor_median _OBJ_KO' : [factor_median_OBJ_KO],
								'factor_mean_OBJ_KO' : [factor_mean_OBJ_KO],
								'factor_total_OBJ_KO' : [factor_total_OBJ_KO],
								'factor_mKO_background_Med' : [factor_mKO_background_Med],
								'factor_mKO_background_Avg' : [factor_mKO_background_Avg],
								'factor_mKO_background_Tot' : [factor_mKO_background_Tot],

								'factor_median_OBJ_mKate' : [factor_median_OBJ_mKate],
								'factor_mean_OBJ_mKate' : [factor_mean_OBJ_mKate],
								'factor_total_OBJ_mKate' : [factor_total_OBJ_mKate],
								'factor_mKate_background_Med' : [factor_mKate_background_Med],
								'factor_mKate_background_Avg' : [factor_mKate_background_Avg],
								'factor_mKate_background_Tot' : [factor_mKate_background_Tot],

								'x80thPercentile_denGFP_RAW':[x80thPercentile_denGFP_RAW],
								'x60thPercentile_denGFP_RAW':[x60thPercentile_denGFP_RAW],
								'x90thPercentile_denGFP_RAW':[x90thPercentile_denGFP_RAW],
								'x95thPercentile_denGFP_RAW':[x95thPercentile_denGFP_RAW],
								'x99thPercentile_denGFP_RAW':[x99thPercentile_denGFP_RAW],
								'max_denGFP_RAW' : [max_denGFP_RAW],
								'x60thPercentile_denmKa_RAW' : [x60thPercentile_denmKa_RAW],
								'x80thPercentile_denmKa_RAW' : [x80thPercentile_denmKa_RAW],
								'x90thPercentile_denmKa_RAW' : [x90thPercentile_denmKa_RAW],
								'x95thPercentile_denmKa_RAW' : [x95thPercentile_denmKa_RAW],
								'x99thPercentile_denmKa_RAW' : [x99thPercentile_denmKa_RAW],
								'max_denmKa_RAW' : [max_denmKa_RAW],
								'x60thPercentile_denmKO_RAW': [x60thPercentile_denmKO_RAW],
								'x80thPercentile_denmKO_RAW': [x80thPercentile_denmKO_RAW],
								'x90thPercentile_denmKO_RAW': [x90thPercentile_denmKO_RAW],
								'x95thPercentile_denmKO_RAW': [x95thPercentile_denmKO_RAW],
								'x99thPercentile_denmKO_RAW': [x99thPercentile_denmKO_RAW],
								'max_denmKO_RAW' : [max_denmKO_RAW]
								}