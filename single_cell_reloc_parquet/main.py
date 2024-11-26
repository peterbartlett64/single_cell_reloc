#, Load the modules
from single_cell_reloc_parquet.global_functions.global_variables import *

#! This is set for easy updating of the source files
#, For STAGE 1: PRE-TRACKING/SEGMENTATION
import single_cell_reloc_parquet.Pre_seg.img_indexer as img_indexer
import single_cell_reloc_parquet.Pre_seg.xml_generator as xml_generator
import single_cell_reloc_parquet.Pre_seg.seg_manger as seg_manager

#, For STAGE 2: TRACKING
import single_cell_reloc_parquet.Pre_track.orgAllmask_er as orgAllmask_er
import single_cell_reloc_parquet.Pre_track.TracX_manager as tracking_manager

#, For STAGE 3: POST-TRACKING/QUANTIFICATION
import single_cell_reloc_parquet.Quantification.Cell_mask_expand_er as mask_expander
import single_cell_reloc_parquet.Quantification.Quantification_cell&bud_wRaw_paraquet_RUSH2 as quantification #. This file needs to be checked as the right version

#, For STAGE 4: POST-QUANTIFICATION
import single_cell_reloc_parquet.Post_quant.post_quant_NOV23 as post_quant
import single_cell_reloc_parquet.Post_quant.best_perc as perc_decider
import single_cell_reloc_parquet.Post_quant.merge_data_FINAL as merge_data

from rich import print #* For better printing

def run_pipeline()
	Global_variables = global_manager() #* Initialize global variables. This will make sure that everything is the right type

	#, STAGE 1: PRE-TRACKING/SEGMENTATION
	img_indexer(Global_variables) #* Index the images
	xml_generator(Global_variables) #* Generate the xml files
	seg_manager(Global_variables) #* Segment the images

	#, STAGE 2: TRACKING
	orgAllmask_er(Global_variables) #* Generate the orgAllmask_er files
	tracking_manager(Global_variables) #* Track the cells

	#, stage 3: POST-TRACKING/QUANTIFICATION
	mask_expander(Global_variables) #* Expand the masks
	quantification(Global_variables) #* Quantify the cells

	#, stage 4: POST-QUANTIFICATION. This is the final stage and the meat of this project
	post_quant(Global_variables) #* Post-quantify the cells
	perc_decider(Global_variables) #* Decide the best percentage
	merge_data(Global_variables)

	#, End of pipeline
	print("Pipeline completed successfully")
	print("Please check the output files for the results")
	print('There are included R scripts for visualization which have not been automated yet')
	print("Thank you for using the pipeline")

if __name__ == 'main'
	run_pipeline()
