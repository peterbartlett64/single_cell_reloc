#%%
import pandas as pd
import os
import single_cell_reloc_parquet.global_functions.global_variables as glv
from icecream import ic
import numpy as np
#%%
#* Dictionaries for renaming reg results
Tkach_dictionary = {'to cyto': ' -> cytoplasm',
	   'to nuc': ' -> nucleus',
	   'nucleus': ' -> nucleus',
	   'nuc foci': ' -> nuclear foci',
	   'cyto foci': ' -> cyto foci',
	   'from budneck/tip': 'budneck ->',
	   'from bud tip': 'budneck ->',
	   'to pm': ' -> plasma membrane',
	   'to bud neck': ' -> budneck',
	   'to nuc/periphery': ' -> nucleus and nuclear periph',
	   'to nuc (from nuc foci)': 'nuclear foci -> nucleus',
	   'nuc foci': ' -> nuclear foci',
	   'from bud neck': 'budneck -> ',
	   'to vac': ' -> vacuole',
	   'to nuc, nuc foci, cyto foci' : ' -> nucleus and nuclear foci and cytoplasm', #* This one is hard to deal with. Almost decided to do a leveled assignment
	   'no cells': np.nan,
	   'to pm foci': ' -> pm foci',
	   'to nucleolus, cyto foci' : ' -> nuceolus and cyto foci',
	   'to vac (from pm)': 'plasma membrane -> vacuole',
	   'to nucleolus': ' -> nucleolus',
	   'to nuc periph (cyto?)': 'cytoplasm -> nuclear periph',
	   'nuc foci (weak)': ' -> nuclear foci',
	   'to cyto (from pm/endosome)': 'plasma membrane -> cytoplasm',
	   'to nuc (from nucleolus)': 'nucleolus -> nucleus',
	   'er foci': ' -> ER foci',
	   'to vac (abund inc)': ' -> valuole',
	   'to vac (focus)': ' -> valuole foci',
	   'to cyto (from pm)': 'plasma membrane -> cytoplasm',
	   'shorter microtubules': ' -> microtubules',
	   'to nuc/nuc periph': ' -> nucleus and nuclear periph',
	   'to cyto (daughter)': ' -> cytoplasm (daughter)',
	   'to diffuse nuc (really, more numerous less intense foci, from foci)': 'nuclear foci -> nucleus and nuclear foci',
	   'to cyto (decreased abund)': 'cytoplasm -> cytoplasm',
	   'from budneck/tip, to pm': 'budneck -> plasma membrane',
	   'to pm foci (endosome)': ' -> pm foci',
	   'to nuc (nucleolus)': ' -> nucleolus',
	   'to nucleolus (weak)': ' -> nucleolus',
	   'not the same strain': np.nan,
	   'to er': ' -> ER',
	   'from budneck': 'budneck -> ',
	   'nuc periph' : ' -> nuclear periph',
	   'to pm (foci)': ' -> pm foci',
	   'to vac (from vac mb)': 'vacuole -> vacuole',
	   'from bud neck (to pm)': 'budneck -> plamsa membrane',
	   'to diffuse nuc (from foci)': 'nuclear foci -> nucleus',
	   'to er foci (weak)': ' -> ER foci',
	   'nuc periph foci': ' -> nuclear foci and nuclear periph',
	   'to nuc (diffuse)': ' -> nucleus',
	   'from bud  tip': 'budneck -> ',
	   'to nuc (nucleolus?)': ' -> nucleus',
	   'to cyto (from vac)': 'vacuole -> cytoplasm',
	   'to er foci(?)': ' -> ER foci',
	   'cyto foci (nuc foci?)': ' -> cyto foci',
	   'to nuc (?)': ' -> nucleus',
	   'to cyto (from cyto foci/vac)' : 'cyto foci and vacuole -> cytoplasm',
	   'to cyto (weak)': ' -> cytoplasm',
	   'vacuole membrane (focus)' : ' -> vacuole foci',
	   'to cyto (from vac, degradation?)': 'vacuole -> cytoplasm',
	   'vac foci?': ' -> vacuole'}



Denervaud_dictionary = {'Cyto->Nuc': 'cytoplasm -> nucleus',
		'Nuc -> Cytoplasm': 'nucleus -> cytoplasm',
		'Nuc -> NucFoci': 'nucleus -> nuclear foci',
		'Nuc Periphery Agg': ' -> nuclear periph',
		'Cyto Agg': 'cytoplasm -> cytoplasm', #!
		'Cyto Disagg': 'cytoplasm -> cytoplasm', #!
		   'From cell Periphery': 'plasma membrane -> cytoplasm',#!
		'To Cell Periphery': 'cytoplasm -> plasma membrane'#!
}

Mazumder_dictionary = {'cytoplasm': ' -> cytoplasm',
		'nucleus': ' -> nucleus',
		'cytoplasm,nucleus': ' -> cytoplasm and nucleus',
		'vacuole': ' -> vacuole',
		'ambiguous': '',
		'cytoplasm,punctate composite': ' -> cyto foci',
		'microtubule': ' -> microtubule',
		'ambiguous,spindle pole': ' -> spindle pole',
		'spindle pole': ' -> spindle pole',
		'mitochondrion': ' -> mitochondrion',
		'peroxisome': ' -> peroxisome',
		'ER,ambiguous,cytoplasm': ' -> ER and cytoplasm',
		'cytoplasm,late Golgi': ' -> cytoplasm and late Golgi',
		'late Golgi': ' -> late Golgi',
		'nucleolus,nucleus': ' -> nucleolus and nucleus',
		'nucleolus': ' -> nucleolus',
		'ER': ' -> ER',
		'actin': ' -> actin',
		'cytoplasm,mitochondrion,nucleus': ' -> cytoplasm and mitochondrion and nucleus',
		'cytoplasm,nucleus,punctate composite': ' -> cytoplasm and nucleus and punctate composite',
		'cytoplasm,nucleolus': ' -> cytoplasm and nucleolus',
		'nuclear periphery': ' -> nuclear periphery',
		'bud neck,cytoplasm,bud': ' -> bud neck and cytoplasm and bud',
		'ER,cell periphery': ' -> ER and cell periphery',
		'punctate composite': ' -> punctate composite',
		'endosome': ' -> endosome',
		'vacuolar membrane': ' -> vacuolar membrane',
		'cell periphery,bud': ' -> cell periphery and bud',
		'bud neck,cytoplasm,mitochondrion,cell periphery': ' -> bud neck and cytoplasm and mitochondrion and cell periphery',
		'cytoplasm,mitochondrion': ' -> cytoplasm and mitochondrion',
		'cytoplasm,nuclear periphery,nucleus': ' -> cytoplasm and nuclear periphery and nucleus',
		'ambiguous,cytoplasm,nucleus': ' -> cytoplasm and nucleus',
		'ambiguous,bud neck,cytoplasm,cell periphery,bud': ' -> bud neck and cytoplasm and cell periphery and bud',
		'cytoplasm,nucleolus,nucleus': ' -> cytoplasm and nucleolus and nucleus',
		'ER,cytoplasm': ' -> ER and cytoplasm',
		'bud neck,cytoplasm': ' -> bud neck and cytoplasm',
		'cytoplasm,nuclear periphery': ' -> cytoplasm and nuclear periphery',
		'ambiguous,nucleus': ' -> nucleus',
		'vacuolar membrane,endosome': ' -> vacuolar membrane and endosome',
		'lipid particle': ' -> lipid particle',
		'bud neck': ' -> bud neck',
		'bud neck,cell periphery': ' -> bud neck and cell periphery',
		'ambiguous,bud neck,cell periphery,bud': ' -> bud neck and cell periphery and bud',
		'cell periphery': ' -> cell periphery',
		'bud neck,cytoplasm,cell periphery,bud': ' -> bud neck and cytoplasm and cell periphery and bud',
		'cytoplasm,actin': ' -> cytoplasm and actin',
		'Golgi': ' -> Golgi',
		'ambiguous,bud neck,cytoplasm,bud': ' -> bud neck and cytoplasm and bud',
		'ER to Golgi': 'ER -> Golgi',
		'ER,ambiguous,bud': ' -> ER and bud',
		'ambiguous,endosome': ' -> endosome',
		'bud neck,cytoplasm,cell periphery': ' -> bud neck and cytoplasm and cell periphery',
		'punctate composite,late Golgi': ' -> punctate composite and late Golgi',
		'cytoplasm,nucleus,bud': ' -> cytoplasm and nucleus and bud',
		'early Golgi': ' -> early Golgi',
		'punctate composite,lipid particle': ' -> punctate composite and lipid particle',
		'bud neck,cell periphery,punctate composite': ' -> bud neck and cell periphery and punctate composite',
		'mitochondrion,nucleus': ' -> mitochondrion and nucleus',
		'ambiguous,cytoplasm,bud': ' -> cytoplasm and bud',
		'cytoplasm,cell periphery': ' -> cytoplasm and cell periphery',
		'ambiguous,cytoplasm': ' -> cytoplasm',
		'cell periphery,vacuole': ' -> cell periphery and vacuole',
		'punctate composite,early Golgi,late Golgi': ' -> punctate composite and early Golgi and late Golgi',
		'nucleus,spindle pole,microtubule': ' -> nucleus and spindle pole and microtubule',
		'cytoplasm,punctate composite,endosome': ' -> cytoplasm and punctate composite and endosome',
		'punctate composite,early Golgi': ' -> punctate composite and early Golgi',
		'bud neck,cytoplasm,nucleus': ' -> bud neck and cytoplasm and nucleus',
		'Golgi,early Golgi': ' -> Golgi and early Golgi',
		'ambiguous,vacuolar membrane': ' -> vacuolar membrane',
		'ambiguous,cell periphery,bud': ' -> cell periphery and bud',
		'cytoplasm,bud': ' -> cytoplasm and bud',
		'ambiguous,late Golgi': ' -> late Golgi',
		'ambiguous,mitochondrion': ' -> mitochondrion',
		'cytoplasm,endosome,lipid particle': ' -> cytoplasm and endosome and lipid particle',
		'bud neck,cytoplasm,cell periphery,punctate composite,bud': 'bud neck and cytoplasm,cell periphery and punctate composite and bud',
		'ER,cytoplasm,nucleus': ' -> ER and cytoplasm and nucleus',
		'ambiguous,bud neck,bud': ' -> bud neck and bud',
		'cytoplasm,vacuole': ' -> cytoplasm and vacuole',
		'mitochondrion,punctate composite': ' -> mitochondrion and punctate composite',
		'early Golgi,late Golgi': ' -> early Golgi and late Golgi',
		'ER,mitochondrion,nuclear periphery': ' -> ER and mitochondrion and nuclear periphery',
		'cytoplasm,spindle pole': ' -> cytoplasm and spindle pole',
		'spindle pole,microtubule': ' -> spindle pole and microtubule',
		'punctate composite,Golgi': ' -> punctate composite and Golgi',
		'vacuole,endosome': ' -> vacuole and endosome',
		'punctate composite,endosome': ' -> punctate composite and endosome',
		'nuclear periphery,nucleus': ' -> nuclear periphery and nucleus',
		'ambiguous,bud neck,cell periphery,vacuole,bud': ' -> bud neck and cell periphery and vacuole and bud',
		'bud neck,cell periphery,bud': ' -> bud neck and cell periphery and bud',
		'cytoplasm,nucleus,spindle pole': ' -> cytoplasm and nucleus and spindle pole',
		'cytoplasm,endosome': ' -> cytoplasm and endosome',
		'ambiguous,cytoplasm,punctate composite': ' -> cytoplasm and punctate composite',
		'cytoplasm,Golgi,early Golgi': ' -> cytoplasm and Golgi and early Golgi',
		'punctate composite,vacuolar membrane,lipid particle': ' -> punctate composite and vacuolar membrane and lipid particle',
		'endosome,lipid particle': ' -> endosome and lipid particle',
		'vacuole,vacuolar membrane': ' -> vacuole and vacuolar membrane',
		'ER,vacuole': ' -> ER and vacuole',
		'ER,cell periphery,bud': ' -> ER,cell periphery and bud',
		'ambiguous,cytoplasm,cell periphery,bud': ' -> cytoplasm and cell periphery and bud',
		'punctate composite,Golgi,early Golgi': ' -> punctate composite and Golgi and early Golgi',
		'cytoplasm,nuclear periphery,punctate composite': ' -> cytoplasm and nuclear periphery and punctate composite',
		'nucleolus,microtubule': ' -> nucleolus and microtubule',
		'ambiguous,punctate composite': ' -> punctate composite',
		'ER,cell periphery,vacuole,bud': ' -> ER and cell periphery and vacuole,bud',
		'punctate composite,spindle pole': ' -> punctate composite and spindle pole',
		'nucleus,spindle pole': ' -> nucleus and spindle pole',
		'ambiguous,bud neck,cytoplasm,vacuole,bud': ' -> bud neck and cytoplasm and vacuole,bud',
		'cytoplasm,cell periphery,bud': ' -> cytoplasm and cell periphery and bud',
		'cytoplasm,early Golgi': ' -> cytoplasm and early Golgi',
		'punctate composite,actin': ' -> punctate composite and actin',
		'bud neck,nucleus': ' -> bud neck and nucleus',
		'ER,nucleus': ' -> ER and nucleus',
		'cytoplasm,punctate composite,spindle pole, microtubule': 'cytoplasm and punctate composite and spindle pole and microtubule',
		'microtubule': ' -> microtubule',
		'ambiguous,Golgi,early Golgi': ' -> Golgi and early Golgi',
		'punctate composite,endosome,early Golgi,late Golgi': ' -> punctate composite and endosome and early Golgi and late Golgi',
		'nucleus,microtubule': ' -> nucleus and microtubule',
		'ambiguous,bud neck,cytoplasm,nucleus,bud': ' -> bud neck and cytoplasm and nucleus and bud',
		'vacuole,Golgi': ' -> vacuole and Golgi',
		'ambiguous,Golgi': ' -> Golgi',
		'ambiguous,bud neck,cell periphery,punctate composite,late Golgi,bud': 'bud neck and cell periphery and punctate composite and late Golgi and bud', 'late Golgi,bud': 'late Golgi and bud',
		'ER,ambiguous': ' -> ER,ambiguous',
		'bud neck,cell periphery,vacuole': ' -> bud neck and cell periphery and vacuole',
		'ambiguous,punctate composite,bud': ' -> punctate composite and bud',
		'ambiguous,nuclear periphery': 'nuclear periphery'}

#%%
Den_ycd_map_dict = {'cytoplasm,punctate': 'cytoplasm and punctate',
'cytoplasm': 'cytoplasm',
'nothing': 'nothing',
'nothing,nucleus': 'nucleus',
'bud,vacuole,punctate': 'bud and vacuole and punctate',
'cytoplasm,bud,punctate': 'cytoplasm and bud and punctate',
'cytoplasm,nucleus': 'cytoplasm and nucleus',
'nucleus': 'nucleus',
'nothing,cytoplasm': 'cytoplasm',
'nothing,cytoplasm,nucleus': 'cytoplasm and nucleus',
'cytoplasm,bud': 'cytoplasm and bud',
'punctate,actin/spindle': 'punctate and actin/spindle',
'cytoplasm,ER,punctate': 'cytoplasm and ER and punctate',
'nucleolus': 'nucleolus',
'nothing,cytoplasm,bud,punctate': 'cytoplasm and bud and punctate',
'nothing,bud': 'bud',
'cytoplasm,nucleus,bud,punctate': 'cytoplasm and nucleus and bud and punctate',
'cytoplasm,nucleus,unclassified': 'cytoplasm and nucleus',
'mitochondrion': 'mitochondrion',
'cytoplasm,unclassified': 'cytoplasm',
'actin/spindle': 'actin/spindle',
'cytoplasm,cell periphery,punctate': 'cytoplasm and cell periphery and punctate',
'nothing,cytoplasm,punctate': 'cytoplasm and punctate',
'cytoplasm,nucleus,punctate,actin/spindle': 'cytoplasm and nucleus and punctate and actin/spindle',
'cytoplasm,bud,punctate,actin/spindle': 'cytoplasm and bud and punctate and actin/spindle',
'nothing,cytoplasm,nuclear periphery': 'cytoplasm and nuclear periphery',
'nothing,cytoplasm,punctate,actin/spindle,unclassified': 'cytoplasm and punctate and actin/spindle',
'nothing,punctate': 'punctate',
'cytoplasm,nucleus,nucleolus,punctate': 'cytoplasm and nucleus and nucleolus and punctate',
'punctate': 'punctate',
'bud': 'bud',
'nothing,cytoplasm,unclassified': 'cytoplasm',
'mitochondrion,punctate,unclassified': 'mitochondrion and punctate',
'cytoplasm,mitochondrion': 'cytoplasm and mitochondrion',
'nucleus,punctate': 'nucleus and punctate',
'bud,unclassified': 'bud',
'cytoplasm,ER': 'cytoplasm and ER',
'cytoplasm,nucleus,vacuole': 'cytoplasm and nucleus and vacuole',
'vacuole': 'vacuole',
'cytoplasm,nucleus,nuclear periphery': 'cytoplasm and nucleus and nuclear periphery',
'cytoplasm,nuclear periphery,vacuole': 'cytoplasm and nuclear periphery and vacuole',
'nuclear periphery': 'nuclear periphery',
'structured,punctate': 'structured and punctate',
'cell periphery,vacuole,punctate,unclassified': 'cell periphery and vacuole and punctate',
'cell periphery': 'cell periphery',
'cytoplasm,vacuole': 'cytoplasm and vacuole',
'cell periphery,unclassified': 'cell periphery',
'nucleolus,structured,punctate': 'nucleolus and structured and punctate',
'cell periphery,ER,punctate': 'cell periphery and ER and punctate',
'cytoplasm,cell periphery': 'cytoplasm and cell periphery',
'nuclear periphery,punctate': 'nuclear periphery and punctate',
'nuclear periphery,cell periphery,ER': 'nuclear periphery and cell periphery and ER',
'cytoplasm,nucleus,bud': 'cytoplasm and nucleus and bud',
'nothing,mitochondrion': 'mitochondrion',
'ER': 'ER',
'nuclear periphery,cell periphery,ER,punctate,unclassified': 'nuclear periphery and cell periphery and ER and punctate',
'cell periphery,punctate,unclassified': 'cell periphery and punctate',
'cytoplasm,nucleus,nucleolus': 'cytoplasm and nucleus and nucleolus',
'nuclear periphery,ER': 'nuclear periphery and ER',
'cytoplasm,cell periphery,structured,punctate': 'cytoplasm and cell periphery and structured and punctate',
'cytoplasm,mitochondrion,ER': 'cytoplasm and mitochondrion and ER',
'nuclear periphery,ER,punctate': 'nuclear periphery and ER and punctate',
'cytoplasm,nuclear periphery,ER': 'cytoplasm and nuclear periphery and ER',
'cell periphery,ER,vacuole': 'cell periphery and ER and vacuole',
'nuclear periphery,vacuole': 'nuclear periphery and vacuole',
'nucleus,nucleolus': 'nucleus and nucleolus',
'cytoplasm,nucleus,cell periphery': 'cytoplasm and nucleus and cell periphery',
'cytoplasm,nuclear periphery': 'cytoplasm and nuclear periphery',
'nothing,cytoplasm,mitochondrion': 'cytoplasm and mitochondrion',
'cytoplasm,nucleolus': 'cytoplasm and nucleolus',
'cytoplasm,structured': 'cytoplasm and structured',
'nothing,cytoplasm,structured': 'cytoplasm and structured',
'cell periphery,bud': 'cell periphery and bud',
'cell periphery,ER': 'cell periphery and ER',
'bud,punctate': 'bud and punctate',
'ER,punctate': 'ER and punctate',
'mitochondrion,structured': 'mitochondrion and structured',
'nucleus,cell periphery,vacuole': 'nucleus and cell periphery and vacuole',
'nucleolus,nuclear periphery,punctate': 'nucleolus and nuclear periphery and punctate',
'cytoplasm,vacuole,actin/spindle': 'cytoplasm and vacuole and actin/spindle',
'ER,unclassified': 'ER',
'cytoplasm,mitochondrion,unclassified': 'cytoplasm and mitochondrion',
'cytoplasm,cell periphery,vacuole': 'cytoplasm and cell periphery and vacuole',
'nucleolus,cell periphery,vacuole,punctate': 'nucleolus and cell periphery and vacuole and punctate',
'cytoplasm,nucleus,cell periphery,punctate,unclassified': 'cytoplasm and nucleus and cell periphery and punctate',
'cytoplasm,nuclear periphery,vacuole,punctate': 'cytoplasm and nuclear periphery and vacuole and punctate',
'nuclear periphery,cell periphery,vacuole,unclassified': 'nuclear periphery and cell periphery and vacuole',
'nucleus,nuclear periphery,vacuole': 'nucleus and nuclear periphery and vacuole',
'vacuole,punctate,unclassified': 'vacuole and punctate',
'cytoplasm,nuclear periphery,punctate': 'cytoplasm and nuclear periphery and punctate',
'cytoplasm,nucleus,nuclear periphery,punctate': 'cytoplasm and nucleus and nuclear periphery and punctate',
'cell periphery,vacuole,punctate': 'cell periphery and vacuole and punctate',
'nucleus,vacuole,punctate': 'nucleus and vacuole and punctate',
'mitochondrion,punctate': 'mitochondrion and punctate',
'nucleus,cell periphery,punctate,unclassified': 'nucleus and cell periphery and punctate',
'ER,vacuole': 'ER and vacuole',
'unclassified': 'unclassified',
'nucleus,nuclear periphery': 'nucleus and nuclear periphery',
'cytoplasm,vacuole,punctate,actin/spindle': 'cytoplasm and vacuole and punctate and actin/spindle',
'cytoplasm,vacuole,punctate': 'cytoplasm and vacuole and punctate',
'mitochondrion,structured,unclassified': 'mitochondrion and structured',
'nucleolus,nuclear periphery': 'nucleolus and nuclear periphery',
'cytoplasm,nucleus,vacuole,punctate': 'cytoplasm and nucleus and vacuole and punctate',
'cytoplasm,cell periphery,mitochondrion,punctate': 'cytoplasm and cell periphery and mitochondrion and punctate',
'nuclear periphery,cell periphery,punctate,unclassified': 'nuclear periphery and cell periphery and punctate',
'cytoplasm,mitochondrion,punctate': 'cytoplasm and mitochondrion and punctate',
'cytoplasm,nuclear periphery,ER,unclassified': 'cytoplasm and nuclear periphery and ER',
'vacuole,punctate,actin/spindle': 'vacuole and punctate and actin/spindle',
'cytoplasm,cell periphery,vacuole,punctate': 'cytoplasm and cell periphery and vacuole and punctate',
'mitochondrion,unclassified': 'mitochondrion',
'nucleus,vacuole,unclassified': 'nucleus and vacuole',
'cytoplasm,nucleus,punctate,unclassified': 'cytoplasm and nucleus and punctate',
'cytoplasm,nucleus,nuclear periphery,unclassified': 'cytoplasm and nucleus and nuclear periphery',
'nuclear periphery,cell periphery,ER,unclassified': 'nuclear periphery and cell periphery and ER',
'nucleus,nucleolus,punctate': 'nucleus and nucleolus and punctate',
'cytoplasm,nucleus,nuclear periphery,ER,unclassified': 'cytoplasm and nucleus and nuclear periphery and ER',
'cytoplasm,mitochondrion,structured': 'cytoplasm and mitochondrion and structured',
'structured': 'structured',
'nucleus,unclassified': 'nucleus',
'cytoplasm,nucleus,cell periphery,vacuole': 'cytoplasm and nucleus and cell periphery and vacuole',
'cell periphery,bud,punctate': 'cell periphery and bud and punctate',
'cytoplasm,punctate,actin/spindle': 'cytoplasm and punctate and actin/spindle',
'nucleus,nuclear periphery,mitochondrion,structured': 'nucleus and nuclear periphery and mitochondrion and structured',
'nothing,cytoplasm,cell periphery': 'cytoplasm and cell periphery',
'cytoplasm,bud,ER,punctate': 'cytoplasm and bud and ER and punctate',
'vacuole,unclassified': 'vacuole',
'vacuole,punctate': 'vacuole and punctate',
'nuclear periphery,ER,unclassified': 'nuclear periphery and ER',
'structured,ER,punctate,unclassified': 'structured and ER and punctate',
'ER,punctate,actin/spindle': 'ER and punctate and actin/spindle',
'cytoplasm,punctate,unclassified': 'cytoplasm and punctate',
'nucleus,nuclear periphery,unclassified': 'nucleus and nuclear periphery',
'cell periphery,bud,ER,punctate': 'cell periphery and bud and ER and punctate',
'nuclear periphery,mitochondrion,unclassified': 'nuclear periphery and mitochondrion',
'nucleolus,structured': 'nucleolus and structured',
'cell periphery,punctate': 'cell periphery and punctate',
'nuclear periphery,cell periphery,vacuole,punctate': 'nuclear periphery and cell periphery and vacuole and punctate',
'nucleolus,punctate': 'nucleolus and punctate',
'nuclear periphery,vacuole,punctate': 'nuclear periphery and vacuole and punctate',
'cell periphery,vacuole': 'cell periphery and vacuole',
'nuclear periphery,structured,punctate': 'nuclear periphery and structured and punctate',
'nucleolus,vacuole,punctate': 'nucleolus and vacuole and punctate',
'cytoplasm,nucleus,nucleolus,mitochondrion,unclassified': 'cytoplasm and nucleus and nucleolus and mitochondrion',
'cytoplasm,mitochondrion,ER,punctate,unclassified': 'cytoplasm and mitochondrion and ER and punctate',
'nucleus,punctate,unclassified': 'nucleus and punctate',
'punctate,unclassified': 'punctate',
'cytoplasm,nucleus,nuclear periphery,structured': 'cytoplasm and nucleus and nuclear periphery and structured',
'cytoplasm,nuclear periphery,unclassified': 'cytoplasm and nuclear periphery',
'nucleus,nuclear periphery,punctate': 'nucleus and nuclear periphery and punctate',
'cytoplasm,structured,bud,punctate': 'cytoplasm and structured and bud and punctate',
'cytoplasm,nucleus,punctate': 'cytoplasm and nucleus and punctate',
'cell periphery,vacuole,unclassified': 'cell periphery and vacuole',
'nucleolus,nuclear periphery,cell periphery,ER,punctate,unclassified': 'nucleolus and nuclear periphery and cell periphery and ER and punctate',
'nucleus,cell periphery,vacuole,punctate': 'nucleus and cell periphery and vacuole and punctate',
'nucleolus,nuclear periphery,cell periphery,vacuole,punctate': 'nucleolus and nuclear periphery and cell periphery and vacuole and punctate'}

Brandon_dict = {'Cytoplasm': 'Cytoplasm',
'Nucleolus Irregular': 'Nucleolus Irregular',
'Cytoplasmic foci': 'Cytoplasmic foci',
'Nucleus': 'Nucleus',
'Vacuole': 'Vacuole',
'Bud Neck': 'Bud Neck',
'ER': 'ER',
'Nuclear Foci': 'Nuclear Foci',
'Cytoplasm irreg.': 'Cytoplasm irreg.',
'PM (Punctate)': 'Plasma Membrane',
'Nucleus P': 'Nucleus Periphery',
'Nucleolus': 'Nucleolus',
'Nucleus Irregular': 'Nucleus Irregular',
'Vacuole Foci': 'Vacuole Foci',
'ER Foci': 'ER Foci',
'Cytoplasmic Foci': 'Cytoplasmic Foci',
'Nuclear / Cyto Foci': 'Nuclear and Cyto Foci',
'Nuclear foci': 'Nuclear foci',
'Nuclear Periphery': 'Nuclear Periphery',
'ER foci': 'ER foci',
'Mitchondria': 'Mitchondria'}

Huh_dict = {
'cytoplasm,nucleus': 'cytoplasm,nucleus',
'endosome':'endosome',
'cytoplasm':'cytoplasm',
'ER':'ER',
'mitochondrion':'mitochondrion',
'nucleus':'nucleus',
'early Golgi,late Golgi':'Golgi',
'ambiguous,bud neck,cytoplasm,cell periphery,bud':'bud neck,cytoplasm,cell periphery,bud',
'spindle pole':'spindle pole',
'ambiguous,bud neck,cell periphery,punctate composite':'bud neck,cell periphery,punctate composite',
'late Golgi,bud': 'Golgi',
'punctate composite':'punctate composite',
'cytoplasm,mitochondrion':'cytoplasm,mitochondrion',
'vacuolar membrane':'vacuolar membrane',
'nuclear periphery':'nuclear periphery',
'ambiguous,spindle pole':'spindle pole',
'ambiguous':'ambiguous',
'vacuole':'vacuole',
'cytoplasm,nucleolus':'cytoplasm,nucleolus',
'cytoplasm,actin':'cytoplasm,actin',
'late Golgi':'Golgi',
'punctate composite,endosome':'punctate composite,endosome',
'nucleolus':'nucleolus',
'nucleolus,nucleus':'nucleolus,nucleus',
'cell periphery':'cell periphery',
'microtubule':'microtubule',
'bud neck,cell periphery':'bud neck,cell periphery',
'cell periphery,bud':'cell periphery,bud',
'bud neck,cytoplasm':'bud neck,cytoplasm',
'bud neck,cytoplasm,nucleus':'bud neck,cytoplasm,nucleus',
'cytoplasm,vacuole':'cytoplasm,vacuole',
'bud neck,cytoplasm,cell periphery,bud':'bud neck,cytoplasm,cell periphery,bud',
'early Golgi':'Golgi',
'ambiguous,bud neck,cytoplasm,bud':'bud neck,cytoplasm,bud',
'cell periphery,vacuole':'cell periphery,vacuole',
'Golgi':'Golgi',
'ER,vacuole':'ER,vacuole',
'lipid particle':'lipid particle',
'bud neck,cytoplasm,cell periphery':'bud neck,cytoplasm,cell periphery',
'Golgi,early Golgi':'Golgi',
'ambiguous,bud neck,cell periphery,bud':'bud neck,cell periphery,bud',
'nucleus,spindle pole':'nucleus,spindle pole',
'bud neck,cell periphery,bud':'bud neck,cell periphery,bud',
'ambiguous,cytoplasm,cell periphery,bud':'cytoplasm,cell periphery,bud',
'nucleus,microtubule':'nucleus,microtubule',
'punctate composite,early Golgi':'punctate composite,Golgi',
'actin':'actin',
'cytoplasm,nucleolus,nucleus':'cytoplasm,nucleolus,nucleus',
'bud neck':'bud neck',
'spindle pole,microtubule':'spindle pole,microtubule',
'cytoplasm,punctate composite':'cytoplasm,punctate composite',
'bud neck,cell periphery,punctate composite':'bud neck,cell periphery,punctate composite',
'peroxisome':'peroxisome',
'cytoplasm,mitochondrion,nucleus':'cytoplasm,mitochondrion,nucleus',
'ambiguous,nuclear periphery':'nuclear periphery',
'mitochondrion,nucleus':'mitochondrion,nucleus',
'punctate composite,early Golgi,late Golgi':'punctate composite, Golgi',
'ER to Golgi':'ER,Golgi',
'ambiguous,cytoplasm':'cytoplasm',
'ambiguous,endosome':'endosome',
'vacuole,vacuolar membrane':'vacuole,vacuolar membrane',
'ambiguous,cell periphery,bud':'cell periphery,bud',
'nuclear periphery,nucleus':'nuclear periphery,nucleus',
'cytoplasm,nucleus,spindle pole':'cytoplasm,nucleus,spindle pole',
'ambiguous,cytoplasm,nucleus':'cytoplasm,nucleus',
'ambiguous,mitochondrion':'mitochondrion',
'ER,cell periphery':'ER,cell periphery',
'ambiguous,bud neck,bud':'bud neck,bud',
'ambiguous,cytoplasm,bud':'cytoplasm,bud',
'cytoplasm,nuclear periphery,nucleus':'cytoplasm,nuclear periphery,nucleus',
'ambiguous,punctate composite':'punctate composite',
'ER,cytoplasm':'ER,cytoplasm',
'cytoplasm,late Golgi':'cytoplasm,Golgi',
'punctate composite,Golgi,early Golgi':'punctate composite,Golgi,Golgi',
'cytoplasm,endosome,lipid particle':'cytoplasm,endosome',
'ambiguous,cytoplasm,punctate composite':'cytoplasm,punctate composite',
'punctate composite,vacuolar membrane,lipid particle':'punctate composite,vacuolar membrane',
'cytoplasm,nucleus,punctate composite':'cytoplasm,nucleus,punctate composite',
'cytoplasm,early Golgi':'cytoplasm,Golgi',
'ambiguous,Golgi':'Golgi',
'ambiguous,Golgi,early Golgi':'Golgi,Golgi',
'cytoplasm,nuclear periphery':'cytoplasm,nuclear periphery',
'bud neck,cytoplasm,cell periphery,punctate composite':'bud neck,cytoplasm,cell periphery,punctate composite,bud',
'cytoplasm,cell periphery':'cytoplasm,cell periphery',
'bud neck,cytoplasm,mitochondrion,cell periphery':'bud neck,cytoplasm,mitochondrion,cell periphery',
'mitochondrion,punctate composite':'mitochondrion,punctate composite',
'ambiguous,vacuolar membrane':'vacuolar membrane',
'nucleus,spindle pole,microtubule':'nucleus,spindle pole,microtubule',
'ER,cell periphery,vacuole,bud':'ER,cell periphery,vacuole,bud',
'punctate composite,late Golgi':'punctate composite,Golgi',
'cytoplasm,bud':'cytoplasm,bud',
'ambiguous,late Golgi':'Golgi',
'vacuole,endosome':'vacuole,endosome',
'ambiguous,nucleus':'nucleus',
'punctate composite,Golgi':'punctate composite,Golgi',
'ambiguous,bud neck,cytoplasm,vacuole,bud':'bud neck,cytoplasm,vacuole,bud',
'cytoplasm,spindle pole':'cytoplasm,spindle pole',
'ambiguous,bud neck,cell periphery,vacuole,bud':'bud neck,cell periphery,vacuole,bud',
'cytoplasm,Golgi,early Golgi':'cytoplasm,Golgi,Golgi',
'vacuolar membrane,endosome':'vacuolar membrane,endosome',
'ER,ambiguous,bud':'ER,bud',
'cytoplasm,punctate composite,endosome':'cytoplasm,punctate composite,endosome',
'ER,cytoplasm,nucleus':'ER,cytoplasm,nucleus',
'punctate composite,endosome,early Golgi,late Golgi':'punctate composite,endosome,Golgi,Golgi',
'ER,cell periphery,bud':'ER,cell periphery,bud',
'cytoplasm,endosome':'cytoplasm,endosome',
'ER,nucleus':'ER,nucleus',
'cytoplasm,nucleus,bud':'cytoplasm,nucleus,bud',
'punctate composite,actin':'punctate composite,actin',
'cytoplasm,nuclear periphery,punctate composite':'cytoplasm,nuclear periphery,punctate composite',
'bud neck,cytoplasm,bud':'bud neck,cytoplasm,bud',
'nucleolus,microtubule':'nucleolus,microtubule',
'punctate composite,spindle pole':'punctate composite,spindle pole',
'cytoplasm,punctate composite,spindle pole,microtubule':'cytoplasm,punctate composite,spindle pole,microtubule',
'vacuole,Golgi':'vacuole,Golgi',
'bud neck,nucleus':'bud neck,nucleus',
'ER,ambiguous,cytoplasm':'ER,cytoplasm',
'endosome,lipid particle':'endosome',
'ambiguous,punctate composite,bud':'punctate composite,bud',
'ER,mitochondrion,nuclear periphery':'ER,mitochondrion,nuclear periphery',
'punctate composite,lipid particle':'punctate composite',
'ER,ambiguous':'ER',
'ambiguous,bud neck,cytoplasm,nucleus,bud':'bud neck,cytoplasm,nucleus,bud',
'bud neck,cell periphery,vacuole':'bud neck,cell periphery,vacuole',
'cytoplasm,cell periphery,bud':'cytoplasm,cell periphery,bud'}

#%%

class dataset_desc:
	def __init__(self, merge_on, information) -> None:
		self.merge_on = merge_on
		self.information = information
		self.func_match = f"{self}_f"

class microcope_info():
	def __init__(self) -> None:
		self.pixel_ratio_microns = 0.1081

# class GO_term:
# 	def __init__(self, go_list) -> None:
# 		self.listed = str.split(go_list, ", ")


def f_GO_map(micro_map): #! This was left out because decided to display based on a network rather than colored scatter
	go = pd.read_excel("C:\\Users\\pcnba\\Grant Brown's Lab Dropbox\\Peter Bartlett\\Peter Bartlett Data\\Code\\Data_copies\\Information_files\\Localization_merging\\GO_Proteins.xlsx")
	terms = go['TERM']

	# def add_GOs(prot:str, go_group_prot:str, go_group_name:str, go_group_list:list): #! This is not in use because deicided to use a network graph to visualize interaction enrichment rather than as scatter with variable grouped GO color
	#. This below code has not been micro_maped
	# 	if prot in go_group_list:
	# 		go_group_prot + "," + go_group_name
	# 	else:
	# 		pass
	# 	return(go_group_prot)

	# micro_map['GO_group_collected'] = ''
	# for r in terms:
	# 	go_matches = go.loc[r, 'ANNOTATED_GENES']
	# 	micro_map['GO_group_collected'] = micro_map.apply(lambda x: add_GOs(x['Protein'], x['GO_group_collected'], r, go_matches), axis = 1)
	return(go) #. , micro_map)

def sgd_map_f():
	# sgd = pd.read_csv("results_best.csv")
	sgd = pd.read_csv("results.tsv", sep = '\t')
	# sgd.rename(columns={'input':'Gene_Standard_Name'}, inplace=True)
	# sgd.rename(columns={'length':'Gene_Length'}, inplace=True)
	sgd.columns = sgd.columns.str.replace(" > ", "_")
	sgd.columns = sgd.columns.str.replace(" ", "_")
	# sgd = sgd.loc[:, ["Gene_Systematic_Name", "Gene_Standard_Name", "Gene_Name", "Gene_Length", "Gene_Phenotype_Summary", "Gene_Length"]]
	sgd['Gene_Standard_Name'] = sgd['Gene_Standard_Name'].fillna(sgd['Gene_Systematic_Name'])
	sgd = sgd.set_index("Gene_Systematic_Name")
	#* global info_sgd
	#* info_sgd = dataset_desc('Standard_Name','information')
	return(sgd)

def tkach_f():
	tkach = pd.read_excel("Tcak_protein_localization.xlsx", sheet_name='Calls')
	tkach.columns = tkach.columns.str.replace(" ", "_")
	tkach['Standard_Name'] = tkach['Standard_Name'].fillna(tkach['Systematic_ORF'])
	tkach = tkach.set_index("Systematic_ORF")

	tkach['EndLOC_Rescreen_MMS_Tcak'] = tkach['EndLOC_Rescreen_MMS_Tcak'].map(Tkach_dictionary).fillna(tkach['EndLOC_Rescreen_MMS_Tcak'])
	tkach['EndLOC_Rescreen_HU_Tcak'] = tkach['EndLOC_Rescreen_HU_Tcak'].map(Tkach_dictionary).fillna(tkach['EndLOC_Rescreen_HU_Tcak'])
	#* global info_tkach
	#* info_tkach = dataset_desc('Standard_Name','localization')
	return(tkach)

def microfluidics_map_f():
	microfluidics_map = pd.read_excel("MicrofluidicsMap_wCol_USE.xlsx", sheet_name="ProteinLocations")
	microfluidics_map.dropna(subset='Protein', inplace=True)
	microfluidics_map.columns = microfluidics_map.columns.str.replace(' ', '_')
	microfluidics_map['Protein'] = microfluidics_map['Protein'].str.upper()
	microfluidics_map = microfluidics_map.set_index('Protein')
	#* global info_microfluidcs
	#* info_microfluidcs = dataset_desc('Standard_Name/Mix','map')
	return(microfluidics_map)

def denervaud_ycd_f():
	denervaud_ycd = pd.read_excel("Den_data_bestgood.xlsx", sheet_name='Sheet1')
	# denervaud_ycd.rename(columns={'denervaud_ycd_Call':'Call'}, inplace=True)
	denervaud_ycd.drop(denervaud_ycd.columns[denervaud_ycd.columns.str.contains('Unnamed',case = False)],axis = 1, inplace = True)
	# denervaud_ycd = denervaud_ycd.add_suffix("_denervaud_ycd")
	# denervaud_ycd['Standard_Name_denervaud_ycd'] = denervaud_ycd['Standard_Name_denervaud_ycd'].str.upper()
	denervaud_ycd.fillna('unclassified', inplace = True)

	denervaud_ycd['initial_localization'] = denervaud_ycd['initial_localization'].map(Den_ycd_map_dict).fillna(denervaud_ycd['initial_localization'])
	denervaud_ycd['end_localization'] = denervaud_ycd['end_localization'].map(Den_ycd_map_dict).fillna(denervaud_ycd['end_localization'])
	denervaud_ycd.sort_values(by = "movieTag", inplace = True) #* Put in order so that the best movie is first before taking the first instance of an ORF label
	denervaud_ycd['geneName'] = denervaud_ycd['geneName'].replace({'-': np.nan})
	denervaud_ycd['geneName'] = denervaud_ycd['geneName'].fillna(denervaud_ycd['yORF'])
	denervaud_ycd = denervaud_ycd.groupby('yORF').aggregate(lambda x: x.iloc[0])
	denervaud_ycd = denervaud_ycd.drop(columns=['geneName', 'exp_cond', 'movieTag'])
	return(denervaud_ycd)

# def denervaud_f():
# 	denervaud = pd.read_excel("Denervaud Calls.xlsx", sheet_name='Sheet1')
# 	denervaud.columns = denervaud.columns.str.replace(' ', '_')
# 	denervaud.rename(columns={'Denervaud_Call':'Call'}, inplace=True)
# 	denervaud.drop(denervaud.columns[denervaud.columns.str.contains('Unnamed',case = False)],axis = 1, inplace = True)
# 	denervaud = denervaud.add_suffix("_Denervaud")
# 	denervaud['Standard_Name_Denervaud'] = denervaud['Standard_Name_Denervaud'].str.upper()

# 	denervaud['Call_Denervaud'] = denervaud['Call_Denervaud'].map(Denervaud_dictionary).fillna(denervaud['Call_Denervaud'])
# 	denervaud = denervaud.set_index('Standard_Name_Denervaud')
# 	#* global info_Denervaud
# 	#* info_Denervaud = dataset_desc('Standard_Name','Localization')
# 	return(denervaud)

def Ho_loc_pen_f():
	Ho_loc = pd.read_excel("Loc_Ho_SuppT5.xlsx", sheet_name="20210809_HUMMS_penetrance")
	Ho_loc.drop(Ho_loc.columns[Ho_loc.columns.str.contains('Unnamed',case = False)],axis = 1, inplace = True)
	Ho_loc.drop(Ho_loc.columns[Ho_loc.columns.str.contains('HU',case = False)],axis = 1, inplace = True)
	Ho_loc.columns = Ho_loc.columns.str.replace(" ", "_")
	Ho_loc = Ho_loc.add_suffix("_Ho")
	Ho_loc = Ho_loc.set_index('Gene_Ho')
	#* global info_Ho
	#* info_Ho = dataset_desc('Standard Name','LocPen')
	return(Ho_loc)

def Mazumder_f():
	mazumder = pd.read_excel("Mazumder_ver2.xlsx", sheet_name="Mod_dest")
	# mazumder.columns = mazumder.columns.str.replace(" ", "_")
	mazumder = mazumder.add_suffix("_Mazumder")
	mazumder['Localization_Mazumder'] = mazumder['Localization_Mazumder'].map(Mazumder_dictionary).fillna(mazumder['Localization_Mazumder'])
	mazumder['CommName_Mazumder'] = mazumder['CommName_Mazumder'].fillna(mazumder['ORF_Mazumder'])
	mazumder = mazumder.set_index('ORF_Mazumder')
	#* global info_Mazumder
	#* info_Mazumder = dataset_desc('ORF','Localization')
	return(mazumder)

def Huh_f():
	Huh = pd.read_csv('Huh2003.txt', sep = '\t').set_index('yORF')
	Huh['localization summary'] = Huh['localization summary'].map(Huh_dict).fillna(Huh['localization summary'])
	return(Huh)


#%%
# def origin_destination(directional:str):
# 	orig_dest = directional.str.split(by = ' -> ')
# 	return(orig_dest) #! must be unpacked

# def merge_manager():
# 	os.chdir(Global_Variables['information_path'])
# 	files = input("What are the files to be micro_map based on protein? [Comma deliminate]").split(', ')
# 	for f in files: #TODO: Add multi-read_functionality
# 		pd.read_csv(f).set_index()
# 	micro_map = pd.DataFrame([])
# 	file_n = 0
# 	for f in files:
# 		exec(f'Global {f}')
# 		exec(f"temp = {f}_f({f})") #! This isn't good form, but will work for now
# 		if file_n > 0:
# 			micro_map = pd.merge(micro_map, temp, left_index = True, right_index = True)
# 			file_n += 1
# 		else:
# 			micro_map = temp.copy()
# 			file_n += 1

#%%
def control_replace(x,search):
	ic(x, search)
	if search.upper() in x.upper():
		return(search.upper())
	else:
		return(x.upper())
#%%
if __name__ == "__main__":
	# Global_Variables = glv.global_manager()
	Global_Variables = {
	# 	"analyze": "F:/Microfluidics/Missing_Analyze2",
		"microfluidics_results": "F:/Microfluidics/RES_N_ULTS",
		"information_path": "C:/Users/pcnba/Grant Brown's Lab Dropbox/Peter Bartlett/Peter Bartlett Data/Code/Data_copies/Information_files/Localization_merging",
		"post_path": "D:/ALL_FINAL"} # * ,
	# 	"subset": False,
	# 	'subset_by': 'range',
	# 	'subset_collection': '',
	# 	"cpu_se": 16,
	# 	"timepoint_gap": 7.5,
	# 	"percentiles": [95, 99],
	# 	"multiplex": True}
	os.chdir(Global_Variables['information_path'])
	# files = input("What are the files to be micro_map based on protein? [Comma deliminate]").split(', ')
	# for f in files: #TODO: Add multi-read_functionality
		# pd.read_csv(f).set_index()

	# micro_map = microfluidics_map_f()
	# # .drop(columns=['Predicted_localization_Change', 'Notes', 'Current_Stage', 'Location', 'Fullmicro_map'])
	# map_drop_columns = ['Predicted_localization_Change', 'Notes', 'Current_Stage', 'Location', 'Fullmicro_map']
	# micro_map.drop(micro_map.columns[micro_map.columns.str.contains('Unnamed',case = False)],axis = 1, inplace = True)
	# micro_map.drop(columns=[col for col in micro_map if col in map_drop_columns], inplace=True)

	sgd = sgd_map_f()
	# micro_map = pd.merge(sgd, micro_map, left_index=True, right_index=True, how = 'left')
	#. Decided that will just use the list of proteins from SGD to merge with the other datasets

	tkach = tkach_f()
	# denervaud = denervaud_f()
	mazumder = Mazumder_f()
	Den_ycd_map_df = denervaud_ycd_f()
	Brandons_map = pd.read_excel("Brandons_Paper.xlsx", sheet_name="Sheet1").set_index('ORF')
	Brandons_map = Brandons_map.drop(columns=['Protein']).rename(columns={'Subcellular Compartment Re-localization': 'Dest_Call'}).add_suffix('_Brandons')
	Huh = Huh_f()

	# loqate = pd.read_excel('proteomesummarylamicro_mapversion.xlsx', sheet_name='Sheet1', usecols=['ORF', 'control Localization']).set_index('ORF').replace('below threshold', np.nan)




	micro_map = sgd.merge(Den_ycd_map_df, left_index= True, right_index= True, how = 'left')
	micro_map = micro_map.merge(tkach, left_index=True, right_index = True, how= "left")

	#Either should work below
	micro_map = micro_map.merge(mazumder, left_index = True, right_index = True, how = "left")
	# micro_map = micro_map.merge(mazumder, left_on = 'Gene_Standard_Name', right_on = 'CommName_Mazumder', how = "left")

	micro_map = pd.merge(micro_map, Brandons_map, right_index=True, left_index=True, how= 'left')
	#Artifact of removed micorfluidics map
	# micro_map = micro_map.sort_values(by = ['Date', 'Run_Number', 'MapID_(Col_Range)'])

	micro_map = micro_map.merge(Huh, left_index=True, right_index=True, how='left')
	# micro_map = micro_map.merge(loqate, left_index=True, right_index=True, how='left')

#%%
#. These are all the protiens which do not have matching localization information
#Check the rows where Na for 'localization_change', 'initial_localization', 'end_localization', 'EndLOC_Rescreen_MMS_Tcak', 'EndLOC_Rescreen_HU_Tcak', Localization_Mazumder', 'Dest_Mazumder', 'Function_Mazumder'
micro_map.loc[micro_map[['localization_change', 'initial_localization', 'end_localization','EndLOC_Rescreen_MMS_Tcak', 'EndLOC_Rescreen_HU_Tcak', 'Localization_Mazumder', 'Dest_Mazumder', 'Dest_Call_Brandons']].isna().all(axis = 1),['Gene_Standard_Name', 'localization_change', 'initial_localization', 'end_localization','EndLOC_Rescreen_MMS_Tcak', 'EndLOC_Rescreen_HU_Tcak', 'Localization_Mazumder', 'Dest_Mazumder']]

#%%
micro_map[['StartLOC_Rescreen_MMS_Tcak', 'EndLOC_Rescreen_MMS_Tcak']] = micro_map['EndLOC_Rescreen_MMS_Tcak'].str.split(' -> ', expand=True).replace('', np.nan, regex=True)
micro_map[['StartLOC_Rescreen_HU_Tcak', 'EndLOC_Rescreen_HU_Tcak']] = micro_map['EndLOC_Rescreen_HU_Tcak'].str.split(' -> ', expand=True).replace('', np.nan, regex=True)

#%%
# micro_map['EndLOC_Rescreen_HU_Tcak'] = micro_map['EndLOC_Rescreen_HU_Tcak'].str.replace(' -> ', '').replace('', np.nan, regex=True)
# micro_map['EndLOC_Rescreen_HU_Tcak'] = micro_map['EndLOC_Rescreen_HU_Tcak'].str.replace('-> ', '').replace('', np.nan, regex=True)

# micro_map['EndLOC_Rescreen_MMS_Tcak'] = micro_map['EndLOC_Rescreen_MMS_Tcak'].str.replace(' -> ', '').replace('', np.nan, regex=True)
# micro_map['EndLOC_Rescreen_MMS_Tcak'] = micro_map['EndLOC_Rescreen_MMS_Tcak'].str.replace('-> ', '').replace('', np.nan, regex=True)

#%%
def simp_pref(*args):
	for arg in args:
		if pd.isna(arg):
			continue
		else:
			return arg
	return np.nan

micro_map['Combined_destination'] = micro_map.apply(lambda x: simp_pref(x['Dest_Call_Brandons'], x['end_localization'], x['EndLOC_Rescreen_MMS_Tcak'], x['EndLOC_Rescreen_HU_Tcak'], x['Dest_Mazumder']), axis=1)
micro_map['Combined_origin'] = micro_map.apply(lambda x: simp_pref(x['StartLOC_Rescreen_MMS_Tcak'], x['StartLOC_Rescreen_HU_Tcak'], x['localization summary']), axis=1)

micro_map['Combined_destination'] = micro_map['Combined_destination'].fillna('unclassified').str.replace(',', ' and ').str.title().str.strip()
micro_map['Combined_origin'] = micro_map['Combined_origin'].fillna('unclassified').str.replace(',', ' and ').str.title().str.strip()
#%%
micro_map[['Single_origin', 'Double_origin', 'Triple_origin', 'Quad_origin']] = micro_map['Combined_origin'].str.split(' And ', expand=True)
micro_map['Double_origin'] = micro_map['Single_origin'] + " and" + micro_map['Double_origin']
micro_map['Triple_origin'] = micro_map['Single_origin'] + " and" + micro_map['Double_origin'] + " and" + micro_map['Triple_origin']
micro_map['Quad_origin'] = micro_map['Single_origin'] + " and" + micro_map['Double_origin'] + " and" + micro_map['Triple_origin'] + " and" + micro_map['Quad_origin']

micro_map[['Single_destination', 'Double_destination', 'Triple_destination', 'Quad_destination']] = micro_map['Combined_destination'].str.split(' And ', expand=True)
micro_map['Double_destination'] = micro_map['Single_destination'] + " and" + micro_map['Double_destination']
micro_map['Triple_destination'] = micro_map['Single_destination'] + " and" + micro_map['Double_destination'] + " and" + micro_map['Triple_destination']
micro_map['Quad_destination'] = micro_map['Single_destination'] + " and" + micro_map['Double_destination'] + " and" + micro_map['Triple_destination'] + " and" + micro_map['Quad_destination']


# reference a pandas column and split into required columns by delimiter



# micro_map[['Primary']] = micro_map['Combined_origin'].str.split(' And ', expand=True).replace('', np.nan, regex=True)

#%%
micro_map.to_parquet('protein_origin_dest.parquet')

#%%
# Read in the percentages and merge with the micro_map
penetrances = pd.read_parquet("D:\ALL_FINAL\Final_combined_comparison.parquet")
penetrances = penetrances.merge(micro_map, left_on = 'Protein', right_on='Gene_Standard_Name', how='left')
yet_percentiles = pd.read_parquet("D:/ALL_FINAL/Combined_by_perc/new_percs.parquet")
updated_yet_perc = yet_percentiles.groupby('Protein').Yet_perc.agg('max').rename('updated_yet_perc')

penetrances = penetrances.merge(updated_yet_perc, left_on = 'Protein', right_on='Protein', how='left') #* This has extra information that is not needed

#%%
Ho = Ho_loc_pen_f()
Ho_agg = Ho.agg('max', axis= 1)*100 #*Convert the decimal to a percentage
Ho_agg = Ho_agg.rename('Ho_max').to_frame()
merged_frame_pens = penetrances.merge(Ho_agg, left_on = 'Protein', right_index = True, how = 'left')
merged_frame_pens.to_parquet('D:/ALL_FINAL/Combined_by_perc/penetrance_updated.parquet')

#%%



#%%
#> This is an older version that was changed to the above
# small = micro_map[['EndLOC_Rescreen_MMS_Tcak','Call_Denervaud','Localization_Mazumder']]
# def combine_logic(tkach, denervaud, mazumder):
# 	if len(tkach) > 0:
# 		t = True

# 	else: t = False

# 	if len(denervaud)>0:
# 		d = True
# 	else: d = False

# 	if len(mazumder) > 0:
# 		m = True
# 	else: m = False

# 	if t == True and d == True and m == True:
# 		if tkach == denervaud == mazumder:
# 			combo = tkach
# 		else:
# 			combo = tkach + denervaud + mazumder

# 	elif t == True:
# 		if d == True and m == False:
# 			if tkach == denervaud:
# 				combo = tkach
# 			else:
# 				combo = tkach + denervau
# 		elif d == False and m == True:
# 			if tkach == mazumder:
# 				combo = tkach
# 			else:
# 				combo = tkach + denervaud
# 		else:
# 			combo = tkach

# 	elif d == True and t == False:
# 		if m == True:
# 			if denervaud == mazumder:
# 				combo = mazumder
# 			else:
# 				combo = denervaud + mazumder
# 		else:
# 			combo = denervaud

# 	elif m == True and t == False and d == False:
# 		if mazumder == denervaud:
# 			combo = mazumder
# 		else:
# 			combo = denervaud + mazumder

# 	if t == False and d == False and m == False:
# 		return(None)

# 	return(combo)


# small['combined'] = small.apply(lambda x: combine_logic(tkach=x['EndLOC_Rescreen_MMS_Tcak'], denervaud= x['Call_Denervaud'], mazumder= x['Localization_Mazumder']), axis = 0)

# def join_strings(s):
	# return ''.join(str(i) for i in s if not pd.isna(i))

# %%
#> This is to match the DDC2 to FLR1 in library
#* The search and join for control should be done last once all the merging has been performed
search = input("What is the name of the control protein? [Esc to skip]")
if search == '':
	pass
else:
	micro_map.reset_index(inplace=True, drop= False)
	micro_map['Protein'] = pd.Series(micro_map['Protein']).apply(lambda x: control_replace(x = x, search = search))
	micro_map.set_index('Protien')
