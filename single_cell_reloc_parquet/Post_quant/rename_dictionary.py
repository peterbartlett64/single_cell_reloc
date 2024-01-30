Tkach_dictionary = {'to cyto': '-> cytoplasm',
       'to nuc': ' -> nucleus',
       'nucleus': ' -> nucleus',
       'nuc foci': ' -> nuclear foci',
       'cyto foci': ' -> cyto foci',
       'from budneck/tip': 'budneck ->',
       'from bud tip': 'budneck ->',
       'to pm': ' -> plasma membrane',
       'to bud neck': ' -> budneck',
       'to nuc/periphery': ' -> nucleus and nuclear periphery',
       'to nuc (from nuc foci)': 'nuclear foci -> nucleus',
       'nuc foci': ' -> nuclear foci',
       'from bud neck': 'budneck -> ',
       'to vac': ' -> vacuole',
       'to nuc, nuc foci, cyto foci' : '-> nucleus and nuclear foci and cytoplasm', #* This one is hard to deal with. Almost decided to do a leveled assignment
       'no cells': '',
       'to pm foci': ' -> pm foci',
       'to nucleolus, cyto foci' : ' -> nuceolus and cyto foci',
       'to vac (from pm)': 'plasma membrane -> vacuole',
       'to nucleolus': ' -> nucleolus',
       'to nuc periph (cyto?)': 'cytoplasm -> nuclear periphery',
       'nuc foci (weak)': ' -> nuclear foci',
       'to cyto (from pm/endosome)': 'plasma membrane -> cytoplasm',
       'to nuc (from nucleolus)': 'nucleolus -> nucleus',
       'er foci': '-> ER foci',
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
       'not the same strain': '',
       'to er': 'ER',
       'from budneck': 'budneck -> ',
       'nuc periph' : ' -> nuclear periph',
       'to pm (foci)': ' -> pm foci',
       'to vac (from vac mb)': 'vacuole -> vacuole',
       'from bud neck (to pm)': 'budneck -> plamsa membrane',
       'to diffuse nuc (from foci)': 'nuclear foci -> nucleus',
       'to er foci (weak)': 'ER foci',
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

Mazumder_aternate_dictionary = {

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
'ambiguous,nuclear periphery': ' -> nuclear periphery'}