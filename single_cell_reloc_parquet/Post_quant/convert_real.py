#
import yaml
#, These are a few function for converting input values to real-world metrics

#* Below is an example. This should still be roughly true with the same lenses and sample depth (in chip)
microscope_info = {'Microscope': 'Ti2',
		'Lens': '60x',
		'Conversion_factor' : 0.1081, #* µm/pixel
		'Unit' : '\N{MICRO SIGN}m'}

def read_microscope_info(path):
	with open(path, 'r') as file:
		microscope_info =  yaml.safeload(file)
	return(microscope_info)

def micro_length(value, analyze_path):
	try:
		microscope_info
	except NameError:
		microscope_info = read_microscope_info
	factor = microscope_info['Conversion_factor']
	micro_length = float(factor) * value
	unit = microscope_info['Unit']
	full_response = f"{micro_length}{unit}"
	os.path.join
	return(micro_length, unit, full_response)

def micro_area(value, analyze_path):
	try:
		microscope_info
	except NameError:
		microscope_info = read_microscope_info
	factor = microscope_info['Conversion_factor']
	micro_length = float(factor) * value
	unit = microscope_info['Unit'] + '\N{Superscript Two}'
	full_response = f"{micro_length}{unit}"
	return(micro_length, unit, full_response)