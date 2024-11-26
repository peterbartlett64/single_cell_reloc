#Reshape list to screen table:

#%%
import pandas as pd
import os
from single_cell_reloc_parquet.global_functions.global_variables import slash_switch

def convert_to_screen_table(file, root, f_type, sheet_name = None):
	os.chdir(root)
	if f_type == 'csv':
		screen_list = pd.read_csv(file)
	elif f_type == 'xlsx':
		screen_list = pd.read_excel(file, sheet_name = sheet_name)
	else:
		print("Invalid file type")
		return

	screen_list  = pd.Series(screen_list.values.flatten()) #* Flatten the list
	padding = len(screen_list) % 5
	if padding != 0:
		screen_list = screen_list.append(pd.Series([""] * (5 - padding)), ignore_index=True) #* Insert empty values to fill the last row
	screen_list = screen_list.values.reshape(-1, 5) #* Reshape the list to 5 columns
	screen_table = pd.DataFrame(screen_list)

	if f_type == 'csv':
		screen_table.to_csv('ScreenTable_fList.csv', index=False, header=False) #* Save the screen table to a csv file
	elif f_type == 'xlsx':
		screen_table.to_excel('ScreenTable_fList.xlsx', index=False, header=False)


if __name__ == '__main__':
	input_path = slash_switch(input("Press Enter to convert the list to a screen table")).replace('"', '') #* These extra functions are just ot deal with the 'Copy as path' feature in Windows and Windows' path format
	name = os.path.basename(input_path)
	root = os.path.dirname(input_path)
	if name.endswith('.csv'):
		f_type = 'csv'
	elif name.endswith('.xlsx'):
		f_type = 'xlsx'
		input_sheet = input("Enter the sheet name: ")
	convert_to_screen_table(file = name, root = root, f_type = f_type, sheet_name= input_sheet)
# %%
