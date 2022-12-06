import customtkinter
import os
import segmentation as segment #* This is my custom package to control the segmentation/pass off to snakemake written by Andreas Cuny

Experimental_info = {
    "timepoint_space": time_point_tk,
	"TracX_version": "TracX",
    "data_subset": False,
    "windows_strings": True
    }


#, Initiate the GUI params and the global functions
customtkinter.set_appearance_mode("dark")
customtkinter.set_default_color_theme("green")

root = customtkinter.CTk()
root.geometry("5000x3500") #* This is the size of frame

def slash_switch(path):
    new = path.replace(os.sep, '/')
    return (new)

#REFER to README.txt for setup protocols
pn = os.cpu_count()# to get the cpu core count
print ("Pipeline will run with", pn, "process nodes (number of system cores)")


def slash_switch() #! This fucntion is found in other files and can be brought into this gui master. The regular master will work now

frame = customtkinter.CTkFrame(master=root)
frame.pack(pady = 200, padx = 600, fill = "both", expand = True)

Windows_paths_checkbox = customtkinter.CTkCheckBox(master=frame, text = "Remember")
Windows_paths_checkbox.pack(pady = 12, padx = 10)

analyze_tk = customtkinter.CTkEntry(master= frame, placeholder_text="Where are the images stored [Full path]? (Windows paths accepted)")
analyze_tk.pack(pady = 12, padx = 10)


if Windows_paths_checkbox == 'True':
	analyze_tk = slash_switch(analyze_tk)


segment_button = customtkinter.CTkButton(master = frame, text="Login", command=segment)
segment_button.pack(pady = 12, padx = 10)



microfluidcs_results_tk = customtkinter.CTkEntry(master= frame, placeholder_text="") 

analyze = slash_switch()
microfluidcs_results = 


