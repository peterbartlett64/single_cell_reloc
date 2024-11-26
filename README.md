This README will walk you through the file standards to use the microfluidics pipeline.
All explanations given should be dated for pipeline protocol changes (DD/MM/YY)

----Comment structure---
Code is commented throughout using the better comments extension for VScode. A copy of the markup types is included in the extra documentation.
Refer to "better-comments.tags": [


	]
#! These are possible issues
#* Description of single line. THis is the most commmon comment type throughout pipeline
#< Start of related loop or logic
#> End of related loop or logic
#%% Delimiter Jupyter notebook cell
#? Things which could be changed/unsur
#// Things that were changed, but could possibly rever
#, Description of function
#todo <-
#- Comment element in list of items
#. This is very rarely used comment type


---Image Aquisition Assumptions ---
Files will be of the format: {"positon""[col 3dig eg. "010"][row 3 dig eg. "200"]"_time"[time 4dig eg. "0006"]}

--- Creating the Cellinformation.(csv/xlsx) file ---
Files should be setup as:
	Date |	Day	| Position ID |	Run Number	|Genotype/Condition element	|	Col	|	Row|
	NOTE: for Genotype/Condition element, there can be multiple elements if doing multiplexing.
	If above is true, you must set the color channels in the "flourescent dictrionary"

--- Python subprocesses ---
In the pipeline file, there both Windos and Unix comand line suprocesses to run Cellx.
Comment out the opposite os based on needs.

***********************************************
The pipeline uses quite a few python extension.
Before running the code, make sure that all programs and python extensions are up to date.
Listed below are the last used versions:

-------------
CellX version
Java [1.6.08_45] (THis is very important. Program likely will not run without the proper version.)

Anaconda 3.8.5 (Python 3):
pandas []
#datetime []
PIL []
numpy []
scipy []
matplotlib []
seaborn []
cv2 []
scipy.io []
PIL []
glob []
os[]
ntpath []
csv []
statistics []
numba []
pymongo []
json []
multiprocessing []
concurrent.futures []
subprocess []
time []

----------- For increased readability, code was written using the better-comments extension in VSCode. These were my settings.

"better-comments.tags": [
        {
            "tag": "!",
            "color" : "#e93f6f",
            "strikethrough": false,
            "underline": false,
            "backgroundColor": "transparent",
            "bold": false,
            "italic": false
        },
        {
            "tag": "<",
            "color" : "#0f9d58",
            "strikethrough": false,
            "underline": true,
            "backgroundColor": "opaque",
            "bold": true,
            "italic": false
        },
        {
            "tag": ">",
            "color" : "#0f9d58",
            "strikethrough": false,
            "underline": true,
            "backgroundColor": "opaque",
            "bold": true,
            "italic": false
        },
        {
            "tag": "%%",
            "color": "#a66ad7",
            "strikethrough": false,
            "underline": true,
            "backgroundColor": "transparent",
            "bold": true,
            "italic": true
        },
        {
            "tag": "?",
            "color": "#3498DB",
            "strikethrough": false,
            "underline": false,
            "backgroundColor": "transparent",
            "bold": false,
            "italic": false
        },
        {
            "tag": ".",
            "color": "#FF2D00",
            "strikethrough": false,
            "underline": false,
            "backgroundColor": "transparent",
            "bold": true,
            "italic": false
        },
        {
            "tag": "//",
            "color": "#f4a462",
            "strikethrough": true,
            "underline": false,
            "backgroundColor": "transparent",
            "bold": false,
            "italic": false
        },
        {
            "tag": ",",
            "color": "#ffffff",
            "strikethrough": false,
            "underline": false,
            "backgroundColor": "transparent",
            "bold": true,
            "italic": true
        },
        {
            "tag": "todo",
            "color": "#FF8C00",
            "strikethrough": false,
            "underline": false,
            "backgroundColor": "transparent",
            "bold": false,
            "italic": false
        },
        {
            "tag": "-",
            "color": "#f4b400",
            "strikethrough": false,
            "underline": false,
            "backgroundColor": "transparent",
            "bold": false,
            "italic": true
        },
        {
            "tag": "*",
            "color": "#da70d6",
            "strikethrough": false,
            "underline": false,
            "backgroundColor": "transparent",
            "bold": false,
            "italic": false
        }
    ]