This README will walk you through the file standards to use the microfluidics pipeline.
All explanations given should be dated for pipeline protocol changes (DD/MM/YY)

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
