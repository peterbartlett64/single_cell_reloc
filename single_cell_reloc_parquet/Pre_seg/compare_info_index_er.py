import os
import pandas as pd
import time

def info_index_er():
    count = 0
    info_index = []
    #This loop should find all the image paths for all the days
    for root, dirs, files, in os.walk(os.getcwd()):
        for name in files:
            if name.endswith(".txt") and name.startswith("TrackingResults", 13):
                info_index.append({'Path': os.path.join(root, name)})
                count = count + 1
                print(count, end="\r")
            elif name.endswith(".txt") and name.startswith("TrackingResults", 14):
                info_index.append({'Path': os.path.join(root, name)})
                count = count + 1
                print(count, end="\r")
            elif name.endswith(".txt") and name.startswith("TrackingResults", 15): #if name.endswith(".txt") and name.startswith("cells"): #changed to 15
                info_index.append({'Path': os.path.join(root, name)})
                count = count + 1
                print(count, end="\r")
            elif name.endswith(".txt") and name.startswith("TrackingResults", 16): #if name.endswith(".txt") and name.startswith("cells"): #changed to 15
                info_index.append({'Path': os.path.join(root, name)})
                count = count + 1
                print(count, end="\r")
            else:
                pass
    del count

    def sync_rem(x):
        if x.startswith(".sync"):
            return(None)
        else:
            return(x)

    def Pos_label(t):
        pstart = t.find("RESULTS")+8 #9 # It is 9 here because "MissingRESULTS2"
        #8 # Must change if "RES_N_ULTS" and to +11
        pend = t.find('GFP_mKO_mKa')-1
        return(t[pstart:pend])

    def Budneck(t):
        bstart = t.find('GFP_mKO_mKa') + 12
        bend = t.find('TrackingResults') - 16
        return(t[bstart:bend])

    def Mod_epoch_time(f):
        return(os.path.getmtime(f))


    def f_non_sync(p):
        substring = ".sync"
        if substring in p:
            return(None)
        else:
            return(p)


    info_index = pd.DataFrame(info_index)
    info_index["Path"] = pd.Series(info_index.iloc[:,0]).apply(f_non_sync)
    info_index.dropna(inplace= True)

    info_index["Pos"] = pd.Series(info_index.iloc[:,0]).apply(Pos_label)
    info_index["Mod_epoch_info"] = pd.Series(info_index.iloc[:,0]).apply(Mod_epoch_time)
    info_index = info_index.loc[info_index.groupby("Pos")["Mod_epoch_info"].idxmax()]
    info_index["Mod_date_info"] = pd.Series(info_index["Mod_epoch_info"]).apply(time.ctime)

    # info_index.to_hdf('MASTER.h5', key= 'info_index', mode='r+')

    info_index.to_csv('info_index.csv', index = False)
    return(info_index)




if __name__ = "__main__":
	singl