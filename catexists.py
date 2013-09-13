import sys,os
import subprocess
import numpy as np
import pylab as plt
import platform

if platform.node() == "bigbang.mit.edu":
    basepath = "/bigbang/data/AnnaGroup/caterpillar/halos"

if platform.node() == "rclogin13.rc.fas.harvard.edu":
    basepath = "/n/home01/bgriffen/data/caterpillar/halos"

levellist = [11,12,13,14,15]
nrvirlist = [3,4,5,6]
haloidlist = []
for filename in os.listdir(basepath):
    if filename[0] == "H":
        haloidlist.append(filename)

fig = plt.figure(figsize=(16,10))
plotinc = 0

for haloid in haloidlist:
    plotinc += 1   
    ax = fig.add_subplot(3,5,plotinc)
    ax.set_title(haloid)
    ax.set_ylim([2,7])
    ax.set_xlim([10,16])
    ax.set_yticks((3,4,5,6))
    ax.set_yticklabels(('3','4','5','6'))
    ax.set_xticks((11,12,13,14,15))
    ax.set_xticklabels(('11','12','13','14','15'))
    print "Searching:",haloid
    for level in levellist:
        for nrvir in nrvirlist:
            ext = haloid + "_BB_Z127_P7_LN7_LX" + str(level) + "_O4_NV" + str(nrvir) + "/"
            corepath =  basepath + "/" + haloid + "/" + ext
            try:
                subdirnames = basepath + "/" + haloid
                for subname in os.listdir(subdirnames)
                    halosubdirlist.append(filename)
                    fileparts =  filename.split("_")
                    print fileparts
                    #levelmax = float(fileparts[5][2:4]
                        
                with open(corepath + "ics.0"):
                    icfilesize = os.path.getsize(corepath + "ics.0")
                    if icfilesize > 0:
                        marker = 'rD'
    
                    if os.path.isdir(corepath + "outputs/snapdir_063"):
                        marker = 'k^'
    
                    if os.path.isdir(corepath + "outputs/group_063"):
                        marker = 'bo'
    
                    ax.plot(int(level),int(nrvir),marker)

            except IOError: pass


ICArtist = ax.plot((0,1),(0,0),'rD', label='IC')
GADGETArtist = ax.plot((0,1),(0,0),'k^', label='Run')
HALOArtist = ax.plot((0,1),(0,0),'bo', label='Halos')
handles, labels = ax.get_legend_handles_labels()
fig.legend(handles, labels, 'upper center',ncol=3,numpoints=1)
        
plt.show()

