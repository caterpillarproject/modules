import sys,os
import subprocess
import numpy as np
import pylab as plt
import platform

if platform.node() == "bigbang.mit.edu":
    basepath = "/bigbang/data/AnnaGroup/caterpillar/halos"

if "harvard.edu" in platform.node():
    basepath = "/n/home01/bgriffen/data/caterpillar/halos"

levellist = [11,12,13,14,15]
nrvirlist = [3,4,5,6]
haloidlist = []
for filename in os.listdir(basepath):
    if filename[0] == "H":
        haloidlist.append(filename)

fig = plt.figure(figsize=(16,10))
plotinc = 0
print "\n"
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
    
    for level in levellist:
        for nrvir in nrvirlist: 
            ext = haloid + "_BB_Z127_P7_LN7_LX" + str(level) + "_O4_NV" + str(nrvir) + "/"
            corepath =  basepath + "/" + haloid + "/" + ext
            marker = 'yD'
            markerface = 'yellow'
            try:
                with open(corepath + "ics.0"):
                    explist = np.loadtxt(corepath+"ExpansionList",delimiter=' ')
		    maxsnap = len(explist)-1
                    icfilesize = os.path.getsize(corepath + "ics.0")

                    if icfilesize > 0:
                        marker = 'rD'
			markerface = 'red'
                    elif icfilesize == 0:
   			 marker = 'gx'
			 markerface = 'green'

                    if os.path.isdir(corepath + "outputs/snapdir_"+str(maxsnap).zfill(3)+"/"):
            	    	marker = 'k^'
			markerface = 'k'

                    elif os.path.isdir(corepath + "outputs/"):
			strprog = ext + " ["
                        subdirnames = basepath + "/" + haloid + "/" + ext + "outputs/"
			snapshotvec = []
                        for subname in os.listdir(subdirnames):
			    if "snapdir" in subname:
				snapshotvec.append(int(subname.replace("snapdir_","")))

			if not snapshotvec:
			    snapshot = -1
			else:
			    snapshot = max(snapshotvec)
			
			if snapshot != -1 and snapshot != maxsnap:
    			    for progsnap in range(0,maxsnap):	
			        if progsnap <= int(snapshot):
				    strprog = strprog + "="
				if progsnap == int(snapshot):
				    strprog = strprog + ">" + str(snapshot)
				if progsnap > int(snapshot):
				    strprog = strprog + "-"
				if progsnap == int(maxsnap)-1:
				    strprog = strprog + "]"
                                
                            print strprog

                        if snapshot != -1: 
                            ax.text(int(level),int(nrvir), str(snapshot), fontsize=9)
                   	    marker = 'k^'
			    markerface = 'white'

                    if os.path.isdir(corepath + "outputs/groups_"+str(maxsnap).zfill(3)+"/"):
                        marker = 'bo'
			markerface = 'b'
    
                    ax.plot(int(level),int(nrvir),marker,markerfacecolor=markerface)
		    

            except IOError: pass


ICArtist = ax.plot((0,1),(0,0),'rD', label='IC')
GADGETArtist = ax.plot((0,1),(0,0),'k^', label='Run')
HALOArtist = ax.plot((0,1),(0,0),'bo', label='Halos')
handles, labels = ax.get_legend_handles_labels()
fig.legend(handles, labels, 'upper center',ncol=3,numpoints=1)
        
plt.show()
