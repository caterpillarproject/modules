import sys,os
import subprocess
import numpy as np
import pylab as plt
import platform
import time
from brendanlib.grifflib import getcurrentjobs

jobids,currentjobs,statusjobs = getcurrentjobs()
tthresh = 6

if "bigbang" in platform.node() or "antares" in platform.node() or "spacebase" in platform.node():
    basepath = "/bigbang/data/AnnaGroup/caterpillar/halos"

if "harvard" in platform.node():
    basepath = "/n/home01/bgriffen/data/caterpillar/halos"

levellist = [11,12,13,14,15]
nrvirlist = [4]
haloidlist = []
for filename in os.listdir(basepath):
    if filename[0] == "H":
        haloidlist.append(filename)

nhalos = len(haloidlist)

nplots_wide = 12
fig,axs = plt.subplots(nplots_wide,10,sharex=True,sharey=True,figsize=(20,12)) 
plt.subplots_adjust( hspace=0 ) 
plt.subplots_adjust( wspace=0 )
#plt.tight_layout()
plotinc = 0

tnow = time.time()
should_time = tnow - (tthresh * 3600)
jobstocancel = []

leveldone = {}
nhalolevels = {}

for level in xrange(11,15):
    leveldone[level] = 0
    nhalolevels[level] = 0

nwide = 0
nhigh = 0
for haloid in haloidlist:

    plotinc += 1
    ax = axs[nwide,nhigh]
    if plotinc % nplots_wide == 0:
	nhigh += 1
        nwide = 0
    else:
        nwide += 1

    #ax.set_title(haloid)
    ax.text(0.5,0.80,haloid,fontsize=14, horizontalalignment='center',verticalalignment='center',transform=ax.transAxes)
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
            if os.path.isfile(corepath + "ics.0"):
                icfilesize = os.path.getsize(corepath + "ics.0")
                if icfilesize > 0:
                    marker = 'rD'
                    markerface = 'red'
                elif icfilesize == 0:
                     marker = 'gx'
                     markerface = 'green'

                if os.path.isdir(corepath + "outputs/"):
                    explist = np.loadtxt(corepath+"ExpansionList",delimiter=' ')
                    maxsnap = len(explist)-1
                    
                    strprog = ext + " ["
                    subdirnames = basepath + "/" + haloid + "/" + ext + "outputs/"
                    snapshotvec = []
                    for subname in os.listdir(subdirnames):
                        if "snapdir" in subname and "BACK" not in subname:
                            snapshotvec.append(int(subname.replace("snapdir_","")))
                        
                    if not snapshotvec:
                        snapshot = -1
  		        strprog = ext + "Completed: 0% -/" + str(maxsnap)
                    else:
                        snapshot = max(snapshotvec)
			strprog = ext + "Completed: %0.2f" % (float(snapshot)*100./float(maxsnap)) + "%, " + str(snapshot)+ "/" + str(maxsnap)

		    try:
		        lastsnapfile = corepath + "/outputs/snapdir_"+str(snapshot).zfill(3)+"/snap_"+str(snapshot).zfill(3)+".0.hdf5"
                        with open(lastsnapfile):
                            file_mod_time = os.stat(lastsnapfile).st_mtime
                            last_time = round((int(tnow) - file_mod_time) / 3600, 2)
			    snaptimestr = "-- wrote out last snapshot {:.2f} hours ago.".format(last_time)
                    except:
			pass

                    try:
                        with open(corepath+"/outputs/cpu.txt"):
                            file_mod_time = os.stat(corepath+"/outputs/cpu.txt").st_mtime
                            last_time = round((int(tnow) - file_mod_time) / 3600, 2)
			    cputimestr = "-- wrote out 'cpu.txt'     {:.2f} hours ago.".format(last_time)
                    except:
                        pass

		    jobname = haloid+"N"+str(nrvir)+"L"+str(level)[-1]

		    if jobname in currentjobs:
			idx = currentjobs.index(jobname)
	                statusi = statusjobs[idx]
  		        if snapshot != maxsnap and "-/" not in strprog and (file_mod_time - should_time) < tthresh and statusi == "R":
			    print
			    print ext
 		            print "--",strprog.replace(ext,"")
			    print snaptimestr
			    print cputimestr			    
			    jobstocancel.append(int(jobids[idx]))

		    if snapshot == maxsnap:
			leveldone[level] += 1

		    nhalolevels[level] += 1

                    if snapshot != -1:
                        marker = 'k^'
                        markerface = 'white'
			if jobname in currentjobs:
  			    idx = currentjobs.index(jobname)
                            statusi = statusjobs[idx]
			    print statusi
		
			    if statusi == 'R':
			        ax.text(int(level),int(nrvir), str(snapshot)+'R', fontsize=9)
			    else:
			        ax.text(int(level),int(nrvir), str(snapshot)+'P', fontsize=9)
			else:
			    ax.text(int(level),int(nrvir), str(snapshot), fontsize=9)

		    if os.path.isdir(corepath + "outputs/snapdir_"+str(maxsnap).zfill(3)+"/"):
                        marker = 'k^'
                        markerface = 'k'

                    if os.path.isdir(corepath + "outputs/groups_"+str(maxsnap).zfill(3)+"/"):
                        marker = 'bo'
                        markerface = 'b'

                    if os.path.isdir(corepath + "halos/halo_"+str(maxsnap).zfill(3)+"/"):
                        marker = 'co'
                        markerface = 'c'

                    if os.path.isdir(corepath + "outputs/groups_"+str(maxsnap).zfill(3)+"/") \
                       and os.path.isdir(corepath + "halos/halos_"+str(maxsnap).zfill(3)+"/"):
                        marker = 'go'
                        markerface = 'g'

            ax.plot(int(level),int(nrvir),marker,markerfacecolor=markerface,markeredgewidth=None)



ICArtist = ax.plot((0,1),(0,0),'rD', label='IC')
GADGETArtist = ax.plot((0,1),(0,0),'k^', label='Gadget Complete')
GADGETncArtist = ax.plot((0,1),(0,0),'k^',markerfacecolor='white',label='Gadget On-going')
SUBFINDArtist = ax.plot((0,1),(0,0),'bo', label='SUBFIND')
ROCKSTARArtist = ax.plot((0,1),(0,0),'mo', label='Rockstar')
HALOSArtist = ax.plot((0,1),(0,0),'go', label='S&R [Complete]')
handles, labels = ax.get_legend_handles_labels()
fig.legend(handles, labels, 'upper center',ncol=6,numpoints=1)

print
print
print "Potential jobs to cancel (been running over: " + str(tthresh) + " hours)"
print jobstocancel

print
print "| Level | % Completed |"
for key,value in leveldone.iteritems():
    if nhalolevels[key] != 0:
	print "   %i        %3.2f" % (key,float(value)/float(nhalolevels[key])*100)

plt.show()
