import readsnapshots.readsnap as rs
import pylab as plt
import sys
import numpy as np
import os
from operator import itemgetter

fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_ylabel(r"$\mathrm{log_{10}\ N}$", fontsize=14)
ax.set_xlabel(r"$\mathrm{Particle\ Type}$", fontsize=14)
ax.set_xlim([-0.5,5.5])
ax.set_xticklabels(('0','1','2','3','4','5'))
ax.set_xticks((0,1,2,3,4,5))
snapshot = 63

pmgridlist = [128,256,512,1024]

def pdist(filename):
    header = rs.snapshot_header(filename)
    title = os.getcwd().split(os.sep)[-1]
    nallnorm = header.nall
    nallnorm[0] = 1.0
    nall = np.log10(header.nall)
    partypes = np.array([0,1,2,3,4,5])
    nlowres = nallnorm[2:5].sum()
    nhighres = nallnorm[1]
    pmgridfloat = nlowres**(1./3.)*2.
    #for pmgrid in pmgridlist:
        #min(enumerate(a), key=itemgetter(1))[0]
        #deltagrid.append(pmgridfloat - float(pmgrid))

    deltagrid = abs(pmgridfloat - np.array(pmgridlist))
    index = min(enumerate(deltagrid), key=itemgetter(1))[0]
    print "\n"
    print "# Run:",title
    print "# Np low-res: ",nlowres
    print "# Np high-res:",nhighres
    print "# Exact PMGRID:",int(pmgridfloat)
    print "# Use PMGRID:",pmgridlist[index]
    print "# Plotting..."
    ax.bar(partypes,nall,width=0.5,align='center')
    ax.set_title(title,fontsize=14)
    plt.show()

if __name__ == '__main__':
    import os
    pdist(os.getcwd() + "/ics.0")
