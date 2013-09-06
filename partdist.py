import readsnapshots.readsnap as rs
import pylab as plt
import sys
import numpy as np
import os

fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_ylabel(r"$\mathrm{log_{10}\ N}$", fontsize=14)
ax.set_xlabel(r"$\mathrm{Particle\ Type}$", fontsize=14)
ax.set_xlim([-0.5,5.5])
ax.set_xticklabels(('0','1','2','3','4','5'))
ax.set_xticks((0,1,2,3,4,5))
snapshot = 63

def pdist(filename):
    header = rs.snapshot_header(filename)
    title = os.getcwd().split(os.sep)[-1]
    nallnorm = header.nall
    nallnorm[0] = 1.0
    nall = np.log10(header.nall)
    partypes = np.array([0,1,2,3,4,5])
    nlowres = nallnorm[2:5].sum()
    nhighres = nallnorm[1]
    print "\n"
    print "# Run:",title
    print "# low np:",nlowres
    print "# high np:",nhighres
    print "# Recommended approx. PMGRID:",nlowres**(1./3.)*2.
    print "# Plotting..."
    ax.bar(partypes,nall,width=0.5,align='center')
    ax.set_title(title,fontsize=14)
    plt.show()

if __name__ == '__main__':
    import os
    pdist(os.getcwd() + "/ics_rewrite.0")
