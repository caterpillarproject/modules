import readsnapshots.readsnap as rs
import pylab as plt
import sys
import numpy as np

fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_ylabel("Number", fontsize=14)
ax.set_xlabel("Particle Type", fontsize=14)
snapshot = 63

def pdist(file):
    print "Reading: ",filename
    header = rs.snapshot_header(filename)
    nall = np.log10(header.nall)
    partypes = np.array([0,1,2,3,4,5])
    ax.plot(partypes,nall,)
    plt.show()

if __name__ == '__main__':
    import os
    pdist(os.getcwd() + "/ics_rewrite.0")