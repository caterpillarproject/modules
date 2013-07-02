import readhalos.readsubf as readsubf
import readsnapshots.readsnap as rs
import numpy as np
import readhalos.AHFDataReader as rah
import readhalos.RSDataReader as rsr
import pandas as pd
import mergertrees.MTCatalogue as MTC

"""
This version of methods is to be used for anything beyond the Subhalo
analysis paper and has elminated most previous functions.
Use the methods.py in ~/SubhaloProj for anything related to that paper.
"""

# distance is periodic
def distance(posA, posB,boxsize=25.):
    dist = abs(posA-posB)
    tmp = dist > boxsize/2.0
    dist[tmp] = boxsize-dist[tmp]
    return np.sqrt(np.sum(dist**2,axis=1))

def getMidpoints(bins):
    """
    Given a range of cut-off values, find the midpoints.
    @return: array of length len(bins)-1
    """
    spacing = bins[1:]-bins[:-1]
    return bins[:-1]+spacing/2.0

def massIndexCutoffs(cutoff, halos):
    """
    cutoff is a list of mass cutoffs to use for partitioning a pandas
    data frame (halos) of halo data.
    returns indexes of data frame corresponding to mass cutoffs,
    with highest masses first.
    halos must be sorted by mass! (the default of RSDataReader.py)
    """
    j = 1
    index_cutoff = []
    for i in range(len(halos)):
        if(halos[i:i+1]['mvir']<cutoff[-j]):
            index_cutoff.append(i)
            j +=1
        if(j>len(cutoff)):
            break
    return index_cutoff
    
