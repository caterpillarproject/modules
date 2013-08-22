import readsnapHDF5 as rsHD
import readsnap as rs
import readsubf
import readids
import numpy as np
#import RSDataReaderv2 as RSDataReader
import pylab as plt
import sys
from random import randint
from matplotlib import *
import GregNewMTCatalogue as GregNewMTCatalogue

# DEFINE PARAMS
singlehalo = True
hubble = 0.6711
nhalo = 2826
snapnum = 63

# SET PATHS

halopath = '/n/home01/bgriffen/data/caterpillar/parent/RockstarData'
gadpath = '/n/scratch2/hernquist_lab/pzukin/gadget/MusicRuns/512Parent/outputs'
filename = '/n/home01/bgriffen/data/caterpillar/halos/candidates.dat'

# LOAD HALOS

#halodata = RSDataReader.RSDataReader(halopath,snapnum,digits=2)
#allhalos = halodata.get_hosts()

# LOAD CANDIDATES

listin = []
for line in open(filename,'r'):
     li=line.strip()
     if not li.startswith("#"):
         line = line.partition('#')[0]
         listin.append(np.array(line.split(' ')[0:6]))

listin = np.array(listin)
listin = listin.astype(np.float)

# SPECIFY CANDIDATE CHOICE


#rvircand = allhalos.ix[candindex]['rvir']

#print "Halo ID:",candindex
#print "Halo Mass:",listin[nhalo,1]
#print "Halo virial radius:",rvircand
#print "Halo z = 0 position:",listin[nhalo,3:6]
#print "Distance from halo enclosed:",Nrvir*rvircand

basepath = '/spacebase/data/AnnaGroup/caterpillarparent/RockstarData/trees/'

fig1 = plt.figure(figsize=(12.0,12.0))
ax1 = fig1.add_subplot(221)
ax2 = fig1.add_subplot(222)
ax3 = fig1.add_subplot(223)
ax4 = fig1.add_subplot(224)

plt.subplots_adjust(hspace=0.05)
xticklabels = ax1.get_xticklabels() + ax2.get_xticklabels()
plt.setp(xticklabels, visible=False)
#len(listin)
if singlehalo == False:
    for i in xrange(0,100):
        mvir = listin[i,1]/hubble
        if mvir >= 1e12 and mvir <= 3e12:
            candindex = idcand[i]
            candindex = int(candindex)
            MTC = GregNewMTCatalogue.NewMTCatalogue(halopath + '/trees/tree.bin',indexfile=halopath + '/trees/treeindex.csv',hostid=candindex)
            tree = MTC.Trees[candindex]
            candMB = tree.getMainBranch(0)
            candMB = np.array(candMB,dtype=MTC.fmttype)
            scalevec = candMB['scale']
            Mvirvec = candMB['mvir']/hubble
            Vmaxvec = candMB['vmax']
    
            normmass = Mvirvec/Mvirvec[0]
            normvMax = Vmaxvec**2/Vmaxvec[0]**2
 
            ax1.plot(scalevec,Mvirvec,linewidth=2)
            ax2.plot(scalevec,normmass,linewidth=2)
            ax3.plot(scalevec,Vmaxvec,linewidth=2)
            ax4.plot(scalevec,normvMax,linewidth=2)

candlist = [1.98301000e5, 80609, 4.72790000e+04, 1.41542000e+05, 3.71610000e+04, 8.94200000e+04, 2.59245000e+05, 3.70390000e+04, 2.78990000e+04, 2.06174000e5, 1.85700000e+03, 1.41522000e+05, 1.98332000e+05, 1.25135000e+05, 1.15657000e+05]
candlist = np.array(candlist)
if singlehalo == True:
    for i in xrange(0,len(candlist)):
        idcand = listin[:,0]
        candindex = candlist[i]
        candindex = int(candindex)
        MTC = GregNewMTCatalogue.NewMTCatalogue(halopath + '/trees/tree.bin',indexfile=halopath + '/trees/treeindex.csv',hostid=candindex)
        tree = MTC.Trees[candindex]
        candMB = tree.getMainBranch(0)
        candMB = np.array(candMB,dtype=MTC.fmttype)
        scalevec = candMB['scale']
        Mvirvec = candMB['mvir']/hubble
        Vmaxvec = candMB['vmax']
    
        normmass = Mvirvec/Mvirvec[0]
        normvMax = Vmaxvec**2/Vmaxvec[0]**2
 
        ax1.plot(scalevec,Mvirvec,linewidth=2)
        ax2.plot(scalevec,normmass,linewidth=2)
        ax3.plot(scalevec,Vmaxvec,linewidth=2)
        ax4.plot(scalevec,normvMax,linewidth=2)

ax1.set_yscale('log')
ax1.set_ylabel(r'$M_v(z)$',size=14)
ax2.set_yscale('log')
ax2.set_ylabel(r'$M_v(z)/M_v(z=0)$',size=14)
ax3.set_xlabel(r'$scale$ $factor$',size=14)
ax3.set_ylabel(r'$V_{max}(z)$',size=14)
ax4.set_xlabel(r'$scale$ $factor$',size=14)
ax4.set_ylabel(r'$V_{max}(z)^2/V_{max}(z=0)^2$',size=14)
fig1.savefig('./mvirvmaxhistory.png')

plt.show()
