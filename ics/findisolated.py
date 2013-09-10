import readsnapHDF5 as rsHD
import readsnap as rs
import readsubf
import readids
import numpy as np
import random
import RSDataReaderv2 as RSDataReader
import pylab as plt
import sys
from random import randint

snapnum = 63    # DO NOT CHANGE
hubble = 0.6711

# SET DATA DIRECTORIES

halopath = '/n/home01/bgriffen/data/caterpillar/parent/RockstarData'
#output_file = '/n/scratch2/hernquist_lab/bgriffen/caterpillar/halos/candidates.dat'
outputfile = '/n/home01/bgriffen/caterpillar/haloselection/candidates.dat'

# SELECT CANDIDATE MASS RANGE
lowermasscut = 7e11
uppermasscut = 7e12

# LOAD ROCKSTAR HALOS
print "READING...",halopath
halodata = RSDataReader.RSDataReader(halopath,snapnum,digits=2)
allhalos = halodata.get_hosts()
MWcand = allhalos[np.logical_and(allhalos['mvir']/hubble>lowermasscut, allhalos['mvir']/hubble<uppermasscut)]

# LOAD HALOS LARGER THAN X MSOL FOR ISOLATION REQUIREMENT

halos12 = allhalos[allhalos['mvir']/hubble>7e12] # HALOS LARGER THAN 7E12
halos13 = allhalos[allhalos['mvir']/hubble>7e13] # HALOS LARGER THAN 7E13
print "Number of MW candidates between 7e11 and 7e12:",len(MWcand)
print "Number of halos larger than 7e12 Msol",len(halos12)
print "Number of halos larger than 7e13 Msol",len(halos13)

xpos12 = np.array(np.float64(halos12['posX']))
ypos12 = np.array(np.float64(halos12['posY']))
zpos12 = np.array(np.float64(halos12['posZ']))

xpos13 = np.array(np.float64(halos13['posX']))
ypos13 = np.array(np.float64(halos13['posY']))
zpos13 = np.array(np.float64(halos13['posZ']))


# OUTPUT FILE
out = open(output_file, 'w')
out.write('# ID M Rvir x y z\n')

Ncandidates = 0
for i in xrange(0,len(MWcand)):
    # CYCLE CANDIDATE i
    xposi = np.array(np.float64(MWcand['posX']))[i]
    yposi = np.array(np.float64(MWcand['posY']))[i]
    zposi = np.array(np.float64(MWcand['posZ']))[i]
    massi = np.array(np.float64(MWcand['mvir']))[i]/hubble
    rviri = np.array(np.float64(MWcand['rvir']))[i]
    idi = np.array(MWcand['id'])[i]

    # CALCULATE DISTANCE TO HALOS LARGER THAN 7E12 and 7E13 MSOL

    R13 = np.sqrt((xposi-xpos13)**2.+(yposi-ypos13)**2.+(zposi-zpos13)**2.)/hubble
    R12 = np.sqrt((xposi-xpos12)**2.+(yposi-ypos12)**2.+(zposi-zpos12)**2.)/hubble

    # SELECT ALL HALOS LARGER THAN HALF THE SIZE OF THE CANDIDATE idi

    largerthanMW = allhalos[allhalos['mvir']/hubble >= 0.5*massi]

    # SINCE 0.5*massi INCLUDES idi, NEED TO REMOVE IT FROM X,Y,Z POS CALC SO MIN(R) != 0.0

    idindex = np.where(largerthanMW['id'] != idi)
    xtmp = np.array(np.float64(largerthanMW['posX']))
    ytmp = np.array(np.float64(largerthanMW['posY']))
    ztmp = np.array(np.float64(largerthanMW['posZ']))

    xposMgtMW = xtmp[idindex[0]]
    yposMgtMW = ytmp[idindex[0]]
    zposMgtMW = ztmp[idindex[0]]

    RMgtMW = np.sqrt((xposi-xposMgtMW)**2.+(yposi-yposMgtMW)**2.+(zposi-zposMgtMW)**2.)/hubble

    # NO HALO LARGER THAN 7e13 CLOSER THAN 4 MPC
    if R13.min() >= 4.:
            # NO HALO LARGER THAN 7e12 CLOSER THAN 3 MPC
        if R12.min() >= 3.:
                # NO HALO HALF THE MASS OF CANDIDATE OR LARGER WITHIN 1.4 MPC
            if RMgtMW.min() >= 1.4:
                Ncandidates += 1
                out.write('%f %e %f %f %f %f \n' % 
                    (int(idi),massi,rviri,xposi,yposi,zposi)) 

out.close()

print "Number of candidates found:",Ncandidates
print "Output written to:",output_file
print "Now use makeics.py to construct lagrangian region."

sys.exit()
