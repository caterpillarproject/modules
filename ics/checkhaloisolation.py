import readsnapHDF5 as rsHD
import readsnap as rs
import readsubf
import readids
import numpy as np
import RSDataReaderv2 as RSDataReader
import pylab as plt
import sys
from grifflib import *
from matplotlib import *

# DEFINE PARAMS
boxwidth = 3. # width of zoomed region
nhalo = 3098
snapnum = 63
titlestr = "halo" + str(nhalo)

# SET PATHS
gadpath = '/n/scratch2/hernquist_lab/pzukin/gadget/MusicRuns/512Parent/outputs'
halopath = '/n/scratch2/hernquist_lab/bgriffen/caterpillar/parent/RockstarData'
filename = '/n/home01/bgriffen/data/caterpillar/halos/candidates.dat'

# LOAD HALOS
halodata = RSDataReader.RSDataReader(halopath,snapnum,digits=2)
allhalos = halodata.get_hosts()

xpos = np.array(allhalos['posX'])
ypos = np.array(allhalos['posY'])
zpos = np.array(allhalos['posZ'])
rvir = np.array(allhalos['rvir'])

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
idcand = listin[:,0]
candindex = idcand[nhalo]
rvircand = allhalos.ix[candindex]['rvir']

print "Halo ID:",candindex
print "Halo Mass:",listin[nhalo,1]
print "Halo virial radius:",rvircand
print "Halo z = 0 position:",listin[nhalo,3:6]
titlestr = "halo" + str(int(candindex))


xmid = listin[nhalo,3]
ymid = listin[nhalo,4]
zmid = listin[nhalo,5]
rvircand = rvircand/1000

# DISPLAY XYZ PROJECTIONS + NEIGHBOURS
fig1 = plt.figure(figsize=(18,5))
ax1 = fig1.add_subplot(1,3,1)
ax2 = fig1.add_subplot(1,3,2)
ax3 = fig1.add_subplot(1,3,3)

xposclose = []
yposclose = []
zposclose = []
rvirclose = []

R = np.sqrt((xmid-xpos)**2 + (ymid-ypos)**2 + (zmid-zpos)**2)
for i in xrange(0,len(xpos)):
    if R[i] < boxwidth:
        xposclose.append(xpos[i])
        yposclose.append(ypos[i])
        zposclose.append(zpos[i])
        rvirclose.append(rvir[i])

xposclose = np.array(xposclose)
yposclose = np.array(yposclose)
zposclose = np.array(zposclose)
rvirclose = np.array(rvirclose)

xcirc,ycirc = drawcircle(xposclose,yposclose,rvirclose/1000,vec=True)
ax1.plot(xcirc,ycirc,'b-',linewidth=2)
xcirc,zcirc = drawcircle(xposclose,zposclose,rvirclose/1000,vec=True)
ax2.plot(xcirc,zcirc,'b-',linewidth=2)
ycirc,zcirc = drawcircle(yposclose,zposclose,rvirclose/1000,vec=True)
ax3.plot(ycirc,zcirc,'b-',linewidth=2)

xcand,ycand = drawcircle(xmid,ymid,rvircand,vec=False)
ax1.plot(xcand,ycand,'r-',linewidth=2)
xcand,zcand = drawcircle(xmid,zmid,rvircand,vec=False)
ax2.plot(xcand,zcand,'r-',linewidth=2)
ycand,zcand = drawcircle(ymid,zmid,rvircand,vec=False)
ax3.plot(ycand,zcand,'r-',linewidth=2)

boxwidth = boxwidth/2
print [xmid - boxwidth, xmid + boxwidth]
print [ymid - boxwidth, ymid + boxwidth]
print [zmid - boxwidth, zmid + boxwidth]

ax1.set_xlim([xmid - boxwidth, xmid + boxwidth])
ax1.set_ylim([ymid - boxwidth, ymid + boxwidth])
ax1.text(0.05, 0.95,titlestr,
    horizontalalignment='left',
    verticalalignment='top',
    color='black',
    transform = ax1.transAxes)
ax2.set_xlim([xmid - boxwidth, xmid + boxwidth])
ax2.set_ylim([zmid - boxwidth, zmid + boxwidth])
ax2.text(0.05, 0.95,titlestr,
    horizontalalignment='left',
    verticalalignment='top',
    color='black',
    transform = ax2.transAxes)
ax3.set_xlim([ymid - boxwidth, ymid + boxwidth])
ax3.set_ylim([zmid - boxwidth, zmid + boxwidth])
ax3.text(0.05, 0.95,titlestr,
    horizontalalignment='left',
    verticalalignment='top',
    color='black',
    transform = ax3.transAxes)

ax1.set_xlabel(r'$\mathrm{x-pos\ [Mpc/h]}$',size=14)
ax2.set_xlabel(r'$\mathrm{x-pos\ [Mpc/h]}$',size=14)
ax3.set_xlabel(r'$\mathrm{y-pos\ [Mpc/h]}$',size=14)

ax1.set_ylabel(r'$\mathrm{y-pos\ [Mpc/h]}$',size=14)
ax2.set_ylabel(r'$\mathrm{z-pos\ [Mpc/h]}$',size=14)
ax3.set_ylabel(r'$\mathrm{z-pos\ [Mpc/h]}$',size=14)

plt.show()
