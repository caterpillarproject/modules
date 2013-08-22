import readsnapHDF5 as rsHD
import readsnap as rs
import readsubf
import readids
import numpy as np
import RSDataReaderv2 as RSDataReader
import pylab as plt
import sys, os, platform
from grifflib import *
from random import randint
from matplotlib import *

# DEFINE PARAMS
rvirlist = [3,4,5,6,7,8,9]  # NUMBER OF TIMES VIRIAL RADIUS TO ENCLOSE AT Z = 0
snapnum = 63 # DO NOT CHANGE
#idhalo=190897
#idhalo=208737
#idhalo=19910
#idhalo=140666
idhalo=28221
hubble=0.6711

gadpath = '/n/scratch2/hernquist_lab/pzukin/gadget/MusicRuns/512Parent/outputs'
ellipsepath = '/n/scratch2/hernquist_lab/bgriffen/caterpillar/halos'
icspath = '/n/home01/bgriffen/caterpillar/ohahn-music-b116f4cb4be2'
icpath = '/n/scratch2/hernquist_lab/pzukin/music/512Parent'
halopath = '/n/scratch2/hernquist_lab/bgriffen/caterpillar/parent/RockstarData'
filename = '/n/home01/bgriffen/data/caterpillar/halos/candidates.dat'
filename = './candidates.dat'

# LOAD HALOS

halodata = RSDataReader.RSDataReader(halopath,snapnum,digits=2)
allhalos = halodata.get_hosts()

# LOAD CANDIDATES

listin = []
for line in open(filename,'r'):
     li=line.strip()
     if not li.startswith("#"):
         line = line.partition('#')[0]
         listin.append(np.array(line.split(' ')[0:6]))

listin = np.array(listin)
listin = listin.astype(np.float)

#Halos away from the edge
#for i in xrange(0,len(listin)):
#    if listin[i][1] < 1.1E12 and listin[i][1] >= 1.0E12:
#        if listin[i][3] > 10.0 and listin[i][3] < 90.0 and \
#           listin[i][4] > 10.0 and listin[i][4] < 90.0 and \
#           listin[i][5] > 10.0 and listin[i][5] < 90.0:
#            print i,listin[i][0:2],listin[i][3:6]

# SPECIFY CANDIDATE CHOICE

idcand = listin[:,0]
del halodata

for index in xrange(0,len(idcand)):
    if idcand[index] == idhalo:
        print "Found index."
        nhalo = index

rvircand = allhalos.ix[idhalo]['rvir']
mvircand = allhalos.ix[idhalo]['mvir']/hubble
posXcand = allhalos.ix[idhalo]['posX']
posYcand = allhalos.ix[idhalo]['posY']
posZcand = allhalos.ix[idhalo]['posZ']

print "Index:",nhalo
print "Halo ID:",idhalo
print "Halo Mass:",allhalos.ix[idhalo]['mvir']/hubble
print "Halo virial radius:",rvircand
print "Halo z = 0 position:",posXcand,posYcand,posZcand

for Nrvir in rvirlist:
    datadir = ellipsepath+'/halo'+ str(int(idhalo)) +'/'
    
    if not os.path.exists(datadir):
        os.makedirs(datadir)

    print "Distance from halo enclosed:",Nrvir*rvircand
    ext = "/snapdir_063/snap_063"
    snapPOS = rsHD.read_block(gadpath+ext,"POS ")
    dx = allhalos.ix[idhalo]['posX'] - snapPOS[:,0]
    dy = allhalos.ix[idhalo]['posY'] - snapPOS[:,1]
    dz = allhalos.ix[idhalo]['posZ'] - snapPOS[:,2]
    R = np.sqrt(dx**2. + dy**2. + dz**2.)
    print rvircand
    Rindex = np.where(R < Nrvir*rvircand/1000)
    currentpos = snapPOS[Rindex[0]]
    del snapPOS, dx, dy, dz

    # IDENTIFY PARTICLES AT Z = 127
    snapIDS=rsHD.read_block(gadpath+ext,"ID  ")
    regionIDS = snapIDS[Rindex[0]]
    extics = "/ics"
    del snapIDS

    snapIDSlagr = rs.read_block(icpath + extics,"ID  ",doubleprec=False)
    mask = np.in1d(snapIDSlagr,regionIDS,assume_unique=True)
    del snapIDSlagr
    snapPOSlagr = rs.read_block(icpath + extics,"POS ",doubleprec=False)
    lagrPos=snapPOSlagr[mask]
    del mask, regionIDS
    
    # OUTPUT REGION IN UNITS OF [X,Y,Z]/BOX WIDTH
    header=rsHD.snapshot_header(gadpath+"/snapdir_063/snap_063.0")
    print "# of particles:",len(Rindex[0])
    #print "Enclosed mass:",'{:.2e}'.format((sum(header.npart)*10**10)*len(Rindex))

    f2=open(ellipsepath+'/halo'+ str(int(idhalo)) +'/halo' + str(int(idhalo)) + 'Nrvir'+str(int(Nrvir)),'w')
    for iv in xrange(0,len(lagrPos[:,0])):
        f2.write(str(lagrPos[iv,0]/header.boxsize)+' '+str(lagrPos[iv,1]/header.boxsize)+' '+str(lagrPos[iv,2]/header.boxsize)+'\n')

    f2.close()

    print "Created lagrangian region here..."
    print ellipsepath+'/halo'+str(int(idhalo))+'/halo'+str(int(idhalo))+'Nrvir'+str(int(Nrvir))

    # CREATE BOX OUTPUTS
    CorrectPos(lagrPos[:,0],header.boxsize)
    CorrectPos(lagrPos[:,1],header.boxsize)
    CorrectPos(lagrPos[:,2],header.boxsize)
    comV=COM(lagrPos[:,0],lagrPos[:,1],lagrPos[:,2])

    xc=comV[0]/header.boxsize
    yc=comV[1]/header.boxsize
    zc=comV[2]/header.boxsize

    print "if [ $nvir ==",Nrvir,"]; then"
    print "xc=",xc
    print "yc=",yc
    print "zc=",zc

    #f.write('Method 2: '+str(comV/header.boxsize) + '\n')
    dx=max(abs(lagrPos[:,0]-comV[0]))
    dy=max(abs(lagrPos[:,1]-comV[1]))
    dz=max(abs(lagrPos[:,2]-comV[2]))

    xe=2.0*dx*1.12/header.boxsize
    ye=2.0*dy*1.12/header.boxsize
    ze=2.0*dz*1.12/header.boxsize

    print "xe=",xe
    print "ye=",ye
    print "ze=",ze
    print "fi"    
    
    #print "Volume:", xe*ye*ze, "(Mpc/h)^3"

#    cond1 = (snapPOSlagr[:,0] >= xc-0.5*xe) & (snapPOSlagr[:,0] <= xc+0.5*xe)
#    cond2 = (snapPOSlagr[:,1] >= yc-0.5*ye) & (snapPOSlagr[:,1] <= yc+0.5*ye)
#    cond3 = (snapPOSlagr[:,2] >= zc-0.5*ze) & (snapPOSlagr[:,2] <= zc+0.5*ze)
    
#    cond = cond1 & cond2 & cond3
#    snapvol = snapPOSlagr[cond]
#    print "Number of particles in minimum cuboid surrounding Lagrangian volume:",np.shape(snapvol)
    #f.write("   Extent: " + str(2.0*dx*1.12/header.boxsize) + ' ' + str(2.0*dy*1.12/header.boxsize) + ' ' + str(2.0*dz*1.12/header.boxsize) + '\n')
    #f.close()
    #del lagrPOS

    minx = np.vstack((lagrPos[:,0],currentpos[:,0])).min()
    maxx = np.vstack((lagrPos[:,0],currentpos[:,0])).max()
    miny = np.vstack((lagrPos[:,1],currentpos[:,1])).min()
    maxy = np.vstack((lagrPos[:,1],currentpos[:,1])).max()
    minz = np.vstack((lagrPos[:,2],currentpos[:,2])).min()
    maxz = np.vstack((lagrPos[:,2],currentpos[:,2])).max()
 
    fig1 = plt.figure(figsize=(20.0,10.0))
    ax1 = fig1.add_subplot(231)
    ax2 = fig1.add_subplot(232)
    ax3 = fig1.add_subplot(233)
    ax4 = fig1.add_subplot(234)
    ax5 = fig1.add_subplot(235)
    ax6 = fig1.add_subplot(236)
    
    header=rs.snapshot_header(icpath + extics + '.0')

    zinit = header.redshift
    zinittitle = '{:.2f}'.format(zinit)

    ax2.set_title('POSITION AT Z = ' + zinittitle)
    ax1.set_xlabel('X-POS [Mpc/h]')
    ax2.set_xlabel('Y-POS [Mpc/h]')
    ax3.set_xlabel('Z-POS [Mpc/h]')

    ax5.set_title('POSITION AT Z = 0')
    ax4.set_xlabel('X-POS [Mpc/h]')
    ax5.set_xlabel('Y-POS [Mpc/h]')
    ax6.set_xlabel('Z-POS [Mpc/h]')

    ax1.hist(lagrPos[:,0])
    ax2.hist(lagrPos[:,1])
    ax3.hist(lagrPos[:,2])
    ax4.hist(currentpos[:,0])
    ax5.hist(currentpos[:,1])
    ax6.hist(currentpos[:,2])

    ax1.set_xlim([minx,maxx])
    ax4.set_xlim([minx,maxx])
    ax2.set_xlim([miny,maxy])
    ax5.set_xlim([miny,maxy])
    ax3.set_xlim([minz,maxz])
    ax6.set_xlim([minz,maxz])

    fig2 = plt.figure(figsize=(20.0,10.0))
    ax7 = fig2.add_subplot(231)
    ax8 = fig2.add_subplot(232)
    ax9 = fig2.add_subplot(233)
    ax10 = fig2.add_subplot(234)
    ax11 = fig2.add_subplot(235)
    ax12 = fig2.add_subplot(236)

    markertype = 'o'
    colselect = 'b'

    ax7.plot(currentpos[:,0],currentpos[:,1],marker=markertype,markerfacecolor=colselect,markeredgecolor=colselect,linestyle='none')
    ax8.plot(currentpos[:,0],currentpos[:,2],marker=markertype,markerfacecolor=colselect,markeredgecolor=colselect,linestyle='none')
    ax9.plot(currentpos[:,1],currentpos[:,2],marker=markertype,markerfacecolor=colselect,markeredgecolor=colselect,linestyle='none')
    ax10.plot(lagrPos[:,0],lagrPos[:,1],marker=markertype,markerfacecolor=colselect,markeredgecolor=colselect,linestyle='none')
    ax11.plot(lagrPos[:,0],lagrPos[:,2],marker=markertype,markerfacecolor=colselect,markeredgecolor=colselect,linestyle='none')
    ax12.plot(lagrPos[:,1],lagrPos[:,2],marker=markertype,markerfacecolor=colselect,markeredgecolor=colselect,linestyle='none')

    ax7.set_xlim([minx,maxx])
    ax7.set_ylim([miny,maxy])
    ax8.set_xlim([minx,maxx])
    ax8.set_ylim([minz,maxz])
    ax9.set_xlim([miny,maxy])
    ax9.set_ylim([minz,maxz])
    ax10.set_xlim([minx,maxx])
    ax10.set_ylim([miny,maxy])
    ax11.set_xlim([minx,maxx])
    ax11.set_ylim([minz,maxz])
    ax12.set_xlim([miny,maxy])
    ax12.set_ylim([minz,maxz])

    ax8.set_title('POSITION AT Z = 0')
    ax11.set_title('POSITION AT Z = ' + zinittitle )

    ax7.set_xlabel('X-POS [Mpc/h]')
    ax7.set_ylabel('Y-POS [Mpc/h]')
    ax8.set_xlabel('X-POS [Mpc/h]')
    ax8.set_ylabel('Z-POS [Mpc/h]')
    ax9.set_xlabel('Y-POS [Mpc/h]')
    ax9.set_ylabel('Z-POS [Mpc/h]')
    ax10.set_xlabel('X-POS [Mpc/h]')
    ax10.set_ylabel('Y-POS [Mpc/h]')
    ax11.set_xlabel('X-POS [Mpc/h]')
    ax11.set_ylabel('Z-POS [Mpc/h]')
    ax12.set_xlabel('Y-POS [Mpc/h]')
    ax12.set_ylabel('Z-POS [Mpc/h]')

    figdir = './figs/halo' + str(int(idhalo)) + '/'

    if not os.path.exists(figdir):
        os.makedirs(figdir)

    fig1.savefig(figdir + 'halo' + str(int(idhalo)) + '_Nrvir' + str(int(Nrvir)) + '_xyzhistpos.png',bbox_inches='tight')
    fig2.savefig(figdir + 'halo' + str(int(idhalo)) + '_Nrvir' + str(int(Nrvir)) + '_xyzscatterpos.png',bbox_inches='tight')

    fig1.clf()
    fig2.clf()

#plt.show()


