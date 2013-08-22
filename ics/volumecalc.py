import numpy as np
import pylab as plt
from grifflib import *
import readsubf
import readgroup
import readsnap as rs
import readsnapHDF5 as rhdf5
import sys
import matplotlib.pyplot as plt
from ellipse import EllipsoidTool
from itertools import product, combinations

basedir = "/n/scratch2/hernquist_lab/bgriffen/caterpillar/ics/halo80609"
snapnum = 63
ET = EllipsoidTool()
includehalos = False
includeparticles = True
nparttypes = 6
padlist = ['p7']
npart = 1000
reslist = ['']
nvirlist= ['nvir4']
partypelist = np.linspace(1,nparttypes-1,nparttypes)
fig1 = plt.figure(figsize=(23.0,16.0))
ax1 = fig1.add_subplot(111, projection='3d')
if includeparticles:
    for res in reslist:
        for pad in padlist:
            for nvir in nvirlist:
                tmppath = basedir + '/' + pad + '/' + nvir + '/outputs'
                filepath = tmppath + "/snapdir_0" + str(snapnum) + "/snap_0" + str(snapnum)
                xuse = []
                yuse = []
                zuse = []
                snapPOS=rs.read_block(filepath, "POS ",parttype=1,doubleprec=False)
                center,radii,rotation = ET.getMinVolEllipse(snapPOS[0:npart,0:3])
                ax1.scatter(snapPOS[0:npart,0], snapPOS[0:npart,1], snapPOS[0:npart,2], zdir='z', s=20, c='b')
                
                ET.plotEllipsoid(center,radii,rotation,ax=ax1,plotAxes=True)
                #draw cube
                cubemin = snapPOS[0:npart,0:3].min()
                cubemax = snapPOS[0:npart,0:3].max()
                  
                r = [cubemin,cubemax]
                print r
                for s, e in combinations(np.array(list(product(r,r,r))), 2):
                    if np.sum(np.abs(s-e)) == r[1]-r[0]:
                        ax1.plot3D(*zip(s,e), color="b")

                vol = ET.getEllipsoidVolume(radii)
                print vol

#draw sphere
#u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
#x=np.cos(u)*np.sin(v)
#y=np.sin(u)*np.sin(v)
#z=np.cos(v)
#ax.plot_wireframe(x, y, z, color="r")
plt.show()
#ax1.plot(xcirc,ycirc,'-',color='k',linewidth=2)
