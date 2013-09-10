import numpy as np
import pylab as plt
from grifflib import *
import readsubf
import readgroup
import readsnap as rs
import readsnapHDF5 as rhdf5
import sys
import matplotlib.pyplot as plt

cm = plt.cm.get_cmap('RdYlBu')
#basedir = "/n/scratch2/hernquist_lab/pzukin/gadget/MusicRuns/H1530_2"
basedir = "/n/scratch2/hernquist_lab/bgriffen/caterpillar/ics/halo80609/ellipsoid"
snapnum = 63
#ET = EllipsoidTool()
includehalos = False
includeparticles = True
nparttypes = 6
padlist = ['p6','p7','p8','p9','p10']
reslist = ['l11']
nvirlist= ['nvir3','nvir4','nvir5','nvir6','nvir7','nvir8','nvir9']
partypelist = np.linspace(1,nparttypes-1,nparttypes)
if includehalos:
    s = readsubf.subfind_catalog(basedir, snapnum)
    sid = readsubf.subf_ids(basedir, snapnum, 0, 0, read_all=1)
    rvir = s.group_r_mean200
    xpos = s.group_pos[:,0]
    ypos = s.group_pos[:,1]

    pccontam = []

    for i in xrange(0,len(xpos)):
        pccontam.append(float(s.group_contamination_count[i])/float(s.group_len[i]))

    pccontam = np.array(pccontam)

nwide = len(nvirlist)
nhigh = len(padlist)
fig1 = plt.figure(figsize=(23.0,16.0))
fig2 = plt.figure(figsize=(23.0,16.0))
fig3 = plt.figure(figsize=(23.0,16.0))
if includeparticles:
    figi = 0
    xmid = 51
    ymid = 48
    zmid = 47
    delta = 6.5
    for res in reslist:
        for pad in padlist:
            for nvir in nvirlist:
                figi += 1
                tmppath = basedir + '/' + res + '/' + pad + '/' + nvir + '/outputs'
                filepath = tmppath + "/snapdir_0" + str(snapnum) + "/snap_0" + str(snapnum)
                #print '/' + pad + '/' + nvir +
                xuse = []
                yuse = []
                zuse = []
                for parttype in partypelist:
                    snapPOS=rs.read_block(filepath, "POS ",parttype=int(parttype),doubleprec=False)
                    xptcle = snapPOS[:,0]
                    yptcle = snapPOS[:,1]
                    zptcle = snapPOS[:,2]
                    #del snapPOS
                    xuse = np.concatenate([xuse,xptcle])
                    yuse = np.concatenate([yuse,yptcle])
                    zuse = np.concatenate([zuse,zptcle])

                    heatmap, xedges, yedges = np.histogram2d(xuse, yuse, bins=1024)
                    extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
                    heatmap = np.flipud(np.rot90(heatmap))
                    ax1 = fig1.add_subplot(nhigh,nwide,figi)
                    ax1.imshow(np.log10(heatmap), extent=extent,cmap = 'jet', origin='lower')
                    
                    ax2 = fig2.add_subplot(nhigh,nwide,figi)
                    heatmap, xedges, yedges = np.histogram2d(xuse, zuse, bins=1024)
                    extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
                    heatmap = np.flipud(np.rot90(heatmap))
                    ax2.imshow(np.log10(heatmap), extent=extent,cmap = 'jet', origin='lower')

                    ax3 = fig3.add_subplot(nhigh,nwide,figi)
                    heatmap, xedges, yedges = np.histogram2d(yuse, zuse, bins=1024)
                    extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
                    heatmap = np.flipud(np.rot90(heatmap))
                    ax3.imshow(np.log10(heatmap), extent=extent,cmap = 'jet', origin='lower')
                
                    #if parttype == 1:
                    #    center,radii,rotation = ET.getMinVolEllipse(snapPOS)
                    #    vol = ET.getEllipsoidVolume(radii)
                    #    print vol

                if includehalos:
                    xcirc,ycirc = drawcircle(xpos,ypos,rvir,vec=True)
                    sc = plt.scatter(xpos,ypos,c=s.group_contamination_count,s=s.group_contamination_count/10, \
                    cmap='gray',facecolor='none',edgecolor='k',linewidth=4)

                ax1.set_xlabel(r'$\mathrm{x-pos\ [Mpc/h}]$',size=14)
                ax1.set_ylabel(r'$\mathrm{y-pos\ [Mpc/h}]$',size=14)
                ax1.set_xlim([xmid - delta, xmid + delta])
                ax1.set_ylim([ymid - delta, ymid + delta])
                titlestr = pad + ", " + nvir
                ax1.text(0.9, 0.9,titlestr,
                    horizontalalignment='right',
                    verticalalignment='top',
                    color='white',
                    transform = ax1.transAxes)

                ax2.set_xlabel(r'$\mathrm{x-pos\ [Mpc/h}]$',size=14)
                ax2.set_ylabel(r'$\mathrm{z-pos\ [Mpc/h}]$',size=14)
                ax2.set_xlim([xmid - delta, xmid + delta])
                ax2.set_ylim([zmid - delta, zmid + delta])
                ax2.text(0.9, 0.9,titlestr,
                    horizontalalignment='right',
                    verticalalignment='top',
                    color='white',
                    transform = ax2.transAxes)

                ax3.set_xlabel(r'$\mathrm{y-pos\ [Mpc/h}]$',size=14)
                ax3.set_ylabel(r'$\mathrm{z-pos\ [Mpc/h}]$',size=14)
                ax3.set_xlim([xmid - delta - 3, xmid + delta - 3])
                ax3.set_ylim([zmid - delta, zmid + delta])
                ax3.text(0.9, 0.9,titlestr,
                    horizontalalignment='right',
                    verticalalignment='top',
                    color='white',
                    transform = ax3.transAxes)
                
              
                if includehalos:
                    cbar = plt.colorbar(sc)
                    cbar.set_label(r'$\mathrm{number\ of\ contamination\ particles}$',size=14)


fig1.savefig('./figs/xy-p6-10nvir3-9.png')
fig2.savefig('./figs/xz-p6-10nvir3-9.png')
fig3.savefig('./figs/yz-p6-10nvir3-9.png')
plt.show()
#ax1.plot(xcirc,ycirc,'-',color='k',linewidth=2)
