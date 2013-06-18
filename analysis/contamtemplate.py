import numpy as np
import pylab as plt
import readsubf
import readgroup
import RSDataReaderv2 as RSDataReader
import readsnap as rs
import readsnapHDF5 as rhdf5

#IMPORTANT STUFF
ictype = "ellipsoid"
nparttypes = 1 # Select only high resolution particles
haloid = 80609

#LESS IMPORTANT STUFF
hubble = 0.6711
snapnum = 63
padlist = ['p6','p7','p8','p9','p10']
reslist = ['l11']
nvirlist= ['nvir3','nvir4','nvir5','nvir6','nvir7','nvir8','nvir9']
partypelist = np.linspace(1,nparttypes-1,nparttypes)

basedir = "/spacebase/data/bgriffen/data/caterpillar/halos/halo" + str(haloid)

#AXES
ymin = 10
ymax = 10**9
xmin = 0
xmax = 100

#SIZE OF FIGURE (MAY NEED TO CHANGE FOR YOUR RES)
figx = 22
figy = 12
fig1 = plt.figure(figsize=(figx,figy))
fig1.subplots_adjust(hspace=0.1)
fig1.subplots_adjust(wspace=0.1)

nwide = len(nvirlist)
nhigh = len(padlist)

figi = 0
for res in reslist:
    for pad in padlist:
        for nvir in nvirlist:
            figi += 1
            tmppath = basedir + '/' + ictype + '/' + res + '/' + pad + '/' + nvir + '/outputs'
            filepath = tmppath + "/snapdir_0" + str(snapnum) + "/snap_0" + str(snapnum)
            
            print "+++++++++++++++++++++++++++++++++++"
            print "+PADDING:",pad,"+N*RVIR:",nvir

            ax1 = fig1.add_subplot(nhigh,nwide,figi)
            titlestr = pad + ", " + nvir

            # HERE YOU CAN LOOP OVER EACH PARTICLE TYPE - CHECK HEADER
            #for parttype in partypelist:
            #    snapPOS=rs.read_block(filepath, "POS ",parttype=int(parttype),doubleprec=False)

            s = readsubf.subfind_catalog(tmppath, snapnum)

            mgroup = s.group_m_mean200*10**10/hubble
            rvirgroup = s.group_r_mean200
            xposgroup = s.group_pos[:,0]
            yposgroup = s.group_pos[:,1]
            zposgroup = s.group_pos[:,2]

            #recall "s." in terminal to bring up catalogue properties

            #ROUGH FOF CENTER - need to generalise but depends on environment in each box :|          
            xcen = 51.24
            ycen = 48.28
            zcen = 47.30

            #SEARCH WINDOW
            deltapos = 1.5 #Mpc/h

            for inx in xrange(0,len(mgroup)):
                if mgroup[inx] > 8E11:
                    if xposgroup[inx] > xcen-deltapos and xposgroup[inx] < xcen+deltapos:
                        if yposgroup[inx] > ycen-deltapos and yposgroup[inx] < ycen+deltapos:
                            if zposgroup[inx] > zcen-deltapos and zposgroup[inx] < zcen+deltapos:
                                posxFOF = xposgroup[inx]
                                posyFOF = yposgroup[inx]
                                poszFOF = zposgroup[inx]
                                rvirFOF = rvirgroup[inx]
                                massFOF = mgroup[inx]
                                print "!FOUND FOF-GROUP CANDIDATE (MMEAN,X,Y,Z,RMEAN)!"
                                print '{0:.2e}'.format(float(massFOF)),'{:.2f}'.format(float(posxFOF)),'{:.2f}'.format(float(posyFOF)),'{:.2f}'.format(float(poszFOF)),'{:.2f}'.format(float(rvirFOF))

            ax1.set_xlim([xmin,xmax])
            ax1.set_ylim([ymin,ymax])
            ax1.set_yscale('log')

            ax1.text(0.95, 0.95,titlestr,
                horizontalalignment='right',
                verticalalignment='top',
                color='black',
                weight='bold',
                transform = ax1.transAxes)

            if figi >= 29:
                ax1.set_xlabel(r'$\mathrm{r\ [kpc/h]}$',size=14)
            else:
                ax1.set_xticklabels([])

            if figi == 1 or figi == 8 or figi == 15 or figi == 22 or figi == 29:
                ax1.set_ylabel(r'$\mathrm{\rho(r)/<\rho>}$',size=14)
            else:
                ax1.set_yticklabels([])

fig1.savefig(ictype + '-NFWprofiles.png',bbox_inches='tight')

plt.show()
