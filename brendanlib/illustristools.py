import pandas as pd
import numpy as np
import glob as glob
import os,time,h5py
# IMPORT ILLUSTRIS MODULES
import readsnapshots.readsnapHDF5 as rsnap
import mergertrees.mergertreeHDF5 as rtree
import readhalos.readsubfHDF5 as rsub
import readsnapshots.readsnapHDF5 as rsnap

def loadsnapspacing():
    filename = '/n/home01/bgriffen/work/illustris-snapshot-spacing.dat'
    snapspacing = np.loadtxt(filename)
    snapshots = snapspacing[:,0]
    expfactors = snapspacing[:,1]
    redshifts = snapspacing[:,2]
    
    return snapshots,expfactors,redshifts

def plothostSFR(treedir,subids,verbose=False):
    snapshots,expfactors,redshifts = loadsnapspacing()
    fig = plt.figure(figsize=(20,5))
    axa = fig.add_subplot(131)
    axb = fig.add_subplot(132)
    axc = fig.add_subplot(133)
    tree = rtree.TreeDB(treedir)
    for subidi in subids_cand:
        branch = tree.get_main_branch(snapshot, subidi)
        if verbose:    
            print "DOING ID:",int(subidi)
            print "SUBID: %i | SUBMASS: %3.2e" % (subidi,branch.SubhaloMassType[:,1][0]*10**10/hubble)
    
        mask = np.in1d(branch.SnapNum,snapshots,assume_unique=True)
        expfact = []
        for snap in branch.SnapNum:
            expfact.append(expfactors[np.where(snap == snapshots)])
    
        axa.plot(np.flipud(expfact),np.cumsum(np.flipud(branch.SubhaloSFR)),'-',linewidth=3)
        axb.plot(np.flipud(expfact),np.cumsum(np.flipud(branch.SubhaloMassType[:,1]*10**10/hubble)),'-',linewidth=3)
    
    axa.set_ylabel(r'$\Sigma$ SFR [Msol/yr]')
    axa.set_xlabel('Expansion Factor')
    axb.set_ylabel(r'$\Sigma$ (DM) Mass [Msol]')
    axb.set_xlabel('Expansion Factor')
    axb.set_yscale('log')
    plt.show()


def loadcandidates(outputpath,filename):
    candidates = np.loadtxt(outputpath + filename,comments='#')
    groupid = candidates[:,0]
    subid = candidates[:,1]
    groupmass = candidates[:,2]*10**10
    subhalomass = candidates[:,3]*10**10
    xposmw = candidates[:,4]
    yposmw = candidates[:,5]
    zposmw = candidates[:,6]
    return groupid,subid,groupmass,subhalomass,xposmw,yposmw,zposmw

def loadstarsfromsnap(outputpath,filename,quant,delimiter='',version='old'):
    fulldata = np.loadtxt(outputpath+filename,delimiter=delimiter)
    fulldata = np.atleast_2d(fulldata)
    if version=='new':
#	print fulldata
#	print len(fulldata)
#	print np.shape(fulldata)
        if quant == "ID":
            subdata = fulldata[...,0]
        if quant == "POS":
            subdata = fulldata[...,1:4]
        if quant == "GAGE":
            subdata = fulldata[...,5]
        if quant == "GIMA":
            subdata = fulldata[...,6]
        if quant == "GMET":
            subdata = fulldata[...,7]
        if quant == "GZ":
            subdata = fulldata[...,8:15]

    if version=='old':
        if quant == "HOSTID":
            subdata = fulldata[:,0]
        if quant == "ID":
            subdata = fulldata[:,1]
        if quant == "POS":
            subdata = fulldata[:,2:5]
        if quant == "GAGE":
            subdata = fulldata[:,6]
        if quant == "GIMA":
            subdata = fulldata[:,7]
        if quant == "GMET":
            subdata = fulldata[:,8]
        if quant == "GZ":
            subdata = fulldata[:,9:16]

    return subdata

def makesubboxinfo(basepath,subbox):
    for sim in xrange(1,4):
        simtype = 'Illustris-'+str(sim)
        arr = []
        print simtype
        print "Subbox:",subbox
        if sim == 3:
            folderlist = glob.glob(basepath + simtype + '/output/subbox' + str(subbox)+'/*.hdf5')
            snappath_subbox = basepath + simtype + '/output/subbox' + str(subbox)
            for headername in folderlist:
                snappath_subbox = headername.rstrip(".hdf5")[:-2]
                header = rsnap.snapshot_header(headername)
                snap = headername.split("_")[-1].replace(".hdf5","")
                arr.append([int(snap),header.redshift,header.time])
        else:  
            folderlist = glob.glob(basepath + simtype + '/output/subbox' + str(subbox)+'/*subbox*')
            snappath_subbox = basepath + simtype + '/output/subbox' + str(subbox)
            for folder in folderlist:
                headername = glob.glob(folder+'/*.hdf5')[0]
                snappath_subbox = headername.rstrip(".hdf5")[:-2]
                header = rsnap.snapshot_header(headername)
                snap = headername.split("_")[-1].replace(".hdf5","").split('.')[0]
                arr.append([int(snap),header.redshift,header.time])
         
        arr = np.array(arr)
        arr = arr[np.argsort(arr[:,0])]
        f = open(simtype+'-subbox-snapshotinfo.txt','w')
        for i in xrange(0,len(arr)):
            f.write(str(arr[i,0])+','+str(arr[i,1])+','+str(arr[i,2])+'\n')

        f.close()

def makesnapinfo(snappath,simtype):
    arr = []
    for snap in range(0,136):
        headername = snappath+"/snapdir_" + str(snap).zfill(3) + "/snap_" + str(snap).zfill(3) + ".0.hdf5"
        header = rsnap.snapshot_header(headername)
        arr.append([int(snap),header.redshift,header.time]) 
    
    arr = np.array(arr)
    arr = arr[np.argsort(arr[:,0])]
    f = open(simtype+'-snapshotinfo.txt','w')
    for i in xrange(0,len(arr)):
       f.write(str(int(arr[i,0]))+','+str(arr[i,1])+','+str(arr[i,2])+'\n')
    
    f.close()

def loadstarparticleblock(snappath,snapshot,block,verbose=False):
    headerpath = snappath + "/snapdir_" + str(snapshot).zfill(3) + "/snap_" + str(snapshot).zfill(3) + ".0.hdf5"
    snappath = headerpath.rstrip(".hdf5")[:-2]

    if verbose:
        print "-- Header:",headerpath

    header = rsnap.snapshot_header(headerpath)
    
    if header.nall[4] > 0:
        if block in ["ID","POS","VEL","GIMA","GMET","GZ"]:
            start = time.time()
            if verbose: print "-- Time to read star particle block..."
            starrawGAGE = rsnap.read_block(snappath,"GAGE",parttype=4)
            starGAGE = starrawGAGE[starrawGAGE >= 0]
            elapsed = time.time() - start
            if verbose: print "---- [Ages] %3.2f seconds" % (elapsed)
            start = time.time()
            starraw  = rsnap.read_block(snappath,block.ljust(4),parttype=4)
            starout = starraw[starrawGAGE >= 0]
            elapsed = time.time() - start
            if verbose: print "---- [" + block + "] %3.2f seconds" % (elapsed)
            del starrawGAGE
        
            if verbose:
                nstars = len(starraw)
                nwindstars = len(starraw) - len(starout)
                nrealstars = len(starout)
                print "---- Total Stars:",nstars
                print "------ Wind Stars: %i (%3.2f)" % (nwindstars,(nstars-nwindstars)*100./nstars)
                print "------ Real Stars: %i (%3.2f)"  % (nrealstars,(nstars-nrealstars)*100./nstars)
	if block in ["GAGE"]:
	    start = time.time()
 	    if verbose: print "-- Time to read star particle block..."
            starrawGAGE = rsnap.read_block(snappath,"GAGE",parttype=4)
            starout = starrawGAGE[starrawGAGE >= 0]
            elapsed = time.time() - start
            if verbose: print "---- [Ages] %3.2f seconds" % (elapsed)
    else:
        print "NO STAR PARTICLES AVAILABLE"
        starout=[]

    return starout

# DEFINING FUNCTIONS, EXECUTION AT BOTTOM
def createcandidatelist(simtype,outputpath,cat,verbose=False):
    hubble = 0.704
    filename = outputpath + simtype +'-mwcandidates.dat'
    if not os.path.isfile(filename):
        f=open(filename,'w')
        f.write('# groupid, subid, mass, position [Mpc/h] (x,y,z)\n')
        jval = []
        xpos = []
        ypos = []
        for j in range(0,cat.ngroups):
            if (cat.Group_M_Crit200[j]/hubble >= 90) & (cat.Group_M_Crit200[j]/hubble<=300):
                dx=(cat.GroupPos[:,0]-cat.GroupPos[j,0])/hubble
                dy=(cat.GroupPos[:,1]-cat.GroupPos[j,1])/hubble
                dz=(cat.GroupPos[:,2]-cat.GroupPos[j,2])/hubble
                dr=np.sqrt(dx**2 + dy**2 + dz**2)
                index = (dr < 1.) & (cat.Group_M_Crit200[j]/hubble > 0.5*cat.Group_M_Crit200[j]/hubble)
                if len(dr[index]) == 1:
                    #if (cat.SubhaloMassType[:,1][j]>0) and (cat.GroupNsubs[j]>0)
                    if (cat.GroupNsubs[j]>0) & (cat.Group_M_Crit200[j]>0) & (cat.Group_M_Mean200[j]>0) & (cat.Group_M_TopHat200[j]>0):
                        if verbose:
                            print "Found candidate:", j
                            #print "FOFID|SUBID|SUBMASS|SUBPOS", cat.SubhaloGrNr[j],cat.GroupFirstSub[cat.SubhaloGrNr[j]],cat.SubhaloMassType[:,1][j]/hubble, cat.SubhaloPos[j]
                            print "GROUPMASS|SUBHALOMASS: %3.2e, %3.2e" % ((cat.Group_M_Crit200[j]/hubble)*10**10,(cat.SubhaloMassType[:,1][cat.GroupFirstSub[j]])*10**10)
                            print "GROUPOS: %3.4f, %3.4f, %3.4f" % (cat.GroupPos[j,0]/1000.,cat.GroupPos[j,1]/1000.,cat.GroupPos[j,2]/1000.)
                            print "SUBPOS:  %3.4f, %3.4f, %3.4f" % (cat.SubhaloPos[cat.GroupFirstSub[j]][0]/1000.,cat.SubhaloPos[cat.GroupFirstSub[j]][1]/1000.,cat.SubhaloPos[cat.GroupFirstSub[j]][2]/1000.)
                            #print "FOFID|SUBID|SUBMASS|SUBPOS", cat.SubhaloGrNr[j],cat.GroupFirstSub[cat.SubhaloGrNr[j]],cat.Group_M_Crit200[:,1][j]/hubble, cat.SubhaloPos[j]
                        
                        f.write(str(j) + ' ' + \
                                str(cat.GroupFirstSub[j]) + ' ' + \
                                str(np.float64(cat.Group_M_Crit200[j])) + ' ' + \
                                str(np.float64(cat.SubhaloMassType[cat.GroupFirstSub[j],1])) + ' ' + \
                                str(np.float64(cat.GroupPos[j,0]/1000.)) + ' ' + \
                                str(np.float64(cat.GroupPos[j,1]/1000.)) + ' ' + \
                                str(np.float64(cat.GroupPos[j,2]/1000.)) + ' ' + \
                str(np.float64(cat.Group_R_Crit200[j])) + '\n')
                
        f.close()
    else:
        print "File Exists:",filename

def findstarparticles(simtype,outputpath,snappath,subid,xposmw,yposmw,zposmw,nrvir,verbose):
    hubble=0.704
    filename = outputpath + simtype +'-stars-snap135.dat'
    if not os.path.isfile(filename):
        enclosedradii = nrvir #Mpc
    
        if verbose:
            print 
            print "> READING SNAPSHOTS"
            print "-- From:",snappath
    
        starttotal = time.time()
        starIDS,starPOS,starGAGE,starIMA,starGMET,starGZ = loadstarparticles(outputpath,snappath,135,verbose=verbose)
        elapsedtotal = time.time() - starttotal

        if verbose:
            print "-- Time to get star particles at z0: %3.2f minutes" % (elapsedtotal/60.)

        f=open(filename,'w')
        f.write('# SUBIDZ0 (0), STARID(1), POSITION(2-4), AGE(5), INITIAL MASS(6), METALICITY(7), METALS(8-16), \n')
        for i in range(0,len(subid)):
            dx = (xposmw[i]-starPOS[:,0])/hubble
            dy = (yposmw[i]-starPOS[:,1])/hubble
            dz = (zposmw[i]-starPOS[:,2])/hubble
            dr = np.sqrt(dx**2 + dy**2 + dz**2)
            mask = (dr <= enclosedradii)

            if mask.any():
                starids_insnap = starIDS[mask]
                #del starIDS
                starpos_insnap = starPOS[mask]
                #del starPOS
                stargage_insnap = starGAGE[mask]
                #del starGAGE
                starima_insnap = starIMA[mask]
                #del starIMA
                stargmet_insnap = starGMET[mask]
                #del starGMET
                stargz_insnap = starGZ[mask]
                #del starGZ
    
                for staridi,starposi,stargagei,starimai,stargmeti,stargzi in \
                    zip(starids_insnap,starpos_insnap,stargage_insnap,starima_insnap,stargmet_insnap,stargz_insnap):
                    f.write(str(np.uint64(subid[i]))+ ' ' + \
                    str(np.uint64(staridi)) + ' ' + \
                    str(np.float64(starposi[0]/1000.)) + ' ' + \
                    str(np.float64(starposi[1]/1000.)) + ' ' + \
                    str(np.float64(starposi[2]/1000.))  + ' ' + \
                    str(np.float64(stargagei)) + ' ' + \
                    str(np.float64(starimai)) + ' ' + \
                    str(np.float64(stargzi)) + ' ' + \
                    str(np.float64(stargmeti[0])) + ' ' + \
                    str(np.float64(stargmeti[1])) + ' ' + \
                    str(np.float64(stargmeti[2])) + ' ' + \
                    str(np.float64(stargmeti[3])) + ' ' + \
                    str(np.float64(stargmeti[4])) + ' ' + \
                    str(np.float64(stargmeti[5])) + ' ' + \
                    str(np.float64(stargmeti[6])) + ' ' + \
                    str(np.float64(stargmeti[7])) + ' ' + \
                    str(np.float64(stargmeti[8])) + '\n')
    
        f.close()
    else:
        print "File Exists:",filename

def tracestars(outputpath,snappath,simtype,verbose):
    filename = simtype +'-stars-snap135.dat'
    #hostIDSz0,starIDSz0,starposXz0,starposYz0,starposZz0 = loadstars(outputpath,filename)
    print outputpath
    print filename
    hostIDSz0 = loadstarsfromsnap(outputpath,filename,quant="HOSTID")
    starIDSz0 = loadstarsfromsnap(outputpath,filename,quant="ID")
    #posdataz0 = loadstarsfromsnap(outputpath,filename,quant="POS")
    #starposXz0 = posdataz0[0]
    #starposYz0 = posdataz0[1]
    #starposZz0 = posdataz0[2]
    for snapshot in range(56,135):
        filename = outputpath + simtype +'-stars-snap'+str(snapshot).zfill(3)+'.dat'
        if not os.path.isfile(filename):
            starIDS,starPOS,starGAGE,starIMA,starGMET,starGZ = loadstarparticles(outputpath,snappath,snapshot,verbose)
            if len(starIDS)>0:
                mask = np.in1d(starIDS,starIDSz0,assume_unique=True)
    
                if mask.any():
                    starids_insnap = starIDS[mask]
                    starpos_insnap = starPOS[mask]
                    stargage_insnap = starGAGE[mask]
                    starima_insnap = starIMA[mask]
                    stargmet_insnap = starGMET[mask]
                    stargz_insnap = starGZ[mask]
                    maskz0 = np.in1d(starIDSz0,starids_insnap,assume_unique=True)
#                    maskz0 = np.in1d(starIDSz0,starids_insnap)
                    #print "Length of z0 mask:",len(maskz0)
                    print "--> length of maskz0:",len(maskz0)
                    print "--> number of 'true' in maskz0:",np.count_nonzero(maskz0)
                    print "--> length of starids_insnap:",len(starids_insnap)
                    print "--> length of starIDSz0:",len(starIDSz0)

                    subhostid = hostIDSz0[maskz0]
                    
                    

                    f=open(filename,'w')
                    f.write('# SUBIDZ0 (0), STARID(1), POSITION(2-4), AGE(5), INITIAL MASS(6), METALICITY(7), METALS(8-16), \n')
                    
                    if verbose:    
                        print
                        print "> DOING: ",filename
                        print "-- # Stars:",len(starIDS),"| # Candidate Stars:",len(starPOS[mask])
                        
                        for staridi,starposi,stargagei,starimai,stargmeti,stargzi,subid in \
                            zip(starids_insnap,starpos_insnap,stargage_insnap,starima_insnap,stargmet_insnap,stargz_insnap,subhostid):
                            f.write(str(np.uint64(subid))+ ' ' + \
                            str(np.uint64(staridi)) + ' ' + \
                            str(np.float64(starposi[0]/1000.)) + ' ' + \
                            str(np.float64(starposi[1]/1000.)) + ' ' + \
                            str(np.float64(starposi[2]/1000.))  + ' ' + \
                            str(np.float64(stargagei)) + ' ' + \
                            str(np.float64(starimai)) + ' ' + \
                            str(np.float64(stargzi)) + ' ' + \
                            str(np.float64(stargmeti[0])) + ' ' + \
                            str(np.float64(stargmeti[1])) + ' ' + \
                            str(np.float64(stargmeti[2])) + ' ' + \
                            str(np.float64(stargmeti[3])) + ' ' + \
                            str(np.float64(stargmeti[4])) + ' ' + \
                            str(np.float64(stargmeti[5])) + ' ' + \
                            str(np.float64(stargmeti[6])) + ' ' + \
                            str(np.float64(stargmeti[7])) + ' ' + \
                            str(np.float64(stargmeti[8])) + '\n')
                  
                    f.close()
        else:
           print "File Exists:",filename

def gettree(fileBase,snapNum,subhaloID,NtreeFiles=4096,verbose=False,getmain=False):
    """ Return the full progenitor tree for a given subhalo at a given snapshot. """
    result = {}
 
    # load the tree number/subindex of this subgroup within the mergertree chunks
    filePath = fileBase + '/postprocessing/offsets/tree_offsets_subgroup_'+str(snapNum)+'_135.hdf5' 

    if verbose: 
	start = time.time()
	
    f = h5py.File(filePath,'r')
 
    result['treeFileNum'] = f["TreeFile"][ subhaloID ]
    result['treeNum']     = f["TreeNum"][ subhaloID ]
    result['treeIndex']   = f["TreeIndex"][ subhaloID ]
 
    f.close()

    if verbose:
	elapsed = time.time() - start
 	print "Reading offset file: %3.2f" % (elapsed)

    if result['treeFileNum'] < 0 or result['treeFileNum'] >= NtreeFiles:
        return[] # subhalo not in tree, and/or general error
 
    # walk tree starting at this subgroup
    filePath = fileBase + 'trees/treedata/trees_sf1_' + str(snapNum).zfill(3)
    filePath += '.' + str(result['treeFileNum']) +'.hdf5'
 
    if verbose:
	start = time.time()

    f = h5py.File(filePath,'r')
    fTree = f['Tree'+str(result['treeNum'])]
 
    masses_all = []
    masses_gas = []
    masses_dm = []
    masses_stars = []
    masses_bh = []
    subhaloid = []
    posx = []
    posy = []
    posz = []
    velx = []
    vely = []
    velz = []
    snapnums = []
    halfmassr = []
    halotype = []
    index = result['treeIndex']
    masses_all.append(np.sum(fTree['SubhaloMassType'][index,:]))
    masses_gas.append(fTree['SubhaloMassType'][index,0])
    masses_dm.append(fTree['SubhaloMassType'][index,1])
    masses_stars.append(fTree['SubhaloMassType'][index,4])
    masses_bh.append(fTree['SubhaloMassType'][index,5])
    subhaloid.append(fTree['SubhaloNumber'][index])
    posx.append(fTree['SubhaloPos'][index,0])
    posy.append(fTree['SubhaloPos'][index,1])
    posz.append(fTree['SubhaloPos'][index,2])
    velx.append(fTree['SubhaloVel'][index,0])
    vely.append(fTree['SubhaloVel'][index,1])
    velz.append(fTree['SubhaloVel'][index,2])
    halfmassr.append(fTree['SubhaloHalfmassRadType'][index,1])
    snapnums.append(fTree['SnapNum'][index])
    halotype.append(0)

    def recProgenitorList( fTree, index ): 
        firstProg = fTree['FirstProgenitor'][index]
        progs = []

        if firstProg == -1:
            return progs
 
        subhaloid.append(fTree['SubhaloNumber'][firstProg])
        masses_all.append(np.sum(fTree['SubhaloMassType'][firstProg,:]))
        masses_gas.append(fTree['SubhaloMassType'][firstProg,0])
        masses_dm.append(fTree['SubhaloMassType'][firstProg,1])
        masses_stars.append(fTree['SubhaloMassType'][firstProg,4])
        masses_bh.append(fTree['SubhaloMassType'][firstProg,5])
        posx.append(fTree['SubhaloPos'][firstProg,0])
        posy.append(fTree['SubhaloPos'][firstProg,1])
        posz.append(fTree['SubhaloPos'][firstProg,2])
        velx.append(fTree['SubhaloVel'][firstProg,0])
        vely.append(fTree['SubhaloVel'][firstProg,1])
        velz.append(fTree['SubhaloVel'][firstProg,2])
        snapnums.append(fTree['SnapNum'][firstProg])
        halfmassr.append(fTree['SubhaloHalfmassRadType'][firstProg,1])
        halotype.append(1)

        recProgenitorList( fTree, firstProg )
        if getmain == False:
            nextProg = fTree['NextProgenitor'][firstProg]
            while nextProg >= 0:
                #print nextProg
                subhaloid.append(fTree['SubhaloNumber'][nextProg])
                masses_all.append(np.sum(fTree['SubhaloMassType'][nextProg,:]))
                masses_gas.append(fTree['SubhaloMassType'][nextProg,0])
                masses_dm.append(fTree['SubhaloMassType'][nextProg,1])
                masses_stars.append(fTree['SubhaloMassType'][nextProg,4])
                masses_bh.append(fTree['SubhaloMassType'][nextProg,5])
                posx.append(fTree['SubhaloPos'][nextProg,0])
                posy.append(fTree['SubhaloPos'][nextProg,1])
                posz.append(fTree['SubhaloPos'][nextProg,2])
                velx.append(fTree['SubhaloVel'][nextProg,0])
                vely.append(fTree['SubhaloVel'][nextProg,1])
                velz.append(fTree['SubhaloVel'][nextProg,2])
                snapnums.append(fTree['SnapNum'][nextProg])
                halfmassr.append(fTree['SubhaloHalfmassRadType'][nextProg,1])
                halotype.append(2)
                recProgenitorList( fTree, nextProg )
                nextProg = fTree['NextProgenitor'][nextProg]

#                recProgenitorList( fTree, nextProg )
#                nextProg = fTree['NextProgenitor'][nextProg]

        return [snapnums,subhaloid,posx,posy,posz, \
            velx,vely,velz,halfmassr,masses_all,masses_gas,masses_dm,masses_stars,masses_bh,halotype]

    if verbose:
        elapsed = time.time() - start
	print "Time to walk tree: %3.2f" % (elapsed)

    varlist = recProgenitorList( fTree, result['treeIndex'] )
    varnames = ['snapnum','subhaloid','posx','posy','posz','velx','vely','velz','halfmassr','mall','mgas','mdm','mstar','mbh','halotype']

    if varlist != []:
        return pd.DataFrame(np.flipud(np.array(varlist).T),columns=varnames)
    else:
        return varlist

def loadSnapSubset(fileBase,snapNum,searchID,getGroup,partType,fields):
    """ Return requested fields for one particle type for all members of group/subgroup. """
 
    groupName = "PartType" + str(partType)
    dsetName = "Subhalo"
    if getGroup > 0:
        dsetName = "Group"
 
    # load the length (by type) of this group/subgroup from the group catalog
    filePath = fileBase + '/postprocessing/offsets/offsets_' + dsetName.lower() + '_'+str(snapNum)+'.npy'
    offsets = np.load(filePath)
 
    offsets = searchID - offsets
    fileNum = np.max( np.where(offsets >= 0) )
    groupOffset = offsets[fileNum]
 
    filePath = fileBase + 'groups_' + str(snapNum).zfill(3) + '/'
    filePath += 'fof_subhalo_tab_' + str(snapNum).zfill(3) + '.' + str(fileNum) + '.hdf5'
 
    f = h5py.File(filePath,'r')
    lenType = f[dsetName][dsetName+"LenType"][groupOffset]
    f.close()
 
    # load the offset (by type) of this group/subgroup within the snapshot chunks
    filePath = fileBase + '/postprocessing/offsets/snap_offsets_' + dsetName.lower() + '_'+str(snapNum)+'.hdf5'
 
    f = h5py.File(filePath,'r')
    offsetType = f["Offsets"][ searchID ]
    f.close()
 
    # load the offsets for the snapshot chunks
    filePath = fileBase + '/postprocessing/offsets/offsets_snap_'+str(snapNum)+'.npy'
    offsets = np.load(filePath)
 
    # determine first snapshot chunk we need to load for this type
    wOffset = 0
    result  = {}
 
    offsetsThisType = offsetType[partType] - offsets[partType,:]
    fileNum         = np.max( np.where(offsetsThisType >= 0) )
    fileOffset      = offsetsThisType[fileNum]
 
    numLeftToRead = lenType[partType]
 
    while numLeftToRead:
 
        # loop over files, for each load the overlapping chunk into the hdf5 file
        curSnapFilePath = fileBase + 'snapdir_' + str(snapNum).zfill(3) + '/'
        curSnapFilePath += 'snap_' + str(snapNum).zfill(3) + '.' + str(fileNum) + '.hdf5'
 
        fSnap = h5py.File(curSnapFilePath,'r')
 
        # set local read length for this file
        readLen = numLeftToRead
 
        if fileNum < offsets.shape[1]-1: # if in last file, assume group is entirely contained in this file
            # if the local length after requested offset is less than the read length, modify read length
            if fileOffset+readLen+offsets[partType,fileNum]-1 >= offsets[partType,fileNum+1]:
                readLen = (offsets[partType,fileNum+1]-offsets[partType,fileNum]) - fileOffset
 
        # loop over each requested field for this particle type
        for fieldName in fields:
            # shape and type
            dtype = fSnap[groupName][fieldName].dtype
            shape = fSnap[groupName][fieldName].shape
 
            # read data local to the current file (allocate dataset if it does not already exist)
            if len(shape) == 1:
                if fieldName not in result:
                    result[fieldName] = np.zeros( (lenType[partType],), dtype=dtype )
                result[fieldName][wOffset:wOffset+readLen] = fSnap[groupName][fieldName][fileOffset:fileOffset+readLen]
            else:
                if fieldName not in result:
                    result[fieldName] = np.zeros( (lenType[partType],shape[1]), dtype=dtype )
                result[fieldName][wOffset:wOffset+readLen,:] = fSnap[groupName][fieldName][fileOffset:fileOffset+readLen,:]
 
        wOffset += readLen
        numLeftToRead -= readLen
        fileNum += 1
        fileOffset = 0
 
        fSnap.close()
 
    # reads across all files done, for all fields of this type
    return result
