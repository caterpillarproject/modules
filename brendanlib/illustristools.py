
def loadsnapspacing():
    filename = './illustris-snapshot-spacing.dat'
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