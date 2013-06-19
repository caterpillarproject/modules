"""
Merger Tree data structure to navigate through halo mergers.
@author: Greg Dooley
@date: 4/4/2012
"""

class MTHalo:
    """
    Merger Tree Halo. Contains links to parents, descendants unlike simple Halo object.
    @param line: String of halo information 
    """
    def __init__(self, line):
        # parse string into array of strings
        split = line.split()
        
        # Halo parameters
        self.scale = float(split[0]) #: cosmological scale factor of halo
        self.ID = int(split[1]) #: ID of halo
        self.desc_scale = float(split[2]) #: Scale of descendant halo, if applicable
        self.descid = int(split[3]) #: ID of descendent halo, if applicable
        self.num_prog = int(split[4]) #: number of pregenitors
        self.pid = int(split[5]) #: Host halo ID (-1 if distinct halo)
        self.upid = int(split[6]) #: most massive host halo ID
        self.desc_pid = int(split[7]) #: Pid of descendant halo
        self.phantom = int(split[8]) #: nonzero for halos interpolated across timesteps
        self.Mvir = float(split[9]) #: Halo mass, smoothed across accretion history; aloways greater than sum of halo masses of contributing progenitors (Msun/h)
        self.orig_Mvir = float(split[10]) #: original halo mass from raw halo catalogs
        self.Rvir = float(split[11]) #: Halo radius (kpc/h comoving)
        self.Rs = float(split[12]) #: scale radius (kpc/h comoving)
        self.Vrms = float(split[13]) #: Velocity dispersion (km/s physical)
        self.mmp = int(split[14]) #: whether the halo is the most massive progenitor or not
        self.scale_of_last_MM = float(split[15]) #: scale factor of the last major merger (Mass ratio > 0.3)
        self.Vmax = float(split[16]) #: Maximum circular velocity (km/s physical)
        self.posn = [float(split[17]), float(split[18]), float(split[19])] #: Halo position (Mpc/h comoving). 3 components
        self.velopcity = [float(split[20]), float(split[21]), float(split[22])] #: Halo velocity (km/s physical). 3 components
        self.L = [float(split[23]), float(split[24]), float(split[25])] #: Halo angular momentum (Msun/h * Mpc/h * km/s physical). 3 components
        self.spin = float(split[26]) #: Halo spin parameter
        self.TreeRootID = int(split[29]) #: ID of the halo at the last timestep in the tree. (redundant info in MT already)
        self.origHaloID = int(split[30]) #: original halo ID from halo finder

        
        self.desc = None #: descendent halo
        self.parent = None #: Most massive parent halo
        self.subhalos = [] #: list of all subhalos
        self.fullParents = [] #: list of all parents

    def getParent(self):
        """
        @return: the most massive parent halo
        """
        return self.parent
    def getDescendant(self):
        """
        @return: the descendant halo. None if no descendant.
        """
        return self.desc

    def containsParticle(self, particle, dir):
        """
        Determine if halo contains particle
        @param particle: ID of particle from Gadget
        @param dir: directory containing Halo Finder files
        @return: True or False
        """
        #need to find halo corresponding to origHaloID
        #self.origHaloID
        #self.scale
        return False
