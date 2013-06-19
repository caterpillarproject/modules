import numpy as np
import pandas as pd

## Mark's reading code
import readsubf
import readgroup

class SUBFDataReader:
    """Creates Subfind Data Reader that reads in halo catalogue. It automatically reads in particles
    for each halo as well, which is an option you can turn off (no_part=1 in the constructor).

    Stores the data as a pandas library DataFrame. http://pandas.pydata.org/

    ++++++THIS SECTION IS NO LONGER TRUE++++++
    We decided to define hosts and subhalos as follows:
    * Subfind groups are considered host halos, and they are put into the halo catalogue.
    * For each group, we do NOT include the LARGEST SUBHALO. This is because the largest
      subhalo is a similar object to the FOF group, but its properties are computed without
      including the particles in its own subhalos.
    * However, we include all other subhalos. So the total number of halos in the catalogue
      should be the total number of subhalos (if every group has at least one subhalo).
    * Groups DO NOT HAVE the following standardized columns:
      rvmax, rs, vmax, vrms, spin, Jx, Jy, Jz
    * Subhalos DO NOT HAVE rvir or rs
    ++++++END SECTION THAT IS NO LONGER TRUE+++++++

    +++++BEGIN CURRENT ACCURATE DESCRIPTION OF HOSTS+++++
    Subfind outputs FOF groups and subhalos of these FOF groups. In each FOF group, there
    is a largest subhalo. We call this largest subhalo the 'host' halo for that FOF group.
    Any other subhalos of the FOF group are considered subhalos of this 'host' halo.
    (This may eventually be amended so that only subhalos within some radius of the host
    halo are considered subhalos.)
    +++++END DESCRIPTION OF HOSTS

    -----------------------------------------------------------------
    Standardized column names:
    id, hostID, mvir, npart, posX, posY, posZ, pecVX, pecVY, pecVZ,
    rvir, rvmax, rs, vmax, vrms, spin, Jx, Jy, Jz
    
    Standardized functions:
    get_hosts(): get a DataFrame with just host halos
    get_subhalos(): get a DataFrame with just subhalos
    get_subhalos_from_halo(id): get a DataFrame with just subhalos corresponding to'id'
    get_particles_from_halo(id): get particle IDs corresponding to 'id'
    -----------------------------------------------------------------

    """
    def __init__(self, basedir, snapnum, no_part=0):
        self.hftype = 'subfind'
        self.no_part = no_part
        s = readsubf.subfind_catalog(basedir, snapnum)
        sid = readsubf.subf_ids(basedir, snapnum, 0, 0, read_all=1)
        
        self.ngroups = s.ngroups
        self.nsubs   = s.nsubs
        self.nhalo = s.nsubs
        self.ncols = 19
        
        ## Information for getting particle IDs
        self.sub_offset   = s.sub_offset
        self.sub_len      = s.sub_len
        self.ids          = sid.SubIDs

        data = np.zeros((self.nhalo,self.ncols))
        data[:,0]          = np.arange(self.nhalo) #id
        data[:,1]          = np.zeros(self.nhalo)-1 ## hostID needs to be computed from the groups
        data[:,2]          = s.sub_mass*10**10 #mvir in Msun/h
        data[:,3]          = s.sub_len #npart
        data[:,[4,5,6]]    = s.sub_cm #posX/Y/Z
        data[:,[7,8,9]]    = s.sub_vel #pecVX/VY/VZ
        ## rvir not given, figure something out...
        data[:,11]         = s.sub_vmaxrad #rvmax
        ## rs not computed
        data[:,13]         = s.sub_vmax #vmax
        data[:,14]         = s.sub_veldisp #vrms
        ## spin: needs to be computed
        data[:,[16,17,18]] = s.sub_spin #Jx, Jy, Jz not normalized yet

        ## use groups to get hostID
        for j in xrange(self.ngroups):
            start = s.group_firstsub[j]
            end   = start+s.group_nsubs[j]
            data[(start+1):end,1] = start #put in the ID of the first subhalo as the host of the others

        self.data = pd.DataFrame(data,index=data[:,0].astype(int),columns=['id','hostID','mvir','npart','posX','posY','posZ','pecVX','pecVY','pecVZ','rvir','rvmax','rs','vmax','vrms','spin','Jx','Jy','Jz'])

    def get_hosts(self):
        return self.data[self.data['hostID']==-1]

    def get_subhalos(self):
        return self.data[self.data['hostID']!=-1]

    def get_subhalos_from_halo(self, id):
        return self.data[self.data['hostID']==id]

    def get_particles_from_halo(self, id):
        start = self.sub_offset[id]
        end   = start + self.sub_len[id]
        return self.ids[start:end]
