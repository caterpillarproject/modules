import numpy as np
import pandas as pd 

import struct

### IMPORTANT
### Must use the ahf_io.c that Alex modified for this all to work! This is in the Amiga src
### Two primary modifications:
### 1) I fixed bugs in the binary halo catalogue output.
###    (it now writes the correct numbers out) This has been fixed in
###    recent versions of AHF.
### 2) I removed the padding bit in the binary particle catalogue output.
###    This is necessary because Python's "struct" assumes the data is aligned when unpacking
###    Unfortunately when words come in chunks of 5 bits, it gets upset (i.e. reading "ixi" needs 12 characters, not 9)
###    So instead of finding out a way to deal with it in python, I just modified the amiga source

class AHFDataReader:
    """
    Creates AHF Data Reader that reads in halo catalogue. It automatically reads in particles
    for each halo as well, which is an option you can turn off (no_part=1 in the constructor).

    Stores the data as a pandas library DataFrame. http://pandas.pydata.org/

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

    Usage examples:
    #reading data
    a = AHFDataReader.AHFDataReader(binary_halo_file,binary_part_file)

    #accessing columns
    all_masses = a.data['mvir'] #column named mvir
    all_positions = np.array(a.data[['posX','posY','posZ']])
    #a.data[0] #error, no column named 0

    #accessing rows
    a.data.ix[0] #pandas Series of the halo with ID 0, NOT the 0th row!
    a.data[[0,2,4]] #DataFrame with three rows, all columns for halos with ID 0, 2, 4. If one of these IDs doesn't exist, will have nan
    a.data[0:2] #DataFrame with first two rows and all columns, NOT 
    a.data[a.data['mvir'] > 1e12] #DataFrame of halos with M>1e12
    """
    def __init__(self, halo_file, particle_file, no_part=0, file=""):
        self.hftype = 'amiga'
        self.no_part = no_part
        f = open(halo_file,'rb')
        c = f.read(8+4)
        self.nhalo, self.ncols = struct.unpack("Qi",c)
        if (self.ncols != 43):
            print "Unexpected number of columns. Expected 43, Receieved " + str(self.ncols)
            print "Output will probably be nonsense..."
        #this variable is taken manually from reading an ascii file output
        #names of columns have been changed to work with RSDataReader
        datatypes = np.dtype([('id', '<i8'), ('hostID', '<i8'),\
                              ('numSubStruct(3)', '<i8'),('mvir', '<f8'),\
                              ('npart', '<i8'),('posX', '<f8'),\
                              ('posY', '<f8'),('posZ', '<f8'),\
                              ('pecVX', '<f8'),('pecVY', '<f8'),\
                              ('pecVZ', '<f8'),('rvir', '<f8'),\
                              ('rvmax', '<f8'),('rs', '<f8'),\
                              ('mbp_offset(15)', '<f8'),('com_offset(16)', '<f8'),\
                              ('vmax', '<f8'),('v_esc(18)', '<f8'),\
                              ('vrms', '<f8'),('lambda(20)', '<f8'),\
                              ('spin', '<f8'),('Jx', '<f8'),\
                              ('Jy', '<f8'),('Jz', '<f8'),\
                              ('b(25)', '<f8'),('c(26)', '<f8'),\
                              ('Eax(27)', '<f8'),('Eay(28)', '<f8'),\
                              ('Eaz(29)', '<f8'),('Ebx(30)', '<f8'),\
                              ('Eby(31)', '<f8'),('Ebz(32)', '<f8'),\
                              ('Ecx(33)', '<f8'),('Ecy(34)', '<f8'),\
                              ('Ecz(35)', '<f8'),('ovdens(36)', '<f8'),\
                              ('nbins(37)', '<i8'),('fMhires(38)', '<f8'),\
                              ('Ekin(39)', '<f8'),('Epot(40)', '<f8'),\
                              ('SurfP(41)', '<f8'),('Phi0(42)', '<f8'),\
                              ('cNFW(43)', '<f8')])
        datatypesstr = "iiifi" + ("f" * 32) + "ifffff"
        numbytes = 43*4
        data = np.zeros(self.nhalo,dtype=datatypes)
        for i in xrange(self.nhalo):
            c = f.read(numbytes)
            data[i] = struct.unpack(datatypesstr,c)
        f.close()
        #Convert positions from kpc to Mpc
        data['posX'] = data['posX']/1000.0
        data['posY'] = data['posY']/1000.0
        data['posZ'] = data['posZ']/1000.0
        self.data = pd.DataFrame(data,index=data['id'])
        print "finished loading " + halo_file
        
        #Read particles unless the no_part flag is set
        if (no_part == 0):
            self.part = []
            f = open(particle_file,'rb')
            #header
            c = f.read(20) 
            nhalo, npart_tot = struct.unpack("QQxxxx",c)
            if (nhalo != self.nhalo):
                print "ERROR: number of halos different between halo file and particle file"
                print "You have likely picked two files associated with different boxes"
                print "Not quitting but program might break"
            # Read in all the particles
            for id in xrange(self.nhalo):
                c = f.read(8)
                npart = struct.unpack("Q",c)[0] 
                c = f.read(npart*4) #skip the actual particle numbers
                self.part.append(np.array(struct.unpack(npart*"i",c)))
            self.part = pd.Series(self.part,index=self.data.index)
            f.close()
            print "finished loading " + particle_file
	
    def get_hosts(self):
        """
        Return a pandas DataFrame that only contains the host halos.
        This is defined as halos that have no parents.
        """
        mask = self.data['hostID'] == -1
        return self.data[mask]
	
    def get_subhalos(self):
        """
        Return a pandas DataFrame that only contains the subhalos. This is
        defined as a halo that has a parent assigned to it.
        """
        mask = self.data['hostID'] != -1
        return self.data[mask]
	
    def get_subhalos_from_halo(self, id):
        """
        Given a single host halo ID, returns a DataFrame of the 
        subhalos of the host halo. Returns an empty array if there
        are no subhalos.
        """
        mask = (id == self.data['hostID'])
        return self.data[mask]
	
    def get_particles_from_halo(self, id):
        """
        Given a single halo id, returns a numpy array of the particle IDs
        associated with that halo. Returns empty numpy array if any errors.
        """
        if (self.no_part != 0):
            raise RuntimeError('AHF Data Reader: Did not read in catalogue particles')
        if (id >= self.nhalo or id < 0):
            return np.array([])
        return np.array(self.part[id])
