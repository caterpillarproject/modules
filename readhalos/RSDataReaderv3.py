"""
By: Brendan Griffen, Alex Ji, Greg Dooley (MIT)
Contact: brendan.f.griffen@gmail.com

Read in all data from Rockstar Halo Finder into a large matrix.
Rows correspond to halos, columns correspond to halo properties.
Much faster access to data than object approach.

RSDataReaderv2 reads in Rockstar version 0.99.9
use RSDataReader for previous Rockstar version.
"""

import numpy as np
import os
import sys
import copy
import operator
import readParentsList as rp
import pandas
import time
import struct

varlist = ['id','posX','posY','posZ','corevelx','corevely','corevelz', \
          'pecVX','pecVY','pecVZ', 'bulkvelX','bulkvelY','bulkvelZ',\
           'mgrav','rvir','child_r','dummy1','mvir','vmax','rvmax','rs','rs_klypin',\
          'vrms','Jx','Jy','Jz','Epot','spin','dummy2','m200c','m500c','m2500c', \
          'Xoff','Voff','b_to_a','c_to_a','A[x]','A[y]','A[z]',
           'b_to_a_2','c_to_a_2','A2[x]','A2[y]','A2[z]',
           'spin_bullock','T/|U|','npart',\
          'num_cp','numstart','desc','flags','n_core',\
           'min_pos_err','min_vel_err','min_bulkvel_err','dummy3', \
          'hostID','offset','particle_offset']

          #'dummy8','dummy9','dummy10','dummy11','dummy12','dummy13','dummy14','dummy15','dummy16','dummy17','dummy18',\
          #'dummy19','dummy20','dummy21']


          #,'npart',\


          #'num_child_p','p_start','desc','flags','n_core','min_pos_err','min_vel_err',\
          #'min_bulkvel_err','hostID','offset','particle_offset']


BinaryHeaderSize = 256 #: in bytes
HeaderInfoSize = 96  #: in bytes
HaloSize = 256
#HaloSize = 176 + 20  #: in bytes
ParticleSize = 8 #: in bytes
KpcToMpc = 1.0/1000

#######################
#numextras = 10
#nfloats = 12+33
#nfloats = len(varlist) - numextras - 3  #33
# hostID, offset, particle_offset added later. not read from file.
#print "nfloats:",nfloats
#datatypesstr = "q"+("f" * nfloats)+"qqqqqqfff" #change
#numbytes = datatypesstr.count('q')*8 + datatypesstr.count('f')*4
datatypesstr = "qfffffffffffffffffffffffffffffffffffffffffffffqqqqqqffff"
numbytes = 256
#numbytes = (nfloats+3)*4+(numq+2)*8+4 #change
#print varlist
#print len(varlist)
#print datatypesstr
#print len(datatypesstr)
#print "numbytes:",numbytes
num_columns = len(varlist) #change
#print "num_columns:",num_columns
#######################

id = 0 
mvir = 17
npart = varlist.index('npart')
hostID = num_columns - 3 #change
offset = num_columns - 2 #change
particle_offset = num_columns - 1 #change


class RSDataReader:    
    """
    Create a Halo Catalogue for a particular snapshot
    """
    def __init__(self, dir, snap_num, base='halos_', digits=2, AllParticles=False,time_me=False):
        start_time = time.clock()
        """
         @param dir: base directory containing binary files from Rockstar Halo Finder. ex: \"/home/gdooley/Rockstar-0.99/Output2\"
         @param snap_num: snap number to be viewed. ex: 50.
         @param digits: number of digits in file name of snapshot. ex: halos_063.0 should have digits = 3.
         """
        self.dir = dir
        self.snap_num = snap_num
        self.AllParticles = AllParticles
        # open first file, test if it exists
        file_num = 0 
        file_name = dir + '/' + base + str(snap_num).zfill(digits) + '/' + base + str(snap_num).zfill(digits) + "." + str(file_num) + ".bin"
        if (file_num == 0 and not os.path.exists(file_name)):
            print "ERROR: file not found", file_name
            sys.exit()

        ## establish total number of halos in order to create matrix of the correct size
        num_rows = 0
        self.total_particles = 0
        while os.path.exists(file_name):
            f = open(file_name)
            # Read in Halo Header
            magic = np.fromfile(f, np.uint64, count = 1)[0]
            snap = np.fromfile(f, np.int64, count = 1)[0]
            chunk = np.fromfile(f, np.int64, count = 1)[0]
            scale = np.fromfile(f, np.float32, count = 1)[0] #: cosmological scale factor of snapshot
            Om = np.fromfile(f, np.float32, count = 1)[0]
            O1 = np.fromfile(f, np.float32, count = 1)[0]
            h0 = np.fromfile(f, np.float32, count = 1)[0]
            bounds = np.fromfile(f, np.float32, count = 6)
            num_halos = np.fromfile(f, np.int64, count = 1)[0]
            num_rows += num_halos
            num_particles = np.fromfile(f, np.int64, count = 1)[0]
            self.total_particles+=num_particles
            #increment to next file
            f.close()
            file_num += 1
            file_name = dir+ '/' + base + str(snap_num).zfill(digits)+'/'+base+str(snap_num).zfill(digits)+"."+str(file_num)+".bin"
        self.num_halos = num_rows
        #reset file name
        file_num = 0
        file_name = dir+ '/' + base + str(snap_num).zfill(digits) +'/'+base+str(snap_num).zfill(digits)+"."+str(file_num)+".bin"
        ### Two important pieces of data storage
        data = np.zeros((num_rows,num_columns)) #: Matrix of all halos and all information
        string_len = len(file_name)*2
        datatype = '|S'+str(string_len)
        files = np.array(['']*num_rows,dtype=datatype)
        self.particles = np.array([])
        ###
        
        i = 0 # counter for halos
        particleID_start2 = 0 # for method of reading all particles at once
        # open up all .bin files and read data.
        while os.path.exists(file_name):
            f = open(file_name)
            # Read in Halo Header
            magic = np.fromfile(f, np.uint64, count = 1)[0]
            self.snap = np.fromfile(f, np.int64, count = 1)[0]
            chunk = np.fromfile(f, np.int64, count = 1)[0]
            self.scale = np.fromfile(f, np.float32, count = 1)[0] #: cosmological scale factor of snapshot
            self.Om = np.fromfile(f, np.float32, count = 1)[0]
            self.O1 = np.fromfile(f, np.float32, count = 1)[0]
            self.h0 = np.fromfile(f, np.float32, count = 1)[0]
            bounds = np.fromfile(f, np.float32, count = 6)
            num_halos = np.fromfile(f, np.int64, count = 1)[0]
            num_particles = np.fromfile(f, np.int64, count = 1)[0]
            self.box_size = np.fromfile(f, np.float32, count = 1)[0]
            self.particle_mass = np.fromfile(f, np.float32, count = 1)[0]
            self.particle_type = np.fromfile(f, np.int64, count = 1)[0]
            # read unused space
            np.fromfile(f,'f',count = (BinaryHeaderSize-HeaderInfoSize)/4)

            particleID_start = BinaryHeaderSize + HaloSize*num_halos

            # Produce array of Halo objects
            for j in range(0,num_halos):
                #line = f.read(numbytes)
                line = f.read(256)
                tmpinx = num_columns -3
                data[i,0:tmpinx] = struct.unpack(datatypesstr, line) #change
                # info to read particle IDs for halo
                data[i,offset] = particleID_start
                #print particleID_start2, data[i,npart], particleID_start2+data[i,npart]
                data[i,particle_offset] = particleID_start2
                files[i] = file_name
                particleID_start += ParticleSize*data[i,npart]
                particleID_start2 += data[i,npart]
                #print particleID_start2, j, ' j listed last'
                i += 1
            # Read in and store particle information (optional)
            if AllParticles:
                self.particles = np.concatenate((self.particles,np.fromfile(f,np.int64)))
            #increment to next file
            f.close()
            file_num += 1
            file_name = dir+ '/' + base + str(snap_num).zfill(digits)+"/halos_"+str(snap_num).zfill(digits)+"."+str(file_num)+".bin"

        sortedIndices = data[:,mvir].argsort()[::-1]
        data = data[sortedIndices]
        files = files[sortedIndices]
        
        np.save('tmp',np.array(data[:,0]))
        self.files = pandas.DataFrame(files,index=data[:,id].astype(int),columns=['file'])

        #print data.shape, 'shape of data'
        #print data[:,id].shape
        
        self.data = pandas.DataFrame(data,index=data[:,id].astype(int),columns=varlist) #change
        #print self.data
        #sys.exit()
        self.particles = self.particles.astype(int)
        ## add column of hostID
        file = 'parents'+'.list'
        parents = rp.readParents(self.dir+'/'+base+str(snap_num).zfill(digits),file, self.num_halos)
        self.data['hostID'].ix[parents[:,0]] = parents[:,1] # fill in hostID column
        #print self.data['id'].ix[0], 'ix method'
        #print self.data['id'][0], 'no ix method'
        if time_me:
            print time.clock()-start_time, ' time'
        
    def get_particles_from_halo(self, haloID):
        """
        @param haloID: id number of halo. Not its row position in matrix
        @return: a list of particle IDs in the Halo

        3/25 change: if AllParticles is on, haloID can be a vector.
        """
        if self.AllParticles:
            if type(haloID) == list or type(haloID) == np.ndarray:
                return np.array([ self.particles[self.data['particle_offset'].ix[ID]: self.data['particle_offset'].ix[ID]+self.data['npart'].ix[ID]] for ID in haloID])
            else:
                start = self.data['particle_offset'].ix[haloID]
                end = start+self.data['npart'].ix[haloID]
                return self.particles[start:end]
        else:
            f = open(self.files['file'].ix[haloID])
            np.fromfile(f,'c',count=int(self.data['offset'].ix[haloID]))
            particleIDs = np.fromfile(f,np.int64,count=int(self.data['npart'].ix[haloID]))
            f.close()
            return particleIDs

    # Retrieve subhalos only one level deep.
    # Does not get sub-sub halos, etc.
    def get_subhalos_from_halo(self,haloID):
        return self.data[self.data['hostID']==haloID]
    # vectorize this at some point. 

    # for multiple halos, returns the subhalos of each halo as an array of arrays.
    def get_subhalos_from_halos(self,haloIDs):
        """
        returns an array of pandas data frames of subhalos. one data frame
        for each host halo. returns only first level of subhalos.
        """
        if type(haloIDs) == list or type(haloIDs) == np.ndarray:
            return np.array([ self.data[self.data['hostID']==ID] for ID in haloIDs])
        else:
            return self.data[self.data['hostID']==haloIDs] 

    def get_subhalos_from_halos_flat(self,haloIDs):
        """
        returns a flattened pandas data frame of all subhalos within
        the hosts given by haloIDs. Returns only first level of subhalos.
        """
        subs = self.get_subhalo_ids_from_halos(haloIDs)
        subIDs = [item for sublist in subs for item in sublist]
        return self.data.ix[subIDs]
    
    def get_hosts(self):
        return self.data[self.data['hostID']==-1]

    def get_subs(self):
        return self.data[self.data['hostID']!=-1]

    def get_subhalo_ids_from_halos(self,haloIDs):
        if type(haloIDs) == list or type(haloIDs) == np.ndarray:
            return np.array([ np.array(self.data[self.data['hostID']==ID]['id']) for ID in haloIDs])
        else:
            return np.array([self.data[self.data['hostID']==haloIDs]['id']])

    # Retrieve all subhalos: sub and sub-sub, etc. 
    def get_all_subs_recurse(self,haloID):
        # just need mask of all subhalos, then return data frame subset
        subs = self.get_subhalo_ids_from_halos(haloID)
        subs = [item for sublist in subs for item in sublist]
        if subs == []:
            return [] 
        else:
            return np.append(subs,self.get_all_subs_recurse(subs))
                                       
    # Retrieve all subhalos: sub and sub-sub, etc.
    # return pandas data fram of subhalos
    def get_all_subhalos_from_halo(self,haloID):
        return self.data.ix[self.get_all_subs_recurse(haloID)]

    # updated 3/26 to include support for array/list input of haloID. Also streamlined the code.
    def get_all_sub_particles_from_halo(self,haloID):
        """
        returns int array of particle IDs belonging to all substructure
        within host of haloID
        """
        if type(haloID) == list or type(haloID)==np.ndarray:
            subids = [self.get_all_subhalos_from_halo(id)['id'] for id in haloID]
            idlist = [[item for s in sid for item in self.get_particles_from_halo(s)] for sid in subids]
            return idlist
        else:
            subids = self.get_all_subhalos_from_halo(haloID)['id']
            return [item for s in subids for item in self.get_particles_from_halo(s)]
            # old method
            # idlist = np.array([])
            #for sid in subids:
            #    idlist = np.append(idlist, self.get_particles_from_halo(sid))
            #return idlist.astype(int)

    def get_all_particles_from_halo(self,haloID):
        """
        returns int array of all particles belonging to haloID
        """
        return np.append(self.get_particles_from_halo(haloID), self.get_all_sub_particles_from_halo(haloID)).astype(int)



#### This section of methods has been re-added to enable RSHaloMatch.py to work
# User not encouraged to use these methods for other purposes.

    def get_host_mask(self):
        """Return mask for all host halos. Example:
        mask = cat.get_host_mask()
        host_masses = cat.get_masses()[mask]
        """
        return self.data['hostID'] == -1

    def get_posns(self):
        """Returns np array of the positions with shape (nhalo, 3), sorted
        from high mass to low mass. Units are Mpc"""
        return np.array(self.data[['posX','posY','posZ']])

    def get_masses(self):
        """Return np array masses for all halos (host and sub), sorted from high to low"""
        return np.array(self.data['mvir'])

    def get_npart(self):
        """Returns np array of number of particles in each halo"""
        return np.array(self.data['npart'])
#### End of section for RSHaloMatch.py functionality
    


 #   def get_all_particle_info_from_halo(self,haloID,gadgetfilepath):
 #       import readsnap as rgadget
 #       particleIDs = np.array(self.get_particles_from_halo(haloID))
#
# #       snapIDS = rgadget.read_block(gadgetfilepath, "ID  ",parttype=1)
# #       index = snapIDS.argsort()
# #       print "Index Complete"
# #       del snapIDS
#
# #       newindex = np.take(index,particleIDs)
# #       del index
#
# #       snapPOS = rgadget.read_block(gadgetfilepath, "POS ",parttype=1)
# #       x = np.take(snapPOS[:,0],newindex)
# #       y = np.take(snapPOS[:,1],newindex)
# #       z = np.take(snapPOS[:,2],newindex)
# #       print "Positions Aquired"
#
# #       snapVEL = rgadget.read_block(gadgetfilepath, "VEL ",parttype=1)
# #       vx = np.take(snapVEL[:,0],newindex)
# #       vy = np.take(snapVEL[:,1],newindex)
# #       vz = np.take(snapVEL[:,2],newindex)
# #       print "Velocities Aquired"
#
 #       return x,y,z,vx,vy,vz
       

## Hard defined constants corresponding to data description and its column
"""
id = 0
posX = 1
posY = 2
posZ = 3
corevelx = 4
corevely = 5
corevelz = 6
pecVX = 7
pecVY = 8
pecVZ = 9
mvir = 10
rvir = 11
child_r = 12
mgrav = 13
vmax = 14
rvmax = 15
rs = 16
vrms = 17
Jx = 18
Jy = 19
Jz = 20
Epot = 21 #need to double check this is really energy
spin = 22
npart = 23
num_child_p = 24
p_start = 25
desc = 26
flags = 27
n_core = 28
min_pos_err = 29
min_vel_err = 30
min_bulkvel_err = 31
"""

"""
single value parameters - 
total_particles
num_halos
snap
scale
Om
O1
h0
box_size
particle_mass
particle_type
data

datatypes = np.dtype([('id', '<i8'), ('posX', '<f8'),\
                      ('posY', '<f8'),('posZ', '<f8'),\
                      ('corevelX', '<f8'),('corevelY', '<f8'),\
                      ('corevelZ', '<f8'),('pecVX', '<f8'),\
                      ('pecVY', '<f8'),('pecVZ', '<f8'),\
                      ('bulkvelX', '<f8'),('bulkvelY', '<f8'),\
                      ('bulkvelZ', '<f8'),('mvir', '<f8'),\
                      ('rvir', '<f8'),('child_r', '<f8'),\
                      ('mgrav', '<f8'),('vmax', '<f8'),\
                      ('rvmax', '<f8'),('rs', '<f8'),\
                      ('vrms', '<f8'),('Jx', '<f8'),\
                      ('Jy', '<f8'),('Jz', '<f8'),\
                      ('Epot', '<f8'),('spin', '<f8'),\
                      ('blank1', '<f4'),('npart', '<i8'),\
                      ('num_child_p', '<i8'),('p_start', '<i8'),\
                      ('desc', '<i8'),('flags', '<i8'),\
                      ('n_core', '<i8'),('min_pos_err', '<f8'),\
                      ('min_vel_err', '<f8'),('min_bulkvel_err', '<f8'),\
                      ('blank2', '<f4')])
"""    

