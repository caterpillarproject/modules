"""
Read in all data from Rockstar Halo Finder into a large matrix.
Rows correspond to halos, columns correspond to halo properties.
Much faster access to data than object approach.
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

BinaryHeaderSize = 256 #: in bytes
HeaderInfoSize = 96  #: in bytes
HaloSize = 176  #: in bytes
ParticleSize = 8 #: in bytes
KpcToMpc = 1.0/1000

datatypesstr = "q"+("f" * 26)+ ("q"*6)+"ffff"
numbytes = 30*4+7*8

id = 0
mvir = 13
npart = 27
hostID = 37
offset = 38
particle_offset = 39

num_columns = 40

class HaloData:    
    """
    Create a Halo Catalogue for a particular snapshot
    """
    def __init__(self, dir, snap_num, base='halos_', digits=2, AllParticles=False):
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
                line = f.read(numbytes)
                data[i,0:37] = struct.unpack(datatypesstr, line)
                # info to read particle IDs for halo
                data[i,offset] = particleID_start
                data[i,particle_offset] = particleID_start2
                files[i] = file_name
                particleID_start += ParticleSize*data[i,npart]
                particleID_start2 += data[i,npart]
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
            
        self.files = pandas.DataFrame(files,index=data[:,id].astype(int),columns=['file'])

        #print data.shape, 'shape of data'
        #print data[:,id].shape

        
        
        self.data = pandas.DataFrame(data,index=data[:,id].astype(int),columns=['id','posX','posY','posZ','corevelx','corevely','corevelz','pecVX','pecVY','pecVZ', 'bulkvelX','bulkvelY','bulkvelZ', 'mvir','rvir','child_r','mgrav','vmax','rvmax','rs','vrms','Jx','Jy','Jz','Epot','spin','blank1','npart','num_child_p','p_start','desc','flags','n_core','min_pos_err','min_vel_err','min_bulkvel_err','blank2','hostID','offset','particle_offset'])
        
        self.particles = self.particles.astype(int)
        ## add column of hostID
        file = 'parents'+'.list'
        parents = rp.readParents(self.dir+'/'+base+str(snap_num).zfill(digits),file, self.num_halos)
        self.data['hostID'].ix[parents[:,0]] = parents[:,1] # fill in hostID column
        #print self.data['id'].ix[0], 'ix method'
        #print self.data['id'][0], 'no ix method'
        print time.clock()-start_time, ' time'
        
    def get_particles_from_halo(self, haloID):
        """
        @param haloID: id number of halo. Not its row position in matrix
        @return: a list of particle IDs in the Halo
        """
        if self.AllParticles:
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
        if type(haloIDs) == list or type(haloIDs) == np.ndarray:
            return np.array([ self.data[self.data['hostID']==ID] for ID in haloIDs])
        else:
            return self.data[self.data['hostID']==haloIDs] 
    
    def get_hosts(self):
        return self.data[self.data['hostID']==-1]

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
    def get_all_subhalos_from_halo(self,haloID):
        
        return self.data.ix[self.get_all_subs_recurse(haloID)]

    def get_all_particles_from_halo(self,id):
        idlist = np.array([])
        subids = self.get_all_subhalos_from_halo(id)
        for sid in subids:
            idlist = np.concatenate(idlist, self.get_particles_from_halo(sid))
        return idlist

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

