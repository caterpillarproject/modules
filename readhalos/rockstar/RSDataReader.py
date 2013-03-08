"""
Read in all data from Rockstar Halo Finder into a large matrix.
Rows correspond to halos, columns correspond to halo properties.
Much faster access to data than object approach.
"""

import numpy as np
import os
import sys
#import copy
#import operator
import readParentsList as rp
import pandas as pd

BinaryHeaderSize = 256 #: in bytes
HeaderInfoSize = 96  #: in bytes
HaloSize = 176  #: in bytes
ParticleSize = 8 #: in bytes
KpcToMpc = 1.0/1000


 ## Hard defined constants corresponding to data description and its column
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
hostID = 32
offset = 33
particle_offset = 34

num_columns = 35


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
"""

class RSDataReader:
    """
	Creates RS Data Reader that reads in halo catalogue. It automatically reads in particles
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
    """
    def __init__(self, dir, snap_num, base='halos_', digits=3, AllParticles=True):
        """
		@param dir: base directory containing binary files from Rockstar Halo Finder. ex: \"/home/gdooley/Rockstar-0.99/Output2\"
		@param snap_num: snap number to be viewed. ex: 50.
		@param digits: number of digits in file name of snapshot. ex: halos_063.0 should have digits = 3.
		"""
        self.hftype = 'rockstar'
        self.dir = dir
        self.snap_num = snap_num
        self.AllParticles = AllParticles
  
        # open first file, test if it exists
        file_num = 0
        file_name = dir + '/' + base + str(snap_num).zfill(digits) + "." + str(file_num) + ".bin"
        if (file_num == 0 and not os.path.exists(file_name)):
            raise IOError('file not found: '+file_name)

        ## establish total number of halos in order to create matrix of the correct size
        nhalo = 0
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
            nhalo += num_halos
            num_particles = np.fromfile(f, np.int64, count = 1)[0]
            self.total_particles+=num_particles
            #increment to next file
            f.close()
            file_num += 1
            file_name = dir+'/'+base+str(snap_num).zfill(digits)+"."+str(file_num)+".bin"
        self.nhalo = nhalo
        #reset file name
        file_num = 0
        file_name = dir+'/'+base+str(snap_num).zfill(digits)+"."+str(file_num)+".bin"
        ### Two important pieces of data storage
        data = np.zeros((nhalo,num_columns)) #: Matrix of all halos and all information
        string_len = len(file_name)*2
        datatype = '|S'+str(string_len)
        files = np.array(['']*nhalo,dtype=datatype)
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
                data[i,id] = np.fromfile(f,np.int64, 1)[0] 
                data[i,posX] = np.fromfile(f,np.float32,count = 1)[0]
                data[i,posY] = np.fromfile(f,np.float32,count = 1)[0]
                data[i,posZ] = np.fromfile(f,np.float32,count = 1)[0]
                data[i,corevelx] = np.fromfile(f,np.float32,count = 1)[0]
                data[i,corevely] = np.fromfile(f,np.float32,count = 1)[0]
                data[i,corevelz] = np.fromfile(f,np.float32,count = 1)[0]
                data[i,pecVX] = np.fromfile(f,np.float32,count = 1)[0]
                data[i,pecVY] = np.fromfile(f,np.float32,count = 1)[0]
                data[i,pecVZ] = np.fromfile(f,np.float32,count = 1)[0]
                bulkvel = np.fromfile(f,np.float32,count = 3) #bulk velocity listed twice in binary output. NO IDEA WHY.
                data[i,mvir] = np.fromfile(f,np.float32,count = 1)[0]
                data[i,rvir] = np.fromfile(f,np.float32,count = 1)[0]
                data[i,child_r] = np.fromfile(f,np.float32,count = 1)[0]
                data[i,mgrav] = np.fromfile(f,np.float32,count = 1)[0]
                data[i,vmax] = np.fromfile(f,np.float32,count = 1)[0]
                data[i,rvmax] = np.fromfile(f,np.float32,count = 1)[0]
                data[i,rs] = np.fromfile(f,np.float32,count = 1)[0]
                data[i,vrms] = np.fromfile(f,np.float32,count = 1)[0]
                data[i,Jx] = np.fromfile(f,np.float32,count = 1)[0]
                data[i,Jy] = np.fromfile(f,np.float32,count = 1)[0]
                data[i,Jz] = np.fromfile(f,np.float32,count = 1)[0]
                data[i,Epot] = np.fromfile(f,np.float32,count = 1)[0]
                data[i,spin] = np.fromfile(f,np.float32,count = 1)[0]
                # weird 4 bytes of blank info
                np.fromfile(f,'f',count = 1) #don't know why
                data[i,npart] = np.fromfile(f,np.int64,count = 1)[0]
                data[i,num_child_p] = np.fromfile(f,np.int64,count = 1)[0]
                data[i,p_start] = np.fromfile(f,np.int64,count = 1)[0]
                data[i,desc] = np.fromfile(f,np.int64,count = 1)[0]
                data[i,flags] = np.fromfile(f,np.int64,count = 1)[0]
                data[i,n_core] = np.fromfile(f,np.int64,count = 1)[0]
                data[i,min_pos_err] = np.fromfile(f,np.float32,count = 1)[0]
                data[i,min_vel_err] = np.fromfile(f,np.float32,count = 1)[0]
                data[i,min_bulkvel_err] = np.fromfile(f,np.float32,count = 1)[0]
                # another 4 byte blank spot
                np.fromfile(f,'f',count = 1)              

                # info to read particle IDs for halo
                data[i,offset] = particleID_start
                data[i,particle_offset] = particleID_start2
                files[i] = file_name
                particleID_start += ParticleSize*data[i,npart]
                particleID_start2 += data[i,npart]
                i += 1
            # Read in and store particle information (optional)
            #if self.AllParticles:
            self.particles = np.concatenate((self.particles,np.fromfile(f,np.int64)))
            #increment to next file
            f.close()
            file_num += 1
            file_name = dir+"/halos_"+str(snap_num).zfill(digits)+"."+str(file_num)+".bin"
        self.files = pd.DataFrame(files,index=data[:,id].astype(int),columns=['file'])

        #print data.shape, 'shape of data'
        #print data[:,id].shape
        
        self.data = pd.DataFrame(data,index=data[:,id].astype(int),columns=['id','posX','posY','posZ','corevelx','corevely','corevelz','pecVX','pecVY','pecVZ','mvir','rvir','child_r','mgrav','vmax','rvmax','rs','vrms','Jx','Jy','Jz','Epot','spin','npart','num_child_p','p_start','desc','flags','n_core','min_pos_err','min_vel_err','min_bulkvel_err','hostID','offset','particle_offset'])
        
        self.particles = self.particles.astype(int)
        ## add column of hostID
        #file = 'parents_'+str(self.snap_num)+'.list'
        file = 'parents.list'
        reload(rp)
        parents = rp.readParents(self.dir,file, self.nhalo)
        self.data['hostID'].ix[parents[:,0]] = parents[:,1] # fill in hostID column
        
    def get_particles_from_halo(self, haloID):
        """
        @param haloID: id number of halo. Not its row position in matrix
        @return: a list of particle IDs in the Halo
        """
##         if self.AllParticles:
        start = self.data['particle_offset'].ix[haloID]
        end = start+self.data['npart'].ix[haloID]
        return self.particles[start:end]
##         else:
##             f = open(self.files['file'].ix[haloID])
##             np.fromfile(f,'c',count=int(self.data['offset'].ix[haloID]))
##             particleIDs = np.fromfile(f,np.int64,count=int(self.data['npart'].ix[haloID]))
##             f.close()
##             return particleIDs

    def get_subhalos_from_halo(self,haloID):
        return self.data[self.data['hostID']==haloID]

    def get_hosts(self):
        return self.data[self.data['hostID']==-1]
            
    def get_subhalos(self):
        return self.data[self.data['hostID']!=-1]

                                       



##     def pos(self, HaloIndex):
##         """
##         @param HaloIndex: index of halo in data matrix
##         @return: a 1x3 array of the x,y,z position in Mpc.
##         """
##         return self.data[HaloIndex, posx:posz+1]

##     def corevel(self, HaloIndex):
##         """
##         @param HaloIndex: index of halo in data matrix
##         @return: a 1x3 array of the core velocity in km/s.
##         """
##         return self.data[HaloIndex, corevelx:corevelz+1]

##     def bulkvel(self, HaloIndex):
##         """
##         @param HaloIndex: index of halo in data matrix
##         @return: a 1x3 array of the bulk velocity in km/s.
##         """
##         return self.data[HaloIndex, bulkvelx:bulkvelz+1]

##     def J(self, HaloIndex):
##         """
##         @param HaloIndex: index of halo in data matrix
##         @return: a 1x3 array of the angular momentum.
##         """
##         return self.data[HaloIndex, Jx:Jz+1]
    
##     def sortMass(self):
##         """
##         change matrix to one sorted by mass. Highest mass first.
##         """
##         sortedIndices = self.data[:,mass].argsort()[::-1]
##         self.data = self.data[sortedIndices]
##         self.files = self.files[sortedIndices]
##         return self.data

##     def sortID(self):
##         """
##         change matrix to one sorted by halo ID. Lowest ID firsts.
##         """
##         sortedIndices = self.data[:,id].argsort()
##         self.data = self.data[sortedIndices]
##         self.files = self.files[sortedIndices]
##         return self.data

##     def getHostIndices(self):
##         file = 'parents_'+str(self.snap_num)+'.list'
##         parents = rp.readParents(self.dir,file, self.nhalo)
##         return np.where(parents[:,1] == -1)[0]
##         #return self.data[hostIndices]

##     def getSubs(self):
##         parents = rp.readParents(self.dir,'parents_'+str(self.snap_num)+'.list', self.nhalo)
##         subIndices = np.where(parents[:,1] != -1)[0]
##         return self.data[subIndices]
