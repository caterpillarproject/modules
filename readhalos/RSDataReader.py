import numpy as np
import readParentsList as rp
import pandas
import struct
import os
import sys

class RSDataReader:
    """
    Alex's 10/4/13 rewrite of RSDataReader, combining v2 and v3 and cleaning everything up
    """
    def __init__(self, dir, snap_num, version=2, sort_by='mvir', base='halos_', digits=2, AllParticles=False):
        self.dir = dir
        self.snap_num = snap_num
        self.AllParticles = AllParticles
        self.version=version
        self.sort_by = sort_by

        self.particlebytes = 8

        def getfilename(file_num):
            return dir+'/'+base+str(snap_num).zfill(digits)+'/'+base+str(snap_num).zfill(digits)+'.'+str(file_num)+'.bin'

        numheaderbytes=256
        if version==1:
            print "VERSION 1 (Greg's cosmboxes version) NOT IMPLEMENTED!"
            datatypesstr = ""
            numbytes = 256
            raise Exception("Not implemented error")
        if version==2:
            headerfmt = "qqqffffffffffqqffq"+"x"*(256-96)
            varlist = np.dtype([('id','<i8'),\
                                ('posX','<f8'),('posY','<f8'),('posZ','<f8'),\
                                ('pecVX','<f8'),('pecVY','<f8'),('pecVZ','<f8'),\
                                ('corevelx','<f8'),('corevely','<f8'),('corevelz','<f8'),\
                                ('bulkvelx','<f8'),('bulkvely','<f8'),('bulkvelz','<f8'),\
                                ('mvir','<f8'),('rvir','<f8'),('child_r','<f8'),('vmax_r','<f8'),\
                                ('mgrav','<f8'),('vmax','<f8'),('rvmax','<f8'),('rs','<f8'),('rs_klypin','<f8'),\
                                ('vrms','<f8'),('Jx','<f8'),('Jy','<f8'),('Jz','<f8'),\
                                ('Epot','<f8'),('spin','<f8'),('altm1','<f8'),('altm2','<f8'),('altm3','<f8'),('altm4','<f8'),\
                                ('Xoff','<f8'),('Voff','<f8'),\
                                ('b_to_a','<f8'),('c_to_a','<f8'),('A[x]','<f8'),('A[y]','<f8'),('A[z]','<f8'),\
                                ('spin_bullock','<f8'),('T/|U|','<f8'),\
                                ('npart','<i8'),('num_cp','<i8'),('numstart','<i8'),\
                                ('desc','<i8'),('flags','<i8'),('n_core','<i8'),\
                                ('min_pos_err','<f8'),('min_vel_err','<f8'),('min_bulkvel_err','<f8'),\
                                ('hostID','<i8'),('offset','<i8'),('particle_offset','<i8')])
            datatypesstr = "qffffffffffffffffffffffffffffffffffffffffqqqqqqxxxxfff"
            numbytes = struct.calcsize(datatypesstr) #232
        if version==3:
            headerfmt = "qqqffffffffffqqffq"+"x"*(256-96)
            varlist = np.dtype([('id','<i8'),\
                                ('posX','<f8'),('posY','<f8'),('posZ','<f8'),\
                                ('pecVX','<f8'),('pecVY','<f8'),('pecVZ','<f8'),\
                                ('corevelx','<f8'),('corevely','<f8'),('corevelz','<f8'),\
                                ('bulkvelx','<f8'),('bulkvely','<f8'),('bulkvelz','<f8'),\
                                ('mvir','<f8'),('rvir','<f8'),('child_r','<f8'),('vmax_r','<f8'),\
                                ('mgrav','<f8'),('vmax','<f8'),('rvmax','<f8'),('rs','<f8'),('rs_klypin','<f8'),\
                                ('vrms','<f8'),('Jx','<f8'),('Jy','<f8'),('Jz','<f8'),\
                                ('Epot','<f8'),('spin','<f8'),('altm1','<f8'),('altm2','<f8'),('altm3','<f8'),('altm4','<f8'),\
                                ('Xoff','<f8'),('Voff','<f8'),\
                                ('b_to_a','<f8'),('c_to_a','<f8'),('A[x]','<f8'),('A[y]','<f8'),('A[z]','<f8'),\
                                ('b_to_a2','<f8'),('c_to_a2','<f8'),('A2[x]','<f8'),('A2[y]','<f8'),('A2[z]','<f8'),\
                                ('spin_bullock','<f8'),('T/|U|','<f8'),\
                                ('npart','<i8'),('num_cp','<i8'),('numstart','<i8'),\
                                ('desc','<i8'),('flags','<i8'),('n_core','<i8'),\
                                ('min_pos_err','<f8'),('min_vel_err','<f8'),('min_bulkvel_err','<f8'),\
                                ('hostID','<i8'),('offset','<i8'),('particle_offset','<i8')])
            datatypesstr = "qfffffffffffffffffffffffffffffffffffffffffffffqqqqqqxxxxfff"
            numbytes = struct.calcsize(datatypesstr) #256
        if version==4:
            headerfmt = "qqqffffffffffqqffq"+"x"*(256-96)
            varlist = np.dtype([('id','<i8'),\
                                ('posX','<f8'),('posY','<f8'),('posZ','<f8'),\
                                ('pecVX','<f8'),('pecVY','<f8'),('pecVZ','<f8'),\
                                ('corevelx','<f8'),('corevely','<f8'),('corevelz','<f8'),\
                                ('bulkvelx','<f8'),('bulkvely','<f8'),('bulkvelz','<f8'),\
                                ('mvir','<f8'),('rvir','<f8'),('child_r','<f8'),('vmax_r','<f8'),\
                                ('mgrav','<f8'),('vmax','<f8'),('rvmax','<f8'),('rs','<f8'),('rs_klypin','<f8'),\
                                ('vrms','<f8'),('Jx','<f8'),('Jy','<f8'),('Jz','<f8'),\
                                ('Epot','<f8'),('spin','<f8'),('altm1','<f8'),('altm2','<f8'),('altm3','<f8'),('altm4','<f8'),\
                                ('Xoff','<f8'),('Voff','<f8'),\
                                ('b_to_a','<f8'),('c_to_a','<f8'),('A[x]','<f8'),('A[y]','<f8'),('A[z]','<f8'),\
                                ('b_to_a2','<f8'),('c_to_a2','<f8'),('A2[x]','<f8'),('A2[y]','<f8'),('A2[z]','<f8'),\
                                ('spin_bullock','<f8'),('T/|U|','<f8'),\
                                ('npart','<i8'),('num_cp','<i8'),('numstart','<i8'),\
                                ('desc','<i8'),('flags','<i8'),('n_core','<i8'),\
                                ('min_pos_err','<f8'),('min_vel_err','<f8'),('min_bulkvel_err','<f8'),\
                                ('num_bound','<i8'),\
                                ('hostID','<i8'),('offset','<i8'),('particle_offset','<i8')])
            datatypesstr = "qfffffffffffffffffffffffffffffffffffffffffffffqqqqqqxxxxfffq"
            numbytes = struct.calcsize(datatypesstr) #256

        if version==5: #modification to include total num bound particles and tidal radius
            # corresponds to rockstar version here: /spacebase/data/gdooley/RockstarSorted/rockstarTidal 
            # with TIDAL defined.
            headerfmt = "qqqffffffffffqqffq"+"x"*(256-96)
            varlist = np.dtype([('id','<i8'),\
                                ('posX','<f8'),('posY','<f8'),('posZ','<f8'),\
                                ('pecVX','<f8'),('pecVY','<f8'),('pecVZ','<f8'),\
                                ('corevelx','<f8'),('corevely','<f8'),('corevelz','<f8'),\
                                ('bulkvelx','<f8'),('bulkvely','<f8'),('bulkvelz','<f8'),\
                                ('mvir','<f8'),('rvir','<f8'),('child_r','<f8'),('vmax_r','<f8'),\
                                ('mgrav','<f8'),('vmax','<f8'),('rvmax','<f8'),('rs','<f8'),('rs_klypin','<f8'),\
                                ('vrms','<f8'),('Jx','<f8'),('Jy','<f8'),('Jz','<f8'),\
                                ('Epot','<f8'),('spin','<f8'),('altm1','<f8'),('altm2','<f8'),('altm3','<f8'),('altm4','<f8'),\
                                ('Xoff','<f8'),('Voff','<f8'),\
                                ('b_to_a','<f8'),('c_to_a','<f8'),('A[x]','<f8'),('A[y]','<f8'),('A[z]','<f8'),\
                                ('b_to_a2','<f8'),('c_to_a2','<f8'),('A2[x]','<f8'),('A2[y]','<f8'),('A2[z]','<f8'),\
                                ('spin_bullock','<f8'),('T/|U|','<f8'),\
                                ('npart','<i8'),('num_cp','<i8'),('numstart','<i8'),\
                                ('desc','<i8'),('flags','<i8'),('n_core','<i8'),\
                                ('min_pos_err','<f8'),('min_vel_err','<f8'),('min_bulkvel_err','<f8'),\
                                ('num_bound','<i8'),('tidal_r','<f8'),\
                                ('hostID','<i8'),('offset','<i8'),('particle_offset','<i8')])
            datatypesstr = "qfffffffffffffffffffffffffffffffffffffffffffffqqqqqqxxxxfffqfxxxx"
            numbytes = struct.calcsize(datatypesstr) #256

        self.datatypesstr = datatypesstr

        file_num = 0
        file_name = getfilename(file_num)
        if (not os.path.exists(file_name)):
            print "ERROR: file not found", file_name
            sys.exit()
            
        ## Count total number of particles/halos in all data blocks
        self.num_halos = 0
        self.total_particles = 0
        while os.path.exists(file_name):
            f = open(file_name)
            h = f.read(numheaderbytes)
            (magic,self.snap_num,chunk,self.scale,self.Om,self.Ol,self.h0,\
             bounds1,bounds2,bounds3,bounds4,bounds5,bounsd6,\
             num_halos,num_particles,\
             self.boxsize,self.particle_mass,self.particle_type) = struct.unpack(headerfmt,h)
            self.num_halos += num_halos
            self.total_particles += num_particles
            f.close()
            file_num += 1
            file_name = getfilename(file_num)
        ## Initialize empty data structure
        data = np.zeros(self.num_halos,dtype=varlist)
        files = np.array(['']*self.num_halos, dtype='|S'+str(len(file_name)*2))
        self.particles = np.array([])

        ## Now, read in the actual data
        file_num = 0 # reset file name
        file_name = getfilename(file_num)
        i = 0
        particleID_start2 = 0 
        while os.path.exists(file_name):
            f = open(file_name)
            h = f.read(numheaderbytes)
            num_halos,num_particles = struct.unpack(("x"*64)+"qq"+("x"*16)+"x"*(256-96),h)            
            # note this is current block's num_halos, not the total self.num_halos
            particleID_start = struct.calcsize(headerfmt)+num_halos*struct.calcsize(datatypesstr)

            for j in xrange(num_halos):
                line = f.read(numbytes)
                data[i] = struct.unpack(datatypesstr,line)+(0,0,0)
                data[i][-2] = particleID_start  # offset
                data[i][-1] = particleID_start2 # particle_offset
                files[i] = file_name
                particleID_start  += self.particlebytes * data['npart'][i]
                particleID_start2 += data['npart'][i]
                i += 1
            if AllParticles:
                line = f.read() # read the rest of the file
                ## DEBUG: this ratio should be 1
                if num_particles != 0:
                    assert(len(line)/float(self.particlebytes))/num_particles == 1
                    self.particles = np.concatenate((self.particles, np.array(struct.unpack("q"*num_particles,line))))
            f.close()
            file_num += 1
            file_name = getfilename(file_num)

        if sort_by != None:
            sortedIndices = data[sort_by].argsort()[::-1]
            data = data[sortedIndices]
            files= files[sortedIndices]

        self.files = pandas.DataFrame(files, index=data['id'].astype(int),columns=['file'])
        self.data = pandas.DataFrame(data,index=data['id'])

        self.particles = self.particles.astype(int)
        parents = rp.readParents(dir+'/'+base+str(snap_num).zfill(digits),'parents.list',self.num_halos)
        self.data['hostID'].ix[parents[:,0]] = parents[:,1]

        self.ix = self.data.ix

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

    def get_subhalos_from_halo(self,haloID):
        """
        # Retrieve subhalos only one level deep.
        # Does not get sub-sub halos, etc.
        """
        return self.data[self.data['hostID']==haloID]

    def get_subhalos_from_halos(self,haloIDs):
        """
        returns an array of pandas data frames of subhalos. one data frame
        for each host halo. returns only first level of subhalos.
        # for multiple halos, returns the subhalos of each halo as an array of arrays.
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

    def get_all_subs_recurse(self,haloID):
        """
        # Retrieve all subhalos: sub and sub-sub, etc. 
        # just need mask of all subhalos, then return data frame subset
        """
        subs = self.get_subhalo_ids_from_halos(haloID)
        subs = [item for sublist in subs for item in sublist]
        if subs == []:
            return [] 
        else:
            return np.append(subs,self.get_all_subs_recurse(subs))
                                       
    def get_all_subhalos_from_halo(self,haloID):
        """
        # Retrieve all subhalos: sub and sub-sub, etc.
        # return pandas data fram of subhalos
        """
        return self.data.ix[self.get_all_subs_recurse(haloID)]

    def get_all_sub_particles_from_halo(self,haloID):
        """
        returns int array of particle IDs belonging to all substructure
        within host of haloID
        # updated 3/26 to include support for array/list input of haloID. Also streamlined the code.
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

    def getversion(self):
        if self.version == 2:
            return "Version 2: Rockstar 0.99.9 RC1"
        if self.version == 3:
            return "Version 3: Rockstar 0.99.9 RC2"
        return "ERROR: Not a valid version number!"

    def __getitem__(self,key):
        return self.data[key]
    def __len__(self):
        return len(self.data)
    def __repr__(self):
        out = "RSDataReader "+self.getversion()+"\n"
        out +="Read from "+self.dir+"\n"
        out +="Snap number "+str(self.snap_num)+"\n"
        out +="Number of halos: "+str(self.num_halos)+"\n"
        return out + "Sorted by "+self.sort_by
    
