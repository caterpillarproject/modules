# code for reading Subfind's fof_subhalo_tab files
# usage e.g.:
#
# import readsubfHDF5
# cat = readsubfHDF5.subfind_catalog("./output/",60)
# print cat.nsubs
# print "largest halo x position = ",cat.sub_pos[0][0] 

import numpy as np
import os
import sys
import tables
 
class subfind_catalog:
  def __init__(self, basedir, snapnum, long_ids = False):
    self.filebase = basedir + "/groups_" + str(snapnum).zfill(3) + "/fof_subhalo_tab_" + str(snapnum).zfill(3) + "."
 
    #print
    #print "reading subfind catalog for snapshot",snapnum,"of",basedir

    if long_ids: self.id_type = np.uint64
    else: self.id_type = np.uint32
 
 
    filenum = 0
    doneflag = False
    skip_gr = 0
    skip_sub = 0
    while not doneflag:
      curfile = self.filebase + str(filenum) + ".hdf5"
      
      if (not os.path.exists(curfile)):
        print "file not found:", curfile
        sys.exit()
      
      f=tables.openFile(curfile)      
              
      ngroups = f.root.Header._v_attrs.Ngroups_ThisFile 
      totngroups = f.root.Header._v_attrs.Ngroups_Total
      nids = f.root.Header._v_attrs.Nids_ThisFile
      totnids = f.root.Header._v_attrs.Nids_Total
      nsubs = f.root.Header._v_attrs.Nsubgroups_ThisFile
      totnsubs = f.root.Header._v_attrs.Nsubgroups_Total
      nfiles = f.root.Header._v_attrs.NumFiles 
  
      if filenum == 0:
        self.ngroups = totngroups
        self.nids = totnids
        self.nsubs = totnsubs

        self.group_len = np.empty(totngroups, dtype=np.uint32)
	#self.group_len_type = np.empty(totngroups, dtype=np.uint32)
        self.group_offset = np.empty(totngroups, dtype=np.uint32)
        self.group_mass = np.empty(totngroups, dtype=np.float32)
        self.group_pos = np.empty(totngroups, dtype=np.dtype((np.float32,3)))
	self.group_vel = np.empty(totngroups, dtype=np.dtype((np.float32,3)))
        self.group_m_mean200 = np.empty(totngroups, dtype=np.float32)
        self.group_r_mean200 = np.empty(totngroups, dtype=np.float32)
        self.group_m_crit200 = np.empty(totngroups, dtype=np.float32)
        self.group_r_crit200 = np.empty(totngroups, dtype=np.float32)
        self.group_m_tophat200 = np.empty(totngroups, dtype=np.float32)
        self.group_r_tophat200 = np.empty(totngroups, dtype=np.float32)
        self.group_contamination_count = np.empty(totngroups, dtype=np.uint32)
        self.group_contamination_mass = np.empty(totngroups, dtype=np.float32)
        self.group_nsubs = np.empty(totngroups, dtype=np.uint32)
        self.group_firstsub = np.empty(totngroups, dtype=np.uint32)
        
        self.sub_len = np.empty(totnsubs, dtype=np.uint32)
        self.sub_offset = np.empty(totnsubs, dtype=np.uint32)
        self.sub_parent = np.empty(totnsubs, dtype=np.uint32)
        self.sub_mass = np.empty(totnsubs, dtype=np.float32)
        self.sub_pos = np.empty(totnsubs, dtype=np.dtype((np.float32,3)))
        self.sub_vel = np.empty(totnsubs, dtype=np.dtype((np.float32,3)))
        self.sub_cm = np.empty(totnsubs, dtype=np.dtype((np.float32,3)))
        self.sub_spin = np.empty(totnsubs, dtype=np.dtype((np.float32,3)))
        self.sub_veldisp = np.empty(totnsubs, dtype=np.float32)
        self.sub_vmax = np.empty(totnsubs, dtype=np.float32)
        self.sub_vmaxrad = np.empty(totnsubs, dtype=np.float32)
        self.sub_halfmassrad = np.empty(totnsubs, dtype=np.float32)
        self.sub_id_mostbound = np.empty(totnsubs, dtype=self.id_type)
        self.sub_grnr = np.empty(totnsubs, dtype=np.uint32)
     
      if ngroups > 0:
        locs = slice(skip_gr, skip_gr + ngroups)
        self.group_len[locs] = f.root.Group.GroupLen
        #self.group_len_type[locs] = f.root.Group.GroupLenType
        self.group_mass[locs] = f.root.Group.GroupMass
	#self.group_mass_type[locs] = f.root.Group.GroupMassType
    	self.group_pos[locs] = f.root.Group.GroupPos
        self.group_vel[locs] = f.root.Group.GroupVel

        self.group_m_mean200[locs] = f.root.Group.Group_M_Mean200 
        self.group_r_mean200[locs] = f.root.Group.Group_R_Mean200
        self.group_m_crit200[locs] = f.root.Group.Group_M_Crit200
        self.group_r_crit200[locs] = f.root.Group.Group_R_Crit200
        self.group_m_tophat200[locs] = f.root.Group.Group_M_TopHat200
        self.group_r_tophat200[locs] = f.root.Group.Group_R_TopHat200
	self.group_nsubs[locs] = f.root.Group.GroupNsubs
	self.group_firstsub[locs] = f.root.Group.GroupFirstSub
#    GroupIDs(skipids:skipids + nlocids -1) = read_dataset(file_id, 'IDs/ID')
        skip_gr += ngroups

        
      if nsubs > 0:
        locs = slice(skip_sub, skip_sub + nsubs)
        self.sub_len[locs] = f.root.Subhalo.SubhaloLen
	#self.sub_len_type[locs] = f.root.Subhalo.SubhaloLenType
        self.sub_mass[locs] = f.root.Subhalo.SubhaloMass
	#self.sub_mass_type[locs] = f.root.Subhalo.SubhaloMassType 
        self.sub_pos[locs] = f.root.Subhalo.SubhaloPos 
        self.sub_vel[locs] = f.root.Subhalo.SubhaloVel
        self.sub_cm[locs] = f.root.Subhalo.SubhaloCM
        self.sub_spin[locs] = f.root.Subhalo.SubhaloSpin
        self.sub_veldisp[locs] = f.root.Subhalo.SubhaloVelDisp
        self.sub_vmax[locs] = f.root.Subhalo.SubhaloVmax
        self.sub_vmaxrad[locs] = f.root.Subhalo.SubhaloVmaxRad
        self.sub_halfmassrad[locs] = f.root.Subhalo.SubhaloHalfmassRad
        self.sub_id_mostbound[locs] = f.root.Subhalo.SubhaloIDMostbound
        self.sub_grnr[locs] = f.root.Subhalo.SubhaloGrNr
	self.sub_parent[locs] = f.root.Subhalo.SubhaloParent

        skip_sub += nsubs
      
      f.close()

      filenum += 1
      if filenum == nfiles: doneflag = True
       
       
    #print
    #print "number of groups =", self.ngroups
    #print "number of subgroups =", self.nsubs
    #if self.nsubs > 0:
    #  print "largest group of length",self.group_len[0],"has",self.group_nsubs[0],"subhalos"
    #  print



