# code for reading Subfind's subhalo_tab files
# usage e.g.:
#
# import readidsHDF5
# ids = readidsHDF5.subid_file("./m_10002_h_94_501_z3_csf/",63)
# print ids.IDs 

import numpy as np
import os
import sys
import tables
 
class subid_file:
	def __init__(self, basedir, snapnum, long_ids = False, name = "fof_subhalo_tab", verbose = False):
		self.filebase = basedir + "/groups_" + str(snapnum).zfill(3) + "/" + name + "_" + str(snapnum).zfill(3) + "."

		if long_ids: self.id_type = np.uint64
		else: self.id_type = np.uint32

		filenum = 0
		doneflag = False
		skip_ids = 0

		while not doneflag:
			curfile = self.filebase + str(filenum) + ".hdf5"

			if (not os.path.exists(curfile)):
				self.filebase = basedir + "/" + name + "_" + str(snapnum).zfill(3)
				curfile = self.filebase + ".hdf5"

			if (not os.path.exists(curfile)):
				print "file not found:", curfile
				sys.exit()

			f=tables.openFile(curfile)
			nfiles = f.root.Header._v_attrs.NumFiles
			nids = f.root.Header._v_attrs.Nids_ThisFile

			idlen = f.root.Header._v_attrs.Nids_Total

			if filenum == 0:
				if (long_ids == False):
					self.IDs = np.empty(idlen, dtype=np.uint32)
				else:
					self.IDs = np.empty(idlen, dtype=np.uint64)
			

			self.IDs[skip_ids:skip_ids+nids] = f.root.IDs.ID

			skip_ids+=nids

			f.close()  
			filenum += 1
			if filenum == self.nfiles: doneflag = True
       
 
