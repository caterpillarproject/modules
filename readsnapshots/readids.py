# code for reading Subfind's subhalo_tab files
# usage e.g.:
#
# import readsubf
# cat = readsubf.subfind_catalog("./m_10002_h_94_501_z3_csf/",63,masstab=True)
# print cat.nsubs
# print "largest halo x position = ",cat.sub_pos[0][0] 

import numpy as np
import os
import sys
 
class subid_file:
  def __init__(self, basedir, snapnum, substart, sublen, long_ids = False, swap = False, verbose = False):
    self.filebase = basedir + "/groups_" + str(snapnum).zfill(3) + "/subhalo_ids_" + str(snapnum).zfill(3) + "."

    if (verbose):	 
	    print
	    print "reading subhalo IDs for snapshot",snapnum,"of",basedir
 
 
    filenum = 0
    doneflag = False
    sublen
    count=substart
    found=0


    while not doneflag:
      curfile = self.filebase + str(filenum)
      
      if (not os.path.exists(curfile)):
        print "file not found:", curfile
        sys.exit()
      
      f = open(curfile,'rb')
      Ngroups = np.fromfile(f, dtype=np.uint32, count=1)[0]
      TotNgroups = np.fromfile(f, dtype=np.uint32, count=1)[0]
      NIds = np.fromfile(f, dtype=np.uint32, count=1)[0]
      TotNids = np.fromfile(f, dtype=np.uint64, count=1)[0]
      NTask = np.fromfile(f, dtype=np.uint32, count=1)[0]      
      Offset = np.fromfile(f, dtype=np.int32, count=1)[0]             #changed from uint32 to int32 to avoid negative overflows
      if swap:     
	      Ngroups = Ngroups.byteswap()
	      TotNgroups = TotNgroups.byteswap()
	      NIds = NIds.byteswap()
	      TotNids = TotNids.byteswap()
	      NTask = NTask.byteswap()
	      Offset = Offset.byteswap()
      if filenum == 0:
	if (verbose):
	        print "Ngroups    = ", Ngroups
	        print "TotNgroups = ", Ngroups	
	        print "NIds       = ", NIds
        	print "TotNids    = ", TotNids
	        print "NTask      = ", NTask	
	        print "Offset     = ", Offset	
	self.nfiles = NTask
	self.SubLen=sublen
	if (long_ids == False):
	        self.SubIDs = np.empty(sublen, dtype=np.uint32)
	else:
                self.SubIDs = np.empty(sublen, dtype=np.uint64)
	
      if count <= Offset+NIds:    
	nskip = count - Offset
	nrem = Offset + NIds - count
	if sublen > nrem:
		n_to_read = nrem
	else:	
		n_to_read = sublen		
	if n_to_read > 0:
		if (verbose):
			print filenum, n_to_read
		if nskip > 0:
			dummy=np.fromfile(f, dtype=np.uint32, count=nskip)
			if (verbose):
				print dummy
	        locs = slice(found, found + n_to_read)
		if (long_ids == False):
		        dummy2 = np.fromfile(f, dtype=np.uint32, count=n_to_read)
		else:
                        dummy2 = np.fromfile(f, dtype=np.uint64, count=n_to_read)
		if (verbose):
			print dummy2
		self.SubIDs[locs]=dummy2
		found += n_to_read
	count += n_to_read
	sublen -= n_to_read
      f.close()  
      filenum += 1
      if filenum == self.nfiles: doneflag = True
       
    if swap:
      self.SubIDs.byteswap(True)
 
