# usage e.g.:
#hsml=readhsml.hsml_file(base, 0)
#range can be given with start and len


import numpy as np
import os
import sys
 
class hsml_file:
  def __init__(self, basedir, snapnum, swap = False, start=-1, len=-1):
    self.filebase = basedir + "/hsmldir_" + str(snapnum).zfill(3) + "/hsml_" + str(snapnum).zfill(3) + "."
 
 
 
    filenum = 0
    doneflag = False
    skip = 0
    ngot = 0
    get_start=start
    get_end=start+len
    get_len=len

    while not doneflag:
      curfile = self.filebase + str(filenum)
      
      if (not os.path.exists(curfile)):
        print "file not found:", curfile
        sys.exit()
      
      f = open(curfile,'rb')
              
      Nhsml = np.fromfile(f, dtype=np.uint32, count=1)[0]
      Nprevious = np.fromfile(f, dtype=np.uint32, count=1)[0]
      Ntotal = np.fromfile(f, dtype=np.uint64, count=1)[0]
      Ntask = np.fromfile(f, dtype=np.uint32, count=1)[0]
      
      if swap:
        Nhsml = Nhsml.byteswap()
        Nprevious = Nprevious.byteswap()
        Ntotal = Ntotal.byteswap()
        Ntask = Ntask.byteswap()

      if filenum == 0:
        self.Nhsml = Ntotal
	self.nfiles = Ntask
	
	if (len>0):
                self.Hsml = np.empty(len, dtype=np.float32)
                self.Density = np.empty(len, dtype=np.float32)
                self.Veldisp = np.empty(len, dtype=np.float32)
		
	else:
	        self.Hsml = np.empty(Ntotal, dtype=np.float32)
        	self.Density = np.empty(Ntotal, dtype=np.float32)
	        self.Veldisp = np.empty(Ntotal, dtype=np.float32)
     
      if Nhsml > 0:
	if (len>0):
                hsml_start=skip
                hsml_end=skip+Nhsml
		#print "FILE :", filenum
		#print "HSML :", hsml_start, hsml_end
		#print "GET  :", get_start, get_end
		if (hsml_end >= get_start) & (hsml_start <= get_end):
			rel_start=max([get_start - hsml_start, 0])
			nget=min([Nhsml-rel_start, get_len])
			locs1=slice(ngot, ngot+nget)
			locs2=slice(rel_start, rel_start+nget)
			#print locs1, Nhsml
			#print locs2, Nhsml
			#print np.fromfile(f, dtype=np.float32, count=Nhsml).shape
	                self.Hsml[locs1] = np.fromfile(f, dtype=np.float32, count=Nhsml)[locs2]
       		        self.Density[locs1] = np.fromfile(f, dtype=np.float32, count=Nhsml)[locs2]
               		self.Veldisp[locs1] = np.fromfile(f, dtype=np.float32, count=Nhsml)[locs2]
			ngot+=nget
			get_start+=nget
			get_len-=nget
			#print "NGOT :", ngot
                skip += Nhsml
			
	else:
	        locs = slice(skip, skip + Nhsml)
	        self.Hsml[locs] = np.fromfile(f, dtype=np.float32, count=Nhsml)
        	self.Density[locs] = np.fromfile(f, dtype=np.float32, count=Nhsml)
	        self.Veldisp[locs] = np.fromfile(f, dtype=np.float32, count=Nhsml)
        	skip += Nhsml
        

      curpos = f.tell()
      f.seek(0,os.SEEK_END)
      if curpos != f.tell(): print "Warning: finished reading before EOF for file",filenum
      f.close()  
      filenum += 1
      if filenum == self.nfiles: doneflag = True
       
    if swap:
      self.Hsml.byteswap(True)
      self.Density.byteswap(True)
      self.Veldisp.byteswap(True)


