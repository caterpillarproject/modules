"""
HDF5toGadget2.py - Convert .hdf5 files to fileformat 1 Gadget2 output files
@author: Greg Dooley gdooley@mit.edu
@date: 06/15/2012. Edited to faster version May 2013

IMPORTANT - only writes out parttype=1 particles. ex: Will not work with
larger contamination particles.
"""
import readsnapshots.readsnapHDF5 as rs
import struct
import readsnapshots.readsnap
import numpy as np
import itertools
import sys

def Convert(base, snapnum, newbase):
    """
    @param base: Full file path containing snapdir folders
    @param snapnum: number of the snapshot
    @param newbase: Full file path to store new snapshot files.
    """
    filebase = base+'/snapdir_'+str(snapnum).zfill(3)+'/snap_'+str(snapnum).zfill(3)
    h = rs.snapshot_header(filebase+'.0')
    headerformat = 'iiiiiiiddddddddiiiiiiiiiidddd'
    
    for k in range(h.filenum):
        file = filebase+'.'+str(k)
        h = rs.snapshot_header(file)
        f = open(newbase+'/snap_'+str(snapnum).zfill(3)+'.'+str(k), 'wb')
        #print 'Starting on Header', snapnum, k
        ## write size of block 1 in bytes
        header_size = 256 #bytes

        #data = struct.pack(headerformat,header_size,h.npart[0],h.npart[1],h.npart[2],h.npart[3],h.npart[4],h.npart[5],h.massarr[0],h.massarr[1],h.massarr[2],h.massarr[3],h.massarr[4],h.massarr[5],h.time,h.redshift,h.sfr,h.feedback,h.nall[0],h.nall[1],h.nall[2],h.nall[3],h.nall[4],h.nall[5],h.cooling,h.filenum,h.boxsize,h.omega0,h.omegaL,h.hubble)
        #f.write(data)
        
        data = struct.pack('i',header_size)
        f.write(data)
        ## Write all the Header Information
        data = struct.pack('i',h.npart[0])
        f.write(data)
        data = struct.pack('i',h.npart[1]) 
        f.write(data)
        data = struct.pack('i',h.npart[2])
        f.write(data)
        data = struct.pack('i',h.npart[3])
        f.write(data)
        data = struct.pack('i',h.npart[4])
        f.write(data)
        data = struct.pack('i',h.npart[5])
        f.write(data)

        data = struct.pack('d',h.massarr[0])
        f.write(data)
        data = struct.pack('d',h.massarr[1])
        f.write(data)
        data = struct.pack('d',h.massarr[2])
        f.write(data)
        data = struct.pack('d',h.massarr[3])
        f.write(data)
        data = struct.pack('d',h.massarr[4])
        f.write(data)
        data = struct.pack('d',h.massarr[5])
        f.write(data)

        data = struct.pack('d',h.time)
        f.write(data)
        data = struct.pack('d',h.redshift)
        f.write(data)
        data = struct.pack('i',h.sfr)
        f.write(data)
        data = struct.pack('i',h.feedback)
        f.write(data)

        data = struct.pack('i',h.nall[0])
        f.write(data)
        data = struct.pack('i',h.nall[1])
        f.write(data)
        data = struct.pack('i',h.nall[2])
        f.write(data)
        data = struct.pack('i',h.nall[3])
        f.write(data)
        data = struct.pack('i',h.nall[4])
        f.write(data)
        data = struct.pack('i',h.nall[5])
        f.write(data)

        data = struct.pack('i',h.cooling)
        f.write(data)
        data = struct.pack('i',h.filenum)
        f.write(data)
        data = struct.pack('d',h.boxsize)
        f.write(data)
        data = struct.pack('d',h.omega0)
        f.write(data)
        data = struct.pack('d',h.omegaL)
        f.write(data)
        data = struct.pack('d',h.hubble)
        f.write(data)
        
        #print f.tell(), 'in bytes', k
        for i in range(header_size-f.tell()+4): # 4 for the size of header_size
            f.write('X')
        #print f.tell(), 'in bytes, before writing 256', k

        ## re-write size of block at the end
        data = struct.pack('i', header_size)
        f.write(data)
        ## end of writing header
        #print 'Finished Header'

        ### Write all the other blocks now
        ## Write positions
        if rs.contains_block(file+'.hdf5',"POS ",parttype=1):
            d1 = rs.read_block(file, "POS ",parttype=1,verbose=False)
            blocksize = struct.pack('I', len(d1)*8*3)
            f.write(blocksize)

            data = struct.pack('ddd'*len(d1), *d1.flatten())
            f.write(data)
            """
            for i in range(len(d1)):
                for j in range(3):
                    data = struct.pack('d',d1[i][j])
                    f.write(data)
            """
            f.write(blocksize)
        else:
            data = struct.pack('I',0)
            f.write(data)
            f.write(data)

        #print 'Finished Positions'
        
        # write velocities
        if rs.contains_block(file+'.hdf5',"VEL ",parttype=1):
            d1 = rs.read_block(file, "VEL ",parttype=1,verbose=False)
            blocksize = struct.pack('I', len(d1)*8*3)
            f.write(blocksize)


            data = struct.pack('ddd'*len(d1), *d1.flatten())
            f.write(data)
            """
            for i in range(len(d1)):
                for j in range(3):
                    data = struct.pack('d',d1[i][j])
                    f.write(data)
            """
            f.write(blocksize)
        else:
            data = struct.pack('I',0)
            f.write(data)
            f.write(data)

        #print 'Finished Velocities'

        # write Ids
        if rs.contains_block(file+'.hdf5',"ID  ",parttype=1):
            d1 = rs.read_block(file, "ID  ",parttype=1,verbose=False)
            blocksize = struct.pack('I', len(d1)*4) # 4 bytes to an id
            f.write(blocksize)
            data = struct.pack('I'*len(d1), *d1)
            f.write(data)
            """
            for i in range(len(d1)):
                data = struct.pack('I',d1[i])
                f.write(data)
            """
            f.write(blocksize)
        else:
            data = struct.pack('I',0)
            f.write(data)
            f.write(data)

        #print 'Finished Ids'

        ## write MASS information now
        if rs.contains_block(file+'.hdf5',"MASS",parttype=1):
            d1 = rs.read_block(file, "MASS",parttype=1,verbose=False)
            blocksize = struct.pack('d', len(d1)*8)
            f.write(blocksize)
            data = struct.pack('I'*len(d1), *d1)
            f.write(data)
            """
            for i in range(len(d1)):
                data = struct.pack('I',d1[i])
                f.write(data)
            """
            f.write(blocksize)
        else:
            data = struct.pack('I',0)
            f.write(data)
            f.write(data)

        #print 'Finished Mass'

        ## write Internal Energy
        if rs.contains_block(file+'.hdf5',"U   ",parttype=1):
            d1 = rs.read_block(file, "U   ",parttype=1,verbose=False)
            blocksize = struct.pack('d', len(d1)*8) # 4 bytes to an id
            f.write(blocksize)

            data = struct.pack('I'*len(d1), *d1)
            f.write(data)
            """
            for i in range(len(d1)):
                data = struct.pack('I',d1[i])
                f.write(data)
            """
            f.write(blocksize)
        else:
            data = struct.pack('I',0)
            f.write(data)
            f.write(data)

        #print 'Finished Internal Energy'

        ## write density
        if rs.contains_block(file+'.hdf5',"RHO ",parttype=1):
            d1 = rs.read_block(file, "RHO ",parttype=1,verbose=False)
            blocksize = struct.pack('d', len(d1)*8) # 4 bytes to an id
            f.write(blocksize)

            data = struct.pack('I'*len(d1), *d1)
            f.write(data)
            """
            for i in range(len(d1)):
                data = struct.pack('I',d1[i])
                f.write(data)
            """
            f.write(blocksize)
        else:
            data = struct.pack('I',0)
            f.write(data)
            f.write(data)

        #print 'Finished Internal Energy'

        ## NE - electron abundance
        if rs.contains_block(file+'.hdf5',"NE  ",parttype=1):
            d1 = rs.read_block(file, "NE  ",parttype=1,verbose=False)
            blocksize = struct.pack('d', len(d1)*8) # 4 bytes to an id
            f.write(blocksize)

            data = struct.pack('I'*len(d1), *d1)
            f.write(data)
            """
            for i in range(len(d1)):
                data = struct.pack('I',d1[i])
                f.write(data)
            """
            f.write(blocksize)
        else:
            data = struct.pack('I',0)
            f.write(data)
            f.write(data)

        #print 'Finished electron abundance'

        ## NH neutral hydrogen abundance
        if rs.contains_block(file+'.hdf5',"NH  ",parttype=1):
            d1 = rs.read_block(file, "NH  ",parttype=1,verbose=False)
            blocksize = struct.pack('d', len(d1)*8) # 4 bytes to an id
            f.write(blocksize)
            
            data = struct.pack('I'*len(d1), *d1)
            f.write(data)
            """
            for i in range(len(d1)):
                data = struct.pack('I',d1[i])
                f.write(data)
            """
            
            f.write(blocksize)
        else:
            data = struct.pack('I',0)
            f.write(data)
            f.write(data)

        #print 'Finished neutral hydrogen abundance'

        ## HSML
        if rs.contains_block(file+'.hdf5',"HSML",parttype=1):
            d1 = rs.read_block(file, "HSML",parttype=1,verbose=False)
            blocksize = struct.pack('d', len(d1)*8) # 4 bytes to an id
            f.write(blocksize)

            data = struct.pack('I'*len(d1), *d1)
            f.write(data)
            """
            for i in range(len(d1)):
                data = struct.pack('I',d1[i])
                f.write(data)
            """
            f.write(blocksize)
        else:
            data = struct.pack('I',0)
            f.write(data)
            f.write(data)

        #print 'Finished HSML'

        ## SFR
        if rs.contains_block(file+'.hdf5',"SFR ",parttype=1):
            d1 = rs.read_block(file, "SFR ",parttype=1,verbose=False)
            blocksize = struct.pack('d', len(d1)*8) # 4 bytes to an id
            f.write(blocksize)

            data = struct.pack('I'*len(d1), *d1)
            f.write(data)
            """
            for i in range(len(d1)):
                data = struct.pack('I',d1[i])
                f.write(data)
            """
            
            f.write(blocksize)
        else:
            data = struct.pack('I',0)
            f.write(data)
            f.write(data)

        #print 'Finished SFR'
   
        ## AGE
        if rs.contains_block(file+'.hdf5',"AGE ",parttype=1):
            d1 = rs.read_block(file, "AGE ",parttype=1,verbose=False)
            blocksize = struct.pack('d', len(d1)*8) # 4 bytes to an id
            f.write(blocksize)
            data = struct.pack('I'*len(d1), *d1)
            f.write(data)
            """
            for i in range(len(d1)):
                data = struct.pack('I',d1[i])
                f.write(data)
            """
            f.write(blocksize)
        else:
            data = struct.pack('I',0)
            f.write(data)
            f.write(data)

        #print 'Finished Age'

        ## Z
        if rs.contains_block(file+'.hdf5',"Z   ",parttype=1):
            d1 = rs.read_block(file, "Z   ",parttype=1,verbose=False)
            blocksize = struct.pack('d', len(d1)*8) # 4 bytes to an id
            f.write(blocksize)
            
            data = struct.pack('I'*len(d1), *d1)
            f.write(data)
            """
            for i in range(len(d1)):
                data = struct.pack('I',d1[i])
                f.write(data)
            """
            f.write(blocksize)
        else:
            data = struct.pack('I',0)
            f.write(data)
            f.write(data)

        #print 'Finished Z'

        #ACCE =?
        #d1 = rs.read_block(file, "ENDT",parttype=1,verbose=False)
        #d1 = rs.read_block(file, "TSTP",parttype=1,verbose=False)
        f.close()

if __name__ == "__main__":
    import sys
    import time
    if len(sys.argv) < 4:
        sys.exit('Usage: %s inbase snapnum outbase' %sys.argv[0])
    start = time.time()
    Convert(sys.argv[1], sys.argv[2], sys.argv[3])
    #print time.time()-start, 'time taken'

## Ordered list of items that were read in
"""
h.npart
h.massarr
h.time
h.redshift
h.sfr
h.feedback
h.nall
h.cooling
h.filenum
h.boxsize
h.omega0
h.omegaL
h.hubble

#h.stellar_age
h.nall_
#h.metals
#h.feedback
#h.double
"""
