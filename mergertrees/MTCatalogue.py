"""
Merger Tree updated version to improve speed


Old Strategy:
1) First read in file line by line searching for subhalos at a=1 and building
a python dictionary matching them to their host.
2) Read in file line by line building trees. Store halo information in a halo
object. Each time a halo is read in, search backwards in the tree (breadth
first search ordering) until its descendant is found. Connect parent and
descendant halo objects. (* note - this procedure is slow and inefficient).
3) Once host halo tree is complete, skip to subhalos in the ascii file and read
in those trees. Connect subhalo trees to the host halo tree.

New Strategy:

First 3 steps are pre-processed.

1) First read in file line by line searching for subhalos at a=1 and building
a python dictionary matching them to their host. Write this data to file in binary format to be easily read in later.

2) Convert ascii file to binary file, re-ordering it into the format of: host halo tree flag + host halo tree line by line, then subhalo tree flag + subhalo tree line by line etc. until no more subhalos. Then to second most massive host halo, etc.

3) Include ascii file with the scale-factors of everything. Make it an ascii file. Open all snapshots, read the header, get the scalefactor.

4) Load scale factor file.

5) Load binary file. Read 

read in ascii file line by line building tree.
Items kept:
Scale
id
desc_id
num_prog
pid
upid
SAM_Mvir
mvir
rvir
rs
vrms
mmp?
scale_of_last_MM
vmax
x
y
z
vx
vy
vz
Jx/Jy/Jz
Spin
Breadth First ID
Depth First ID
Orig_halo_ID
Xoff
Voff
Spin_Bullock



Items discarded:
desc_scale
desc_pid
Tree_root_ID
Next_copogenitor_depthfirst_ID
M200c
M500c
M2500c
Rs_klypin
"""

import time
import struct
import csv
import sys
import numpy as np

def storeSubInfo(file):
    """
    search linearly, get line number and host halo ID for all a=0 halos.
    then put information in dictionary. Each hostID has a corresponding
    list of line numbers of its hosts.
    """
    lines = []
    hostIDs = []
    f = open(file,'r')
    line = f.readline()
    while line!='':
        if line[0:5] == "#tree":
            num = int(line[6::])
            loc = f.tell()
            line = f.readline()
            upid = int(line.split()[6]) # most massive host
            if upid != -1: #if subhalo, store hostID and location
                hostIDs.append(upid)
                lines.append(loc)
        line = f.readline()
    f.close()
    # Now create python dictionary (hash table) of hosts and sub locations
    host2sub = {}
    for i in range(len(hostIDs)):
        if hostIDs[i] in host2sub:
            host2sub[hostIDs[i]].append(lines[i])
        else:
            host2sub[hostIDs[i]]=[lines[i]]
    return host2sub

def writeline(line,fout,fmt):
    s = line.split()
    data = struct.pack(fmt,
                       float(s[0]), #scale
                       int(s[1]),   #id
                       int(s[3]),   #desc_id
                       int(s[4]),   #num_prog
                       int(s[5]),   #pid
                       int(s[6]),   #upid
                       float(s[9]), #sam_mvir
                       float(s[10]),#mvir
                       float(s[11]),#rvir
                       float(s[12]),#rs
                       float(s[13]),#vrms
                       int(s[14]),  #mmp
                       float(s[15]),#scale_of_last_MM
                       float(s[16]),#vmax
                       float(s[17]),#x
                       float(s[18]),#y
                       float(s[19]),#z
                       float(s[20]),#vx
                       float(s[21]),#vy
                       float(s[22]),#vz
                       float(s[23]),#jx
                       float(s[24]),#jy
                       float(s[25]),#jz
                       float(s[26]),#spin
                       int(s[27]),  #Breadth_first_ID
                       int(s[28]),  #Depth_first_ID
                       int(s[30]),  #Orig_halo_ID
                       int(s[33]),  #Last_progenitor_depthfirst_ID
                       float(s[35]),#M200c_all
                       float(s[36]),#M200b
                       float(s[40]),#Xoff
                       float(s[41]),#Voff
                       float(s[42]),#Spin_Bullock
                       float(s[43]),#b_to_a
                       float(s[44]),#c_to_a
                       float(s[45]),#A[x]
                       float(s[46]),#A[y]
                       float(s[47]),#A[z]
                       float(s[48]))#T/|U|
    fout.write(data)

def writeline_old(line,fout,fmt):
    s = line.split()
    data = struct.pack(fmt,
                       float(s[0]), #scale
                       int(s[1]),   #id
                       int(s[3]),   #desc_id
                       int(s[4]),   #num_prog
                       int(s[5]),   #pid
                       int(s[6]),   #upid
                       float(s[9]), #mvir
                       float(s[10]),#orig_mvir
                       float(s[11]),#rvir
                       float(s[12]),#rs
                       float(s[13]),#vrms
                       int(s[14]),  #mmp
                       float(s[15]),#scale_of_last_MM
                       float(s[16]),#vmax
                       float(s[17]),#x
                       float(s[18]),#y
                       float(s[19]),#z
                       float(s[20]),#vx
                       float(s[21]),#vy
                       float(s[22]),#vz
                       float(s[23]),#jx
                       float(s[24]),#jy
                       float(s[25]),#jz
                       float(s[26]),#spin
                       int(s[27]),  #Breadth_first_ID
                       int(s[28]),  #Depth_first_ID
                       int(s[30]))  #Orig_halo_ID
    fout.write(data)


def convertmt(dir,time_me=False,oldversion=False,verbose=False):
    """
    
    """
    filenamein = dir+"/tree_0_0_0.dat"
    filenameout = dir+"/tree.bin"
    fileindexname = dir+"/treeindex.csv"

    if time_me:
        print "Reading subhalo positions"
        start = time.time()
    host2sub = storeSubInfo(filenamein)
    if time_me:
        print time.time()-start, 'Time to get subhalo positions'

    fin = open(filenamein,'r')
    fout= open(filenameout,'wb')
    findex=open(fileindexname,'w')
    writer = csv.writer(findex)

    if oldversion:
        fmt = "fiiiiifffffiffffffffffffiii"
    else:
        fmt = "fiiiiifffffiffffffffffffiiiifffffffffff"
    fmtsize = struct.calcsize(fmt)
    if verbose==True:
        print "Bytes per line:",fmtsize

    line = fin.readline()
    i = 0
    ##TODO save scale factor <-> snap number
    while line != '':
        if line[0:5] == "#tree":
            line = fin.readline()
            hostsplit = line.split()
            ## if host halo, readwrite this halo and its subhalos
            if int(hostsplit[5]) == -1:
                if time_me:
                    print "Reading in host halo"
                    start = time.time()
                #save host ID and file location to index
                writer.writerow([hostsplit[30],fout.tell()])
                #write header tag and dummy numlines
                fout.write(struct.pack("ii",0,-1))
                numlines = 0
                #write host tree
                start2 = time.time()
                while line != '' and line[0:5] != "#tree":
                    if oldversion:
                        writeline_old(line,fout,fmt)
                    else:
                        writeline(line,fout,fmt)
                    numlines = numlines+1
                    line = fin.readline()
                #seek back and fill in numlines
                fout.seek(-1*(numlines*fmtsize + struct.calcsize("i")),1) #from current location backwards
                fout.write(struct.pack("i",numlines))
                fout.seek(numlines*fmtsize,1) #from current location forwards
                #remember host location to come back to after writing subs
                host_loc = fin.tell()

                ## read/write in all the subhalos
                try:
                    if time_me:
                        print time.time()-start, 'Time to read host halo',i
                        print "Reading in subhalos"
                        start = time.time()
                    for subhalo_loc in host2sub[int(hostsplit[1])]: #key is host halo id
                        #write header tag and dummy numlines
                        fout.write(struct.pack("ii",1,-1))
                        numlines = 0
                        #write sub tree
                        fin.seek(subhalo_loc) #put file head at first line of subhalo
                        line = fin.readline()
                        while line != '' and line[0:5] != "#tree":
                            numlines = numlines+1
                            if oldversion:
                                writeline_old(line,fout,fmt)
                            else:
                                writeline(line,fout,fmt)
                            line = fin.readline()
                        #seek back and fill in num lines
                        fout.seek(-1*(numlines*fmtsize + struct.calcsize("i")),1) #from current location backwards
                        fout.write(struct.pack("i",numlines))
                        fout.seek(numlines*fmtsize,1) #from current location forwards
                    fin.seek(host_loc)
                    if time_me:
                        print time.time()-start, 'Time to read all subhalos'
                except KeyError: #Host with no subs, do nothing
                    pass
            ## else, is a subhalo, which means it will have already been read/written
            else: 
                #skip down to next tree
                while line != '' and line[0:5] != "#tree":
                    line = fin.readline()
            i = i+1
            if i % 5000 == 0 and verbose==True:
                print i,"host halos completed"
        else:
            line = fin.readline()
    fin.close()
    fout.close()
    findex.close()

class NewMTCatalogueTree:
    """
    Creates a Tree
    Either from an input datatable (e.g. used internally to make subtrees)
    Or from a file that is already pointing to the right place (e.g. from NewMTCatalogue)
    """
    def __init__(self,
                 datatable=np.dtype([]),
                 f=None,halotype=-1,nrow=-1,fmt="",fmttype=np.dtype([])):
        if datatable != np.dtype([]):
            self.fileloc = -1
            self.halotype = 2
            self.nrow = len(datatable)
            self.data = datatable
            self.rockstar_id = self.data[0]['origid']
        else:
            if f==None or halotype==-1 or nrow==-1 or fmt=="" or fmttype==np.dtype([]):
                print "ERROR: must specify all of these variables:"
                print "f,halotype,nrow,fmt,fmttype"
                raise RuntimeException("Didn't specify all variables")
            self.fileloc = f.tell()
            self.halotype=halotype #0 or 1
            self.nrow = nrow
            self.data = np.zeros(self.nrow,dtype=fmttype)
            fmtsize = struct.calcsize(fmt)
            for i in xrange(nrow):
                self.data[i] = struct.unpack(fmt,f.read(fmtsize))
            self.rockstar_id = self.data[0]['origid']

    def getMainBranch(self, row):
        """
        repeatedly call getMMP until end is reached.
        @return: all halos that are in the main branch of the halo
        specified by row.
        """
        rows = [row]
        next = self.getMMP(row)
        while next != None:
            rows.append(next)
            next = self.getMMP(next)
        return self.data[rows]

    def getMMP(self, row):
        """
        @param row: row number (int) of halo considered
        @ return: row number of the most massive parent.
        """
        parent = np.where(np.logical_and(self.data[row]['id']==self.data['desc_id'], self.data['mmp']==1))[0]
        if len(parent) == 0:
            return None # if it has no parent
        else:
            return parent[0]

    def getNonMMPprogenitors(self,row):
        # return row index of all progenitors that are not the most massive
        # These are all the subhalos that were destroyed in the interval
        return np.where(np.logical_and(self.data[row]['id']==self.data['desc_id'], self.data['mmp']!=1))[0]

    def __getitem__(self,key):
        return self.data[key]
    
class MTCatalogue:
    """
    Read in the binary file created by convertmt

    Two ways of creating this object.
    The first way is to read in the first N hosts (and their subs).
    Specify numHosts = N (default reads in everything) to read N hosts.
    Also specify byID=False to access trees in mass order (i.e. self.trees[0]
    self.trees[1], self.trees[2] for 3 most massive hosts)
    By default trees must be accessed by their rockstar halo ID.
    
    The second way is to specify the indexfile (csv file) and a haloid.
    This will use indexfile to pinpoint the location in datafile corresponding
    to haloid and read in that host halo with all the subhalos as the catalogue.
    """
    # need array pointing to location of hosts.
    def __init__(self,file,
                 indexfile=None,hostid=-1,oldversion=False,numHosts=np.infty,byID=True,verbose=False):
        self.file = file
        self.Trees = {} #key: rockstar halo ID; value: MT file
	if oldversion:
	    self.fmt = "fiiiiifffffiffffffffffffiii"
	else:	        
	    self.fmt = "fiiiiifffffiffffffffffffiiiifffffffffff" 
        self.fmtsize = struct.calcsize(self.fmt)	
        ## Create structured dtype

        if oldversion:
             self.fmttype = np.dtype([('scale','<f8'),('id','<i8'),('desc_id','<i8'),
                                     ('num_prog','<i8'),('pid','<i8'),('upid','<i8'),
                                     ('mvir','<f8'),('orig_mvir','<f8'),
                                     ('rvir','<f8'),('rs','<f8'),
                                     ('vrms','<f8'),('mmp','<i8'),
                                     ('scale_of_last_MM','<f8'),('vmax','<f8'),
                                     ('posX','<f8'),('posY','<f8'),('posZ','<f8'),
                                     ('pecVX','<f8'),('pecVY','<f8'),('pecVZ','<f8'),
                                     ('Jx','<f8'),('Jy','<f8'),('Jz','<f8'),('spin','<f8'),
                                     ('bfid','<i8'), #breadth first ID
                                     ('dfid','<i8'), #depth first ID
                                     ('origid','<i8')]) #rockstar cat ID
        else:
	    self.fmttype = np.dtype([('scale','<f8'),('id','<i8'),('desc_id','<i8'),
                                     ('num_prog','<i8'),('pid','<i8'),('upid','<i8'),
                                     ('sam_mvir','<f8'),('mvir','<f8'),
                                     ('rvir','<f8'),('rs','<f8'),
                                     ('vrms','<f8'),('mmp','<i8'),
                                     ('scale_of_last_MM','<f8'),('vmax','<f8'),
                                     ('posX','<f8'),('posY','<f8'),('posZ','<f8'),
                                     ('pecVX','<f8'),('pecVY','<f8'),('pecVZ','<f8'),
                                     ('Jx','<f8'),('Jy','<f8'),('Jz','<f8'),('spin','<f8'),
                                     ('bfid','<i8'), #breadth first ID
                                     ('dfid','<i8'), #depth first ID
                                     ('origid','<i8'), #rockstar cat ID
                                     ('lastprog_dfid','<i8'), #depth first ID last progenitor
                                     ('m200c_all','<f8'),('m200b','<f8'),
                                     ('xoff','<f8'),('voff','<f8'),
                                     ('spin_bullock','<f8'),
                                     ('b_to_a','<f8'),('c_to_a','<f8'),
                                     ('A[x]','<f8'),('A[y]','<f8'),('A[z]','<f8'),
                                     ('T/|U|','<f8')])

        f = open(file,'rb')
        if indexfile==None and hostid==-1:
            #Read everything!
            if verbose:
                print "Reading whole catalogue"
            start = time.time()
            tag = f.read(8)
            counter = 0
            nhosts = 0
            self.HostLocs = []
            while tag != '' and nhosts <= numHosts:
                halotype,nrow = struct.unpack("ii",tag)
                thistree = NewMTCatalogueTree(f=f,halotype=halotype,nrow=nrow,fmt=self.fmt,fmttype=self.fmttype)
                if byID: #Good for retrieving known trees.
                    # bad for computing averages of trees.
                    rsid = thistree.rockstar_id
                    self.Trees[rsid] = thistree
                else:
                    self.Trees[counter] = thistree
                tag = f.read(8)
                if halotype == 0:
                    self.HostLocs.append(counter)
                    nhosts+=1
                counter+=1
            if verbose:
                print "Time to finish reading:",time.time()-start
        else:
            #Read just the host halo and its subs
            reader = csv.reader(open(indexfile,'r'))
            index = dict(x for x in reader)
            try:
                file_loc = int(index[str(hostid)]) #raises KeyError if problem
                #Read in tree for host halo
                f.seek(file_loc)
                tag = f.read(8)
                halotype,nrow = struct.unpack("ii",tag)
                if halotype != 0:
                    raise ValueError
                hosttree = NewMTCatalogueTree(f=f,halotype=halotype,nrow=nrow,fmt=self.fmt,fmttype=self.fmttype)
                rsid = hosttree.rockstar_id
                if rsid != hostid:
                    raise ValueError
                self.Trees[hostid]=hosttree
                tag = f.read(8)
                halotype,nrow = struct.unpack("ii",tag)
                #Read in trees for all subhalos
                while halotype==1:
                    thistree = NewMTCatalogueTree(f=f,halotype=halotype,nrow=nrow,fmt=self.fmt,fmttype=self.fmttype)
                    rsid = thistree.rockstar_id
                    self.Trees[rsid] = thistree
                    tag = f.read(8)
                    halotype,nrow = struct.unpack("ii",tag)                    
            except KeyError:
                print "ERROR: No halo with this ID in the listed index!"
                print "Did you put in the ID of a host halo? The index Alex made only has host halos"
                print "(Catalogue object is still created but empty)"
            except ValueError:
                print "ERROR: Problem with the index. Did not point to a host halo"
                print "or did not point to a host halo that corresponded to the input halo"
                print "(Catalogue object is still created but empty)"
        f.close()


    def getSubTrees(self,treenum):
        """
        returns list of indices of subhalo trees of tree given by treenum
        """
        if self.Trees[treenum].halotype !=0:
            print 'ERROR: input halo not a host'
            return []
        row = treenum+1
        while self.Trees[row].halotype == 1:
            row+=1
        return np.arange(treenum+1, row)


## if __name__ == "__main__":
##     convertmt(time_me=False)

## Read file once to get subhalo position information in dictionary.

## treefile ='/spacebase/data/AnnaGroup/caterpillarparent/RockstarData/trees/tree_0_0_0.dat'

## print 'Getting subhalo positions'
## start = time.time()
## host2sub = storeSubInfo(treefile)
## print time.time()-start, 'Time to get subhalo positions'

## f = open(treefile)
## [f.seek(posn) for posn in host2sub[host] for host in host2sub.keys()[0:1]]
