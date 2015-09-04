"""
Merger Tree updated version to improve speed

Includes MTCatalogueTree, which has simple visualization built in

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
Scale,id,desc_id,num_prog,pid,upid,phantom,sam_mvir,mvir,rvir,rs,vrms,mmp,scale_of_last_MM,
vmax,posX,posY,posZ,pecVX,pecVY,pecVZ,Jx,Jy,Jz,spin,bfid,dfid,origid,lastprog_dfid,
m200c_all,m200b,xoff,voff,spin_bullock,b_to_a,c_to_a,A[X/Y/Z],T/|U|

Items discarded:
desc_scale,desc_pid,Tree_root_ID,snapnum,Next_coprogenitor_depthfirst_ID,M200c,M500c,M2500c,Rs_klypin
"""

import time
import struct
import csv
import sys
import numpy as np
import itertools
import os
import pydot
import subprocess
import glob
import functools
import asciitable

def getScaleFactors(path,minsnap=0,digits=3,returnminsnap=False):
    scalefile = path+'/outputs/scales.txt'

    f = open(scalefile)
    line = f.readline()
    linesplit = line.split()
    # pad beginning with extra -1's 
    minsnap = int(linesplit[0])
    scale_list = [-1 for x in xrange(minsnap)]
    # read scales
    while line != '':
        scale_list.append(float(line.split()[1]))
        line = f.readline()
    f.close()
    if returnminsnap: return minsnap,np.array(scale_list)
    return np.array(scale_list)

def writeline_v4(line,fout,fmt,snapoffset):
    s = line.split()
    data = struct.pack(fmt,
                       float(s[0]), #scale
                       int(s[1]),   #id
                       int(s[3]),   #desc_id
                       int(s[4]),   #num_prog
                       int(s[5]),   #pid
                       int(s[6]),   #upid
                       int(s[8]),   #phantom
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
                       int(s[27]),  #bfid [Breadth_first_ID]
                       int(s[28]),  #dfid [Depth_first_ID]
                       int(s[30]),  #origid [Orig_halo_ID]
                       int(s[33]),  #lastprog_dfid [Last_progenitor_depthfirst_ID]
                       float(s[35]),#m200c_all
                       float(s[36]),#m200b
                       float(s[40]),#xoff
                       float(s[41]),#voff
                       float(s[42]),#spin_bullock
                       float(s[48]),#b_to_a(500c)
                       float(s[49]),#c_to_a(500c)
                       float(s[50]),#A[x](500c)
                       float(s[51]),#A[y](500c)
                       float(s[52]),#A[z](500c)
                       float(s[53]),#T/|U|
                       int(s[31])+snapoffset)  #snap
    fout.write(data)

def writeline_v3(line,fout,fmt):
    s = line.split()
    data = struct.pack(fmt,
                       float(s[0]), #scale
                       int(s[1]),   #id
                       int(s[3]),   #desc_id
                       int(s[4]),   #num_prog
                       int(s[5]),   #pid
                       int(s[6]),   #upid
                       int(s[8]),   #phantom
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
                       int(s[27]),  #bfid [Breadth_first_ID]
                       int(s[28]),  #dfid [Depth_first_ID]
                       int(s[30]),  #origid [Orig_halo_ID]
                       int(s[33]),  #lastprog_dfid [Last_progenitor_depthfirst_ID]
                       float(s[35]),#m200c_all
                       float(s[36]),#m200b
                       float(s[40]),#xoff
                       float(s[41]),#voff
                       float(s[42]),#spin_bullock
                       float(s[48]),#b_to_a(500c)
                       float(s[49]),#c_to_a(500c)
                       float(s[50]),#A[x](500c)
                       float(s[51]),#A[y](500c)
                       float(s[52]),#A[z](500c)
                       float(s[53]))#T/|U|
    fout.write(data)

def writeline_v2(line,fout,fmt):
    s = line.split()
    data = struct.pack(fmt,
                       float(s[0]), #scale
                       int(s[1]),   #id
                       int(s[3]),   #desc_id
                       int(s[4]),   #num_prog
                       int(s[5]),   #pid
                       int(s[6]),   #upid
                       int(s[8]),   #phantom
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
                       int(s[27]),  #bfid [Breadth_first_ID]
                       int(s[28]),  #dfid [Depth_first_ID]
                       int(s[30]),  #origid [Orig_halo_ID]
                       int(s[33]),  #lastprog_dfid [Last_progenitor_depthfirst_ID]
                       float(s[35]),#m200c_all
                       float(s[36]),#m200b
                       float(s[40]),#xoff
                       float(s[41]),#voff
                       float(s[42]),#spin_bullock
                       float(s[43]),#b_to_a
                       float(s[44]),#c_to_a
                       float(s[45]),#A[x]
                       float(s[46]),#A[y]
                       float(s[47]),#A[z]
                       float(s[48]))#T/|U|
    fout.write(data)

def writeline_v1(line,fout,fmt):
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

def _read_host(fin, hostline, mywriteline, fmt, fmtsize, fout):
    fout.write(struct.pack("ii",0,-1))
    mywriteline(hostline,fout,fmt)
    numlines = 1
    line = fin.readline()
    while line != "" and line[0:5] != "#tree":
        mywriteline(line,fout,fmt) #error used to be here
        numlines += 1
        lastloc = fin.tell()
        line = fin.readline()
    #seek back and fill in numlines
    fout.seek(-1*(numlines*fmtsize + struct.calcsize("i")),1) #from current location backwards
    fout.write(struct.pack("i",numlines))
    fout.seek(numlines*fmtsize,1) #from current location forwards
    return line,lastloc

def _read_subhalos(fsubfiles,subloclist,mywriteline,fmt,fmtsize,fout):
    if len(subloclist) > 0:
        for fileloc in subloclist:
            #pick the right file to look at
            ifile,subhalo_loc = fileloc
            fin = fsubfiles[ifile]
            #write header tag and dummy numlines
            fout.write(struct.pack("ii",1,-1))
            numlines = 0
            #write sub tree
            fin.seek(subhalo_loc) #put file head at first line of subhalo
            line = fin.readline()
            while line != '' and line[0:5] != "#tree":
                mywriteline(line,fout,fmt)
                numlines += 1
                line = fin.readline()
            #seek back and fill in num lines
            fout.seek(-1*(numlines*fmtsize + struct.calcsize("i")),1) #from current location backwards
            fout.write(struct.pack("i",numlines))
            fout.seek(numlines*fmtsize,1) #from current location forwards
        return len(subloclist)
    else:
        return 0

def _skip_tree(fin):
    line = fin.readline()
    while line != "" and line[0:5] != "#tree":
        line = fin.readline()
    return line

def storeSubInfo(fpaths):
    """
    search linearly, get line number and host halo ID for all a=1 halos .
    then put information in dictionary. Each hostID has a corresponding
    list of line numbers of its subs.
    """
    def get_corrected_upid(treeid,treeid2upid):
        # recursively go up the upid chain until you find a host halo
        upid = treeid2upid[treeid]
        if upid == -1:
            return treeid
        else:
            return get_corrected_upid(upid,treeid2upid)

    host2sub = {}
    treeid2upid = {}

    for ifile,fpath in enumerate(fpaths):
        f = open(fpath,'r')
        line = f.readline()
        while line!='':
            if line[0:5] == "#tree":
                treeid = int(line[6::]) #tree mtid
                loc = f.tell()
                
                line = f.readline(); linesplit = line.split()
                upid = int(linesplit[6]) # ultimate parent
                #pid = int(linesplit[5]) # immediate parent
                scale = float(linesplit[0]) #scale
                if scale != 1:
                    print "%i not rooted at a=1 (at a=%3.2f)" % (treeid,scale)

                treeid2upid[treeid] = upid
                if upid != -1: #if subhalo, store hostID and location
                    #assume the trees are sorted so host appears before sub (e.g. by mass)
                    #corrected_upid = get_corrected_upid(treeid,treeid2upid)
                    #IMPORTANT: there is an edge case where upid is not a host halo
                    #For example, if a subsubhalo is inside a subhalo, but not in the host halo.
                    corrected_upid = upid
                    if corrected_upid in host2sub:
                        host2sub[corrected_upid].append((ifile,loc))
                    else:
                        host2sub[corrected_upid] = [(ifile,loc)]
            line = f.readline()
        f.close()
    return host2sub

def convertmt(dir,version=2,verbose=True):
    """
    Converts dir/tree_*.dat into dir/tree.bin (with dir/treeindex.csv also output)
    Should correctly deal with multiple tree_x_x_x.dat files now.

    aka turns ASCII output into binary output, while also resorting the data to allow
      fast retrieval of all trees (host and subs) from a single host halo

    IMPORTANT: there is an edge case where upid is not a host halo
      For example, if a subsubhalo is inside a subhalo, but not in the host halo,
      then the tree for that subsubhalo may NOT be written out!
      There are potentially other edge cases.
      
    version=1: very old Rockstar (e.g. Greg's paper)
    version=2: RC1 before the shape parameters were fixed (default)
    version=3: RC2, probably correct for RC3 (also works for RC3+)
    version=4: same as version=3 but with snap written out
    """

    fpaths = glob.glob(dir+'/tree_*.dat')
    fnames = [os.path.basename(os.path.normpath(fpath)) for fpath in fpaths]
    print "acting over these "+str(len(fnames))+" files"
    print fnames
    filenameout = dir+'/tree.bin'
    fileindexname = dir+'/treeindex.csv'

    print "Storing Sub Info (takes a while)..."
    start = time.time()
    host2sub = storeSubInfo(fpaths)
    print "Took %.2f sec" % (time.time()-start)

    fsubfiles = [open(fpath,'r') for fpath in fpaths]
    fout= open(filenameout,'wb')
    findex=open(fileindexname,'w')
    writer = csv.writer(findex)

    if version==1:
        fmt = "fiiiiifffffiffffffffffffiii"
        mywriteline = writeline_v1
    elif version==2:
        fmt = "fiiiiiifffffiffffffffffffiiiifffffffffff"
        mywriteline = writeline_v2
    elif version==3:
        fmt = "fiiiiiifffffiffffffffffffiiiifffffffffff"
        mywriteline = writeline_v3
    elif version==4:
        fmt = "fiiiiiifffffiffffffffffffiiiifffffffffffi"
        snapoffset,scalelist = getScaleFactors(dir+'/..',returnminsnap=True)
        mywriteline = functools.partial(writeline_v4,snapoffset=snapoffset)
    else:
        raise ValueError("ERROR, version must be 1, 2, 3, or 4")
    fmtsize = struct.calcsize(fmt)
    if verbose: print "Bytes per line:",fmtsize

    host_counter = 0
    skip_counter = 0
    comment_counter = 0
    written_subhalos_counter = 0
    host_with_no_subs_counter = 0
    numtrees=0

    for fpath in fpaths:
        fin = open(fpath,'r')
        line = fin.readline()
        nlineskip = 1 #skip the line with the number of trees in this file
        start = time.time()
        while True: #loop through this file
            if line == "":
                break
            elif line[0:5] == "#tree":
                hostline = fin.readline(); hostsplit = hostline.split()
                if int(hostsplit[5]) == -1: #is host
                    hostid = int(hostsplit[30])
                    writer.writerow([hostid,fout.tell()])
                    line,lastloc = _read_host(fin,hostline,mywriteline,fmt,fmtsize,fout)
                    assert line[0:5] == "#tree" or line == ""
                    host_mtid = int(hostsplit[1])
                    if host_mtid in host2sub:
                        subloclist = host2sub[host_mtid]
                        written_subhalos_counter += _read_subhalos(fsubfiles,subloclist,mywriteline,fmt,fmtsize,fout)
                    else:
                        host_with_no_subs_counter += 1
                    fin.seek(lastloc) #reset back to the end of the host
                    line = fin.readline()
                    assert line[0:5] == "#tree" or line == ""
                    host_counter += 1
                else: #is subhalo
                    line = _skip_tree(fin)
                    skip_counter += 1
            elif line[0] == "#":
                line = fin.readline()
                comment_counter += 1
            elif nlineskip > 0:
                nlineskip -= 1
                numtrees += int(line)
                #print nlineskip,line
                line = fin.readline()
            else:
                print line
                raise RuntimeError("line is not a #tree, comment, or empty string indicating end of file")
        fin.close()

        print "finished "+os.path.basename(os.path.normpath(fpath))
        print "%.2f" % (time.time()-start)

    fout.close()
    findex.close()
    for f in fsubfiles:
        f.close()
    if verbose:
        print "Should have",numtrees,"trees"
        print "Found",comment_counter,"comment lines"
        print "Found",host_counter,"hosts"
        print "Found",host_with_no_subs_counter,"hosts with no subs"
        print "Found",skip_counter,"total trees with no parent (subs)"
        print "Wrote",written_subhalos_counter,"subs (found using host2sub)"

class MTCatalogueTree:
    """
    Creates a Tree
    Either from an input datatable (e.g. used internally to make subtrees)
    Or from a file that is already pointing to the right place (e.g. from MTCatalogue)
    """
    def __init__(self,scale_list=[],datatable=np.array([]),
                 f=None,halotype=-1,nrow=-1,fmt="",fmttype=np.dtype([]),
                 readtxt=False,txtheader=None):
        self.scale_list = scale_list
        if datatable != np.array([]):
            self.fileloc = -1
            self.halotype = 2
            self.nrow = len(datatable)
            self.data = datatable
            self.rockstar_id = self.data[0]['origid']
        else:
            if readtxt:
                assert txtheader != None
                self.fileloc = f.tell()
                line = f.readline() #first row of the MT
                mtlist = []
                while line != "" and line[0:5] != "#tree":
                    mtlist.append(line.split())
                    line = f.readline()
                self.data = asciitable.read(mtlist,Reader=asciitable.Memory,names=txtheader)
                self.data = self.data.astype(fmttype)
                self.halotype = halotype
                self.nrow = len(self.data)
                self.rockstar_id = self.data[0]['origid']
            else:
                if f==None or halotype==-1 or nrow==-1 or fmt=="" or fmttype==np.dtype([]):
                    print "ERROR: must specify all of these variables:"
                    print "f,halotype,nrow,fmt,fmttype"
                    raise RuntimeError("Didn't specify all variables")
                self.fileloc = f.tell()
                self.halotype=halotype #0 or 1
                self.nrow = nrow
                self.data = np.zeros(self.nrow,dtype=fmttype)
                fmtsize = struct.calcsize(fmt)
                for i in xrange(nrow):
                    self.data[i] = struct.unpack(fmt,f.read(fmtsize))
                self.rockstar_id = self.data[0]['origid']

    def plottree(self,filename='mytree',makepdf=True,savedot=False,mask=None):
        """
        http://www.graphviz.org/doc/info/attrs.html
        @param filename: the file name you want to save to. Don't put on a file extension,
                         e.g. if you want 'mytree.pdf' then filename='mytree'
        @param makepdf: defaults to True, writes a .ps2 and converts to .pdf using the
                        shell's convert program
        @param savedot: writes out the .dot file. May be useful if your graph is ginormous
                        and takes too long to compute.
        @param mask: UNDER DEVELOPMENT an optional boolean array to mask the tree's data. 
                     Intended for pruning small objects out of gigantic trees when plotting.
                     WILL NOT WORK if you prune a halo from the middle of a branch!
                     (It'll add at least one more halo after the one you want,
                     it won't reconnect halos if you splice one in the middle,
                     and no guarantees that it will work properly even after that.)
                     The user is responsible for doing something that won't break the
                     plotting (i.e., don't get rid of the base halo), and for keeping
                     track of what mask is applied.
                     Note the mask is just for making the plot, i.e. it will not be saved
                     into the tree's internal data structures.
        TODO allow some flexibility in coloring/text? I don't think so, it may be best
             to write more functions that can modify the returned graph object.

        @return graph: the pydot.Dot graph object
        """
        if makepdf and filename[-4:] == '.pdf': filename = filename[:-4]
        if mask != None:
            plotdata = self.data[mask]
        else:
            plotdata = self.data
        print "Generating graph"
        maxvalue = np.max(plotdata['rvir'])

        graph = pydot.Dot(graph_type='graph',size="8, 8")
        for row in reversed(xrange(len(plotdata))):
            nodesize = (plotdata[row]['rvir']/maxvalue)
            if nodesize < 0.01: continue #skip over things too small to plot
            #plotdata[row]['scale'],
            label = "%i \n %i" % (plotdata[row]['id'],plotdata[row]['lastprog_dfid'])
            graph.add_node(pydot.Node(plotdata[row]['id'],
                                      shape='circle',fixedsize="true",
                                      width=nodesize,#label=str(plotdata[row]['id'])
                                      #,#height=nodesize,
                                      #label=" "
                                      label=label #,fontsize=str(nodesize)
                                      ))
            graph.add_edge(pydot.Edge(plotdata[row]['id'],plotdata[row]['desc_id']))
        #Delete the extra last node
        graph.del_edge(plotdata[0]['id'],plotdata[0]['desc_id'])
        graph.del_node(plotdata[0]['desc_id'])
        if savedot:
            print "Writing dot file"
            graph.write_dot(filename+'.dot')
        if makepdf:
            print "Writing "+str(len(graph.get_nodes()))+" nodes and "+str(len(graph.get_edges()))+" edges to "+filename+'.ps2'
            print "(may freeze/take hours if graph is too big)"
            graph.write_ps2(filename+'.ps2')
            print "Converting to "+filename+'.pdf'
            subprocess.call(["convert",filename+'.ps2',filename+'.pdf'])
            print "Removing "+filename+'.ps2'
            subprocess.call(["rm",filename+'.ps2'])
        return graph

    def getSubTree(self,row,dfid=None,lastprog_dfid=None):
        """
        Returns a MTCatalogueTree object that is the subtree of the halo specified by row
        Uses depth-first ID to get the subtree
        (and may be a shallow copy of the data, so don't delete the original tree)
        """
        if dfid != None and lastprog_dfid != None:
            dfid_base = dfid
            dfid_last = lastprog_dfid
        else:
            dfid_base = self.data[row]['dfid']
            dfid_last = self.data[row]['lastprog_dfid']
        mask = np.logical_and(self.data['dfid']>=dfid_base, self.data['dfid']<=dfid_last)
        return MTCatalogueTree(scale_list=self.scale_list,datatable=self.data[mask])




    def get_mmp_map(self):
        mask = self.data['mmp']==1
        return dict(zip(self.data[mask]['desc_id'], np.arange(len(self.data))[mask]))

    def get_non_mmp_map(self):
        non_mmp_map={}
        mask2 = self.data['mmp']==0
        for key, val in zip(self.data[mask2]['desc_id'], np.arange(len(self.data))[mask2]):
            non_mmp_map.setdefault(key, []).append(val)

    def getMainBranch(self, row=0, mmp_map=None):
        """
        @param row: row of the halo you want the main branch for. Defaults to row 0
        Uses getSubTree, then finds the smallest dfid that has no progenitors
        @return: all halos that are in the main branch of the halo specified by row (in a np structured array)
        """
        if mmp_map is None:
            subtree = self.getSubTree(row)
            dfid_base = self.data[row]['dfid'] #==subtree[0]['dfid']
            # Get the smallest dfid of the tree "leaves" (halos with no progenitors)
            dfid_last = np.min(subtree[subtree['num_prog']==0]['dfid'])
            return subtree[np.logical_and(subtree['dfid']>=dfid_base,subtree['dfid']<=dfid_last)]
        else:
            rows=[]
            while not row is None:
                rows.append(row)
                row = self.getMMP(row,mmp_map)
            return self.data[rows]
            
    def getMMP(self, row, mmp_map=None):
        """
        @param row: row number (int) of halo considered
        @return: row number of the most massive parent, or None if no parent
        """
        if mmp_map is None:
            parent = np.where(np.logical_and(self.data[row]['id']==self.data['desc_id'], self.data['mmp']==1))[0]
            try:
                return parent[0]
            except:
                return None # if it has no parent
        else:
            try:
                return mmp_map[self.data[row]['id']]
            except:
                return None

    def getNonMMPprogenitors(self,row, non_mmp_map=None):
        """
        return row index of all progenitors that are not the most massive
        These are all the subhalos that were destroyed in the interval
        """
        if non_mmp_map is None:
            return np.where(np.logical_and(self.data[row]['id']==self.data['desc_id'], self.data['mmp']!=1))[0]
        else:
            try:
                return np.array(non_mmp_map[self.data[row]['id']])
            except:
                return np.array([])

    def __getitem__(self,key):
        return self.data[key]

    def __repr__(self):
        return "MTCatalogueTree for RSid: "+str(self.rockstar_id)+", "+str(self.nrow)+" rows."
    
class MTCatalogue:
    """
    Read in the binary file created by convertmt

    Three ways of creating this object.
    The first way is to read in the first N hosts (and their subs).
    Specify numHosts = N (default reads in everything) to read N hosts.
    Also specify indexbyrsid=True to access trees by rsid instead of mass order
    By default trees are accessed by their mass order
    
    The second way is to specify the indexfile (csv file) and a haloid.
    This will use indexfile to pinpoint the location in datafile corresponding
    to haloid and read in that host halo with all the subhalos as the catalogue.

    The third way is to use version=5, which reads the txt files and parses trees 
    using the locations.dat files to read in the haloids you give it.

    """

    def __init__(self,dir,haloids=[],version=2,numHosts=np.infty,indexbyrsid=False,verbose=False):
        self.dir = dir
        self.Trees = {} #key: rockstar halo ID; value: MT file
        self.indexbyrsid = indexbyrsid
        self.scale_list = getScaleFactors(dir+'/../') #assumes /trees is current
                                                      #folder in dir
        if version==1:
            self.fmt = "fiiiiifffffiffffffffffffiii"
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
        elif version==2:
            self.fmt = "fiiiiiifffffiffffffffffffiiiifffffffffff"
            self.fmttype = np.dtype([('scale','<f8'),('id','<i8'),('desc_id','<i8'),
                                     ('num_prog','<i8'),('pid','<i8'),('upid','<i8'),
                                     ('phantom','<i8'),
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
        elif version==3:
            self.fmt = "fiiiiiifffffiffffffffffffiiiifffffffffff"
            self.fmttype = np.dtype([('scale','<f8'),('id','<i8'),('desc_id','<i8'),
                                     ('num_prog','<i8'),('pid','<i8'),('upid','<i8'),
                                     ('phantom','<i8'),
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
                                     ('b_to_a(500c)','<f8'),('c_to_a(500c)','<f8'),
                                     ('A[x](500c)','<f8'),('A[y](500c)','<f8'),('A[z](500c)','<f8'),
                                     ('T/|U|','<f8')])
        elif version==4:
            self.fmt = "fiiiiiifffffiffffffffffffiiiifffffffffffi"
            self.fmttype = np.dtype([('scale','<f8'),('id','<i8'),('desc_id','<i8'),
                                     ('num_prog','<i8'),('pid','<i8'),('upid','<i8'),
                                     ('phantom','<i8'),
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
                                     ('b_to_a(500c)','<f8'),('c_to_a(500c)','<f8'),
                                     ('A[x](500c)','<f8'),('A[y](500c)','<f8'),('A[z](500c)','<f8'),
                                     ('T/|U|','<f8'),('snap','<i8')])
        elif version==5:
            self._readtxt(dir,haloids,indexbyrsid,verbose)
            return #don't do the rest of init since it's already done!
        else:
            print "ERROR, version must be 1, 2, 3, 4, or 5"
            sys.exit()
        self.fmtsize = struct.calcsize(self.fmt)
        ## NOTE: the reason I have separated the while loops is to speed up
        ## the reading (don't have to call if every single time)
        f = open(dir+"/tree.bin",'rb')
        if haloids==[]:
            if verbose: print "Reading whole catalogue"
            start = time.time()
            tag = f.read(8)
            nhosts = 0
            if indexbyrsid: #index by rsid
                while tag != '' and nhosts <= numHosts:
                    halotype,nrow = struct.unpack("ii",tag)
                    if halotype == 0:
                        nhosts+=1
                        if  nhosts>numHosts:
                            continue
                    thistree = MTCatalogueTree(f=f,scale_list=self.scale_list,halotype=halotype,nrow=nrow,fmt=self.fmt,fmttype=self.fmttype)
                    rsid = thistree.rockstar_id ##use rsid
                    self.Trees[rsid] = thistree 
                    tag = f.read(8)
            else: #index by mass order
                counter = 0
                self.HostLocs = []
                while tag != '' and nhosts <= numHosts:
                    halotype,nrow = struct.unpack("ii",tag)
                    if halotype == 0:
                        nhosts+=1
                        if nhosts>numHosts:
                            continue
                        self.HostLocs.append(counter)
                    thistree = MTCatalogueTree(f=f,scale_list=self.scale_list,halotype=halotype,nrow=nrow,fmt=self.fmt,fmttype=self.fmttype)
                    self.Trees[counter] = thistree ##use counter instead of rsid
                    tag = f.read(8)
                    counter+=1
            if verbose: print "Time to finish reading:",time.time()-start
        else:
            #Read just the host halos and their subs
            reader = csv.reader(open(dir+"/treeindex.csv",'r'))
            index = dict(x for x in reader)
            #print "Reading these IDs:",haloids
            if numHosts!=np.infty: print "  Warning: ignoring numHosts variable"
            try:
                file_locs = [int(index[str(x)]) for x in haloids] #raises KeyError if problem
                #Read in tree for host halos
                if ~indexbyrsid:
                    counter=0
                    self.HostLocs = []
                for file_loc,haloid in itertools.izip(file_locs,haloids):
#                    print haloid
                    f.seek(file_loc)
                    tag = f.read(8)
                    halotype,nrow = struct.unpack("ii",tag)
                    if halotype != 0:
                        print "halotype != 0, instead",halotype
                        raise ValueError
                    hosttree = MTCatalogueTree(f=f,scale_list=self.scale_list,halotype=halotype,nrow=nrow,fmt=self.fmt,fmttype=self.fmttype)
                    rsid = hosttree.rockstar_id
                    if rsid != haloid:
                        print "rsid != haloid"
                        print "rsid",rsid,"; haloid",haloid
                        raise ValueError
                    if indexbyrsid:
                        self.Trees[haloid]=hosttree
                    else:
                        self.Trees[counter]=hosttree
                        self.HostLocs.append(counter)
                        counter+=1
                    tag = f.read(8)
                    if len(tag)==0: continue #last host halo in the list, no subs
                    halotype,nrow = struct.unpack("ii",tag)
                    #Read in trees for all subhalos
                    if indexbyrsid: #index by rsid                        
                        while halotype==1:
                            thistree = MTCatalogueTree(f=f,scale_list=self.scale_list,halotype=halotype,nrow=nrow,fmt=self.fmt,fmttype=self.fmttype)
                            rsid = thistree.rockstar_id ##use rsid
                            self.Trees[rsid] = thistree
                            tag = f.read(8)
                            halotype,nrow = struct.unpack("ii",tag)
                    else: #index by mass order
                        while halotype==1:
                            thistree = MTCatalogueTree(f=f,scale_list=self.scale_list,halotype=halotype,nrow=nrow,fmt=self.fmt,fmttype=self.fmttype)
                            self.Trees[counter] = thistree ##use counter instead of rsid
                            counter+=1
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

    def _readtxt(self,tpath,haloids,indexbyrsid,verbose):
        """
        Used with version=5 to read in text files
        """
        assert len(haloids)>0,'version=5 Requires nonzero number of haloids to create catalogue'
        if verbose: print "Using text files to read %i halos" % len(haloids)
        if not os.path.exists(tpath+'/origtreeindex.dat'):
            self._create_origtreeindex(tpath)
        origtreetab = asciitable.read(tpath+'/origtreeindex.dat')
        self.origtreeindex = dict(zip(origtreetab['origid'],origtreetab['treeid']))
        loctab = asciitable.read(tpath+'/locations.dat')
        self.locations = dict(zip(loctab['TreeRootID'],zip(loctab['FileID'],loctab['Offset'])))
        
        # create allfilenames based on consistenttrees mapping
        # open files beforehand so you don't have to open a file for every tree
        numfiles = len(glob.glob(tpath+'/tree_*.dat'))
        box_divisions = numfiles**(1./3)
        assert box_divisions == int(box_divisions)
        box_divisions = int(box_divisions)
        allfilenames = [tpath+'/tree_%i_%i_%i.dat' % (i/(box_divisions*box_divisions),
                                                      (i % (box_divisions*box_divisions))/box_divisions,
                                                      i % box_divisions) for i in range(numfiles)]
        filelist = [open(f,'r') for f in allfilenames]

        ## TODO make fmttype and headers
        header = ['scale','id','desc_scale','desc_id','num_prog','pid','upid','desc_pid','phantom',
                  'sam_mvir','mvir','rvir','rs','vrms','mmp','scale_of_last_MM','vmax','posX','posY','posZ',
                  'pecVX','pecVY','pecVZ','Jx','Jy','Jz','spin','bfid','dfid','treerootid','origid',
                  'snapnum','nextcoprog_dfid','lastprog_dfid','rs_klypin',
                  'm200c_all','m200b','m200c','m500c','m2500c','xoff','voff','spin_bullock',
                  'b_to_a','c_to_a','A[x]','A[y]','A[z]',
                  'b_to_a(500c)','c_to_a(500c)','A[x](500c)','A[y](500c)','A[z](500c)',
                  'T/|U|','m_pe_b','m_pe_d','halfmassrad']
        types = ['<f8','<i8','<f8','<i8','<i8','<i8','<i8','<i8','<i8',
                 '<f8','<f8','<f8','<f8','<f8','<i8','<f8','<f8','<f8','<f8','<f8',
                 '<f8','<f8','<f8','<f8','<f8','<f8','<f8',
                 '<i8','<i8','<i8','<i8','<i8','<i8','<i8','<f8',
                 '<f8','<f8','<f8','<f8','<f8','<f8','<f8','<f8',
                 '<f8','<f8','<f8','<f8','<f8','<f8','<f8','<f8','<f8','<f8',
                 '<f8','<f8','<f8','<f8']
        self.fmttype = np.dtype({'names':header,'formats':types})
        line = filelist[0].readline()
        # kludge for the new last_mainleaf_dfid
        if len(line.split())==len(header)+1:
            header = ['scale','id','desc_scale','desc_id','num_prog','pid','upid','desc_pid','phantom',
                      'sam_mvir','mvir','rvir','rs','vrms','mmp','scale_of_last_MM','vmax','posX','posY','posZ',
                      'pecVX','pecVY','pecVZ','Jx','Jy','Jz','spin','bfid','dfid','treerootid','origid',
                      'snapnum','nextcoprog_dfid','lastprog_dfid','lastmainleaf_dfid','rs_klypin',
                      'm200c_all','m200b','m200c','m500c','m2500c','xoff','voff','spin_bullock',
                      'b_to_a','c_to_a','A[x]','A[y]','A[z]',
                      'b_to_a(500c)','c_to_a(500c)','A[x](500c)','A[y](500c)','A[z](500c)',
                      'T/|U|','m_pe_b','m_pe_d','halfmassrad']
            types = ['<f8','<i8','<f8','<i8','<i8','<i8','<i8','<i8','<i8',
                     '<f8','<f8','<f8','<f8','<f8','<i8','<f8','<f8','<f8','<f8','<f8',
                     '<f8','<f8','<f8','<f8','<f8','<f8','<f8',
                     '<i8','<i8','<i8','<i8','<i8','<i8','<i8','<i8','<f8',
                     '<f8','<f8','<f8','<f8','<f8','<f8','<f8','<f8',
                     '<f8','<f8','<f8','<f8','<f8','<f8','<f8','<f8','<f8','<f8',
                     '<f8','<f8','<f8','<f8']
            self.fmttype = np.dtype({'names':header,'formats':types})
        assert len(line.split())==len(header),'number of columns in file {} is not same as in MTCatalogue {}'.format(len(line.split()),len(header))

        counter = 0
        for haloid in haloids:
            treeid = self.origtreeindex[haloid]
            fileid,loc = self.locations[treeid]
            f = filelist[fileid]
            f.seek(loc)
            thistree = MTCatalogueTree(f=f,scale_list=self.scale_list,fmttype=self.fmttype,readtxt=True,txtheader=header)
            if indexbyrsid:
                rsid = thistree.rockstar_id
                self.Trees[rsid] = thistree
            else:
                self.Trees[counter] = thistree
                counter+=1
        for f in filelist:
            f.close()

    def _create_origtreeindex(self,tpath,stoponrooted=True):
        start = time.time()
        print "creating origtreeindex (may take some time)..."
        outfile = tpath+'/origtreeindex.dat'
        infiles = glob.glob(tpath+'/tree_*.dat')
        origids = []
        treeids = []
        skipcount=len(infiles)
        for filename in infiles:
            treecount=0
            with open(filename,'r') as f:
                line = f.readline()
                while True:
                    if line=="": break
                    elif line[0:5]=="#tree":
                        line = f.readline()
                        s = line.split()
                        if float(s[0]) != 1.0: 
                            if stoponrooted: raise IOError("Tree not rooted at a=1")
                            else: continue
                        # get ids
                        origids.append(int(s[30])); treeids.append(int(s[1]))
                        treecount += 1
                        # skip to next tree
                        while line != "" and line[0:5] != "#tree":
                            line = f.readline()
                    elif line[0]=="#": line = f.readline()
                    elif skipcount>0: skipcount -= 1; numtrees = int(line); line = f.readline()
                    else: raise IOError("Didn't read the tree files properly, {}".format(line))
                assert treecount == numtrees
        asciitable.write({'origid':origids,'treeid':treeids},outfile,names=['origid','treeid'])
        print "done! Time: {0:.1f}".format(time.time()-start)

    def getSubTrees(self,treenum):
        """
        returns list of indices of subhalo trees of tree given by treenum
        """
        if self.indexbyrsid:
            # TODO
            print 'getSubTrees not supported for indexbyrsid yet'
            return 0
        else:
            if self.Trees[treenum].halotype !=0:
                print 'ERROR: input halo not a host'
                return []
            row = treenum+1
            while row<len(self.Trees) and self.Trees[row].halotype == 1:
                row+=1
            return np.arange(treenum+1, row)

    def __getitem__(self,key):
        return self.Trees[key]

    def __repr__(self):
        out = "MTCatalogue from directory: "+self.dir+"\n"
        out = out + "Number of trees: "+str(len(self.Trees))+"\n"
        if self.indexbyrsid:
            return out+"Indexed by Rockstar halo ids"
        else:
            return out+"Indexed by sequential integers in order of mass"

    def __iter__(self):
        """
        Allows loops such as this..
        for mt in mtc: print mt
        """
        return self.Trees.itervalues()
