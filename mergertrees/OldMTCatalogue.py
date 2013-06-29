"""
Read Merger Tree data from ascii file and produce merger trees
@author: Greg Dooley

6/28/2013 - added verbose option to print to stdout. default prints nothing.
"""
import numpy as np
import MTHalo as MTH
import MergerTree as MT
import os
import operator
import time

class OldMTCatalogue:
    """
    Create a Merger Tree Catalogue for Rockstar output
    @param file: ascii file from Rockstar Consistent Tree code ex: \"/home/gdooley/Rockstar-0.99/Output2\tree\tree_0_0_0.dat\"
    @param numTrees: Specifies number of trees to read. Leave blank to read in all trees.


    NOTE** In this version, halo.fullParents gives all parents EXCEPT most massive parent. Should be renamed halo.otherParents, but MTHalo is read only for me (Greg).
    """
    
    def __init__(self, file, numTrees=-1, verbose=False):
        self.file = file #: input file name
        self.Trees = [] #: List of all Merger Trees in data file
        start = time.time()
        if verbose==True:
            print 'finding z=0 host-sub match ups'
        host2sub = storeSubInfo(file)
        if verbose==True:
            print 'done with match ups. Time = ', time.time()-start
        
        f = open(file, 'r')
        g = open(file, 'r')
        # Read Merger Tree information
        i = 0 #used only for printing progress updates
        line = f.readline()
        while line != '' and (i < numTrees or numTrees == -1):
            if line[0:5] == "#tree":
                # found a new tree. Build it.
                i +=1
                tree = MT.MergerTree(file, int(line[6::]))
                line = f.readline()
                index = 0 # line number within a tree
                while line != '' and line[0:5] != '#tree':
                    #Create Halo for input line
                    halo = MTH.MTHalo(line)
                    tree.haloList.append(halo)
                    # If it has no parents, add to progenitors list
                    if halo.num_prog == 0:
                        tree.progenitors.append(halo)
                    # Update parent-descendant links
                    updateLinks(tree.haloList, index)
                    line = f.readline()
                    index += 1
                    if index%10000==0:
                        if verbose==True:
                            print 'finished item # ', index, ' in tree ', i
                # Now find all subhalos of z=0 root halo if its a strict host.
                if tree.haloList[0].pid == -1:
                    #print 'about to find subhalos'
                    updateSubhalos(tree.haloList[0],file,host2sub)
                self.Trees.append(tree)
                if i < 50 or i%100 == 0:
                    if verbose==True:
                        print 'finished tree ', i, ' at ', time.time()-start, ' seconds'
            else:
                line = f.readline()
        f.close()

         

def updateLinks(haloList, index):
    """
    Assign parent and descendant links for last halo in haloList
    @param haloList: current list of halos in MT.
    @param index: index location in MT haloList to start search for descendant
    @return 
    """
    found = False
    while not found and index != -1:
        if haloList[-1].descid == haloList[index].ID:
            found = True
            haloList[-1].desc = haloList[index]
            if haloList[-1].mmp == 1:
                haloList[index].parent = haloList[-1]
            else:
                haloList[index].fullParents.append(haloList[-1])
        else:
            index -= 1
    return haloList

## Want to re-write the updateSubhalos routine.
## Need to read in entire file at once, store links
def storeSubInfo(file):
    #search linearly, get line number and host halo ID for all a=0 halos.
    #then sort host halo IDs, and for each has a corresponding list of line
    # numbers of its hosts.

    #file = '/spacebase/data/AnnaGroup/cosm1/rockstar/trees/tree_0_0_0.dat'
    lines = []
    hostIDs = []
    f = open(file,'r')
    line = f.readline()
    while line!='':
        if line[0:5] == "#tree":
            num = int(line[6::])
            loc = f.tell()
            line = f.readline()
            sub = MTH.MTHalo(line)
            if sub.pid != -1:
                #print sub.pid, sub.scale
                hostIDs.append(sub.pid)
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
    #print lines[0:10]
    #print hostIDs[0:10]
    return host2sub


def updateSubhalos(host,file, host2sub):
    """
    Use the host2sub dictionary to find subhalos of host
    in Merger Tree text output.
    Run tree creation on all subhalos to get their information.

    @param host: host halo (MTHalo object)
    @param g: file descriptor of already opened file.
    """
    if not (host.ID in host2sub):
        return
    g = open(file,'r')
    for posn in host2sub[host.ID]:
        g.seek(posn)
        line = g.readline()
        sub = MTH.MTHalo(line)
        if sub.pid != host.ID:
            print 'WARNING: ERROR: halo not sub of host! Proceeding anyway'
        tree = MT.MergerTree(file,sub.ID)
        tree.haloList.append(sub)
        if sub.num_prog==0:
            tree.progenitors.append(sub)
        # Now deal with all other halos in the tree
        index = 1
        line = g.readline()
        while line !='' and line[0:5] != '#tree':
            halo = MTH.MTHalo(line)
            tree.haloList.append(halo)
            if halo.num_prog ==0:
                tree.progenitors.append(halo)
            updateLinks(tree.haloList, index)
            line = g.readline()
            index += 1
        host.subhalos.append(sub)
    g.close()

def updateSubhalos_old(host, file):
    """
    Search entire Merger Tree text output for subhalos of host.
    Add subhalos MTHalo subhalos list.
    Run tree creation on all subhalos to get their information.
    
    @param host:  host halo (MTHalo object)
    @param file: name of file to be opened
    Assumes that all subhalos come listed after host halo.
    """
    f = open(file, 'r')
    line = f.readline()
    i = 0
    while line != '':
        if line[0:5] == "#tree":
            #if i%10000 == 0:
                #print 'subhalo finder scanned ', i, ' trees'
            i+=1
            num = int(line[6::])
            # Deal with a=0 halo independently
            line = f.readline()
            sub = MTH.MTHalo(line)
            if sub.pid == host.ID: # not upid. only subhalos, not subsub etc.
                #build tree, add to subhalo list of host
                tree = MT.MergerTree(file, num)
                tree.haloList.append(sub)
                if sub.num_prog ==0:
                    tree.progenitors.append(sub)

                # Now deal with all other halos in the tree
                index = 1
                line = f.readline()
                while line !='' and line[0:5] != '#tree':
                    halo = MTH.MTHalo(line)
                    tree.haloList.append(halo)
                    if halo.num_prog ==0:
                        tree.progenitors.append(halo)
                    updateLinks(tree.haloList, index)
                    line = f.readline()
                    index +=1
                # add a=1 subhalo to subhalo list of host (maybe should add tree?)
                host.subhalos.append(sub)
            else:
                line = f.readline()
        else:
            line = f.readline()
    f.close()

def particleLastHalo(particle, halo):
    """
    Find which halo particle was in in the previous timestep.
    @param particle: gadget ID of particle to trace
    @param halo: MTHalo where particle is found
    @return: halo in previous timestep hosting particle. -1 if particle not contained in any parent halo
    """
    for i in range(0, len(halo.fullParents)):
        if halo.fullParents[i].containsParticle(particle):
            return halo.fullParents[i]
    return -1
