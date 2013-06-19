"""
Read Merger Tree data from ascii file and produce merger trees
@author: Greg Dooley
"""
import numpy as np
import MTHalo as MTH
import MergerTree as MT
import os
import operator


class MTCatalogue:
    """
    Create a Merger Tree Catalogue for Rockstar output
    @param file: ascii file from Rockstar Consistent Tree code ex: \"/home/gdooley/Rockstar-0.99/Output2\tree\tree_0_0_0.dat\"
    @param numTrees: Specifies number of trees to read. Leave blank to read in all trees.
    """
    
    def __init__(self, file, numTrees=-1):
        self.file = file #: input file name
        self.Trees = [] #: List of all Merger Trees in data file
        f = open(file, 'r')

        # Read Merger Tree information
        i = 0
        line = f.readline()
        while line != '' and (i < numTrees or numTrees == -1):
            if line[0:5] == "#tree":
                i +=1
                tree = MT.MergerTree(file, int(line[6::]))
                line = f.readline()
                index = 0
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
                self.Trees.append(tree)
                if i%1000==0:
                    print 'read %d trees'%i
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
            haloList[-1].desc = haloList[index] #or just index?
            if haloList[-1].mmp == 1:
                haloList[index].parent = haloList[-1]
            haloList[index].fullParents.append(haloList[-1])
        else:
            index -= 1
    return haloList


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
