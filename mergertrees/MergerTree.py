"""
Merger Tree object for data produced from Rockstar Merger Tree code
@author: Greg Dooley
"""

import numpy as np

class MergerTree:
    """
    Merger Tree data structure. Just a list of all halos in tree.
    Links to other halos contained in the MTHalo objects themselves.
    @param file: full path of ascii output file from Rockstar Consistent Tree code
    @param num: ID of tree given in ascii output file
    """
    def __init__(self, file, num):
        self.file = file #: full path to file w/ information
        self.haloList = [] #: list of all halos in tree
        self.progenitors = [] #: list of all progenitors of tree
        self.tree_num = num #: ID of tree given in ascii output file
        
    
