# Python HDF5 merger tree reader 
# (requires util/hdf5lib.py)
#
# see example_X.py for usage
#
#
# Mark Vogelsberger (mvogelsb@cfa.harvard.edu)

import numpy as np
import os
import sys
import hdf5lib

mergertree_datablocks = {"Descendant":                 ["int32",   1, True],
                         "FirstProgenitor":            ["int32",   1, True],
                         "NextProgenitor":             ["int32",   1, True],
                         "FirstHaloInFOFGroup":        ["int32",   1, True],
                         "NextHaloInFOFGroup":         ["int32",   1, True],
                         "SubhaloLen":                 ["int32",   1, True],
                         "Group_M_Mean200":            ["float32", 1, True],
                         "Group_M_Crit200":            ["float32", 1, True], 
                         "Group_M_TopHat200":          ["float32", 1, True],
                         "SubhaloPos":                 ["float32", 3, True],
                         "SubhaloVel":                 ["float32", 3, True], 
                         "SubhaloVelDisp":             ["float32", 1, True],
                         "SubhaloVMax":                ["float32", 1, True], 
                         "SubhaloSpin":                ["float32", 3, True], 
                         "SubhaloIDMostBound":         ["int64",   1, True],
                         "SnapNum":                    ["int32",   1, True],
                         "FileNr":                     ["int32",   1, True], 
                         "SubhaloGrNr":                ["int32",   1, True],
                         "SubhaloNumber":              ["int32",   1, True],
                         "SubhaloSFR":                 ["float32", 1, True],
                         "SubhaloGasMetallicity":      ["float32", 1, True],
                         "SubhaloGasMetallicitySfr":   ["float32", 1, True],
                         "SubhaloStarMetallicity":     ["float32", 1, True],
                         "SubhaloOffsetType":          ["int64",   6, True],
                         "SubhaloLenType":             ["int32",   6, True],
	                 "SubhaloMassType":            ["float32", 6, True],
                         "SubhaloMassInRadType":       ["float32", 6, True],
		         "SubhaloHalfmassRadType":     ["float32", 6, True],
                         "SubhaloBHMass":              ["float32", 1, True],
                         "SubhaloBHMdot":              ["float32", 1, True], 
                         "SubhaloSFRinRad":            ["float32", 1, True],
                         "SubhaloStellarPhotometrics": ["float32", 8, True]}


class merger_tree:
        def __init__(self, basedir, skipfac, snapnum, filenum = 0, tree_start = -1, tree_num = -1, keysel = None):

		self.filebase = basedir + "trees_sf"+str(skipfac)+"_"+str(snapnum).zfill(3)
		self.basedir = basedir
		self.filenum = filenum
		filename = self.filebase + "." + str(filenum) + ".hdf5"
		f=hdf5lib.OpenFile(filename)
		self.NtreesPerFile = hdf5lib.GetAttr(f, "Header", "NtreesPerFile") 
		self.NumberOfOutputFiles = hdf5lib.GetAttr(f, "Header", "NumberOfOutputFiles") 
		self.ParticleMass = hdf5lib.GetAttr(f, "Header", "ParticleMass") 
		if (self.ParticleMass == 0):
			print "WARNING: ParticleMass = 0, needed for merger rate calculation"
		self.TreeNHalos = hdf5lib.GetData(f, "Header/TreeNHalos")[:] 
		self.TotNsubhalos = hdf5lib.GetData(f, "Header/TotNsubhalos")[:] 
		self.Redshifts = hdf5lib.GetData(f, "Header/Redshifts")[:] 
		if (tree_start == -1 ) | (tree_num == -1):
			tree_start = 0
			tree_num = self.NtreesPerFile
		self.trees = np.empty(tree_num - tree_start, dtype='object')
		self.tree_start = tree_start
		self.tree_num = tree_num
		for ntree in range(tree_start, tree_start + tree_num):
			list = []
			if (keysel==None):
				for datablock in mergertree_datablocks.keys():
					data = hdf5lib.GetData(f, "Tree"+str(ntree)+"/"+datablock)[:] 
					list.append((datablock,data))
			else:
				for datablock in keysel:
					if hdf5lib.Contains(f, "Tree"+str(ntree), datablock):
						data = hdf5lib.GetData(f, "Tree"+str(ntree)+"/"+datablock)[:]
						list.append((datablock,data))
			self.trees[ntree - tree_start] = dict(list)
		f.close()

        def __count_unique(self, keys):
                uniq_keys = np.unique(keys)
                bins = uniq_keys.searchsorted(keys)
                return uniq_keys, np.bincount(bins)

	def getNumberOfMergers(self, snapnum, bins_halo = 10, bins_ratio = 10, halo_min = 8, halo_max = 13, ratio_min = 0, ratio_max = 1):
		htot = np.zeros([bins_halo, bins_ratio])
		xtot = 0
		ytot = 0
		for ntree in range(0,self.tree_num):
			idx = (self.trees[ntree]["SnapNum"][:] == snapnum) & (self.trees[ntree]["Descendant"] >= 0)
			if (idx.any()):
				halos = np.arange(0, self.TreeNHalos[ntree])[idx]
				descs = self.trees[ntree]["Descendant"][idx]
				d_tmp, n_tmp = self.__count_unique(descs)
				merger_descs = d_tmp[n_tmp > 1]
				if (len(descs) > 0):
					for md in merger_descs:
						len_desc = self.trees[ntree]["Len"][md]
						len_halos = self.trees[ntree]["Len"][md == descs] 
						ratio = 1.0 * len_halos / len_desc		
						x = np.log10(len_halos * self.ParticleMass * 1e10)
						y = ratio
						h, xtot, ytot = np.histogram2d(x, y, bins=(bins_halo, bins_ratio), range = [[halo_min,halo_max], [ratio_min, ratio_max]])
						htot += h
		return [xtot, ytot, htot]

	def getAllProgenitors(self, ntree, nhalo):
		list_next = []
		list_first = []
		list_first.append(self.trees[ntree]["FirstProgenitor"][nhalo])
		while (len(list_first) > 0):
			next = list_first.pop()
			while (next >= 0):
				list_next.append(next)
				new_next = self.trees[ntree]["NextProgenitor"][next]
				new_first = self.trees[ntree]["FirstProgenitor"][next]
				list_first.append(new_first)
				next = new_next
		return list_next

	def getProgenitors(self, ntree, nhalo):
		list = []
		next = self.trees[ntree]["FirstProgenitor"][nhalo]
		while (next >= 0):
			list.append(next)
			next = self.trees[ntree]["NextProgenitor"][next]
		return list

        def getFirstProgenitors(self, ntree, nhalo):
                list = []
                next = nhalo 
                while (next >= 0):
                        list.append(next)
                        next = self.trees[ntree]["FirstProgenitor"][next]
                return list

        def getHalosInFOFGroup(self, ntree, nhalo):
                list = []
                next = self.trees[ntree]["FirstHaloInFOFGroup"][nhalo]
                while (next >= 0):
                        list.append(next)
                        next =self.trees[ntree]["NextHaloInFOFGroup"][next]
                return list

	def getDescendants(self, ntree, nhalo):
		list = []
		next = self.trees[ntree]["Descendant"][nhalo]
		while (next >=0 ):
			list.append(next)
			next = self.trees[ntree]["Descendant"][next]
		return list

	def constructSubhaloLookup(self, snapnum):
		self.SubhaloLookupTable = np.zeros([self.TotNsubhalos[snapnum],3], dtype='int32') - 1
		for ntree in range(0,self.tree_num):
			idx = (self.trees[ntree]["SnapNum"][:] == snapnum) 
			halos = np.arange(0, self.TreeNHalos[ntree], dtype='int32')[idx]
			subnums = self.trees[ntree]["SubhaloNumber"][idx]
			self.SubhaloLookupTable[subnums,0] = self.filenum
			self.SubhaloLookupTable[subnums,1] = ntree
			self.SubhaloLookupTable[subnums,2] = halos
		f=open(self.basedir+"/SubhaloLookup_"+str(snapnum).zfill(3)+"."+str(self.filenum)+".dat","wb")
		self.SubhaloLookupTable.astype("int32").tofile(f)
		f.close()
	
	def combineSubhaloLookup(self, snapnum):
		self.SubhaloLookupTable = np.zeros([self.TotNsubhalos[snapnum],3], dtype='int32') - 1
		for filenum in range(0, self.NumberOfOutputFiles):
			f=open(self.basedir+"/SubhaloLookup_"+str(snapnum).zfill(3)+"."+str(filenum)+".dat","rb")			
			tmp = np.fromfile(f, dtype="int32", count=3 * self.TotNsubhalos[snapnum]).reshape([self.TotNsubhalos[snapnum],3])
			f.close()
			idx = tmp != -1
			self.SubhaloLookupTable[idx] = tmp[idx]
	
	def saveSubhaloLookup(self, base, snapnum):
		f=open(base+"/SubhaloLookup_"+str(snapnum).zfill(3)+".dat","wb")
		self.SubhaloLookupTable.astype("int32").tofile(f)
		f.close()

	def loadSubhaloLookup(self, base, snapnum):
		f=open(base+"/SubhaloLookup_"+str(snapnum).zfill(3)+".dat","rb")
		self.SubhaloLookupTable = np.fromfile(f, dtype="int32", count=3 * self.TotNsubhalos[snapnum]).reshape([self.TotNsubhalos[snapnum],3])
		f.close()
	
	def SubhaloLookupTable(self):
		return self.SubhaloLookupTable
	
	def lookupSubhalo(self, subhalo_num):
		filenum = self.SubhaloLookupTable[subhalo_num,0]
		ntree = self.SubhaloLookupTable[subhalo_num,1]
		nhalo = self.SubhaloLookupTable[subhalo_num,2]
		return [filenum, ntree, nhalo]
