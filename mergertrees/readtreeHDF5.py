import numpy as np
import tables
import sys
import os

"""
Simple Python script for reading merger tree HDF5 files
in "database mode," which is optimized for extracting
quantities along the main branch of a subhalo. This code
requires PyTables (http://www.pytables.org).

There is also a "linked-list mode" which allows for more flexibility.
However, the linked-list approach is only feasible with C++,
at least for the larger simulations.

Vicente Rodriguez-Gomez (vrodriguez-gomez@cfa.harvard.edu)

--------------- USAGE EXAMPLE: PRINT STELLAR MASS HISTORY ---------------

import readtreeHDF5
treedir = '/n/ghernquist/vrodrigu/MergerTrees/output/Subhalos/Illustris/L75n1820FP'
tree = readtreeHDF5.TreeDB(treedir)
snapnum = 135; subfind_id = 0
branch = tree.get_main_branch(snapnum, subfind_id)
print branch.SubhaloMassType[:, 4]

-------------------------------------------------------------------------

"""

##########################################################################
# --------------------- General purpose functions ------------------------
##########################################################################

def _get_column_info_extended(treefilename):
    """
    INPUT:  Filename of the HDF5 merger tree.
    OUTPUT: Dictionary with the following information:
            Key: field_name
            Value: [dimension, numpy.dtype]
    """
    # Read column names and numpy dtypes from HDF5 file
    f = tables.openFile(treefilename)
    coldtypes = f.root.Tree.coldtypes
    f.close()

    # Create dictionary
    columns = {}
    for k in coldtypes:
        if len(k.split('_')) > 1:  # e.g. 'SubhaloVel_0', 'SubhaloVel_1', ...
            field_name = k.split('_')[0]
            cur_dim = int(k.split('_')[1])
            if columns.has_key(field_name):
                if cur_dim+1 > columns[field_name][0]:
                    columns[field_name] = [cur_dim+1, coldtypes[k]]
            else:
                columns[field_name] = [cur_dim+1, coldtypes[k]]
        else:
            field_name = k
            columns[field_name] = [1, coldtypes[k]]

    return columns

def _get_column_info_minimal():
    """
    Minimal column information. This should be the same
    for all simulations.
    """
    columns_minimal = ["SubhaloID",
                       "SubhaloIDRaw",
                       "LastProgenitorID",
                       "MainLeafProgenitorID",
                       "RootDescendantID",
                       "TreeID",
                       "SnapNum",
                       "FirstProgenitorID",
                       "NextProgenitorID",
                       "DescendantID",
                       "FirstSubhaloInFOFGroupID",
                       "NextSubhaloInFOFGroupID",
                       "NumParticles",
                       "Mass",
                       "MassHistory",
                       "SubfindID"]

    return columns_minimal


##########################################################################
# --------------------------- DATABASE MODE ------------------------------
##########################################################################

class _Subset(object):
    """
    Used to represent any subset of _AdjacentRows.
    Initialized with an integer or boolean array.
    """
    def __init__(self, adj_rows, indices):

        # Minimal columns
        for field_name in adj_rows._columns_minimal:
            setattr(self, field_name, getattr(adj_rows, field_name)[indices])
        
        # Other columns
        if adj_rows._keysel == None:
            for field_name in adj_rows._columns_extended.keys():
                setattr(self, field_name, getattr(adj_rows, field_name)[indices])
        else:
            for field_name in adj_rows._keysel:
                setattr(self, field_name, getattr(adj_rows, field_name)[indices])

class _AdjacentRows(object):
    """
    Used by the TreeDB class. Consists of 
    a set of adjacent rows from the merger tree file.
    Since subhalo IDs are assigned in a depth-first fashion,
    a "chunk" of adjacent rows can represent, e.g., a main branch
    or a "subtree."
    For a given file number and range of row numbers,
    create arrays containing information from the merger tree
    for the specified rows.
    """
    def __init__(self, treedir, name, columns_minimal, columns_extended,
            row_start, row_end, row_original=None, filenum=-1, keysel=None):

        # Private attributes
        self._columns_minimal = columns_minimal
        self._columns_extended = columns_extended
        self._keysel = keysel
        
        # Public attributes
        self.row_start = row_start
        self.row_end = row_end
        if row_original == None:
            self.index_given_sub = 0
        else:
            self.index_given_sub = row_original - row_start

        # Only interested in these row numbers:
        self.nrows = row_end - row_start + 1
        locs = slice(row_start, row_end+1)

        # Tree filename
        if filenum == -1:
            filename = '%s/%s.hdf5' % (treedir, name)
        else:
            filename = '%s/%s.%d.hdf5' % (treedir, name, filenum)

        if keysel == None:
            # Read all columns at once
            f = tables.openFile(filename)
            data_extended = f.root.Tree[locs]
            f.close()

            # First add the "minimal" tree columns
            for field_name in columns_minimal:
                data_tmp = data_extended[field_name][:]
                setattr(self, field_name, data_tmp)

            # Then add the other columns
            for field_name, value in columns_extended.iteritems():
                field_size = value[0]
                field_dtype = value[1]
                if field_size == 1:
                    data_tmp = data_extended[field_name][:]
                elif field_size > 1:
                    data_tmp = np.empty((self.nrows, field_size), dtype=field_dtype)
                    for k in range(field_size):
                        data_tmp[:,k] = data_extended['%s_%d' % (field_name, k)][:]
                else:
                    print 'Invalid dimension:', field_size

                # Create attribute
                setattr(self, field_name, data_tmp)

        else:
            # In this case we read the columns one by one
            f = tables.openFile(filename)
            
            # First add the "minimal" tree columns
            for field_name in columns_minimal:
                data_tmp = f.root.Tree.cols._f_col(field_name)[locs]
                setattr(self, field_name, data_tmp)
            
            # Then add the other columns
            for field_name in keysel:
                field_size = columns_extended[field_name][0]
                field_dtype = columns_extended[field_name][1]
                if field_size == 1:
                    data_tmp = f.root.Tree.cols._f_col(field_name)[locs]
                elif field_size > 1:
                    data_tmp = np.empty((self.nrows, field_size), dtype=field_dtype)
                    for k in range(field_size):
                        data_tmp[:,k] = f.root.Tree.cols._f_col('%s_%d' % (field_name, k))[locs]
                else:
                    print 'Invalid dimension:', field_size

                # Create attribute
                setattr(self, field_name, data_tmp)
            
            # Close (and flush) HDF5 file
            f.close()

    def get_subset(self, indices):
        return _Subset(self, indices)

class TreeDB(object):
    """
    Python class to extract information from merger tree files
    in "database mode."
    
    --------------------------- USAGE EXAMPLE -----------------------------
    import readtreeHDF5
    treedir = '/n/ghernquist/vrodrigu/MergerTrees/output/Subhalos/Illustris/L75n1820FP'
    tree = readtreeHDF5.TreeDB(treedir)
    snapnum = 135; subfind_id = 0
    branch = tree.get_main_branch(snapnum, subfind_id)
    print branch.SubhaloMassType[:, 4]
    -----------------------------------------------------------------------
    """

    def __init__(self, treedir, name='tree_extended'):
        self._treedir = treedir
        self._name = name

        # Check that a few files exist
        for rel_path in ['/tree_index.hdf5']:
            if not os.path.exists(self._treedir + rel_path):
                print 'File not found: ' + self._treedir + rel_path
                sys.exit()

        # Get info about columns
        self._columns_extended = _get_column_info_extended(self._treedir + '/' + self._name + '.0.hdf5')
        self._columns_minimal = _get_column_info_minimal()

    def get_main_branch(self, snapnum, subfind_id, keysel=None):
        """
        For a subhalo specified by its snapshot number and Subfind ID,
        return the progenitors along its main branch, i.e. all subhalos
        with IDs between SubhaloID and MainLeafProgenitorID.
        
        Parameters
        ----------
        snapnum : int
        subfind_id : int
        keysel: list of strings or None, optional
                This argument specifies which fields from the Subfind catalog
                should be read. 
        
        """
        # Create raw ID
        raw_id = snapnum*10**12 + subfind_id
        
        # Branch is defined by the subhalos between SubhaloID and
        # MainLeafProgenitorID. Open minimal tree file to get this info.
        f = tables.openFile('%s/tree.hdf5' % (self._treedir))
        result = [[row["SubhaloID"], row["MainLeafProgenitorID"], row.nrow] for row in f.root.Tree.where(
                  """SubhaloIDRaw == %d""" % (raw_id))]
        f.close()
        subhalo_id, main_leaf_progenitor_id, rownum = result[0]
        
        # Create branch instance
        row_start = rownum
        row_end = rownum + (main_leaf_progenitor_id - subhalo_id)
        branch = _AdjacentRows(self._treedir, self._name, self._columns_minimal, self._columns_extended,
                row_start, row_end, keysel=keysel)
        
        return branch

    def get_all_progenitors(self, snapnum, subfind_id, keysel=None):
        """
        For a subhalo specified by its snapshot number and Subfind ID,
        return all the objects in the subtree which is rooted on the
        subhalo of interest, i.e. all subhalos with IDs between SubhaloID
        and LastProgenitorID. Note that this includes the given subhalo itself.
        
        Parameters
        ----------
        snapnum : int
        subfind_id : int
        keysel: list of strings or None, optional
                This argument specifies which fields from the Subfind catalog
                should be read. 
        
        """
        # Create raw ID
        raw_id = snapnum*10**12 + subfind_id
        
        # Subtree is defined by the subhalos between SubhaloID and
        # LastProgenitorID. Open minimal tree file to get this info.
        f = tables.openFile('%s/tree.hdf5' % (self._treedir))
        result = [[row["SubhaloID"], row["LastProgenitorID"], row.nrow] for row in f.root.Tree.where(
                  """SubhaloIDRaw == %d""" % (raw_id))]
        f.close()
        subhalo_id, last_progenitor_id, rownum = result[0]
        
        # Create branch instance
        row_start = rownum
        row_end = rownum + (last_progenitor_id - subhalo_id)
        subtree = _AdjacentRows(self._treedir, self._name, self._columns_minimal, self._columns_extended,
                row_start, row_end, keysel=keysel)
        
        return subtree

    def get_all_progenitors_of_root_descendant(self, snapnum, subfind_id, keysel=None):
        """
        Return the subtree rooted on the root descendant of the given subhalo,
        i.e. all subhalos with IDs between RootDescendantID and
        RootDescendant->LastProgenitorID. Note that this includes
        the given subhalo itself.
        
        Parameters
        ----------
        snapnum : int
        subfind_id : int
        keysel: list of strings or None, optional
                This argument specifies which fields from the Subfind catalog
                should be read. 
        
        """
        # Create raw ID
        raw_id = snapnum*10**12 + subfind_id
        
        # Subtree is defined by the subhalos between RootDescendantID and
        # RootDescendant->LastProgenitorID. Open minimal tree file to get this info.
        f = tables.openFile('%s/tree.hdf5' % (self._treedir))

        # First locate the given subhalo and get its root descendant ID.
        result = [[row["SubhaloID"], row["RootDescendantID"], row.nrow] for row in f.root.Tree.where(
                  """SubhaloIDRaw == %d""" % (raw_id))]
        subhalo_id, root_descendant_id, rownum = result[0]
        # We know the row number of the root descendant without searching for it
        row_start = rownum - (subhalo_id - root_descendant_id)
        row_end = f.root.Tree.cols.LastProgenitorID[row_start]
        # Close file
        f.close()
        
        # Create branch instance
        branch = _AdjacentRows(self._treedir, self._name, self._columns_minimal, self._columns_extended,
                row_start, row_end, row_original=rownum, keysel=keysel)
        
        return branch

    def get_direct_progenitors(self, snapnum, subfind_id, **kwargs):
        """
        Return the subhalos for which DescendantID corresponds to the
        current subhalo.
        """
        subtree = self.get_all_progenitors(snapnum, subfind_id, **kwargs)
        subhalo_id = subtree.SubhaloID[0]  # unique ID of given subhalo
        indices = subtree.DescendantID == subhalo_id
        return subtree.get_subset(indices)

    def get_fellow_progenitors(self, snapnum, subfind_id, **kwargs):
        """
        Return all the subhalos that will merge into the same object
        during the next snapshot (or two), i.e. those subhalos for which
        DescendantID equals DescendantID of the current subhalo.
        """
        subtree = self.get_all_progenitors_of_root_descendant(snapnum, subfind_id, **kwargs)
        desc_id = subtree.DescendantID[subtree.index_given_sub]  # unique ID of descendant
        indices = subtree.DescendantID == desc_id
        return subtree.get_subset(indices)

    def get_all_fellow_progenitors(self, snapnum, subfind_id, **kwargs):
        """
        Return all the subhalos that will merge into the same object
        at any point in the future, i.e. those subhalos for which
        RootDescendantID equals RootDescendantID of the current subhalo.
        """
        subtree = self.get_all_progenitors_of_root_descendant(snapnum, subfind_id, **kwargs)
        root_desc_id = subtree.RootDescendantID[subtree.index_given_sub]  # unique ID of root descendant
        indices = subtree.RootDescendantID == root_desc_id
        return subtree.get_subset(indices)

    def get_future_branch(self, snapnum, subfind_id, **kwargs):
        """
        Return the subhalos found in a sort of "forward" branch between
        SubhaloID and RootDescendantID. Note that these subhalos are not
        necessarily stored in adjacent rows, as is the case
        with a main branch (following FirstProgenitor links).
        """
        subtree = self.get_all_progenitors_of_root_descendant(snapnum, subfind_id, **kwargs)
        # Unfortunately, there are no shortcuts in this case and we must
        # proceed iteratively. This is almost at the limit of what one
        # can do when reading trees in "database mode."
        desc_id = subtree.DescendantID[subtree.index_given_sub]
        root_desc_id = subtree.RootDescendantID[subtree.index_given_sub]
        indices = [subtree.index_given_sub]
        while desc_id != root_desc_id:
            cur_index = np.where(subtree.SubhaloID == desc_id)[0][0]
            desc_id = subtree.DescendantID[cur_index]
            indices.append(desc_id)
        indices = indices[::-1]  # reverse
        return subtree.get_subset(indices)
        
