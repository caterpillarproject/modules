import numpy as np
import operator

op = {}
op['lt'] = operator.lt
op['le'] = operator.le
op['eq'] = operator.eq
op['ne'] = operator.ne
op['ge'] = operator.ge
op['gt'] = operator.gt

class read:
    '''
    A class for reading in data from the Millennium database

    Makes use of dictionaries so that an arbitrary number of columns with 
    arbitrary names can be read in. Each set of data can then be accessed with
    a call class_name['column_name']
    To find what data is available, use the function keys()
    '''

    def __init__(self):
        self._registry={} # names of available data types
        self._data={}     # the data
        self._N=0         # length of data arrays

    def keys(self):
        '''Returns a list of available data types'''
        return self._registry.keys()

    def readfile(self,filename,first=0,last=0):
        '''
        Read in the data
        
        first and last specify which columns to take data from (using pythonic list notation)
        If last==-1 all columns after first are used
        '''
        column_headings_read = False
        listin = []
        for line in open(filename,'r'):
            li=line.strip()
        
            if not li.startswith("#"):
                line = line.partition('#')[0] 

                # first line of data after commented section is a list of 
                # column headings. This is used for the dictionary names
                if not column_headings_read:
                    column_headings_read = True
                    if last==0: last=len(line.rstrip().split(','))
                    column_headings = line.rstrip().split(',')[first:last]
                    continue
                   
                listin.append([token for token in line.split(',')][first:last])

        if not column_headings_read:
            raise RuntimeError, 'Column headings not read from file'

        if len(listin)==0:
            print 'Warning: No data in file',filename
            return False

        listin = np.asarray(listin)  

        for i in xrange(len(column_headings)):
            self._registry[column_headings[i]] = {}
            if 'ID' in column_headings[i] or 'Id' in column_headings[i] or column_headings[i]=='np' or 'snapnum' in column_headings[i].lower():
              self._data[column_headings[i]] = np.array(listin[:,i]).astype(np.int64)
            else:
              self._data[column_headings[i]] = np.array(listin[:,i]).astype(np.float64)

        self._N = len(self._data[self._registry.keys()[0]])

        return True

    def where(self, name, cmpr, value): 
        '''
        Finds all items where classobjectname[name] cmpr value
        and cmpr is in the set of comparisons lt, le, eq, ne, ge, gt
        Output is a new readmodule.read() class
        '''
        try:
            idx = np.where(op[cmpr](self._data[name],value))[0] # where returns tuple
        except:
            return None
        newobject = read()
        for item in self._registry.keys():
            newobject[item] = self._data[item][idx]
        return newobject

    def getindex(self, name, value):
        '''
        Finds the indices of items where classobjectname[name]==value
        Returns an array
        '''
        try:
            idx = np.where(self._data[name]==value)[0]
        except:
            return None
        return idx

    def write(self,filename):
        f = open(filename, 'w')
        f.write( ','.join(self._registry.keys())+'\n' ) # headings
        for i in xrange(self._N):
            line = []
            for name in self._registry.keys():
                line.append(str(self._data[name][i]))
            f.write(','.join(line)+'\n')
        f.close()

    def sort(self,name,direction='asc'):
        '''
        Sort all arrays accoring to classojectname[name] in the given direction
        Possible sorting directions are 'asc'(ascending) and 'desc'(descending)
        '''
        if direction.lower()=='asc': reverse=False
        elif direction.lower()=='desc': reverse=True
        else: 
            print 'Error: sort failed, direction not recognised'
            return False

        arr = zip(self._data[name], np.arange(self._N))
        arr = sorted(arr, key=lambda arr: arr[0], reverse=reverse)
        arr, sort_key = zip(*arr)
        sort_key = np.array(sort_key)
        del arr

        for key in self._registry.keys():
            self._data[key] = self._data[key][sort_key]

        return True

    #
    # private functions
    #

    def __setitem__(self,name,item):
        self._data[name] = np.array(item)
        self._registry[name] = {}
        if self._N==0:
            self._N = len(self._data[name])

    def __getitem__(self,name):
        if name in self._data:
            return self._data[name]
        else:
            raise KeyError, "Data type "+name+" is not valid, check keys() for defined data"


def readMBPs(filename):
  mbp = {}
  mbp['id'] = []
  mbp['x'] = []
  mbp['y'] = []
  mbp['z'] = []
  mbp['index'] = {}
  i=0
  for line in open(filename,'r'):
    li = line.split()
    mbp['x'].append(float(li[0]))
    mbp['y'].append(float(li[1]))
    mbp['z'].append(float(li[2]))
    mbp['id'].append(int(li[6]))
    mbp['index'][int(li[6])] = i
    i+=1
  mbp['id'] = np.array(mbp['id'])
  mbp['x'] = np.array(mbp['x'])
  mbp['y'] = np.array(mbp['y'])
  mbp['z'] = np.array(mbp['z'])
  return mbp
