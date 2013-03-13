import numpy as np
import os

def readParents(dir, file, numhalos):
    """
    @ param dir: directory of file (string)
    @ param file: name of .list file produced by ./find_parents
    @ numhalos: number of halos in the list.
    @ return:  
    """
    file_name = dir+'/'+file
    if(not os.path.exists(file_name)):
		raise IOError('file not found: '+file_name)
    f = open(file_name)
    data = np.zeros((numhalos,2))
    f.readline() # read the header
    i = 0
    for line in f:
        a = line.split()
        try:
            data[i][0] = int(a[0])
            data[i][1] = int(a[-1])
        except:
            print "ERROR: AHHHH"
            print a
            print i
            print "numhalos: ",numhalos
            print data[i][0], data[i][1]
            print data[i-1][0], data[i-1][1]
        i += 1
    if i != numhalos:
        print 'numhalos does not match number of lines in readParentsList.py'
    #sortedIndices = data[:,0].argsort()
    #data = data[sortedIndices]
    return data