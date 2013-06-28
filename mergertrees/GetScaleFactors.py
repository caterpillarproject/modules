import numpy as np
import os
import sys

"""
produce list of scale factors of halos
"""

def getScaleFactors(path):
    snap_num = 7
    file = path + '/' + 'halos_' + str(snap_num) + '/' + 'halos_' + str(snap_num) + ".0.bin"
    f = open(file)
    f.close()
    
    if (not os.path.exists(file)):
            print "ERROR: file not found", file
            sys.exit()

    scale_list = []
    while os.path.exists(file):
        f = open(file)
        magic = np.fromfile(f, np.uint64, count = 1)[0]
        snap = np.fromfile(f, np.int64, count = 1)[0]
        chunk = np.fromfile(f, np.int64, count = 1)[0]
        scale = np.fromfile(f, np.float32, count = 1)[0]
        scale_list.append(scale)
        snap_num+=1
        file = path + '/' + 'halos_' + str(snap_num) + '/' + 'halos_' + str(snap_num) + ".0.bin"

    return scale_list


#path = '/spacebase/data/AnnaGroup/cosm1/rockstar'
#scales = getScaleFactors(path)
#print scales
