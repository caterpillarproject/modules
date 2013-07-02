import numpy as np
import os
import sys
import scipy.interpolate as interpolate

"""
produce list of scale factors of halos
"""

def getScaleFactors(path,minsnap=0):
    snap_num = minsnap
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

class getsnap:
    """
    Class that allows you to easily turn scale factor into snap number.
    Usage:
    from mergertrees.GetScaleFactors import getsnap
    getsnap = getsnap() #create object named getsnap. Note this destroys your import!
    # also works with parentheses instead of brackets if you so desire
    getsnap[1.0]
    getsnap[[1.0,0.9,0.8]]
    

    Implementation: spline the snapnums against the scale factors, then int(round(spline))
    (So note that you can give it scale factors that aren't close to a spline)
    """
    def __init__(self,path='/spacebase/data/AnnaGroup/caterpillar/parent/RockstarData',
                 minsnap=0,maxsnap=63):
        self.minsnap = minsnap
        self.maxsnap = maxsnap
        self.scale_list = getScaleFactors(path)
        self.snap_list = range(minsnap,maxsnap+1) #minsnap to maxsnap inclusive
        self.spl = interpolate.UnivariateSpline(self.scale_list,self.snap_list,s=0)
    def getsnap(self,x):
        snap=int(round(self.spl(x)))
        if snap>max(self.snap_list):
            print "WARNING: snap is "+str(snap)+", larger than largest snap!"
        if snap<min(self.snap_list):
            print "WARNING: snap is "+str(snap)+", less than smallest snap!"
        return snap
    def __getitem__(self,key):
        try:
            return [self.getsnap(x) for x in key]
        except:
            return self.getsnap(key)
    def __call__(self,arg):
        return self.__getitem__(arg)

#path = '/spacebase/data/AnnaGroup/cosm1/rockstar'
#scales = getScaleFactors(path)
#print scales
