import asciitable
import scipy.interpolate as interpolate
import numpy as np

class sd93cool:
    """
    Sutherland and Dopita 1993 Cooling Functions
    http://www.mso.anu.edu.au/~ralph/data/cool/ABOUT.txt
    """

    def __init__(self,NEQ=False,zf=False,filedir="/spacebase/data/alexji/lib/python/caterpillarmodules/alexlib/croton06/sd93cool/"):
        """
        Load in a cooling function grid

        Only implemented for CIE models right now
        """
        self.filedir = filedir
        self.logZarr = [-99,-3.0,-2.0,-1.5,-1.0,-0.5,0.0]
        self.logTarr = np.linspace(4,8.5,91)
        self.logLgrid = ['' for x in self.logZarr]
        for i,file in enumerate([self.filedir+x+'.cie' for x in ['mzero','m-30','m-20','m-15','m-10','m-05','m-00']]):
            data = asciitable.read(file,data_start=2)
            self.logLgrid[i] = data['col5']
        self.logLgrid = np.array(self.logLgrid)
        self.logLfns = ['' for x in self.logTarr]
        for i in xrange(len(self.logTarr)):
            self.logLfns[i] = interpolate.InterpolatedUnivariateSpline(self.logZarr[1:],self.logLgrid[1:,i],k=1)

    def get_cooling_function(self,logZ=None):
        """
        Returns the cooling function logLambda(logT) for gas with metallicity 10^logZ
        If logz = None, get the cooling function for metal-free (Z=0) gas

        Note that I'm saying 'Z' for metallicity, but actually SD93 use [Fe/H].
        Since almost all the cooling functions are for solar abundance ratios,
        this should be the same. But for CIE
        """
        
        if logZ==None or logZ < -3.0:
            ## metal free gas
            if logZ != None and logZ < -3.0:
                print "IMPORTANT: entered logZ=",logZ
                print "But rounding down to Z=0 since the cooling functions are almost the same"
            return interpolate.InterpolatedUnivariateSpline(self.logTarr,self.logLgrid[0])
        else:
            Larr = [f(logZ) for f in self.logLfns]
            return interpolate.InterpolatedUnivariateSpline(self.logTarr,Larr)

    def get_Tarr(self):
        return self.logTarr

if __name__ == "__main__":
    print "Plotting default cooling curves"
    import pylab as plt
    s = sd93cool()
    for logZ in [None, -3.0, -2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0]:
        f = s.get_cooling_function(logZ=logZ)
        plt.plot(s.logTarr,f(s.logTarr))
    plt.ylim((-24,-20))
    plt.xlabel('log T [K]')
    plt.ylabel(r'log $\Lambda$ [ergs/cm^3]')
    plt.show()
          

##         ##Load Z=0 separately
##         data = asciitable.read(self.filedir+"mzero.cie",data_start=2)
##         self.logTarr = data['col1']
##         self.nometLnet = data['col5']

##         plt.plot(self.logTarr,self.nometLnet)

##         self.logTarr = np.arange(4,7.55,.05)

##         ##Pick which files to open
##         if NEQ:
##             self.logZgrid = [-3.0,-2.0,-1.5,-1.0,-0.5]
##             if zf:
##                 self.filelist = ['pk6zf75'+x+'.neq' for x in ['m-30','m-20','m-15','m-10','m-05']]
##             else:
##                 self.filelist = ['pk6ff75'+x+'.neq' for x in ['m-30','m-20','m-15','m-10','m-05']]
##         else:
##             self.logZgrid = [-3.0,-2.0,-1.5,-1.0,-0.5,0.0]
##             self.filelist = [x+'.cie' for x in ['m-30','m-20','m-15','m-10','m-05','m-00']]

##         ##Load Z != 0 grid
##         self.Lnetarr = ['' for x in self.logZgrid]
##         for idx,file in enumerate([self.filedir+x for x in self.filelist]):
##             data = asciitable.read(file,data_start=2)
##             logTarr = (data['col1'])[::-1]
##             print logTarr
##             logLnet = (data['col5'])[::-1]
##             print logLnet
##             f = interpolate.InterpolatedUnivariateSpline(logTarr,logLnet,k=1)
##             self.Lnetarr[idx] = f(self.logTarr)
##         print self.logTarr
##         print f(self.logTarr)

##         plt.plot(self.logTarr,f(self.logTarr))
##         plt.show()
        
##         #print self.Lnetarr
##         #self.Lnetarr = np.array(self.Lnetarr)
##         #print self.Lnetarr.shape

##         #self.LnetFns = ['' for x in self.logTarr]
##         #for i in xrange(len(self.logTarr)):
##         #    self.LnetFns[i] = interpolate.InterpolatedUnivariateSpline(self.logZgrid,

