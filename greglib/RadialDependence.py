import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import haloutils as htils
from scipy import stats
import methods
from PlotParams import *
from scipy.optimize import curve_fit
from scipy import interpolate

from scipy.integrate import dblquad


def plot_radial_fits():
    radius2 = np.array([0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2])
    radius = np.array([0.0, 0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2])

    Khost = np.array([0.0, 0.000008, 0.000047,0.000129,0.00023,0.000362,0.000504,0.000655,0.000811,0.000971,0.001122,0.001269,0.001421,0.001549,0.001681,0.001802,0.001905,0.002003,0.002095,0.00218,0.002257])

    Kfield = np.array([0.0, 0.000008, 0.00006,0.000144,0.00025,0.000386,0.000511,0.00078,0.000961,0.001139,0.001314,0.001526,0.001677,0.001805,0.001972,0.002118,0.002277,0.002452,0.002603,0.002761,0.003])

    Khost2 = np.array([0, 0.000008, 0.000045, 0.000139, 0.000245, 0.000384, 0.000521,0.000672,0.000815,0.000992,0.001133,0.001282,0.001427,0.001545,0.001662,0.001784,0.001879,0.001962,0.002068,0.002148,0.002228])

    volume = 4/3.*np.pi*radius**3

    # fit a line to the data
    slope, intercept, r_value, p_value, std_err = stats.linregress(radius,Khost)
    slopeV, interceptV, r_valueV, p_valueV, std_errV = stats.linregress(np.log10(volume),np.log10(Khost))
    slopeF, interceptF, r_valueF, p_valueF, std_errF = stats.linregress(radius,Kfield)
    slopeFV, interceptFV, r_valueFV, p_valueFV, std_errFV = stats.linregress(np.log10(volume),np.log10(Kfield))
    #print slope, intercept, r_value, p_value, std_err, 'host slope, intercept, r, p, std error'
    #print slopeF, interceptF, r_valueF, p_valueF, std_errF, 'field slope, intercept, r, p, std error'
    #print slopeV, interceptV, r_valueV, p_valueV, std_errV, 'host by volume slope, intercept, r, p, std error'

    #slope, a,b,c = np.linalg.lstsq(np.vstack([radius, np.ones(len(radius))]).T, Kfield)
    slopeF2, _,_,_ = np.linalg.lstsq(radius[:,np.newaxis],Kfield)

    (c1,c2,c3,c4) = np.polyfit(radius,Khost,3)
    (d1,d2,d3,d4) = np.polyfit(radius2,Khost2,3)

    # try a logistics fit
    def logistic(x, a,b,k):
        return a/(1+ b* np.e**(-k * (x)))
    # end line fitting

    popt, pcov = curve_fit(logistic, radius2[1:], Khost2[1:])
    print logistic(1,popt[0],popt[1],popt[2]), 'logistic at 1. should be .0011'    




    #plt.plot(radius, interceptF+slopeF*radius,color='red',ls='--')
    #plt.plot(radius, slopeF2*radius,color='red',ls='--')
    
    norm = c4+c3+c2+c1
    print norm, 'should be close to K0 = 0.001133'
    #print (c4+c3+c2+c1)/norm, 'value at rvir is 1 rvir'
    #print c4,c3,c2,c1, 'c4,c3,c2,c1'
    print c4/norm, c3/norm, c2/norm, c1/norm, 'normed constants'


    norm = d4+d3+d2+d1
    #print (d4+d3+d2+d1)/norm, 'sum of d values. value at rvir is 1 rvir'
    #print d4,d3,d2,d1, 'd4,d3,d2,d1'
    print d4/norm, d3/norm, d2/norm, d1/norm, 'normed constants with no rmax cut'


    rr = np.arange(0.05,2,.01)
    plt.plot(rr, logistic(rr, popt[0],popt[1],popt[2]), color='magenta',ls='--',label='logistic fit')
    plt.plot(rr,d4+d3*rr+d2*rr**2+d1*rr**3,color='cyan',ls='--')
    plt.plot(rr,c4+c3*rr+c2*rr**2+c1*rr**3,color='green',ls='--')
    #plt.plot(radius,-.00015+.00071*radius+.00086*radius**2-.00031*radius**3,color='purple',ls='--')
    plt.plot(radius,Khost,color='blue',label='host')
    plt.plot(radius,Kfield,color='orange',label='field')
    plt.plot(radius2,Khost2,color='red',label='host no cut')
    plt.legend(loc='upper left')
    plt.xlabel('radius [kpc]')
    plt.xlim((0.1,2))
    plt.yscale('log')
    plt.xscale('log')
    plt.savefig('radialfit_logscale')
    
    plt.close()



    rr = np.arange(0.0,2,.01)
    plt.plot(rr, logistic(rr, popt[0],popt[1],popt[2]), color='magenta',ls='--',label='logistic fit')
    plt.plot(rr,d4+d3*rr+d2*rr**2+d1*rr**3,color='cyan',ls='--')
    plt.plot(rr,c4+c3*rr+c2*rr**2+c1*rr**3,color='green',ls='--')
    print c4+c3*.1+c2*.1**2+c1*.1**3, 'value of green fit at 0.1'
    #plt.plot(radius,-.00015+.00071*radius+.00086*radius**2-.00031*radius**3,color='purple',ls='--')
    plt.plot(radius,Khost,color='blue',label='host')
    plt.plot(radius,Kfield,color='orange',label='field')
    plt.plot(radius2,Khost2,color='red',label='host no cut')
    plt.legend(loc='upper left')
    plt.xlabel('radius [kpc]')
    plt.xlim((0.0,2))
    plt.savefig('radialfit')
    plt.close()

    
    """
    plt.plot(volume,Khost,color='blue',label='host')
    plt.plot(volume,Kfield,color='orange',label='field')
    plt.xlabel('volume $[kpc^3]$')
    plt.xscale('log')
    plt.yscale('log')
    plt.legend(loc='upper left')
    plt.savefig('radialfit_volume')
    plt.xlim((0,2))
    plt.close()
    """


# values correspond to host halo values and the commented normed values k1=-.1366
def getK_not_normed(r):
    return -0.00015148503612 + 0.000708840718105*r + 0.000861092068406*r**2 - 0.00030949197861*r**3
    #return -.00015+.00071*r+.00086*r**2-.00031*r**3

# input the radial ratio r/rvir. normalized to 1 at r=1
def getK(r):
    return -0.0439735566925 + 0.39131605932*r + 0.996495770479*r**2 - 0.343838273107*r**3  # fit to rmaxcut = False. best so far

# normed values of k
#k1 = -0.136601512859; k2=0.639196563237; k3=0.77648909931; k4= -0.279084149688  # fit to rmaxcut=True above 2.0
k1 = -0.0439735566925; k2=0.39131605932; k3=0.996495770479; k4=-0.343838273107  # fit to rmaxcut = False. best so far
#k1=-0.0362134075824; k2=0.297178749054; k3=1.11253297655; k4=-0.37349831802  # fit to rmaxcut = True


def get_normed_factor(r):
    # constatns are a,b,c,d in the order written to correspond tomy math notes
    return k1 + k2*r + k3*r**2 + k4*r**3


def plot_normed_fit():
    radius = np.arange(0.2,2.1,.1)
    plt.plot(radius, get_normed_factor(radius))
    plt.savefig('normed_radial_fit')
    plt.close()


def get_normed_los(R,Z):
    g = np.arctan2(R,Z)
    term1 = Z*k2/2. * np.log(1+R**2/Z**2)
    term2 = k3*Z**2 * (np.sqrt(1+R**2/Z**2) - 1)
    term3 = k4*R**2*Z/2.

    term4 = k2*R*(np.pi/2 - g)
    term5 = k3*R**2 * np.log( Z/R *(1+np.sqrt(1+R**2/Z**2)))
    term6 = k4*R**2*Z

    return term1+term2+term3+term4+term5+term6 + k1  # pretty sure I need that constant in front!!



# could add numerical integration to the plot to show how well it fits
# 
def get_numerical_normed_los(R,Z, tck):
    def dkdr(r):
        return interpolate.splev(r,tck,der=1) / 0.001133

    t0 = np.arctan2(R,Z)
    term1 = dblquad(lambda r, theta: dkdr(r)*np.sin(theta), 0, t0, lambda x: 0, lambda x: Z/np.cos(x))
    term2 = dblquad(lambda r, theta: dkdr(r)*np.sin(theta), t0, np.pi/2, lambda x: 0, lambda x: R/np.sin(x))
    print 'finished an iteration of integrating', R
    return term1[0] + term2[0]


# this is for khost2
#[ 0.01189194  0.03825338  0.1338373   0.40731802  0.69055662  0.94740314  1.16469735 
# 1.33790466  1.4688462   1.52461371] values 3 for [0.05,.1,.2,.4,.6,.8,1,1.2,1.4,1.5]

# this is for khost
#[ 0.01172879  0.03766099  0.12990579  0.39270223  0.67768313  0.93840293
#  1.15882717  1.34039798  1.48181768  1.53956582] values 3 for [0.05,.1,.2,.4,.6,.8,1,1.2,1.4,1.5]

# this is z=0.5 for khost
#[ 0.01007687  0.03110473  0.10446132  0.30095863  0.49895946  0.66964709
#  0.80700626  0.91751962  1.00130305  1.0353223 ] values 4 for [0.05,.1,.2,.4,.6,.8,1,1.2,1.4,1.5]

# this is z=0.5 for khost2
#[ 0.01027001  0.03181417  0.10881709  0.31673135  0.51349503  0.68128896
#  0.81923617  0.92638784  1.00304656  1.03662456] values 4 for [0.05,.1,.2,.4,.6,.8,1,1.2,1.4,1.5]


def plot_los_frac():
    radius = np.array([0.0, 0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2])
    Khost2 = np.array([0, 0.000008, 0.000045, 0.000139, 0.000245, 0.000384, 0.000521,0.000672,0.000815,0.000992,0.001133,0.001282,0.001427,0.001545,0.001662,0.001784,0.001879,0.001962,0.002068,0.002148,0.002228])
    Khost = np.array([0.0, 0.000008, 0.000047,0.000129,0.00023,0.000362,0.000504,0.000655,0.000811,0.000971,0.001122,0.001269,0.001421,0.001549,0.001681,0.001802,0.001905,0.002003,0.002095,0.00218,0.002257])

    # numerical integration
    tck =  interpolate.splrep(radius, Khost2, k=1)

    Z = 1
    r = np.arange(0.1,1.91,.05)
    values = get_normed_los(r,Z)
    values2 = get_normed_los(r,0.5)
    r3 = [0.05,.1,.2,.4,.6,.8,1,1.2,1.4,1.5]
    values3 = [0.01189194, 0.03825338, 0.1338373, 0.40731802, 0.69055662, 0.94740314, 1.16469735, 1.33790466, 1.4688462, 1.52461371]
    values4 = np.array([get_numerical_normed_los(R,0.5,tck) for R in r3])
    print values4, 'values 4 for [0.05,.1,.2,.4,.6,.8,1,1.2,1.4,1.5]'

    plt.plot(r3,values3, lw=linewidth, label='$Z=1.0$ raw data',color='red',ls='--')
    plt.plot(r3,values4, lw=linewidth, label='$Z=0.5$ raw data',color='orange',ls='--')

    plt.plot(r,values, lw=linewidth, label='$Z=1.0$',color='blue')
    plt.plot(r,values2, lw=linewidth,label='$Z=0.5$',color='green')
    plt.ylabel('K(R)',fontsize=label_font)
    plt.xlabel('$R \equiv R_{\mathrm{fov}}/ R_{\mathrm{vir}}$', fontsize=label_font)
    plt.tick_params(axis='both', which='major', labelsize=tick_size)
    plt.gcf().subplots_adjust(bottom=0.15)
    plt.xlim((0,1.5))
    plt.ylim((0,1.8))
    plt.grid()
    plt.legend(frameon=False,loc='upper left',fontsize=legend_size)
    plt.savefig('FieldOfView_MultiplicativeFactor')
    plt.close()


#plot_radial_fits()
#plot_los_frac()


"""
#plot_normed_fit()
K = get_normed_los(1./np.sqrt(2),1./np.sqrt(2))
print K, 'should be less than 1 and greater than 0.26'
K = get_normed_los(1.,1.)
print K,'should be less than 1.5 and greater than 1'
K = get_normed_los(0.001,0.001)
print  K, 'should get close to 0'
"""

#K = get_normed_los(0.5,1)
#print K, 'value for half of the virial radius'



