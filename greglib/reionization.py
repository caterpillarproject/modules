import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate
import statsmodels.nonparametric.smoothers_lowess as smoother
import DwarfMethods as dm
from PlotParams import *


# z = 9.33, 10.73, 12.1, 13.57
h0 = 0.6711



def get_Barber_data():
    path = '/nfs/blank/h4231/gdooley/DwarfsOfDwarfs/code/for_Dooley/'
    classical = np.loadtxt(path+'classical.dat')
    luminous = np.loadtxt(path+'luminous.dat')
    allsats = np.loadtxt(path+'allSats.dat')
    # convert masses to log10 in Msun, not log10(msun/h)
    classical = np.log10(10**classical/.73)  # hubble value from the Barber paper is 0.73
    luminous = np.log10(10**luminous/.73)
    allsats = np.log10(10**allsats/.73)
    return classical, luminous,allsats
    

# generate log_mass array and a dark fraction array.
# then create an interpolation (linear) between the points
def get_lum_fraction_Barber():
    classical, luminous,allsats = get_Barber_data()
    log_mass = np.arange(6.7,10.5,.1) # play with binnings
    hist_all, rvals = np.histogram(allsats,bins=log_mass)
    hist_lum, rvals = np.histogram(luminous,bins=log_mass)
    x_axis = dm.getMidpoints(rvals)

    fracs = np.array(hist_lum,dtype=float)/np.array(hist_all)

    smoothed_fracs= smoother.lowess(fracs,x_axis,frac=0.25)[:,-1]
    tck = interpolate.splrep(x_axis, smoothed_fracs,k=1)
    return tck, x_axis, fracs # smoothed_fracs

    """ # for plotting the luminous fraction
    plt.bar(rvals[:-1], hist_all,width=0.1)
    plt.bar(rvals[:-1], hist_lum,width=0.1,color='purple')
    plt.ylabel('Total Number of halos')
    plt.xlabel('$\log_{10}(M_{200}) M_\odot$')
    plt.yscale('log')
    plt.savefig('lum_frac_barber_hist')
    plt.close()


    plt.plot(x_axis, smoothed_fracs,color='pink',ls='-',lw=2)
    plt.plot(x_axis, interpolate.splev(x_axis,tck),color='red',ls='--',lw=2)
    plt.plot(x_axis, fracs)
    plt.ylabel('fraction luminous')
    plt.xlabel('$\log_{10}(M_{200}) M_\odot$')
    plt.savefig('lum_frac_barber')   
    plt.close()
    """


# given halo mass, what fraction are luminous? Halo mass must be in terms of M200 infall.
# For all but moster, must therefore make a conversion from Mhalo to M200 infall.
def frac_luminous(M,tck):  
    if len(M)==0:
        return np.array([])
    return interpolate.splev(np.log10(M),tck)
    # if >1 or < 0, doesn't matter. always hosts stars or never does when compared to random int
    #mask0=values<0
    #values[mask0]=0
    #mask1=values>1
    #values[mask1]=1
    #return values


def get_lum_fraction_caterpillar(mdef,z=13):
    if mdef=='m200':
        mstring='infall_mass200'
    if mdef=='mpeak':
        mstring='max_mass'
    if mdef=='m350':
        mstring='max_mass350NFW'
    h0=.6711
    # want to collect infall 200 masses for all halos
    # and for just the luminous ones
    m200_all = np.array([])
    m200_lum9 = np.array([])
    m200_lum11 = np.array([])
    m200_lum12 = np.array([])
    m200_lum13 = np.array([])
    hpaths = dm.get_hpaths(False)
    for hpath in hpaths[0:25]:
        data = dm.get_extant_data(hpath,False)
        m200_all = np.append(m200_all, np.array(data[mstring]/0.6711,dtype=np.float64))

        # now do the reionization checks for z=11
        #mask9 = (data['vmax_9'] > 10)
        #mask11 = (data['vmax_11'] > 10)
        #mask12 = (data['vmax_12'] > 10)
        #mask13 = (data['vmax_13'] > 10)

        mask9 = (data['tvir_9'] > 3000)
        mask11 = (data['tvir_11'] > 3000)
        mask12 = (data['tvir_12'] > 3000)
        mask13 = (data['tvir_13'] > 3000)

        mask9 = (data['m200_9']/h0 > 10**6.8)
        mask11 = (data['m200_11']/h0 > 10**6.8)
        mask12 = (data['m200_12']/h0 > 10**6.8)
        mask13 = (data['m200_13']/h0 > 10**6.8)


        #print np.sum(mask), 'num bigger than 15'
        mask2 = data['peak_vmax'] >= 25  #25
        #print np.sum(mask&mask2), 'num also peaking higher than 50'
        m200_lum9 = np.append(m200_lum9,  np.array(data[mstring][mask9 | mask2]/0.6711,dtype=np.float64))
        m200_lum11 = np.append(m200_lum11,  np.array(data[mstring][mask11 | mask2]/0.6711,dtype=np.float64))
        m200_lum12 = np.append(m200_lum12,  np.array(data[mstring][mask12 | mask2]/0.6711,dtype=np.float64))    
        m200_lum13 = np.append(m200_lum13,  np.array(data[mstring][mask13 | mask2]/0.6711,dtype=np.float64))

    log_mass = np.arange(7.5,10.5,.1) # play with binnings
    hist_all, rvals = np.histogram(np.log10(m200_all),bins=log_mass)
    hist_lum9, rvals = np.histogram(np.log10(m200_lum9),bins=log_mass)
    hist_lum11, rvals = np.histogram(np.log10(m200_lum11),bins=log_mass)
    hist_lum12, rvals = np.histogram(np.log10(m200_lum12),bins=log_mass)
    hist_lum13, rvals = np.histogram(np.log10(m200_lum13),bins=log_mass)

    x_axis = dm.getMidpoints(rvals)
    fracs9 = np.array(hist_lum9,dtype=float)/np.array(hist_all)
    fracs11 = np.array(hist_lum11,dtype=float)/np.array(hist_all)
    fracs12 = np.array(hist_lum12,dtype=float)/np.array(hist_all)
    fracs13 = np.array(hist_lum13,dtype=float)/np.array(hist_all)

    """
    # for plotting the luminous fraction
    plt.bar(rvals[:-1], hist_all,width=0.1)
    plt.bar(rvals[:-1], hist_lum11,width=0.1,color='purple')
    plt.ylabel('Total Number of halos')
    plt.xlabel('$\log_{10}(M_{200}) M_\odot$')
    plt.yscale('log')
    plt.savefig('lum_frac_caterpillar_hist')
    plt.close()
    """

    if z==13:
        yvals=fracs13
        print 'getting reionization at 13.57'
    elif z==12:
        yvals=fracs12
    elif z==11:
        yvals=fracs11
    elif z==9:
        yvals=fracs9
    else:
        print "ERROR: Z=13,12,11 or 9. No other values can be specified now"

    smoothed_fracs= smoother.lowess(yvals,x_axis,frac=0.25)[:,-1]
    tck = interpolate.splrep(x_axis, smoothed_fracs,k=1)
    return tck, x_axis, fracs9, fracs11, fracs12, fracs13

def plot_frac_lum_mdef():
    tck_barb, xbarber, frac_barber = get_lum_fraction_Barber()
    tck_m350, xm350, frac_m350 = get_lum_fraction_m350()
    tck_mpeak, xmpeak, frac_mpeak = get_lum_fraction_mpeak()
    plt.plot(xbarber, interpolate.splev(xbarber,tck_barb),ls='-',lw=linewidth,label='Barber Data ($M_{\mathrm{200}}^{\mathrm{infall}}$)')
    plt.plot(xm350, interpolate.splev(xm350,tck_m350),ls='-',lw=linewidth,label='$M_{\mathrm{350}}^{\mathrm{peak}}$')
    plt.plot(xmpeak, interpolate.splev(xmpeak,tck_mpeak),ls='-',lw=linewidth,label='$M_{\mathrm{vir}}^{\mathrm{peak}}$')
    plt.legend(loc='upper left', fontsize=legend_size, frameon=False)

    plt.ylim((0,1))
    plt.ylabel('Fraction Luminous',fontsize=label_font )
    plt.xlabel('$\log_{10}(M_{\mathrm{halo}}) \mathrm{M_\odot}$',fontsize=label_font)
    plt.gcf().subplots_adjust(bottom=0.15)
    plt.savefig('lum_frac_mdef')
    plt.close()


def plot_frac_luminous():
    tck_barb, xbarber, frac_barber = get_lum_fraction_Barber()
    tck_cat, xcat, fracs9,fracs11,fracs12,fracs13 = get_lum_fraction_caterpillar(mdef='m200')
    tck_cat_350, xcat_350, fracs9_350,fracs11_350,fracs12_350,fracs13_350 = get_lum_fraction_caterpillar(mdef='m350')
    tck_cat_peak, xcat_peak, peak_fracs9,peak_fracs11,peak_fracs12,peak_fracs13 = get_lum_fraction_caterpillar(mdef='mpeak')

    plt.plot(xbarber, interpolate.splev(xbarber,tck_barb),color='red',ls='--',lw=3)
    plt.plot(xbarber, frac_barber, label='Barber',lw=3)

    plt.plot(xcat, fracs13,lw=3,label='m200')
    plt.plot(xcat_350, fracs13_350,lw=3,label='m350')
    plt.plot(xcat_peak, peak_fracs13,lw=3,label='mpeak')

    #plt.plot(xcat, fracs9,lw=3,label='$z_{reion} = 9.33$')
    #plt.plot(xcat, fracs11,lw=3,label='$z_{reion} = 10.73$')
    #plt.plot(xcat, fracs12,lw=3,label='$z_{reion} = 12.1$')
    #plt.plot(xcat, fracs13,lw=3,label='$z_{reion} = 13.57$')
    plt.legend(loc='upper left', fontsize=legend_size, frameon=False)

    plt.ylabel('Fraction Luminous', fontsize=label_font)
    plt.xlabel('$\log_{10}(M_{\mathrm{halo}}) M_\odot$',fontsize=label_font)

    plt.gcf().subplots_adjust(bottom=0.15)
    plt.savefig('lum_frac_comparison')
    plt.close()




# can take lum_frac_barber applied to all the m200 masses and generate random yes or no on them.
# get data on yes or no of the m350 values and plot the fractions.
# re-generate it 100 times and combine the values.
def barber_to_other():
    tck_barb, xbarber, frac_barber = get_lum_fraction_Barber()    

    lum_m350=[]; all_m350 = []
    lum_mpeak = []; all_mpeak = []
    hpaths = dm.get_hpaths(False)

    for hpath in hpaths:
        data = dm.get_extant_data(hpath,False)
        m200 = data['infall_mass200']/h0
        lum_chances = frac_luminous(m200, tck_barb)

        for i in range(20):
            randnums = np.random.rand(len(lum_chances))
            mask = randnums < lum_chances # chooses which halos will form stars
        
            lum_m350.append( np.array(data['max_mass350NFW'][mask]/h0))
            all_m350.append( np.array(data['max_mass350NFW']/h0))
            
            lum_mpeak.append( np.array(data['max_mass'][mask]/h0))
            all_mpeak.append( np.array(data['max_mass']/h0))

    lum_m350 = np.array([item for arr in lum_m350 for item in arr])
    all_m350 = np.array([item for arr in all_m350 for item in arr])
    lum_mpeak = np.array([item for arr in lum_mpeak for item in arr])
    all_mpeak = np.array([item for arr in all_mpeak for item in arr])

    np.save('lum_m350', lum_m350)
    np.save('all_m350', all_m350)
    np.save('lum_mpeak', lum_mpeak)
    np.save('all_mpeak', all_mpeak)



def get_m350_data():
    path = '/nfs/blank/h4231/gdooley/DwarfsOfDwarfs/code/'
    return np.load('lum_m350.npy'), np.load('dark_m350.npy')
    

# generate log_mass array and a dark fraction array.
# then create an interpolation (linear) between the points
def get_lum_fraction_m350():
    luminous, allsats = np.load('lum_m350.npy'), np.load('all_m350.npy')
    log_mass = np.arange(6.7,10.5,.1) # play with binnings
    hist_all, rvals = np.histogram(np.log10(allsats),bins=log_mass)
    hist_lum, rvals = np.histogram(np.log10(luminous),bins=log_mass)
    x_axis = dm.getMidpoints(rvals)
    fracs = np.array(hist_lum,dtype=float)/np.array(hist_all)

    smoothed_fracs= smoother.lowess(fracs,x_axis,frac=0.25)
    tck = interpolate.splrep(smoothed_fracs[:,0], smoothed_fracs[:,-1],k=1)
    return tck, x_axis, fracs # smoothed_fracs


def get_lum_fraction_mpeak():
    luminous, allsats = np.load('lum_mpeak.npy'), np.load('all_mpeak.npy')
    log_mass = np.arange(7.5,10.5,.1) # play with binnings
    hist_all, rvals = np.histogram(np.log10(allsats),bins=log_mass)
    hist_lum, rvals = np.histogram(np.log10(luminous),bins=log_mass)
    x_axis = dm.getMidpoints(rvals)
    fracs = np.array(hist_lum,dtype=float)/np.array(hist_all)

    smoothed_fracs= smoother.lowess(fracs,x_axis,frac=0.25)
    tck = interpolate.splrep(smoothed_fracs[:,0], smoothed_fracs[:,-1],k=1)
    return tck, x_axis, fracs # smoothed_fracs
