#import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
#import DwarfMethods as dm
#from PlotParams import *
#import haloutils as htils
#import MTanalysis3 as mta
#import abundance_matching as am

# SHMF for MW size for peak mass
alpha = 1.889667
K = .003906

# SHMF for field halos

def Nsubs_bigger(msub, Mhost,alpha,K):
    return K * Mhost /(alpha-1) * (msub**(1-alpha) - Mhost**(1-alpha))


"""
def peak_ngreater(hpath,cat,haloid,target_halo_mass=1.6e12):
    AE = mta.AllExtantData()
    dataE = AE.read(hpath)
    peak_masses = dataE['max_mass']/cat.h0
    bins_per_dex = 5
    min_mass=7.5; max_mass=np.log10(cat.ix[haloid]['mgrav']/cat.h0)-1.5 # peak mass must be 7.5, or else dataE doesn't hav enough values
    histrange = np.arange(min_mass,max_mass+.01,1./bins_per_dex)
    y_sub = np.array([np.sum(peak_masses > lowmass) for lowmass in 10**histrange])
    y_sub = y_sub*(target_halo_mass)/(cat.ix[haloid]['mgrav']/cat.h0)
    #print y_sub, 'normed to host mass values'
    return 10**histrange, y_sub


def plot_ngreater():
    target_halo_mass = 1.6e12
    hpaths = dm.get_hpaths(field=False)[0:20]  #htils.get_all_halo_paths_lx(lx)
    y_values_full = []
    for hpath, color in zip(hpaths[0:],colorlist[0:len(hpaths)]):
        snap_z0 = htils.get_numsnaps(hpath)-1
        cat=htils.load_rscat(hpath,snap_z0,rmaxcut=False)
        hostID = htils.load_zoomid(hpath)
        x_axis, y_axis = peak_ngreater(hpath,cat,hostID,target_halo_mass)
        y_values_full.append(y_axis)
        plt.plot(x_axis, y_axis, linewidth=1.0,color=color)

"""



"""
mhost = 1.6e12
dm_masses = 10**np.linspace(7,12,20)
ngreater = Nsubs_bigger(dm_masses, mhost, alpha,K)
model = am.GarrisonKimmel()
halo_cd = model.stellar_to_halo_mass(10**6)
mthresh = model.stellar_to_halo_mass(4.5*10**5)
print Nsubs_bigger(halo_cd,mhost,alpha,K), 'nsubs bigger GK14 10**6 Mstar in 1.6e12 host'
print 3.85*( halo_cd / .01 / mhost)**(-.9), 'nsubs bigger GK14 from his paper'


print Nsubs_bigger(mthresh,mhost,alpha,K), 'nsubs bigger GK14 4.5e10**5 Mstar in 1.6e12 host'
print 3.85*( mthresh / .01 / mhost)**(-.9), 'nsubs bigger GK14 from his paper 4.5e5'

ngreater_GK = 3.85*( dm_masses / .01 / mhost)**(-.9)


#plot_ngreater()

plt.vlines(halo_cd, 1,1000)
plt.plot(dm_masses, ngreater,label='Caterpillar',lw=2.5)
plt.plot(dm_masses, ngreater_GK,label='GK',lw=2.5)
#plt.plot(dm_masses, ngreater_GK*1.5,label='GK*1.5')
plt.legend()
plt.xscale('log')
plt.yscale('log')
plt.ylim((.8, 1000 ))
plt.xlim((5*10**7, 8*10**11))
plt.grid()
plt.savefig('compare_to_fig3_elvis')
"""


# for m peak
alpha = 1.889667
K = .003906

dm_masses = 10**np.linspace(7,12,20)
for M in [1.6]:
    ngreater_cater = Nsubs_bigger(dm_masses, 10**12*M, alpha,K)
    ngreater_GK = 3.85*( dm_masses / .01 / (M* 10**12))**(-.9)

    plt.plot(dm_masses, ngreater_cater,label='Caterpillar %.1f' %M,lw=2.5)
    plt.plot(dm_masses, ngreater_GK,label='Elvis %.1f' %M,lw=2.5)
    plt.plot(dm_masses, ngreater_GK/.71,label='Elvis/h %.1f' %M,lw=2.5, ls= '-.')
    
plt.legend()
plt.xscale('log')
plt.yscale('log')
plt.ylim((.8, 1000 ))
plt.xlim((5*10**7, 8*10**11))
plt.ylabel('$N(>M_{peak})$')
plt.xlabel('$ M_{bound}^{z=0} \, (\mathrm{M_\odot})$')
plt.grid()
plt.savefig('compare_to_eq1_elvis')
plt.close()





# for z=0 bound mass
alpha = 1.9
K = .001133
dm_masses = 10**np.linspace(7,12,20)
for M in [1.0,1.6,2.2]:
    ngreater_cater = Nsubs_bigger(dm_masses, 10**12*M, alpha,K)
    ngreater_GK = 1.11*( dm_masses / .01 / (M*10**12))**(-.95)

    plt.plot(dm_masses, ngreater_cater,label='Caterpillar %.1f' %M,lw=2.5)
    plt.plot(dm_masses, ngreater_GK,label='Elvis %.1f' %M,lw=2.5)
    

plt.legend()
plt.xscale('log')
plt.yscale('log')
plt.ylim((.8, 1000 ))
plt.xlim((5*10**7, 10**11))
plt.grid()
plt.ylabel('$N(>M_{bound}^{z=0})$')
plt.xlabel('$ M_{bound}^{z=0} \, (\mathrm{M_\odot})$')
plt.savefig('compare_to_eq2_elvis')
plt.close()


