import numpy as np
import reionization as reion

ML_ratio = 1.0
N_iter = 1000 # make 20,000 later


# subhalo mass function next. Make alpha positive
def Nsubs_bigger(msub, Mhost,alpha,K):
    return K * Mhost /(alpha-1) * (msub**(1-alpha) - Mhost**(1-alpha))

# generate a list of halo masses from mlow=10**7.5 to mhigh = Mhost that follows the mean of the SHMF
def generate_shm_sample(Mhost,alpha,K): 
    mlow = 10**7.4 # lower limit of where stars can form above 10**3 Lsun, and not dark from reionization
    # as in sawala paper. Could go to 10**7 to be safer, but the luminous fraction is almost 0 there.
    mean = Nsubs_bigger(mlow,Mhost,alpha,K)
    nsubs = np.random.poisson(mean)
    randnums = np.random.rand(nsubs) # a number between 0 and 1
    masses = (mlow**(1-alpha) - (alpha-1)/(K*Mhost) * randnums*mean)**(1/(1-alpha))
    return masses


class AMbase(object):
    """
    methods that are standard to all AM models.
    methods that are unique should be over written by particular AM models
    """
    def __init__(self):
	return

    # given halo mass, what is stellar mass?
    def halo_to_stellar_mass(self,mvir,a=1):
	if type(mvir)==list:
            mvir = np.array(mvir)
        return self.getStellarMass(mvir,a)

    # given stellar mass/luminosity what is DM halo mass?
    def stellar_to_halo_mass(self, mstar, a=1):
	#USE INTERPOLATION     
 	from scipy import interpolate
 	mass_array = 10**np.arange(7,12,.07)
	mstar_array = self.halo_to_stellar_mass(mass_array,a)
        tck = interpolate.splrep(np.log10(mstar_array),np.log10(mass_array), k=1)
        return 10**interpolate.splev(np.log10(mstar),tck)

    def P_at_least_one(self, halo_masses, min_lum):
        if self.isPoisson:
            lowest_mass = self.stellar_to_halo_mass(ML_ratio*min_lum)
            mean = Nsubs_bigger(lowest_mass,halo_masses,self.alpha,self.K)
            Pgt1 = 1 - np.e**-mean
            return Pgt1
        else:
            frac_gte_one = []
            for halo_mass in halo_masses:
                samples = self.generate_stellar_mass_samples(halo_mass,N=N_iter) # generate list of lists of stellar masses
                lums = samples/ML_ratio
                Ngreater = np.array([np.sum(lum>min_lum) for lum in lums])
                frac_gte_one.append(np.sum(Ngreater>0) / float(len(Ngreater)))
            Pgt1 = np.array(frac_gte_one)
            return Pgt1


    def ngreater(self,halo_masses,min_lum):
        if self.isPoisson:
            lowest_mass = self.stellar_to_halo_mass(ML_ratio*min_lum)
            N_visible = Nsubs_bigger(lowest_mass, halo_masses,self.alpha,self.K)
            return N_visible, np.sqrt(N_visible)
        else:
            if type(halo_masses)==list or type(halo_masses)==np.ndarray:
                num_gte_min = []; spread = []
                for halo_mass in halo_masses:
                    samples = self.generate_stellar_mass_samples(halo_mass,N=N_iter) # generate list of lists of stellar masses
                    lums = samples/ML_ratio
                    Ngreater = np.array([np.sum(lum>min_lum) for lum in lums])
                    num_gte_min.append( np.mean(Ngreater))
                    spread.append( np.std(Ngreater))
                N_visible = np.array(num_gte_min); N_std = np.array(spread)
            else:
                samples = self.generate_stellar_mass_samples(halo_masses,N=N_iter) # generate list of lists of stellar masses
                lums = samples/ML_ratio
                Ngreater = np.array([np.sum(lum>min_lum) for lum in lums])
                N_visible = np.mean(Ngreater)
                N_std = np.std(Ngreater)
            return N_visible, N_std



    # generate list of lists of stellar masses
    # must have getStellarMass_random() implemented in subclass
    def generate_stellar_mass_samples(self,halo_mass,N):
        stellar_mass_lists=[]
    	for i in range(N):
            dm_masses = generate_shm_sample(halo_mass,self.alpha,self.K)
            if self.reionization:
                lum_chances = reion.frac_luminous(dm_masses, self.lum_tck) # hopefully a list of fractions between 0 and 1
                randnums = np.random.rand(len(lum_chances))
                mask = randnums < lum_chances # chooses which halos will form stars
                dm_masses = dm_masses[mask] # dm masses that host stars

            stellar_mass_lists.append(self.getStellarMass_random(dm_masses,a=1.0))
   	return np.array(stellar_mass_lists)



    def get_field_total(self,dm_masses,min_lum):
        N, _ = self.ngreater(dm_masses,min_lum)
        return np.sum(N)

    def get_field_distr(self,dm_masses,min_lum,N=10000):
        distr=np.array([0]*N)
        for halo_mass in dm_masses:
            samples = self.generate_stellar_mass_samples(halo_mass,N)
            lums = samples/ML_ratio
            Ngreater = np.array([np.sum(lum>min_lum) for lum in lums]) #should be 10,000 long
            distr+=Ngreater
        return distr






