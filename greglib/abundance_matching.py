import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PlotParams import *
from AMbase import *
import infall_times
import random
#import DwarfMethods as dm
import reionization as reion

##### HOST SHMF DEFAULT VALUES MUST BE SET TO THE CORRECT VALUES AFTER THEY ARE COMPUTED!!!!!!######

ML_ratio = 1
N_iter = 1000  # 10000 for final plots, but for speed of testing things, make 2000

# assume a slope of -1.90 for the host at z=0. 
# r_ratio is r/r_vir, where r_vir depends on the host mass Mvir.
def radial_dependence(r_ratio):
    return -.00015+.00071*r_ratio + .00086*r_ratio**2 - .00031*r_ratio**3


class Sawala(AMbase):
    def __init__(self, hostSHMF=False, hostAlpha=1.90,hostK=.001122): # Host SHMF from the z=0 radial fit.
	self.slope1 = 2.2#2.00
	self.slope2 = 0.41#.546
	self.M1 = 10**10.347 / (10**(7.6143/self.slope1))  #34.7e5
	self.M2 = 10**8.0214 / (10**(4.9/self.slope2)) #0.11144
        # these values correspond to equation 1 in the Sawala abundance of not just dark matter haloes paper
	self.a=0.69
	self.Mt = 10**11.6
	self.w = 0.79
        self.b = 0.98
        self.ls = '-'

        # SHMF parameters. Slope alpha and constant, K
        self.alpha = 1.931428  #1.870
        self.K = .002578    #.000833
        if hostSHMF:
            self.alpha = hostAlpha
            self.K = hostK

        self.color = 'cyan'
        self.label = 'Sawala'

        self.sigma = 0.5  ## the lognormal scatter in the SHMF relationship. Estimated from plot
        self.isPoisson=False

    # Infered from figure 4 lower panel of Bent by Baryons
    # I should add in some scatter too
    # high mass end needs to be adjusted - force it to converge to moster for instance?????
    # yes, at 10**11, make it switch over to Moster et al.
    def getStellarMass_single(self,M,a=1):
        if M < 10**11:
            return (M/(self.M2))**self.slope2 + (M/(self.M1))**self.slope1
        else:
            model=Moster()
            return model.getStellarMass(M,a)

    def getStellarMass(self,M,a=1):
        if type(M)==list or type(M)==np.ndarray:
            stellar=[-1]*len(M)
            for i,mass in enumerate(M):
                stellar[i]=self.getStellarMass_single(mass)
            return np.array(stellar)
        else:
            return self.getStellarMass_single(M)

    def getStellarMass_up1sig(self,M,a=1):
        return 10**( np.log10(self.getStellarMass(M)) + self.sigma)
    
    def getStellarMass_down1sig(self,M,a=1):
        return 10**( np.log10(self.getStellarMass(M)) - self.sigma)

    # this is for subhalos, not centrals
    # from equation 2 in "The abundance of not just dark matter haloes"
    def reduceDMmassEQ2(self,M):
	return M* (.65+(M/10**11.4)**.51)/(1+ (M/10**11.4)**.51)

    # from equation 1 in "The abundance of not just dark matter haloes"
    def reduceDMmass(self,M):
        return M* self.b * (self.a/self.b + (M/self.Mt)**self.w ) / (1+ (M/self.Mt)**self.w)

    # Taken from The Chosen Few Figure 3
    # given halo mass, what fraction are luminous?
    def frac_luminous(self,M,a=1):
        from scipy import interpolate
        frac_array = np.array([0, 0.01389365137212395, 0.03227229587377156,  0.0402324326871204, 0.09594819785818098,   0.17944972724833808,0.3697539417495146, 0.7910273213219392, 0.8458747813611525, 1  ])
        log_mass = np.array([7.283067472982055, 7.998634558928878,8.247742935525217,  8.49887928058564, 8.747863613591717, 8.996755635088757,9.247286185871147, 9.497049397234692, 9.746036614975427, 10])
        tck = interpolate.splrep(log_mass,frac_array, k=1)
        if len(M)==0:
            return np.array([])
        return interpolate.splev(np.log10(M),tck)

    # should only input reduced mass into this function
    def getStellarMass_random(self,M,a=1):
        # get random scatter, perturbation
        if type(M)==list or type(M)==np.ndarray:
            perturbation = np.array([random.gauss(0,self.sigma) for i in range(len(M))])
        else:
            perturbation = random.gauss(0,self.sigma)

        return 10**( np.log10(((M/(self.M2))**self.slope2 + (M/(self.M1))**self.slope1)) +perturbation )

    # !!! This function includes all of the new reionization and DM mass suppression stuff!! 
    def generate_stellar_mass_samples(self,halo_mass,N):
        stellar_mass_lists=[]
    	for i in range(N):
            # generate dm_masses
            dm_masses = generate_shm_sample(halo_mass,self.alpha,self.K)
            # reduce dm_masses by certain amount
            dm_masses = self.reduceDMmass(dm_masses)
            lum_chances = self.frac_luminous(dm_masses) # hopefully a list of fractions between 0 and 1
            randnums = np.random.rand(len(lum_chances))
            mask = randnums < lum_chances # chooses which halos will form stars
            dm_masses = dm_masses[mask] # dm masses that host stars

            stellar_mass_lists.append(self.getStellarMass_random(dm_masses,a=1.0))
   	return np.array(stellar_mass_lists)


    # must implement this differently
    def P_at_least_one(self, halo_masses, min_lum):
        frac_gte_one = []
        for halo_mass in halo_masses:
            samples = self.generate_stellar_mass_samples(halo_mass,N=N_iter) # generate list of lists of stellar masses
            lums = samples/ML_ratio
            Ngreater = np.array([np.sum(lum>min_lum) for lum in lums])
            frac_gte_one.append(np.sum(Ngreater>0) / float(len(Ngreater)))
        Pgt1 = np.array(frac_gte_one)
        return Pgt1


    def ngreater(self,halo_masses,min_lum):
        if type(halo_masses)==list or type(halo_masses)==np.ndarray: 
            num_gte_min = []
            spread = []
            for halo_mass in halo_masses:
                samples = self.generate_stellar_mass_samples(halo_mass,N=N_iter) # generate list of lists of stellar masses
                lums = samples/ML_ratio
                Ngreater = np.array([np.sum(lum>min_lum) for lum in lums])
                num_gte_min.append( np.mean(Ngreater) ) # this is taking the average
                spread.append( np.std(Ngreater))
            N_visible = np.array(num_gte_min); N_std = np.array(spread)
        else:
            samples = self.generate_stellar_mass_samples(halo_masses,N=N_iter) # generate list of lists of stellar masses
            lums = samples/ML_ratio
            Ngreater = np.array([np.sum(lum>min_lum) for lum in lums])
            N_visible = np.mean(Ngreater)
            N_std = np.std(Ngreater)
        return N_visible, N_std




"""        
# data points for interpolating
10, 1
9.746036614975427, 0.8458747813611525
9.497049397234692, 0.7910273213219392
9.247286185871147, 0.3697539417495146
8.996755635088757, 0.17944972724833808
8.747863613591717, 0.09594819785818098
8.49887928058564, 0.0402324326871204
8.247742935525217, 0.03227229587377156
7.998634558928878, 0.01389365137212395
7.283067472982055, -0.00041347863421270503
"""
	

class Brook(AMbase):
    def __init__(self,hostSHMF=False,hostAlpha=1.895953, hostK=.004976,reionization=False,z=13): # host SHMF should match plot in the hostfigs folder, and value written in paper
        super(Brook,self).__init__()

        self.alphaAM = 3.1
        self.M0 = 79.6 #63.1

        # SHMF parameters. Slope alpha and constant, K
        self.alpha = 1.802206    #1.818
        self.K = .000778     #.001073
        self.color = 'red'
        self.ls = '--'
        self.label = 'Brook'
        self.isPoisson=True
        self.z=z
        if hostSHMF:
            self.alpha = hostAlpha
            self.K = hostK

        self.reionization = reionization
        if reionization:
            #self.color = 'midnightblue'
            self.label = 'Brook + Reion'
            #self.lum_tck, _, _ = reion.get_lum_fraction_Barber()
            #self.lum_tck, _, _,_,_,_ = reion.get_lum_fraction_caterpillar(mdef='m350',z=self.z)
            self.lum_tck, _, _ = reion.get_lum_fraction_m350()

            self.isPoisson=False
            self.ls = '-'

    def getStellarMass_single(self,M,a=1):
        if M < 10**10.48:
            return (M/(self.M0*10**6))**self.alphaAM
        else:
            model=Moster()
            return model.getStellarMass(M,a)

    def getStellarMass(self,M,a=1):
        if type(M)==list or type(M)==np.ndarray:
            stellar=[-1]*len(M)
            for i,mass in enumerate(M):
                stellar[i]=self.getStellarMass_single(mass)
            return np.array(stellar)
        else:
            return self.getStellarMass_single(M)

    def getStellarMass_random(self,M,a=1):
        return self.getStellarMass(M,a)

                
        

class GarrisonKimmel(AMbase):
    def __init__(self,hostSHMF=False,hostAlpha=1.889667, hostK=0.003906,reionization=False,z=13): # host SHMF should match plot in the hostfigs folder, and value written in paper
        self.epsilon0=-1.777
        self.M10 = 11.514
        self.alphaAM0=-1.92
        self.delta0=3.508
        self.gamma0=0.316
        self.ls = '--'
        self.z=z
        # SHMF parameters. Slope alpha and constant, K
        #self.alpha = 1.867  # these values at z=0 mass
        #self.K = .000703
        self.alpha = 1.181523   #1.828   # these values for Mpeak, bryan & norman
        self.K = .000985      #.001169
        if hostSHMF:
            self.alpha= hostAlpha
            self.K = hostK

        self.color = 'green'
        self.label = 'GK14'
        self.isPoisson=True

        self.reionization = reionization
        if reionization:
            #self.color = 'midnightblue'
            self.label = 'GK14 + Reion'
            #self.lum_tck, _, _ = reion.get_lum_fraction_Barber()
            #self.lum_tck, _, _,_,_,_ = reion.get_lum_fraction_caterpillar(mdef='mpeak',z=self.z)
            self.lum_tck, _, _ = reion.get_lum_fraction_mpeak()

            self.isPoisson=False
            self.ls = '-'

    def f(self,x):
        return -np.log10(10**(x*self.alphaAM())+1)+ self.delta()*(np.log10(1+np.e**x))**self.gamma() / (1+np.e**(10**-x))

    def getStellarMass(self,M,a=1):
        M = M*.71 #MODIFY ELVIS PLOT HERE
        return 10**(self.epsilon()+self.M1() + self.f(np.log10(M) -self.M1()) - self.f(0))

    def getStellarMass_random(self,M,a=1):
        return self.getStellarMass(M,a)

    def M1(self):
        return self.M10
    def epsilon(self):
        return self.epsilon0
    def alphaAM(self):
        return self.alphaAM0
    def delta(self):
        return self.delta0
    def gamma(self):
        return self.gamma0


class GarrisonKimmel16(AMbase):
    # default values for hostK and hostAlpha are not corrrect for peak values!!!
    def __init__(self,sigma=0.8,panel='satellites',variable=False,hostSHMF=False,hostAlpha=1.889667,hostK=0.003906,reionization=False,z=13):
        if panel=='satellites':
            self.alphaAM0 = -(0.14* sigma**2 + .14*sigma+1.79)
        if panel=='field':
            self.alphaAM0 = -(0.24*sigma**2+0.16*sigma+1.99)        

        self.epsilon0=-1.777
        self.M10 = 11.514
        self.delta0=3.508
        self.gamma0=0.316
        self.sigma = sigma
        self.ls = '--'
        self.z=z
        # SHMF parameters. Slope alpha and constant, K
        self.alpha = 1.818523  #1.828
        self.K = .000985 #.001169
        if hostSHMF:
            self.alpha= hostAlpha
            self.K = hostK

        self.color = 'purple'
        self.label = 'GK16'
        self.isPoisson=False

        self.reionization = reionization
        if reionization:
            #self.color = 'midnightblue'
            self.label = 'GK16 + Reion'
            #self.lum_tck, _, _ = reion.get_lum_fraction_Barber()
            #self.lum_tck, _, _,_,_,_ = reion.get_lum_fraction_caterpillar(mdef='mpeak',z=self.z)
            self.lum_tck, _, _ = reion.get_lum_fraction_mpeak()
            self.isPoisson=False
            self.ls = '-'

    def f(self,x):
        return -np.log10(10**(x*self.alphaAM())+1)+ self.delta()*(np.log10(1+np.e**x))**self.gamma() / (1+np.e**(10**-x))

    def getStellarMass(self,M,a=1):
        M = M*.71 #MODIFY ELVIS PLOT HERE
        return 10**(self.epsilon()+self.M1() + self.f(np.log10(M) -self.M1()) - self.f(0))

    def getStellarMass_up1sig(self,M,a=1):
        M = M*.71 #MODIFY ELVIS PLOT HERE
        return 10**( (self.epsilon()+self.M1() + self.f(np.log10(M) -self.M1()) - self.f(0)) +self.sigma)
    
    def getStellarMass_down1sig(self,M,a=1):
        M = M*.71 #MODIFY ELVIS PLOT HERE
        return 10**( (self.epsilon()+self.M1() + self.f(np.log10(M) -self.M1()) - self.f(0)) -self.sigma)

    def getStellarMass_random(self,M,a=1):
        M = M*.71 #MODIFY ELVIS PLOT HERE
        if type(M)==list or type(M)==np.ndarray:
            perturbation = np.array([random.gauss(0,self.sigma) for i in range(len(M))])
        else:
            perturbation = random.gauss(0,self.sigma)
        return 10**( (self.epsilon()+self.M1() + self.f(np.log10(M) -self.M1()) - self.f(0)) +perturbation )

    def M1(self):
        return self.M10
    def epsilon(self):
        return self.epsilon0
    def alphaAM(self):
        return self.alphaAM0
    def delta(self):
        return self.delta0
    def gamma(self):
        return self.gamma0


class GK16_grow(AMbase):
    def __init__(self,plotgamma=-0.2,panel='satellites',variable=False,hostSHMF=False,hostAlpha=1.889667, hostK=0.003906,reionization=False,z=13):
        # plotgamma is from x axis of figure 4 in the paper
        if panel=='satellites':
            self.alphaAM0 = -(0.25*plotgamma**2 - 1.37*plotgamma + 1.69)
            print self.alphaAM0, 'should be close to 1.8'
        if panel=='field':
            self.alphaAM0 = -(0.47*plotgamma**2 - 1.48*plotgamma + 1.81)

        self.plotgamma = plotgamma
        self.epsilon0=-1.777
        self.M10 = 11.514
        self.delta0=3.508
        self.gamma0=0.316
        self.ls = '--'
        self.z=z
        #self.sigma = sigma
        
        # SHMF parameters. Slope alpha and constant, K
        self.alpha = 1.818523  # 1.828
        self.K = .000985  #.001169
        self.color = 'orange'
        self.label = 'GK16_grow'
        self.isPoisson=False
        if hostSHMF:
            self.alpha= hostAlpha
            self.K = hostK
            self.ls = '-'

        self.reionization = reionization
        if reionization:
            #self.color = 'midnightblue'
            self.label = 'GK16_Grow + Reion'
            #self.lum_tck, _, _ = reion.get_lum_fraction_Barber()
            #self.lum_tck, _, _,_,_,_ = reion.get_lum_fraction_caterpillar(mdef='mpeak',z=self.z)
            self.lum_tck, _, _ = reion.get_lum_fraction_mpeak()

    def sigma(self,M):
        M = M*.71 #MODIFY ELVIS PLOT HERE

        # above Mpeak, sigma = 0.2, then it grows linearly.
        # plotgamma negative, so this should really grow!!!
        if type(M)==list or type(M)==np.ndarray:
            sigs = np.array([0.2]*len(M))                                    
            sigs[M<10**11]+=self.plotgamma * (np.log10(M[M<10**11]) - self.M1())
        else:
            if M>10**11:
                return 0.2
            else:
                return self.plotgamma * (np.log10(M) - self.M1())
        return np.array(sigs)


    def f(self,x):
        return -np.log10(10**(x*self.alphaAM())+1)+ self.delta()*(np.log10(1+np.e**x))**self.gamma() / (1+np.e**(10**-x))

    def getStellarMass(self,M,a=1):
        M = M*.71 #MODIFY ELVIS PLOT HERE
        return 10**(self.epsilon()+self.M1() + self.f(np.log10(M) -self.M1()) - self.f(0))

    def getStellarMass_up1sig(self,M,a=1):
        M = M*.71 #MODIFY ELVIS PLOT HERE
        print M, 'values of M'
        return 10**( (self.epsilon()+self.M1() + self.f(np.log10(M) -self.M1()) - self.f(0)) +self.sigma(M))
    
    def getStellarMass_down1sig(self,M,a=1):
        M = M*.71 #MODIFY ELVIS PLOT HERE
        return 10**( (self.epsilon()+self.M1() + self.f(np.log10(M) -self.M1()) - self.f(0)) -self.sigma(M))

    def getStellarMass_random(self,M,a=1):
        M = M*.71 #MODIFY ELVIS PLOT HERE
        if type(M)==list or type(M)==np.ndarray:
            perturbation = np.array([random.gauss(0,self.sigma(mpeak)) for i,mpeak in zip(range(len(M)), M) ])
        else:
            perturbation = random.gauss(0,self.sigma(M))
        return 10**( (self.epsilon()+self.M1() + self.f(np.log10(M) -self.M1()) - self.f(0)) +perturbation )

    def M1(self):
        return self.M10
    def epsilon(self):
        return self.epsilon0
    def alphaAM(self):
        return self.alphaAM0
    def delta(self):
        return self.delta0
    def gamma(self):        
        self.gamma0=0.316
        self.gammaa=1.319
        self.gammaz=0.279

        # SHMF parameters. Slope alpha and constant, K
        self.alpha = 1.818523   # 1.828
        self.K = .000985   #.001169
        self.color = 'pink'
        self.label = 'Behroozi'
        self.isPoisson=True
        self.ls = '--'
        self.z=z

        if hostSHMF:
            self.alpha= hostAlpha
            self.K = hostK

        self.reionization = reionization
        if reionization:
            #self.color = 'midnightblue'
            self.label = 'Behroozi + Reion'
            #self.lum_tck, _, _ = reion.get_lum_fraction_Barber()
            #self.lum_tck, _, _,_,_,_ = reion.get_lum_fraction_caterpillar(mdef='mpeak',z=self.z)
            self.lum_tck, _, _ = reion.get_lum_fraction_mpeak()
            self.isPoisson=False
            self.ls = '-'


    def f(self,x,a):
        return -np.log10(10**(x*self.alphaAM(a))+1)+ self.delta(a)*(np.log10(1+np.e**x))**self.gamma(a) / (1+np.e**(10**-x))

    def getStellarMass(self,M,a):
        return 10**(self.epsilon(a)+self.M1(a) + self.f(np.log10(M) -self.M1(a),a) - self.f(0,a))

    def getStellarMass_random(self,M,a=1):
        return self.getStellarMass(M,a)

    def nu(self,a):
        return np.e**(-4*a**2)
    def M1(self,a):
        z=1/float(a)-1
        return self.M10 + (self.M1a*(a-1)+self.M1z*z) * self.nu(a)
    def epsilon(self,a):
        z=1/float(a)-1
        return self.epsilon0 + (self.epsilona*(a-1)+self.epsilonz*z)*self.nu(a)+self.epsilona2*(a-1)

    def alphaAM(self,a):
        return self.alphaAM0+(self.alphaAMa*(a-1))*self.nu(a)
    def delta(self,a):
        z=1/float(a)-1
        return self.delta0+(self.deltaa*(a-1)+self.deltaz*z)*self.nu(a)
    def gamma(self,a):
        z=1/float(a)-1
        return self.gamma0+(self.gammaa*(a-1)+self.gammaz*z)*self.nu(a)
    

class Moster(AMbase):
    def __init__(self,hostSHMF=False,hostAlpha=1.851473,hostK=0.001445,reionization=False,z=13,t_offset=0):
        self.M10 = 11.590
        self.M11 = 1.195
        self.N10 = .0351
        self.N11 = -0.0247
        self.B10 = 1.376
        self.B11 = -0.826
        self.G10 = 0.608
        self.G11 = 0.329

        # SHMF parameters. Slope alpha and constant, K
        self.alpha = 1.805410   # 1.812
        self.K = .000674     #.000761
        self.color = 'blue'
        self.label = 'Moster'
        self.isPoisson=False
        self.ls = '--'
        self.z=z
        self.t_offset = t_offset

        if hostSHMF:
            self.alpha= hostAlpha
            self.K = hostK

        self.infall_tck = infall_times.get_infall_tck()
        self.lum_tck = None

        self.reionization = reionization
        if reionization:
            #self.color = 'midnightblue'
            self.label = 'Moster + Reion'
            self.lum_tck, _, _ = reion.get_lum_fraction_Barber()
            #self.lum_tck, _, _,_,_,_ = reion.get_lum_fraction_caterpillar(mdef='m200',z=self.z)
            self.isPoisson=False
            self.ls = '-'

    def getStellarMass(self,M, a):
        return M*(2*self.N(a)*( (M/self.M1(a) )**-self.beta(a) + ( M/self.M1(a) )**self.gamma(a) )**-1)

    def M1(self,a):
        return 10**(self.M10+self.M11*(1-a))

    def N(self,a):
        return self.N10 + self.N11*(1-a)

    def beta(self,a):
        return self.B10 + self.B11*(1-a)

    def gamma(self,a):
        return self.G10 + self.G11*(1-a)


    # should only input reduced mass into this function
    def getStellarMass_random(self,M,a=1):
        # get random scatter, perturbation
        if type(M)==list or type(M)==np.ndarray:
            perturbation = np.array([random.gauss(0,self.sigma) for i in range(len(M))])
        else:
            perturbation = random.gauss(0,self.sigma)
        return 10**( ((M/(self.M2))**self.slope2 + (M/(self.M1))**self.slope1) + perturbation)


    ### NEW FUNCTION TO GENERATE ###
    # why is it so slow????
    # is it the append?
    # is it get stellar mass function?
    # generate dm_masses?
    # scale_factors??
    def generate_stellar_mass_samples(self,halo_mass,N):
        #import time
        #ts=time.time()
        #tlum=0; tinfall=0
        stellar_mass_lists=[]
    	for i in range(N):
            # generate dm_masses
            dm_masses = generate_shm_sample(halo_mass,self.alpha,self.K)
            if self.reionization:  # reduce the sample
                #if i%100==0:
                #    print i
                #    print tlum, tinfall, 'lum and infall time'
                #t0=time.time()
                lum_chances = reion.frac_luminous(dm_masses, self.lum_tck) # hopefully a list of fractions between 0 and 1
                #t1=time.time()
                #tlum+=t1-t0
                randnums = np.random.rand(len(lum_chances))
                mask = randnums < lum_chances # chooses which halos will form stars
                dm_masses = dm_masses[mask] # dm masses that host stars
            
            #t0=time.time()
            scale_factors = infall_times.random_infall_scale(self.infall_tck, len(dm_masses), offset=self.t_offset) ## the luminous ones should be shifted towards earlier times
            #t1=time.time()
            #tinfall +=t1-t0
            
            stellar_mass_lists.append(self.getStellarMass(dm_masses,scale_factors))
        #tf=time.time()
        #print tf-ts, 'total time'
        #print tlum, tinfall, 'lum and infall time'
   	return np.array(stellar_mass_lists)


def plotMstar_v_Mhalo():
    a = 1.0
    M = 10**np.arange(7.5,12.5,.25)
    AMmodels=[Moster(), Behroozi(), GarrisonKimmel(), Brook(), GarrisonKimmel16(), Sawala()]
    for model in AMmodels:
        mstar = model.getStellarMass(M,a)
        #print mstar, model_label
        plt.plot(M, mstar, label=model.label,lw=linewidth, color=model.color)
        if isinstance(model, GarrisonKimmel16) or isinstance(model, Sawala):
            mstar_up = model.getStellarMass_up1sig(M,a)
            mstar_down = model.getStellarMass_down1sig(M,a)
            plt.fill_between(M,mstar_up,mstar_down, facecolor=model.color, alpha=0.2)
            #mstar_rand = model.getStellarMass_random(M,a)
            #plt.scatter(M,mstar_rand,color=model.color)

    plt.legend(loc='lower right',fontsize=legend_size,frameon=False)
    plt.yscale('log')
    plt.xscale('log')
    plt.xlabel('$M_{\mathrm{halo}} \mathrm{[M_\odot]}$',fontsize=label_font)
    plt.ylabel('$M_* \mathrm{[M_\odot]}$',fontsize=label_font)
    plt.ylim((10**3,10**11))
    plt.xlim((10**7.5,10**12.5))
    plt.tick_params(axis='both', which='major', labelsize=tick_size)
    plt.gcf().subplots_adjust(bottom=0.15)
    plt.savefig('smhm_relation') # compare to figure 5 in paper
    plt.close()


def plotBentbyBaryons():
    #f,ax = plt.subplots(ncols=1)
    a = 1.0
    M = 10**np.arange(7.5,np.log10(5e10),.15)
    AMmodels= [Moster(), Sawala()]
    model_labels = ['Moster','Sawala']
    colors = ['blue','red']
    for model,model_label,color in zip(AMmodels,model_labels,colors):
        mstar = model.getStellarMass(M,a)
        #print mstar, model_label
        plt.plot(M, mstar, label=model_label,linewidth=5, color=color)

    plt.legend(loc='lower right',frameon=False,fontsize=26)
    plt.yscale('log')
    plt.xscale('log')
    plt.xlabel('$M_{halo} \mathrm{[M_\odot]}$',fontsize=30)
    plt.ylabel('$M_* \mathrm{[M_\odot]}$',fontsize=30)
    plt.ylim((10**3,10**8))
    plt.xlim((10**8,5e10))
    plt.tick_params(axis='both', which='major', labelsize=21)
    #plt.gca().tight_layout()	
    plt.gcf().subplots_adjust(bottom=0.17)
    plt.gcf().subplots_adjust(left=0.17)
    plt.savefig('BentByBaryons') # compare to figure 4 in bent by baryons
    plt.close()



#plotMstar_v_Mhalo()
