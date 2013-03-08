import numpy as np
from scipy.integrate import quad
import matplotlib.pyplot as plt
from scipy import interpolate

class MassFunction:
    
    def __init__(self, h=0.7, Omega_m=.27, sigma8=.809, n=0.95, rho_c=1.35972365653e11, delta_c=1.686):
        """
        @param h: hubble parameter.
        @param omega_m: current matter fraction of universe.
        @param sigma8: Value of variance of density field with spherical smoothing at 8 Mpc/h
        @param n: spectral scalar index for primordial power spectrum
        @param delta_c: amplitude of perturbation at collapse using linear theory.
        @param rho_c: critical density of universe in Msun/Mpc
        """
        self.h = h
        self.Omega_m = Omega_m
        self.n = n
        self.delta_c = delta_c
        self.rho_m = rho_c*Omega_m/self.h**2
        self.sigma8 = sigma8
        
    # Initial Power Spectrum
    def Po(self,k):
        return k**self.n

    ### Transfer Function from Bond and Efstathiou
    def T(self,k,a=6.4,b=3.,c=1.7, nu=1.13):
        Gamma = self.Omega_m*self.h
        q = k/Gamma
        return (1+(a*q+(b*q)**1.5+(c*q)**2)**nu)**(-1./nu)

    def P(self,k):
        return self.T(k)**2*self.Po(k)
    
    def b(self,z):
        return 1./(1.+z)

    def W(self,k,M):
        R = (3*M/(4*np.pi*self.rho_m))**(1/3.)
        u = k*R
        return 3*(np.sin(u)-u*np.cos(u))/u**3

    def sigma_func(self,k,M):
        return k**2*self.P(k)*self.W(k,M)**2

    # Returns normalized sigma, not sigma^2
    def Sigma(self,M, z=0):
        if type(M) is np.ndarray:
            results = np.array([0.0]*M.size,dtype='float')
            for i in range(0,M.size):
                results[i] = quad(self.sigma_func, 0.0, np.inf, args=(M[i]),limit=2000)[0]
        else:
            results = quad(self.sigma_func, 0.0, np.inf, args=(M),limit=2000)[0]
        M8 = 4./3*np.pi*8**3*self.rho_m
        s8 = self.b(z)/(2*np.pi**2)*quad(self.sigma_func, 0.0, np.inf, args=(M8),limit=2000)[0]
        results *= self.sigma8**2/s8
        return np.sqrt(self.b(z)*results/(2*np.pi**2))

    def N(self,M,z=0):
        s = self.Sigma(M,z)
        f = np.sqrt(2/np.pi)*self.delta_c/s*np.exp(-self.delta_c**2/(2*s**2))
        return f*self.rho_m/M

    def PSNofM(self,M, z=0):
        """
        Use Press Schecter Formalism to compute the mass function N(m),
        fractional number density of halos in mass range (m,m+dm). Should
        integrate to total number of halos/volume.
        @param M: x-range of masses in Msun/h
        @param z: red shift at which to evaluate mass function
        @return: [x-axis in log10(mass), y-axis in log10(N(m)), xlabel, ylabel]
        """
        s = self.Sigma(M,z)
        tck = interpolate.splrep(M,s)
        dsdm = -interpolate.splev(M, tck, der=1)
        return self.rho_m*self.delta_c/M*np.sqrt(2/np.pi)*np.exp(-self.delta_c**2/(2*s**2))/s**2*dsdm

    def PSdNdLogM(self,M, z=0):
        """
        Use Press Schecter Formalism to compute the mass function
        dn(m)/dlog(M).
        @param M: x-range of masses in Msun/h
        @param z: red shift at which to evaluate mass function
        @return: [x-axis in log10(mass), y-axis in log10(N(m)), xlabel, ylabel]
        """
        s = self.Sigma(M,z)
        tck = interpolate.splrep(np.log10(M),s)
        dsdlogm = -interpolate.splev(np.log10(M), tck, der=1)
        return self.rho_m*self.delta_c/M*np.sqrt(2/np.pi)*np.exp(-self.delta_c**2/(2*s**2))/s**2*dsdlogm
    
    def STNofM(self,M, z=0):
        """
        Use Sheth & Tormen Formalism to compute the mass function N(m),
        fractional number density of halos in mass range (m,m+dm).
        @param M: x-range of masses in Msun/h
        @param z: red shift at which to evaluate mass function
        @return: [x-axis in log10(mass), y-axis in log10(N(m)), xlabel, ylabel]
        """
        s = self.Sigma(M,z)
        tck = interpolate.splrep(M,s)
        dsdm = -interpolate.splev(M, tck, der=1)
        A = 0.3222
        a = 0.707
        p = 0.3
        f = A*np.sqrt(2*a/np.pi)*(1+(s**2/a/self.delta_c**2)**p)*self.delta_c/s*np.exp(-.5*a*self.delta_c**2/s**2)*dsdm
        
        return f/s*self.rho_m/M

    def STdNdLogM(self,M, z=0):
        """
        Use Sheth & Tormen Formalism to compute the mass function
        dn(m)/dlog(M).
        @param M: x-range of masses in Msun/h
        @param z: red shift at which to evaluate mass function
        @return: [x-axis in log10(mass), y-axis in log10(N(m)), xlabel, ylabel]
        """
        s = self.Sigma(M,z)
        tck = interpolate.splrep(np.log10(M),s)
        dsdlogm = -interpolate.splev(np.log10(M), tck, der=1)
        A = 0.3222
        a = 0.707
        p = 0.3
        f = A*np.sqrt(2*a/np.pi)*(1+(s**2/a/self.delta_c**2)**p)*self.delta_c/s*np.exp(-.5*a*self.delta_c**2/s**2)
        
        return f/s*self.rho_m/M*dsdlogm

    def NofM(self,masses, numbins, boxsize):
        """
        Produce mass function data for N(m). Fractional number density of halos in
        mass range (m,m+dm). Integrates to total number of halos/volume.
        @param masses: list of all masses in Msun/h
        @param numbins: number of bins to use in histogram
        @param boxsize: Size of box in MPc/h.
        @return: [x-axis in log10(mass), y-axis in log10(N(m)), xlabel, ylabel]
        """
        logmasses = np.log10(masses)
        hist, r_array = np.histogram(logmasses, numbins)
        dlogM = r_array[1]-r_array[0]
        x_array = r_array[1:] - .5*dlogM
        dM = 10.**r_array[1:]-10.**r_array[0:numbins] #Mass size of bins in non-log space.
        volume = np.float(boxsize**3) # in MPc^3
        return [x_array, np.log10(hist/volume/dM)]
    
    def dNdLogM(self,masses, numbins, boxsize):
        """
        Produce mass function data for dn(m)/dlog10(M). n(m) is number of
        halos with mass < m, divided by box volume.
        mass range (m,m+dm). Integrates to total number of halos.
        @param masses: list of all masses in Msun/h
        @param numbins: number of bins to use in histogram
        @param boxsize: Size of box in MPc
        @return: [x-axis in log10(mass), y-axis in log10(dn/dlog10(M)), xlabel, ylabel]
        """
        logmasses = np.log10(masses)
        ### Make n(M) number of halos with mass < M.
        hist, r_array = np.histogram(logmasses, numbins)
        dlogM = r_array[1]-r_array[0]
        x_array = r_array[1:] - .5*dlogM
        dM = 10.**r_array[1:]-10.**r_array[0:numbins] #Mass size of bins in non-log space.
        dn_dlogM = np.log(10)*hist*10**x_array/dM
        volume = np.float(boxsize**3) # in MPc^3
        return [x_array, np.log10(dn_dlogM/volume)]
