import numpy as np
import abundance_matching as am
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PlotParams import *

# How does the mean number of satellites above 10^4 change for Moster, GK16, and Brook for each of the following changes:


"""
1) reionization model: at 4 different redshifts
2) Moster - shift infall distribution by 0.5 Gyrs
3) GK16 - increase the sigma by a lot
4) change M* of field halos by 10% up and down
5) change slope of SHMF - force it to -1.9, -1.85 and get best fit K
6) magnitude of SHMF (no, why would it be uncertain?)
7) Change from Field SHMF to Host SHMF  (must use rmaxcut=False)
8) take SHMF of half of the halos and the other half. like jackknife?
9) host halo mass +/- 10%

"""



N_iter = 5000
StellarMasses = 10**6*np.array([0.14,0.54,0.77,0.82,1.3,1.4,1.6,2.1,2.7,3.5,3.5,4.5,6.0,6.4,6.6,7.8,8.3,16,16,17,19,37,44,47,51,52,62,76,100,100,190,270])
# my metric is the sum of the 5 largest systems with Moster
star_masses = StellarMasses[-5:]


def get_sum_fields(model, min_lum, star_masses):
    dm_masses = model.stellar_to_halo_mass(star_masses,a=1.0)
    if model.isPoisson:
        mean = model.get_field_total(dm_masses,min_lum)
    else:
        samples = model.get_field_distr(dm_masses,min_lum,N=N_iter)
        mean = np.mean(samples)
    return mean



# 1) reionization model: at 4 different redshifts
def change_reion():
    model=am.Moster(reionization=True,z=13)
    z13 = get_sum_fields(model,10**4,star_masses)
    model=am.Moster(reionization=True,z=12)
    z12 = get_sum_fields(model,10**4,star_masses)
    model=am.Moster(reionization=True,z=11)
    z11 = get_sum_fields(model,10**4,star_masses)
    model=am.Moster(reionization=True,z=9)
    z9 = get_sum_fields(model,10**4,star_masses)

    del12 = z12/z13 -1
    del11 = z11/z13 - 1
    del9 = z9/z13 - 1
    return del12*100, del11*100, del9*100
# 10.069396252602369, 18.065579458709213, 31.179736294240112
"""
10.5652 z=13
11.4308 z=12
11.9286 z=11
12.9832 z=9
2.4768   z=13 with 15km/s and 50km/s cut-offs instead of 10 and 25
"""


#2) Moster - shift infall distribution by 0.5 Gyrs
def change_infall():
    model=am.Moster(reionization=True, t_offset = 0.0)
    norm = get_sum_fields(model,10**4,star_masses)

    model=am.Moster(reionization=True, t_offset = 0.5)
    earlier = get_sum_fields(model,10**4,star_masses)

    model=am.Moster(reionization=True, t_offset = -0.5)
    later = get_sum_fields(model,10**4,star_masses)

    model=am.Moster(reionization=True, t_offset = 1.0)
    earlier1 = get_sum_fields(model,10**4,star_masses)

    model=am.Moster(reionization=True, t_offset = -1)
    later1 = get_sum_fields(model,10**4,star_masses)

    del_earlier= earlier/norm -1
    del_later= later/norm -1
    del_earlier1= earlier1/norm -1
    del_later1= later1/norm -1
    return del_earlier*100, del_later*100, del_earlier1*100, del_later1*100
# 1.9029528578829114, -1.804524261785545, 3.9060611293386227, -3.6746675876359935



#3) GK16 - increase the sigma by a lot    
def change_scatter():
    model = am.GarrisonKimmel16(reionization=True) 
    norm = get_sum_fields(model,10**4,star_masses)
    
    model = am.GarrisonKimmel16(reionization=True, sigma=0.4)
    sigdown = get_sum_fields(model,10**4,star_masses)
    
    model = am.GarrisonKimmel16(reionization=True, sigma=1.2)
    sigup = get_sum_fields(model,10**4,star_masses)
    
    model = am.GK16_grow(reionization=True, plotgamma=-.2)
    gk162 = get_sum_fields(model,10**4,star_masses)
    
    model = am.GK16_grow(reionization=True, plotgamma=-.5)
    gk165 = get_sum_fields(model,10**4,star_masses)

    del_sigdown = sigdown/norm - 1
    del_sigup = sigup/norm - 1
    del_gk162 = gk162/norm - 1
    del_gk165 = gk165/norm - 1
    return del_sigdown*100, del_sigup*100, del_gk162*100, del_gk165*100
    # 8.3766425593837965, -8.2617591387025726, 5.0373644649167826, -18.566254909898571


# 4) change M* of field halos by +/- 10%
def change_mstar():
    model=am.Moster(reionization=True)
    base = get_sum_fields(model,10**4,star_masses)
    down = get_sum_fields(model,10**4,star_masses*.9)
    up = get_sum_fields(model,10**4,star_masses*1.1)
    delup = up/base - 1.
    deldown = down/base - 1
    return delup*100, deldown*100
# delup: 4.4271149938343468 deldown: - 4.7953175747260102




## std error = 0.02 for alpha = -1.8054104588, K = 0.000673927822796


## first 16
# alpha = -1.91225193546
# K = 0.00481748488821
# stderr = 0.0888785222922


## second 17
# alpha = -1.77700887695
# K = 0.000398018832465
# 0.0287554270974 std err




# 5) use the uncertainty on the fit on alpha and K and get the extremes
def shmf_uncertainty():
    model=am.Moster(reionization=True) # self.alpha = 1.812, sel
    norm = get_sum_fields(model,10**4,star_masses)

    # first16
    model=am.Moster(reionization=True,hostSHMF=True, hostAlpha= 1.91225193546, hostK= 0.00481748488821 ) 
    first16 = get_sum_fields(model,10**4,star_masses)

    # second17
    model=am.Moster(reionization=True,hostSHMF=True, hostAlpha=1.77700887695, hostK= 0.000398018832465  ) 
    second17 = get_sum_fields(model,10**4,star_masses)

    del_first16= first16/norm -1
    del_second17= second17/norm -1
    return del_first16*100, del_second17*100
    # (-25.844947584753818, 10.982159819007631)





#7) Change from Field SHMF to Host SHMF  (must use rmaxcut=False)
def field_v_host():
    model=am.Moster(reionization=True) # self.alpha = 1.812, self.K = .000761
    field = get_sum_fields(model,10**4,star_masses)
    model=am.Moster(reionization=True,hostSHMF=True)
    host = get_sum_fields(model,10**4,star_masses)

    # 7.0 as lower limit as all 33 halos
    model=am.Moster(reionization=True,hostSHMF=True, hostAlpha=1.811243 , hostK= 0.000753) #1.811243 & 0.000753
    field2 = get_sum_fields(model,10**4,star_masses)

    # 7.5 as lower limit   1.805410 & 0.000674
    model=am.Moster(reionization=True,hostSHMF=True, hostAlpha=1.805410, hostK= 0.000674) 
    field3 = get_sum_fields(model,10**4,star_masses)

    del_host= host/field -1
    del_field2= field2/field -1
    del_field3= field3/field -1
    return del_host*100, del_field2*100, del_field3*100
    # -18.3101763625654, 0.97363223727522374, 2.2608364228543731
    
    
    
#9) host halo mass +/- 10%
def change_mhost():
    model=am.Moster(reionization=True) # self.alpha = 1.812, self.K = .000761
    base = get_sum_fields(model,10**4,star_masses)

    dm_masses =  model.stellar_to_halo_mass(star_masses,a=1.0)
    star_masses_up = model.getStellarMass(dm_masses*1.1,a=1)
    star_masses_down = model.getStellarMass(dm_masses*0.9,a=1)

    down = get_sum_fields(model,10**4,star_masses_down)
    up = get_sum_fields(model,10**4,star_masses_up)
    delup = up/base - 1.
    deldown = down/base - 1
    return delup*100, deldown*100
  # (11.196523485517407, -10.275454907460302)  # need higher N_iter. shouldn't it be 10% exactly?


