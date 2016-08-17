import numpy as np
import abundance_matching as am
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PlotParams import *


min_lum = 10**4
DwarfNames = ['Leo T', 'KKR 3', 'Phoenix', 'KKH 86','Antlia','KKR 25','Aquarius','DDO 113','ESO 294- G 010','Sagittarius dIrr',
'ESO 410- G 005','KKH 98','Leo A','GR 8','Pegasus dIrr','UGC 9128','UGC 4879','UGCA 86','DDO 99','UKS 2323-326',
'UGC 8508','NGC 4163','Sextans A','DDO 125','DDO 190','Sextans B','IC 3104','NGC 3109','NGC 6822','IC 1613',
'IC 4662','IC 5152']

StellarMasses = 10**6*np.array([0.14,0.54,0.77,0.82,1.3,1.4,1.6,2.1,2.7,3.5,3.5,4.5,6.0,6.4,6.6,7.8,8.3,16,16,17,19,37,44,47,51,52,62,76,100,100,190,270])

N_iter = 5000  # make 20000 for the real thing. 5000 min acceptable

def get_ngreater(mstar,min_lum,model):
        dm_masses = model.stellar_to_halo_mass(mstar,a=1.0)
        return model.ngreater(dm_masses,min_lum)
        
def get_P_at_least_one(mstar,min_lum,model):
        dm_masses = model.stellar_to_halo_mass(mstar,a=1.0)
        return model.P_at_least_one(dm_masses,min_lum)


def generate_latex_table(Names,Mstar,reion=True):
    with open('/nfs/blank/h4231/gdooley/Dropbox/DwarfsOfDwarfs/Paper/Table1.tex' ,'w') as f:
        f.write("\\begin{table*}\n")
        f.write("\\tablewidth{0.97\\textwidth}\n")
        f.write("\\centering\n")
        f.write("\\caption{Local Group dIrrs}\n")
        f.write("\\label{table:results}\n")
        f.write("\\begin{tabular}{lccccccc} \n")
        f.write("\\hline\n")
        f.write("\\hline\n")
        f.write("Name &  $M_* [10^6 \\, \\mathrm{M_\odot}]$ & N Moster & N Behroozi & N GK14 & N Brook & N GK16 & Sawala  \\\\\n" )
        Ngreater = []; Pgt1=[]
        for model in [am.Moster(reionization=reion),am.Behroozi(reionization=reion), am.GarrisonKimmel(reionization=reion), am.Brook(reionization=reion), am.GarrisonKimmel16(reionization=reion), am.Sawala()]:
            print model.label, 'on this model'
            Ngreater.append(get_ngreater(StellarMasses,min_lum,model))
            Pgt1.append(get_P_at_least_one(StellarMasses,min_lum,model))

        Nsubs = np.array(Ngreater)[:,0]
        Nstd = np.array(Ngreater)[:,1]

        for name, mstar, nsubs, pgt1, nstd in zip(Names,Mstar, Nsubs.T, np.array(Pgt1).T, Nstd.T):
            f.write("\\hline\n")
            f.write("%s & %.2f & %.2f $\pm$ %.2f & %.2f $\pm$ %.2f & %.2f $\pm$ %.2f & %.2f $\pm$ %.2f & %.2f $\pm$ %.2f & %.2f $\pm$ %.2f \\\\\n" %(name,mstar*1e-6,nsubs[0],nstd[0],  nsubs[1],nstd[1], nsubs[2],nstd[2], nsubs[3],nstd[3], nsubs[4],nstd[4], nsubs[5],nstd[5] ))
        f.write("\\hline\n")
        f.write("\\hline\n")
        f.write("\\end{tabular}\n")
        f.write("\\tablecomments{Expected number of satellites with $M^* > 10^4$\msun to be found around known local group Dwarf Irregular and Dwarf Spheroidal galaxies for various AM models and minimum detection thresholds. Also shown is the probability of finding at least one satellite around each galaxy for various AM models.}\n")
        f.write("\\end{table*}\n")



def add_sum_distr_ax(ax,min_lum,subset=True,fixedDM=False, re=False):
    StellarMasses = 10**6*np.array([0.14,0.54,0.77,0.82,1.3,1.4,1.6,2.1,2.7,3.5,3.5,4.5,6.0,6.4,6.6,7.8,8.3,16,16,17,19,37,44,47,51,52,62,76,100,100,190,270])
    if subset:
        StellarMasses = StellarMasses[-5:]
    import scipy.misc
    # formerly had reionization together with no reionization am.Moster(reionization=True)
    for model in [am.Sawala(),am.GarrisonKimmel16(reionization=re), am.Moster(reionization=re), am.GarrisonKimmel(reionization=re), am.Brook(reionization=re)]:
        if fixedDM:  # fix the dm mass to just moster
            tmpmodel = am.Moster()
            dm_masses = tmpmodel.stellar_to_halo_mass(StellarMasses,a=1.0)
            if subset:
                print np.log10(dm_masses), 'moster masses'
                print model.stellar_to_halo_mass(StellarMasses,a=1.0), model.label, 'masses'
        else:
            dm_masses = model.stellar_to_halo_mass(StellarMasses,a=1.0)
        print model.label
        if model.isPoisson:
            mean = model.get_field_total(dm_masses,min_lum)
            nocc = np.arange(0,mean*3+1)
            prob = (mean**nocc * np.e**-mean /scipy.misc.factorial(nocc))
            tmp=np.repeat(nocc,2)
        else:
        # for monte carlo methods, need to do 10,000 instances, looping over each halo and collecting all 10,000 numbers
        # then loop over all halos, and keep adding to the running total list.
        # then make histogram of the final 10,000 long distribution, normalized to 1
        # plot each in their own color.
            samples = model.get_field_distr(dm_masses,min_lum,N=N_iter)
            distr,nocc = np.histogram(samples,bins=np.arange(min(samples),max(samples)+2))
            prob = distr/float(len(samples))
            tmp=np.repeat(nocc[0:-1],2)  # last value of bins is not inclusive
        ax.plot(np.append(tmp[1:],tmp[-1]+1), np.repeat(prob,2),linewidth=linewidth,color=model.color,label=model.label)

# want to generate mean and variance of the sum of all distributions
#for non monte carlo ones, just get the sum total, and plot poisson distribution
def plot_sum_distr(subset=True):
    plt.figure()
    ax = plt.subplot(111)
    add_sum_distr_ax(ax,min_lum=10**3,subset=True)

    plt.ylabel('Probability',fontsize=label_font)
    plt.xlabel('N satellites total',fontsize=label_font)
    plt.legend(loc='upper right',frameon=False,fontsize=legend_size)
    plt.tick_params(axis='both', which='major', labelsize=tick_size)
    plt.gcf().subplots_adjust(bottom=0.15)
    
    if subset:
        plt.savefig('largest10_field_sats')
    else:
        plt.savefig('total_field_sats')



# determine dark matter mass instead of stellar mass from just Moster.
# what are those masses? 
# refine field halos used for SHMF to those field halos that most closely mirror these 10 largest
# make a version of the plot below with the mass fixed to just that of Moster

def plot_sum_distr_3panel(subset=True,fixedDM=False,re=False):
    #w,h = plt.figaspect(1.5)
    f=plt.figure()
    h=f.get_figheight()
    w=f.get_figwidth()
    plt.close()

    fig = plt.figure(figsize=(w,h*1.5))
    ax1 =fig.add_subplot(3,1,1)
    ax2 =fig.add_subplot(3,1,2)
    ax3 =fig.add_subplot(3,1,3)
    #f,(ax1,ax2,ax3) = plt.subplots(nrows=3)

    add_sum_distr_ax(ax1,min_lum=10**3,subset=subset,fixedDM=fixedDM,re=re)
    add_sum_distr_ax(ax2,min_lum=10**4,subset=subset,fixedDM=fixedDM,re=re)
    add_sum_distr_ax(ax3,min_lum=10**5,subset=subset,fixedDM=fixedDM,re=re)
    
    ax1.set_ylabel('Probability',fontsize=label_font)
    ax2.set_ylabel('Probability',fontsize=label_font)
    ax3.set_ylabel('Probability',fontsize=label_font)
    ax3.set_xlabel('N satellites total',fontsize=label_font)
    ax1.legend(loc='upper right',frameon=False,fontsize=legend_size-3)
    ax1.tick_params(axis='both', which='major', labelsize=tick_size)
    ax2.tick_params(axis='both', which='major', labelsize=tick_size)
    ax3.tick_params(axis='both', which='major', labelsize=tick_size)
    #plt.gcf().subplots_adjust(bottom=0.15)
    plt.gcf().subplots_adjust(bottom=0.10)
    #plt.gcf().subplots_adjust(top=0.05)

    ax1.text(.25,.8,'$M_* > 10^3 \, \mathrm{M_\odot}$',transform=ax1.transAxes,fontsize=legend_size)
    ax2.text(.25,.8,'$M_* > 10^4 \, \mathrm{M_\odot}$',transform=ax2.transAxes,fontsize=legend_size)
    ax3.text(.25,.8,'$M_* > 10^5 \, \mathrm{M_\odot}$',transform=ax3.transAxes,fontsize=legend_size)
    
    extra=''
    if fixedDM:
        extra = 'fixed_DM_mass'
    if re:
        extra+='reionization'
    if subset:
        plt.savefig('largest10_field_sats_3panel'+extra)
    else:
        plt.savefig('total_field_sats_3panel'+extra)



def plot_sum_distr_2panel(subset=True,fixedDM=False):
    f,(ax1,ax2) = plt.subplots(nrows=2)

    add_sum_distr_ax(ax1,min_lum=10**3,subset=subset,fixedDM=fixedDM)
    add_sum_distr_ax(ax2,min_lum=10**4,subset=subset,fixedDM=fixedDM)
    
    ax1.set_ylabel('Probability',fontsize=label_font)
    ax2.set_ylabel('Probability',fontsize=label_font)
    ax2.set_xlabel('N satellites total',fontsize=label_font)
    ax1.legend(loc='upper right',frameon=False,fontsize=legend_size)
    ax1.tick_params(axis='both', which='major', labelsize=tick_size)
    ax2.tick_params(axis='both', which='major', labelsize=tick_size)
    plt.gcf().subplots_adjust(bottom=0.15)

    ax1.text(.45,.85,'$M_* > 10^3 \, \mathrm{L_\odot}$',transform=ax1.transAxes,fontsize=legend_size)
    ax2.text(.45,.85,'$M_^* > 10^4 \, \mathrm{L_\odot}$',transform=ax2.transAxes,fontsize=legend_size)
    
    extra=''
    if fixedDM:
        extra = 'fixed_DM_mass'
    if subset:
        plt.savefig('largest10_field_sats_2panel'+extra)
    else:
        plt.savefig('total_field_sats_2panel'+extra)


#plot_sum_distr()
#plot_sum_distr_2panel(subset=True)
plot_sum_distr_3panel(subset=True,re=True)
#plot_sum_distr_2panel(subset=False)
#plot_sum_distr_2panel(subset=True,fixedDM=True)  # fixed DM doesn't really work, it must be changed to 350 or 200
#plot_sum_distr_2panel(subset=False,fixedDM=True)

#generate_latex_table(DwarfNames, StellarMasses)




"""
NGC 6822 & 100.00 & 3.88 & 60.94 & 3.24 & 2.27 & 11.41 \\
\hline
IC 1613 & 100.00 & 3.88 & 60.94 & 3.24 & 2.27 & 11.43 \\
\hline
IC 4662 & 190.00 & 5.12 & 97.71 & 4.43 & 3.01 & 15.40 \\
\hline
IC 5152 & 270.00 & 5.96 & 125.09 & 5.18 & 3.51 & 17.86 \\
"""
