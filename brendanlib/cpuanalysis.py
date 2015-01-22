import string
import numpy as np
import matplotlib.pylab as plt
import matplotlib.colors
import matplotlib.font_manager
from matplotlib.ticker import MaxNLocator
from matplotlib.font_manager import FontProperties
import haloutils as htils

cm = plt.cm.get_cmap('RdYlBu')

snapnum = 63
fontsize = 12
TimeStepMod=100
subtags = ["total", "treegrav", "pmgrav", "voronoi", "hydro", "domain", "sph", "i/o"]
label_subtags={"total": "total", 
               "treegrav": "tree gravity", 
               "pmgrav": "PM gravity", 
               "voronoi": "mesh", 
               "hydro": "hydro fluxes", 
               "domain": "domain decomp.", 
               "sph": "SPH", 
               "i/o": "I/O"}

tcolors = ["Azure", "Blue", "Green", "Red", "SpringGreen", "DarkTurquoise", "Red", "Black", "FireBrick", "Purple", "OrangeRed", "SeaGreen", "Navy", "Yellow", "Tomato", "Violet", "Tan", "Olive", "Gold"] 

fig1 = plt.figure(figsize=(16,8))
fig2 = plt.figure(figsize=(16,8))
figi = 0

hpaths = htils.get_paper_paths_lx(14)
#for res in reslist:
#    for pad in padlist:
#        for nvir in nvirlist:
for hpath in hpaths:
            figi += 1
            ax1 = fig1.add_subplot(3,4,figi)
            ax2 = fig2.add_subplot(3,4,figi)
            #tmppath = basedir + '/' + res + '/' + pad + '/' + nvir + '/outputs'
            tmppath = hpath + '/outputs'
            print tmppath
            titlestr = hpath.split("/")[-2]
            #titlestr = pad + ", " + nvir
            cputxt = tmppath + '/cpu.txt'
            timingstxt = tmppath + '/timings.txt'
            energytxt = tmppath + '/energy.txt'
            infotxt = tmppath + '/info.txt'
            
            try:
                with open(cputxt): pass
                value_total = [] 
                value_rel = [] 
                found_subtags = []
                TimeStep=-1
                logfile = open(cputxt, "r").readlines()
                maintag = "Step"
                #loop over file
                for line in logfile:
                    parsed_line = line.split()	
                    if (parsed_line!=[]):
                        if (parsed_line[0].find(maintag)>-1):
                            TimeStep=TimeStep+1
                            if TimeStep%TimeStepMod == 0:
                                ReadFlag=1
                            else:
                                ReadFlag=0
                        else:
                            if (ReadFlag):
                                s = parsed_line[1]
                                try:
                                    i = float(s)
                                except ValueError:
                                    parsed_line[0]=parsed_line[0]+" "+parsed_line[1]
                                    parsed_line[1]=parsed_line[2]
                                    parsed_line[2]=parsed_line[3]    
                                for word in parsed_line:
                                    if (word in subtags):
                                        if (TimeStep==TimeStepMod):
                                            found_subtags.append(word)
                                        value_total.append(float(parsed_line[1]))
                                        value_rel.append(float(string.translate(parsed_line[2],None,"%")))
                
                
                #construct arrays
                number_of_found_subtags=len(found_subtags)
                steps=len(value_total)/number_of_found_subtags
                All_value_total=np.array(value_total).reshape(steps,number_of_found_subtags)
                All_value_rel=np.array(value_rel).reshape(steps,number_of_found_subtags)
                All_value_cum_rel=np.zeros([steps,number_of_found_subtags])
                for i in range(1,number_of_found_subtags):
                    All_value_cum_rel[:,i]=100.0*All_value_total[:,i].cumsum()/All_value_total[:,0].cumsum()

                # relative cpu usage over time
                
                data1=np.zeros(steps)
                data2=np.zeros(steps)
                times=np.arange(steps)
                for i in range(1,number_of_found_subtags):
                    data2=data2+All_value_rel[:,i]
                    x=TimeStepMod*np.concatenate((times[:],[times[steps-1]],times[::-1],[times[0]]))/1000 
                    y=np.concatenate((data1[:],[data2[steps-1]],data2[::-1],[data1[0]]))
                    ax1.fill(x,y,label=label_subtags[found_subtags[i]],color=tcolors[subtags.index(found_subtags[i])])
                    data1=data2

                ax1.axis([x.min(),0.99*x.max(),0,100.])
                if figi == 1:
                    handles, labels = ax1.get_legend_handles_labels()
                # reverse the order
                    prop = matplotlib.font_manager.FontProperties(size=fontsize)
                    #ax1.legend(handles[::-1], labels[::-1],loc="lower right",prop=prop, bbox_to_anchor=(0., -5.0, 1., .102)).draggable()
                    fig1.legend(handles[::-1], labels[::-1], 'upper center',prop=prop,ncol=number_of_found_subtags-1)
                #plt.savefig(tmppath+"/times_rel_diff.pdf")
                #plt.clf()

                # cumulative cpu usage over time
                #fig2=plt.figure(2, figsize=(20,10))
                data1=np.zeros(steps)
                data2=np.zeros(steps)
                times=np.arange(steps)
                for i in range(1,number_of_found_subtags):
                    data2=data2+All_value_cum_rel[:,i]
                    x=TimeStepMod*np.concatenate((times[:],[times[steps-1]],times[::-1],[times[0]]))/1000 
                    y=np.concatenate((data1[:],[data2[steps-1]],data2[::-1],[data1[0]]))
                    ax2.fill(x,y,label=label_subtags[found_subtags[i]],color=tcolors[subtags.index(found_subtags[i])])
                    ax2.text(0.4, 0.935, "miscellaneous", transform = ax2.transAxes, size=fontsize, color="brown")
                    data1=data2
                ax2.axis([x.min(),0.99*x.max(),0,100.])
                if figi == 1:
                # reverse the order
                    prop = matplotlib.font_manager.FontProperties(size=fontsize)
                    #ax2.legend(handles[::-1], labels[::-1],loc="lower right",prop=prop, bbox_to_anchor=(0., -1.0, 1., .102)).draggable()
                    fig2.legend(handles[::-1], labels[::-1], 'upper center',prop=prop,ncol=number_of_found_subtags-1)
                #plt.savefig(tmppath+"/times_rel_cum.pdf", bbox_inches='tight')
            
            except IOError:
                print 'cpu.txt does not exist for:',titlestr
                ax1.text(0.5,0.5,"IOError", transform=ax1.transAxes,horizontalalignment='center',verticalalignment='center')
                ax2.text(0.5,0.5,"IOError", transform=ax2.transAxes,horizontalalignment='center',verticalalignment='center')
                ax1.set_yticks([])
                ax1.set_xticks([])
                ax2.set_yticks([])
                ax2.set_xticks([])

            if figi == 15:
                ax1.set_ylabel("relative (%)", fontsize=fontsize)
                ax2.set_ylabel("cumulative (%)", fontsize=fontsize)
            if figi >= 29:
                ax1.set_xlabel("time-step (10$^3$)", fontsize=fontsize)
                ax2.set_xlabel("time-step (10$^3$)", fontsize=fontsize)

            ax1.text(0.5,0.05,titlestr,
                    horizontalalignment='center',
                    verticalalignment='bottom',
                    color='white',
                    transform = ax1.transAxes)
            ax2.text(0.5, 0.05,titlestr,
                    horizontalalignment='center',
                    verticalalignment='bottom',
                    color='white',
                    transform = ax2.transAxes)

            if figi == 5:
                ax1.set_ylabel('relative cpu usage')
                ax2.set_ylabel('cumulative cpu usage')

        #if figi == 
            #plt.clf()

fig1.savefig("/bigbang/data/bgriffen/modules/brendanlib/cpu_relative.png")
fig2.savefig("/bigbang/data/bgriffen/modules/brendanlib/cpu_cumulative.png")

plt.show()
