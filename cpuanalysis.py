import string
import numpy as np
import matplotlib.pylab as plt
import matplotlib.colors
import matplotlib.font_manager
from matplotlib.ticker import MaxNLocator
from matplotlib.font_manager import FontProperties
cm = plt.cm.get_cmap('RdYlBu')
import os
import sys

def cpu(basedir):
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
    
    fig = plt.figure(figsize=(16,7))
    
    ax1 = fig.add_subplot(1,2,1)
    ax2 = fig.add_subplot(1,2,2)
    
    cputxt = basedir + '/cpu.txt'
    timingstxt = basedir + '/timings.txt'
    energytxt = basedir + '/energy.txt'
    infotxt = basedir + '/info.txt'
    
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
        handles, labels = ax1.get_legend_handles_labels()
        # reverse the order
        prop = matplotlib.font_manager.FontProperties(size=fontsize)
        #ax1.legend(handles[::-1], labels[::-1],loc="lower right",prop=prop, bbox_to_anchor=(0., -5.0, 1., .102)).draggable()
        fig.legend(handles[::-1], labels[::-1], 'upper center',prop=prop,ncol=number_of_found_subtags-1)
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
            ax1.text(0.4, 0.935, "miscellaneous", transform = ax1.transAxes, size=fontsize, color="brown")
            data1=data2
        ax2.axis([x.min(),0.99*x.max(),0,100.])
        #if figi == 1:
        ## reverse the order
        #    prop = matplotlib.font_manager.FontProperties(size=fontsize)
        #    #ax2.legend(handles[::-1], labels[::-1],loc="lower right",prop=prop, bbox_to_anchor=(0., -1.0, 1., .102)).draggable()
        #    fig.legend(handles[::-1], labels[::-1], 'upper center',prop=prop,ncol=number_of_found_subtags-1)
        ##plt.savefig(tmppath+"/times_rel_cum.pdf", bbox_inches='tight')
    
    except IOError:
        print 'cpu.txt does not exist for:',nvir,pad
    
    ax1.set_ylabel("relative (%)", fontsize=fontsize)
    ax2.set_ylabel("cumulative (%)", fontsize=fontsize)
    ax1.set_xlabel("time-step (10$^3$)", fontsize=fontsize)
    ax2.set_xlabel("time-step (10$^3$)", fontsize=fontsize)
    #ax1.text(0.5,0.05,"",
    #        horizontalalignment='center',
    #        verticalalignment='bottom',
    #        color='white',
    #        transform = ax1.transAxes)
    #ax2.text(0.5, 0.05,titlestr,
    #        horizontalalignment='center',
    #        verticalalignment='bottom',
    #        color='white',
    #        transform = ax2.transAxes)
    #plt.clf()
    #fig1.savefig('./figs/' + regionlist[0] + 'cpu_relative.png')
    #fig2.savefig('./figs/' + regionlist[0] + 'cpu_cumulative.png')
    
    plt.show()

if __name__ == '__main__':
    import os
    cpu(os.getcwd())
