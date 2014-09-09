import numpy as np
import glob,sys
import pylab as plt

base_path = "/Users/griffen/Desktop/H706754_BB_Z127_P7_LN7_LX13_O4_NV4/"
filename = "rockstar.e"

full_path = base_path + filename

f = open(full_path,'r')
ti = 0.
info = {"read":[],"transfer_boundary":[],"transfer_particles":[],"load":[],"constructing_mt":[],"snap":[],"fof":[],"halos":[]}

tcolors = ["Azure", "Blue", "Green", "Red", "SpringGreen", "DarkTurquoise", "Red", "Black", "FireBrick", "Purple", "OrangeRed", "SeaGreen", "Navy", "Yellow", "Tomato", "Violet", "Tan", "Olive", "Gold"] 

for line in f:
    if "[" in line and "Success" not in line and "[Finished" not in line:
        description = line.split("]")[1]
        tf = int(line.split()[1][:-2])
        dt = tf - ti
        #print tf,dt,description.strip()

        if "snap" in description:
           snap = description.split("snapshot")[1]
           snap = int(snap.replace("...",""))
           info['snap'].append(snap)

        if "Reading" in description:
           info['read'].append(dt)

        if "FoF" in description:
            info['fof'].append(dt)

        if "halos" in description:
            info['halos'].append(dt)

        if "Transferring particles" in description:
            info['transfer_particles'].append(dt)

        if "Transferring boundary" in description:
            info['transfer_boundary'].append(dt)

        if "Loading" in description:
            info['load'].append(dt)
    
        if "Constructing" in description:
            info['constructing_mt'].append(dt)

        ti = tf

print "----------------------------------"
print "Key:      Length / Time [min.]"
print "----------------------------------"
print "Snap:        %3i / %3.2f   "         % (len(info['snap']),np.sum(np.array(info['snap']))/60.)
print "Read:        %3i / %3.2f         "   % (len(info['read']),np.sum(np.array(info['read']))/60.)
print "FoFs:        %3i / %3.2f          "  % (len(info['fof']),np.sum(np.array(info['fof']))/60.)
print "Halos:       %3i / %3.2f        "   % (len(info['halos']),np.sum(np.array(info['halos']))/60.)
print "T Boundary:  %3i / %3.2f   "   % (len(info['transfer_boundary']),np.sum(np.array(info['transfer_boundary']))/60.)
print "T Particles: %3i / %3.2f   "  % (len(info['transfer_particles']),np.sum(np.array(info['transfer_particles']))/60.)
print "Load:        %3i / %3.2f         "   % (len(info['load']),np.sum(np.array(info['load']))/60.)
print "Construct MT:%3i / %3.2f    "   % (len(info['constructing_mt']),np.sum(np.array(info['constructing_mt']))/60.)
print "----------------------------------"
fig1 = plt.figure(figsize=(16,6))
fig2 = plt.figure()
ax2 = fig2.add_subplot(1,1,1)
figi = 1
for item in info.keys():
    if "snap" not in item:
        len_vec = len(info[item])
        ax1 = fig1.add_subplot(2,4,figi)
        ax1.plot(info['snap'][:len_vec],info[item][:len_vec],'b-',linewidth=3)
        ax2.plot(info['snap'][:len_vec],np.cumsum(np.array(info[item][:len_vec]))*100./float(tf),label=item,color=tcolors[figi],linewidth=3)
        ax1.set_ylabel(item + " [s]")
        
        if figi > 4:
            ax1.set_xlabel('snapshot')

        figi += 1

ax2.set_xlim([0,len_vec])
ax2.set_xlabel("snapshot")
ax2.set_ylabel("cumulative fraction [%]")
ax2.legend(fancybox=True,fontsize=10,loc=2)
plt.show()
