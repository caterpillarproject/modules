import sys,os
import subprocess
import numpy as np
import pylab as plt
import platform
import time,glob
from brendanlib.grifflib import getcurrentjobs,get_folder_size,get_last_modified

currentjobs,jobids,statusjobs = getcurrentjobs()

if platform.node() == "bigbang.mit.edu":
    basepath = "/bigbang/data/AnnaGroup/caterpillar/halos/oldhalos"

if "bigbang" in platform.node() or "antares" in platform.node() or "spacebase" in platform.node():
    basepath = "/bigbang/data/AnnaGroup/caterpillar/halos"

if "harvard" in platform.node():
    basepath = "/n/home01/bgriffen/data/caterpillar/halos"

halo_list = glob.glob("/bigbang/data/AnnaGroup/caterpillar/halos/low_mass_halos/H*")
halo_list += glob.glob("/bigbang/data/AnnaGroup/caterpillar/halos/high_mass_halos/H*")
halo_list += glob.glob("/bigbang/data/AnnaGroup/caterpillar/halos/middle_mass_halos/H*")

haloidlist = []
for folder_name in halo_list:
    halo_name = folder_name.split("/")[-1]
    haloidlist.append(halo_name)

print "------------------------------------------------------------------------------------------------------"
print "  HALO   | NV | ICG | LX | ICS | SNAP | %PC | S | R | CPU.TXT | OUTPUT | IC SIZE | SNAP SIZE | STATUS"
print "------------------------------------------------------------------------------------------------------"
for haloid in haloidlist:
    run_list = glob.glob("/bigbang/data/AnnaGroup/caterpillar/halos/*_mass_halos/"+haloid+"/H*")
    for run in run_list:
        level = int(run.split("LX")[1][:2])
        nrvir = int(run.split("NV")[1][0])
        geometry = run.split(haloid + "_")[1][:2]

        if os.path.isfile(run + "/ics.0"):
            ic_file_size = os.path.getsize(run + "/ics.0")
            if ic_file_size > 0:
                ic_marker = '+'
            elif ic_file_size == 0:
                 ic_marker = "x"
        else:
            ic_marker = "o"

        subdirnames = glob.glob(run + "/outputs/snapdir_*")
        snapshotvec = []
        for subname in subdirnames:
            snapshotvec.append(int(subname.split("snapdir_")[-1]))

        if snapshotvec:
            snapshot = max(snapshotvec)
            pc_complete = "%3i" % (float(snapshot)*100/255.)

            snap_folder = run + "/outputs/snapdir_"+str(snapshot).zfill(3)
            snap_file_size = "%5i" % (get_folder_size(snap_folder)/1e6)
            #pc_complete = .zfill(3)
        else:
            snap_file_size = "    -"
            snapshot = " - "
            pc_complete = "  -"
            
        if os.path.isdir(run + "/outputs/groups_255"):
            subfind = "x"
        else:
            subfind = "-"

        if os.path.isdir(run + "/halos/"):
            rockstar = "x"
        else:
            rockstar = "-"

        jobname = haloid+"N"+str(nrvir)+"L"+str(level)[-1]

        if jobname in currentjobs:
            idx = currentjobs.index(jobname)
            statusi = statusjobs[idx]
        else:
            statusi = "NR"
            #print statusi

        cpu_last_modified = get_last_modified(run + "/outputs/cpu.txt")
        output_last_modified = get_last_modified(run + "/OUTPUT")

        print haloid.ljust(10) + "  " + str(nrvir) + "    " + geometry + "   " + str(level) + "    " + ic_marker + "    " + str(snapshot).zfill(3) \
            + "    " + str(pc_complete) + "   " + subfind + "   " + rockstar + "   "  + str(cpu_last_modified).rjust(6) + "    "  + str(output_last_modified).rjust(6) \
            + "     " + str("%5i" % (ic_file_size/1e6)) + "       " + snap_file_size + "     " + statusi
    if len(run_list) > 0:
        print "------------------------------------------------------------------------------------------------------"
    
#print 
