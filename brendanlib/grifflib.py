import pylab as plt
import platform
from array import array
import numpy as np
from scipy.optimize import curve_fit
from matplotlib import rc
from scipy.optimize import leastsq
from math import *
import numpy.random as nprnd
import matplotlib.colors as col
import math,glob,subprocess,shlex,os,asciitable
import readsnapshots.readsnapHDF5_greg as rs
import mergertrees.MTCatalogue as MT
import asciitable
import sys,os,time
import pickle
import haloutils as htils
import readsnapshots.readhsml as readhsml
#import viz.SphMap as SphMap
import calc.CalcHsml as ca
from matplotlib.colors import LogNorm
import readhalos.RSDataReader as RS
import viz.SphMap as SphMap

fig_for_talks = {  'lines.color': 'white',
            'patch.edgecolor': 'white',
            'text.color': 'white',
            'axes.facecolor': 'black',
            'axes.edgecolor': 'white',
            'axes.labelcolor': 'white',
            'xtick.color': 'white',
            'ytick.color': 'white',
            'grid.color': 'white',
            'figure.facecolor':  'white',
            'figure.edgecolor':  'white',
            'savefig.facecolor': 'white',
            'savefig.edgecolor': 'white',
	    'savefig.dpi':300,            
            'axes.linewidth': 4,
            'font.size' : 24,
            'xtick.labelsize': 20,
            'ytick.labelsize': 20,

            'ytick.major.width':3,
            'ytick.minor.width':3,
            'xtick.major.width':3,
            'xtick.minor.width':3,

            'xtick.major.size':10,
            'xtick.minor.size':5,
            'ytick.major.size':10,
            'ytick.minor.size':5,
            'xtick.major.pad':10,
            'xtick.minor.pad':10,
            'ytick.major.pad':10,
            'ytick.minor.pad':10,
            'figure.figsize': (12.5,7.5),
            'legend.framealpha':1,
	    'legend.frameon': False,
	    'legend.borderpad':1,
            'font.weight': 'bold',
            'axes.labelweight':'bold',
            'savefig.transparent':True,
            'savefig.bbox': 'tight'}

fig_for_papers = {  
            
            'axes.linewidth': 4,
            'font.size' : 28,
            'xtick.labelsize': 22,
            'ytick.labelsize': 22,

            'xtick.major.size':12,
            'xtick.minor.size':8,
            'ytick.major.size':12,
            'ytick.minor.size':8,

            'ytick.major.width':4,
            'xtick.major.width':4,
            'ytick.minor.width':3,
            'xtick.minor.width':3,
            'xtick.major.pad':10,
            'xtick.minor.pad':10,
            'ytick.major.pad':10,
            'ytick.minor.pad':10,
            'figure.figsize': (10,10),
            'legend.fontsize': 16,
            'legend.frameon': False,
            #'font.weight': 'bold',
            #'axes.labelweight':'bold',
            'legend.borderpad':1,
            'savefig.bbox': 'tight',
            'savefig.dpi':300,
            'text.latex.preamble': r"\usepackage{amsmath}"}

fig_for_normal = {
	    'axes.linewidth': 3,
            'font.size' : 18,
	    'xtick.major.size':12,
            'xtick.minor.size':8,
            'ytick.major.size':12,
            'ytick.minor.size':8,

            'ytick.major.width':3,
            'xtick.major.width':3,
            'ytick.minor.width':2,
            'xtick.minor.width':2,

            'xtick.major.pad':10,
            'xtick.minor.pad':10,
            'ytick.major.pad':10,
            'ytick.minor.pad':10,
	    
            'figure.figsize': (8,6),
            'legend.fontsize': 16,
	    'legend.borderpad':1,
            'legend.frameon': True,
            'savefig.bbox': 'tight',
            'savefig.dpi':300,
            'text.latex.preamble': r"\usepackage{amsmath}"}

def savefig(fname, quant, pixel, vmin = None, vmax = None, cmap = None, hires = None, text = None, xtext = None, ytext = None, textsize = None, textcolor = None,nrvir=None):
    fig,ax = plt.subplots(1,1,frameon=False,figsize=(5,5))
#    ax.set_axis_off()

#    fig.add_axes(ax)
    im=ax.imshow(quant, vmin=vmin, vmax=vmax, origin='lower', interpolation='nearest', aspect='auto', cmap=cmap)                   
    if (text != None):
        print "adding text:", text
        ax.text(xtext, ytext, text, horizontalalignment='left', verticalalignment='center', transform = ax.transAxes, size=textsize, color=textcolor,weight='bold') #, bbox=dict(facecolor='gray', alpha=0.8, boxstyle="round, pad=0.25", edgecolor="pink"))

    if (nrvir!=None):
	xcirc,ycirc = drawcircle(float(pixel)/2.,float(pixel)/2.,(1/(nrvir*2.))*pixel)
	ax.plot(xcirc,ycirc,'w-',linewidth=0.5)

    ax.set_xlim([0,pixel])
    ax.set_ylim([0,pixel])
    ax.set_axis_off()
    ax.xaxis.set_major_locator(plt.NullLocator())
    ax.yaxis.set_major_locator(plt.NullLocator())

    
    
    if hires: 
        if ".png" in fname: 
            fig.savefig(fname.replace(".png","_hires.png"), bbox_inches='tight', dpi=pixel, pad_inches=0)
        elif ".pdf" in fname:
            fig.savefig(fname.replace(".pdf","_hires.pdf"), bbox_inches='tight', dpi=pixel, pad_inches=0)
    else:
        fig.savefig(fname, bbox_inches='tight', dpi=pixel/8., pad_inches=0)
        
    plt.close(fig)

def getRSCat(label, lx, hid):
    snap_num = 255
    hpath = htils.get_hpath_lx(hid,lx)
    if label=='M':
        version=6
    elif label=='B':
        version=8
    elif label =='K' or label=='L':
        version=9
    elif label in ['O','O2','P']:
        version=10

    print "READING: halos_%s/" % (label)

    cat = RS.RSDataReader(hpath+'/halos_'+label,snap_num,version=version,digits=1,unboundfrac=None,minboundpart=None)
    return cat

def distance(posA, posB,boxsize=100.):
    dist = abs(posA-posB)
    tmp = dist > boxsize/2.0
    dist[tmp] = boxsize-dist[tmp]
    if dist.shape == (3,):
        return np.sqrt(np.sum(dist**2))
    else:
        return np.sqrt(np.sum(dist**2,axis=1))
        
def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)

def angle_between(v1, v2):

    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)

    angle = np.arccos(np.dot(v1_u, v2_u))    
    if np.isnan(angle):
        if (np.abs(v1_u-v2_u)<1e-7).all():
            return 0.0
        else:
            return np.pi
    
    return angle

def get_size_of_files(file_list):
    total_size = 0

    for filei in file_list:
        total_size += os.path.getsize(filei)
    
    return total_size

def get_size_of_folder(folder):
    total_size = os.path.getsize(folder)
    for item in os.listdir(folder):
        itempath = os.path.join(folder, item)
        if os.path.isfile(itempath):
            total_size += os.path.getsize(itempath)
        elif os.path.isdir(itempath):
            total_size += get_size_of_folder(itempath)
    return total_size

def get_size_of_folders(folders):
    total_size = 0

    for filei in folders:
        total_size += get_size_of_folder(filei)
    
    return total_size

def first_twelve():
    #first_twelve = ["H1631506", "H1195448", "H1725139", "H447649", "H5320", "H581141", "H94687", "H1130025", "H1387186", "H581180", "H1725372", "H1354437"]
    first_twelve = ["H1631506","H264569","H1725139","H447649","H5320","H581141","H94687","H1130025","H1387186","H581180","H1725372","H1354437"]

    return first_twelve

def load_dict(file_name):
    dict = pickle.load( open( file_name, "rb" ) )
    return dict

def save_dict(file_name, dict):
    pickle.dump( dict, open( file_name, "wb" ) )

def load_halo_alias():
    base_path = "/bigbang/data/AnnaGroup/caterpillar/halos/"
    file_name = base_path+"public_halo_names.p"
    dict = load_dict(file_name)

    return dict

def create_slurm_job(submit_line,job_name,ncores,queue,snap,memory,hours):
    f = open("sjob.slurm",'w')
    f.write("#!/bin/bash\n")
    f.write("#SBATCH -n " + str(ncores) + "\n")
    f.write("#SBATCH -o job.o" + str(snap) + "\n")
    f.write("#SBATCH -e job.e" + str(snap) + "\n")
    f.write("#SBATCH -J " + job_name + "\n")
    f.write("#SBATCH -p " + str(queue) + "\n")
    f.write("#SBATCH -t " + str(hours) + ":00:00\n")
    f.write("#SBATCH --mem=" + str(memory) + "gb\n")
    f.write("\n")
    f.write("source new-modules.sh; module load python\n")
    f.write("python construct_index.py " + str(snap) + "\n")
    f.close()

#def get_quant_zoom(halo_path,quant):
#    import alexlib.haloutils as haloutils
#    htable = haloutils.get_parent_zoom_index()
#    halo_split = halo_path.split("_")
#    haloid = int(halo_split[0].split("H")[1])
#    geom,lx,nrvir = haloutils.get_zoom_params(halo_path)

#    mask = (haloid == htable['parentid']) & \
#           (geom == htable['ictype']) & \
#           (int(lx) == htable['LX']) & \
#           (int(nrvir) == htable['NV']) 

#    print np.sum(geom == htable['ictype'])
#    print np.sum(int(lx) == htable['LX'])
#    print np.sum(int(nrvir) == htable['NV'])
#    print np.sum(haloid == htable['parentid'])
#    print np.sum(mask)
#    print
#    return htable[mask][quant]
    
def get_folder_size(folder):
    total_size = os.path.getsize(folder)
    for item in os.listdir(folder):
        itempath = os.path.join(folder, item)
        if os.path.isfile(itempath):
            total_size += os.path.getsize(itempath)
        elif os.path.isdir(itempath):
            total_size += getFolderSize(itempath)
    return total_size
    
def get_last_modified(file_name):
    if os.path.isfile(file_name):
        tnow = time.time()
        #should_time = tnow - (tthresh * 3600)
        file_mod_time = os.stat(file_name).st_mtime
        last_modified = round((int(tnow) - file_mod_time) / 3600, 2)
        last_modified = "%3.2f" % (last_modified)
    else:
        last_modified = ""
        
    return last_modified
    
def convert_pid_zid(pid_in,lx_in):
    htable = asciitable.read("/bigbang/data/AnnaGroup/caterpillar/halos/parent_zoom_index.txt",Reader=asciitable.FixedWidth)
    key = [str(pid)+'_'+str(lx) for pid,lx in zip(htable['parentid'],htable['LX'])]
    hindex = dict(zip(key,htable['zoomid']))
    zoomid = hindex[str(pid_in)+'_'+str(lx_in)]
    return zoomid

def get_completed_list(suite_paths,verbose=True):
    gadget_done = []
    subfind_done = []
    ic_done = []
    haloids = []

    suite_paths = sorted(suite_paths)
    for suite in suite_paths:
        if os.path.isdir(suite + "/outputs/snapdir_255"):
            gadget_done.append(1)
        else:
            gadget_done.append(0)

        if os.path.isfile(suite + "/ics.0") and os.path.getsize(suite + "/ics.0") > 0:
            ic_done.append(1)
        else:
            ic_done.append(0)

        if os.path.isdir(suite + "/outputs/groups_255"):
            subfind_done.append(1)
        else:
            subfind_done.append(0)

    if verbose:
        print "-----------------------------------"
        print "     RUN                  I G S"
        print "-----------------------------------"
        #print "H299792_A            "
        for suite,gadget,ic,subfind in zip(suite_paths,gadget_done,ic_done,subfind_done):
            print suite.split("/")[-1].replace("Z127_P7_LN7_","").replace("O4_","").ljust(25),ic,gadget,subfind 

    return gadget_done,ic_done,subfind_done



def make_subfind_submission_script(runpath,job_name):
    f = open(runpath + "/ssubfind",'w')
    f.write('#!/bin/bash\n')
    f.write('#SBATCH -o gadget.o%j\n')
    f.write('#SBATCH -e gadget.e%j\n')
    f.write('#SBATCH -J '+ job_name + '\n')

    if "LX11" in runpath or "LX12" in runpath:
        f.write('#SBATCH -p RegNodes\n')
        ncores = 8

    if "LX13" in runpath:
        f.write('#SBATCH -p AMD64\n')
        ncores = 64

    if "LX14" in runpath:
        f.write('#SBATCH -p AMD64\n')
        ncores = 256
    
    f.write('#SBATCH -n ' + str(ncores) + '\n')
    f.write("\n")
    f.write("cd " + runpath + "\n")
    f.write("\n")
    f.write('mpirun -np ' + str(ncores) +  ' ./P-Gadget3_sub param_sub.txt --bind-to-core 1>OUTPUTsub 2>ERRORsub\n')
    f.close()

def make_music_submission_script(runpath,cfgname,job_name):
    f = open(runpath + "/smusic",'w')
    f.write("#!/bin/bash \n")
    f.write("#SBATCH -o music.o%j \n")
    f.write("#SBATCH -e music.e%j \n")
    f.write("#SBATCH -N 1\n")
    f.write("#SBATCH --exclusive\n")
    
    #if "LX11" in runpath:
    #    f.write('#SBATCH -p AMD64\n')
    #    ncores = 24
	    
    if "LX11" in runpath or "LX12" in runpath or "LX13" in runpath or "LX14" in runpath:
        f.write('#SBATCH -p AMD64\n')
        ncores = 64

    f.write("#SBATCH -J " + job_name + "\n")
    f.write("\n")
    f.write("export OMP_NUM_THREADS=" + str(ncores) + "\n")
    f.write("\n")
    f.write("cd " + runpath + "\n")
    f.write("\n")
    f.write("./MUSIC ./" + cfgname + ".conf 1>OUTPUTmusic 2>ERRORmusic\n")
    f.close()

def get_sim_info_from_lx(lx,ncores,queue):
    if lx == 1:
        pmgrid = 256
    if lx == 2:
        pmgrid = 512
    if lx == 3:
        pmgrid = 512
        queue = "HyperNodes"
        ncores = 256
    if lx == 4:
        pmgrid = 512
        queue = "AMD64"
        ncores = 512

    return pmgrid,queue,ncores	

def make_gadget_submission_script(runpath,job_name,restart_flag,gadget3):
    f = open(runpath + "/sgadget",'w')
    f.write('#!/bin/bash\n')
    
    if gadget3:
        f.write('#SBATCH -o gadget3.o%j\n')
        f.write('#SBATCH -e gadget3.e%j\n')
        
    if not gadget3:
        f.write('#SBATCH -o gadget4.o%j\n')
        f.write('#SBATCH -e gadget4.e%j\n')

    f.write('#SBATCH -J '+ job_name + '\n')

    if "LX11" in runpath:
        queue = "HyperNodes"
        ncores = 24
        core_setup = '-n ' + str(ncores)
        #queue = "AMD64"
        #ncores = 64
        #core_setup = '-n ' + str(ncores)

    if "LX12" in runpath:
        queue = "HyperNodes"
        ncores = 24
        core_setup = '-n ' + str(ncores)
        #queue = "AMD64"
        #ncores = 64
        #core_setup = '-n ' + str(ncores)
        
    if "LX13" in runpath:
        queue = "AMD64"
        ncores = 64
        core_setup = '-n ' + str(ncores)
        #queue = "HyperNodes"
        #ncores = 144
        #core_setup = '-N 6 -n ' + str(ncores)

    if "LX14" in runpath:
        queue = "AMD64"
        ncores = 256
        core_setup = '-n ' + str(ncores)

    f.write('#SBATCH -p ' + queue + '\n')
    f.write('#SBATCH ' + core_setup + '\n')

    f.write("\n")
    f.write("cd " + runpath + "\n")
    f.write("\n")
    
    if not restart_flag:
        if gadget3:
            if queue != "HyperNodes":
                f.write('mpirun --bind-to-core -np ' + str(ncores) +  ' ./P-Gadget3 param.txt 1>OUTPUT 2>ERROR\n')
            else:
                f.write('mpirun -np ' + str(ncores) +  ' ./P-Gadget3 param.txt 1>OUTPUT 2>ERROR\n')
        else:
            if queue != "HyperNodes":
                f.write('mpirun --bind-to-core -np ' + str(ncores) +  ' ./Gadget4 param.txt 1>OUTPUT 2>ERROR\n')
            else:
                f.write('mpirun -np ' + str(ncores) +  ' ./Gadget4 param.txt 1>OUTPUT 2>ERROR\n')
    else:
        if gadget3:
            if queue != "HyperNodes":
                f.write('mpirun --bind-to-core -np ' + str(ncores) +  ' ./P-Gadget3 param.txt 1 1>>OUTPUT 2>ERROR\n')
            else:
                f.write('mpirun -np ' + str(ncores) +  ' ./P-Gadget3 param.txt 1 1>>OUTPUT 2>ERROR\n')
        else:
            if queue != "HyperNodes":
                f.write('mpirun --bind-to-core -np ' + str(ncores) +  ' ./Gadget4 param.txt 1 1>>OUTPUT 2>ERROR\n')
            else:
                f.write('mpirun -np ' + str(ncores) +  ' ./Gadget4 param.txt 1 1>>OUTPUT 2>ERROR\n')
    f.close()

def run_gadget(suite_paths,gadget_file_path,lx_list,submit=True):
    job_name_list = []
    current_jobs,jobids,jobstatus = getcurrentjobs()
    for folder in suite_paths:
        topfolder = "/".join(folder.split("/")[:8])
        #print topfolder
        folder_single = folder.split("/")[-1]
        halo_label = folder_single.split("_")[0]
        nrvir_label = folder_single.split("NV")[1][0]
        lx_label = folder_single.split("LX")[1][1:2]
        job_name = halo_label+"X"+lx_label+"N"+nrvir_label+folder_single.split(halo_label+"_")[1][:2]
        
        if os.path.isfile(topfolder+"/.gadget3"):
            gadget3 = True
            lastsnap = '255'
            job_name += 'G3'
        elif os.path.isfile(topfolder+"/.gadget4"):
            gadget3 = False
            lastsnap = '319'
            job_name += 'G4'

      
        if os.path.isfile(folder + "/ics.0") and os.path.getsize(folder + "/ics.0") > 0 and \
            job_name not in current_jobs and job_name not in job_name_list:
            if not os.path.isdir(folder + "/outputs/snapdir_"+lastsnap) and lx_label in lx_list:
                print
                print "COPYING GADGET FILES...",folder_single
                mkdir_outputs = "mkdir -p " + folder + "/outputs/"
                
                if "LX11" in folder:
                    ncores = 24
                    pmgrid = 256
                    queue = "HyperNodes"
                    #ncores = 64
                    #pmgrid = 256
                    #queue = "AMD64"

                if "LX12" in folder:
                    ncores = 24
                    pmgrid = 256
                    queue = "HyperNodes"
                    #ncores = 64
                    #pmgrid = 256
                    #queue = "AMD64"

                if "LX13" in folder:
                    ncores = 64
                    pmgrid = 512
                    queue = "AMD64"
                    #queue = "HyperNodes"
                    #ncores = 144
                    #pmgrid = 512
                    
                if "LX14" in folder:
                    queue = "AMD64"
                    ncores = 512
                    pmgrid = 512
                
                if gadget3:
                    if "contamination" in folder:
                        file_times  = "cp " + gadget_file_path + "/gadget3/ExpansionList_contam " + folder + "/ExpansionList"
                    else:
                        file_times  = "cp " + gadget_file_path + "/gadget3/ExpansionList " + folder + "/ExpansionList"
                else:
                    if "contamination" in folder:
                        file_times  = "cp " + gadget_file_path + "/gadget4/ExpansionList_contam_hi " + folder + "/ExpansionList"
                    else:
                        file_times  = "cp " + gadget_file_path + "/gadget4/ExpansionList_hi " + folder + "/ExpansionList"

                if len(glob.glob(folder+"/outputs/snapdir*")) > 0 and len(glob.glob(folder+"/outputs/restartfiles/*.bak")) == ncores:
                    restart_flag = True
                else:
                    restart_flag = False

                if len(glob.glob(folder+"/outputs/snapdir*")) == 0 and len(glob.glob(folder+"/outputs/restartfiles/*.bak")) == 0:
                    print "OUTPUTS/ IS CLEAN > STARTING AT START!"
                    skip_flag = False
                elif len(glob.glob(folder+"/outputs/snapdir*")) != 0 and len(glob.glob(folder+"/outputs/restartfiles/*.bak")) == ncores:
                    print "RESTARTING FROM RESTART FILES > " + glob.glob(folder+"/outputs/snapdir*")[-1].split("/")[-1]
                    skip_flag = False
                elif len(glob.glob(folder+"/outputs/snapdir*")) != 0 and len(glob.glob(folder+"/outputs/restartfiles/*.bak")) != ncores:
                    print "WANTED TO RESTART FROM RESTART FILES BUT # RESTART FILES != ncores > " + glob.glob(folder+"/outputs/snapdir*")[-1].split("/")[-1]
                    skip_flag = True
                    print "# restart files:",len(glob.glob(folder+"/restartfiles/*.bak")),"# cores:",ncores

                if not skip_flag: 
                    if gadget3:
                        file_exe    = "cp " + gadget_file_path + "/gadget3/P-Gadget3_"+str(pmgrid)+ " "  + folder + "/P-Gadget3"
                        file_param  = "cp " + gadget_file_path + "/gadget3/param_1"+lx_label+"_" + queue + ".txt " + folder + "/param.txt"
                        file_config = "cp " + gadget_file_path + "/gadget3/Config_"+str(pmgrid)+".sh " + folder + "/Config.sh"
                    else:
                        file_exe    = "cp " + gadget_file_path + "/gadget4/Gadget4_"+str(pmgrid)+ " "  + folder + "/Gadget4"
                        file_param  = "cp " + gadget_file_path + "/gadget4/param_1"+lx_label+"_" + queue + ".txt " + folder + "/param.txt"
                        file_config = "cp " + gadget_file_path + "/gadget4/Config_"+str(pmgrid)+".sh " + folder + "/Config.sh"

                    cmd_copy_all_files = [mkdir_outputs,file_times,file_exe,file_param,file_config]
                    subprocess.call([";".join(cmd_copy_all_files)],shell=True)
                    make_gadget_submission_script(folder,job_name,restart_flag,gadget3)
                    cd_folder = "cd " + folder
                    cmd_submit_gadget = "sbatch sgadget"

                    if submit == True:
                        print "SUBMITTING GADGET..."
                        subprocess.call([cd_folder+"; "+cmd_submit_gadget],shell=True)
                        job_name_list.append(job_name)

                else:
                    print "NOT SUBMITTING!"

def run_subfind(suite_paths,gadget_file_path):
    job_name_list = []
    current_jobs,jobids,jobstatus = getcurrentjobs()
    for folder in suite_paths:
        folder_single = folder.split("/")[-1]
        halo_label = folder_single.split("_")[0]
        job_name = "G"+halo_label[1:]+"LX11"+folder_single.split(halo_label+"_")[1][:2]
        if os.path.isfile(folder + "/ics.0") and os.path.getsize(folder + "/ics.0") > 0 and not os.path.isdir(folder + "/outputs/") and \
            job_name not in current_jobs and job_name not in job_name_list:
            if not os.path.isdir(folder + "/outputs/snapdir_255") and not os.path.isdir(folder + "/outputs/groups_255"):
                print "COPYING SUBFIND FILES..."
                file_exe    = "cp " + gadget_file_path + "P-Gadget3_256sub " + folder + "/P-Gadget3_sub"
                file_param  = "cp " + gadget_file_path + "paramsub_11.txt " + folder + "/param_sub.txt"
                file_config = "cp " + gadget_file_path + "Configsub_256.sh " + folder + "/Config_sub.sh"
                cmd_copy_all_files = [file_times,file_exe,file_param,file_config]
                subprocess.call([";".join(cmd_copy_all_files)],shell=True)
                print "SUBMITTING SUBFIND..."
                make_subfind_submission_script(folder,folder_single,job_name)
                cd_folder = "cd " + folder
                cmd_submit_gadget = "sbatch ssubfind"
                subprocess.call([cd_folder+"; "+cmd_submit_subfind],shell=True)
                job_name_list.append(job_name)


def run_music(suite_paths,music_path,lagr_path,lx_list):
     for folder in suite_paths:
         folder_single = folder.split("/")[-1]
         halo_label = folder_single.split("_")[0]
         nrvir_label = folder_single.split("NV")[1][0]
         lx_label = folder_single.split("LX")[1][1:2]
         #job_name = halo_label[:4]+"X"+lx_label+"N"+nrvir_label+folder_single.split("_")[-1]
         job_name = "I"+halo_label[1:]+"X"+lx_label+"N"+nrvir_label+folder_single.split(halo_label+"_")[1][:2]
         job_name_list = []
         current_jobs,jobids,jobstatus = getcurrentjobs()
#         if os.path.isfile(folder + "/ics.0") and os.path.getsize(folder + "/ics.0") == 0:
#             print "NEED TO RERUN ICS:",job_name

         if not os.path.isfile(folder + "/ics.0") and job_name not in current_jobs and job_name not in job_name_list and lx_label in lx_list:
             print
             print "RUNNING:",folder_single
             print "MAKING MUSIC SUBMISSION SCRIPT..."

             make_music_submission_script(folder,folder_single,job_name)
             
             print "COPYING MUSIC FILES..."
             cmd_cp_music = "cp " + music_path + " " + folder
             subprocess.call([cmd_cp_music],shell=True)
         
             print "CONSTRUCTING MUSIC CONFIGURATION FILES..."
             lagr_file = lagr_path + folder_single.split("_")[0]
             master_music_cfg_dest = folder + "/" + folder_single + ".conf"
             
             region_point_file = lagr_path+"H"+halo_label[1:]+"NRVIR"+nrvir_label
             seed = int(halo_label[1:])
             
             if os.path.isfile(region_point_file + ".head"):
                 if halo_label+"_B" in folder_single:
                    ictype = 'box'
                    
                 if halo_label+"_E" in folder_single:
                    ictype = "ellipsoid"

                 if halo_label+"_C" in folder_single:
                    ictype = "convex_hull"

             #print master_music_cfg_dest,folder_single,ictype,seed,region_point_file,nrvir_label
             #sys.exit()
             make_music_file(master_music_cfg_dest,ictype,seed,region_point_file,nrvir_label)

             #np.loadtxt(lagr_path+"H"+halo_label[1:]+".head")
             #with open(lagr_path+"H"+halo_label[1:]+".head") as myfile:

             #make_LX11_musicfile(master_music_cfg_dest,ictype,seed,region_point_file)
             print "SUBMITTING INITIAL CONDITIONS..."
             cd_folder = "cd " + folder
             cmd_submit_ics = "sbatch smusic"
             subprocess.call([cd_folder+"; "+cmd_submit_ics],shell=True)
             job_name_list.append(job_name)


def run_music_higher_levels(halo_geometries,base_path,music_path,lagr_path,lx_list):
    #print halo_geometries
    #special_attention = glob.glob(base_path+"special_attention/*")
    major_mergers = glob.glob(base_path+"massive_mergers/*")

    skip_list = []

    #for halo in special_attention:
    #    skip_list.append(halo.split("/")[-1])

    for halo in major_mergers:
        skip_list.append(halo.split("/")[-1])

    #print skip_list
    for halo_name,ic_info in halo_geometries.iteritems():
        if halo_name not in ["H861036","H1599902","H1232127","H706710"]:
            geometry = ic_info.split("_")[0]
            nvir = ic_info.split("_")[1]
            for LX in ["11","12","13","14"]:
                if LX[1] in lx_list and halo_name not in skip_list:
                    folder = base_path + halo_name + "/" + halo_name +  "_"+geometry+"_Z127_P7_LN7_LX"+LX+"_O4_NV"+nvir
                    if not os.path.isdir(folder):
                        #print halo_name
                        cmd_make_next_level = "mkdir -p " + folder
                        print cmd_make_next_level
                        subprocess.call([cmd_make_next_level],shell=True)
            
            sub_suite_paths = glob.glob(base_path + "/"+halo_name+"/H*")
            run_music(sub_suite_paths,music_path,lagr_path,lx_list)

def constructresimconf(confname,boxlength,zstart,lmin,lTF,lmax,padding,overlap,refcentx,refcenty,refcentz, \
                        refextx,refexty,refextz,align,baryons,use2LPT,useLLA,omegam,omegal,omegab,hubble, \
                        sigma8,nspec,transfer,seednum,seedlevel,outformat,icfilename,fftfine,accuracy,presmooth,postsmooth, \
                        smoother,laplaceorder,gradorder,boxlevel,periodicTFstr,pointfile,boxtype,noutput,haloid,nrvir):

    f = open(confname,'w')
    f.write('[setup]' + '\n')
    f.write('boxlength            = ' + str(boxlength) + '\n')
    f.write('zstart               = ' + str(zstart) + '\n')
    f.write('levelmin             = ' + str(lmin) + '\n')
    f.write('levelmin_TF          = ' + str(lTF) + '\n')
    f.write('levelmax             = ' + str(lmax) + '\n')
    f.write('padding              = ' + str(padding) + '\n')
    f.write('overlap              = ' + str(overlap) + '\n')
    
    if boxtype == 'box':
        f.write('ref_center           = ' + str(refcentx) + ',' + str(refcenty) + ',' + str(refcentz) + '\n')
        f.write('ref_extent           = ' + str(refextx) + ',' + str(refexty) + ',' + str(refextz) + '\n')

    if boxtype == 'convex':
      f.write('region               = convex_hull' + '\n')
    else:
      f.write('region               = ' + str(boxtype) + '\n')

    #    f.write('ref_center           = ' + str(refcentx) + ',' + str(refcenty) + ',' + str(refcentz) + '\n')
    #    f.write('ref_extent           = ' + str(refextx) + ',' + str(refexty) + ',' + str(refextz) + '\n')
   
    f.write('region_point_file    = ' + str(pointfile) + '\n')
    f.write('align_top            = ' + str(align) + '\n')
    f.write('baryons              = ' + str(baryons) + '\n')
    f.write('use_2LPT             = ' + str(use2LPT) + '\n')
    f.write('use_2LLA             = ' + str(useLLA) + '\n')
    f.write('periodic_TF          = ' + str(periodicTFstr) + '\n')
    f.write('\n')
    f.write('[cosmology]'+ '\n')
    f.write('Omega_m              = ' + str(omegam) + '\n')
    f.write('Omega_L              = ' + str(omegal) + '\n')
    f.write('Omega_b              = ' + str(omegab) + '\n')
    f.write('H0                   = ' + str(hubble) + '\n')
    f.write('sigma_8              = ' + str(sigma8) + '\n')
    f.write('nspec                = ' + str(nspec) + '\n')
    f.write('transfer             = ' + str(transfer) + '\n')
    f.write('\n')
    f.write('[random]' + '\n')
    delta = int(nrvir)
    for level in range(seedlevel,int(lmax)+1):
        delta += delta*2 + 1
        if level != seedlevel:
            seeduse = int(haloid) + int(delta)
            #seeduse = seednumnew[seedi]
            f.write('seed[' + str(level) + ']              = ' + str(seeduse) + '\n')
        elif level == seedlevel:
            f.write('seed[' + str(seedlevel) + ']              = ' + str(seednum) + '\n')
            #seedi += 1

    f.write('\n')
    f.write('[output]' + '\n')
    f.write('format               = ' + str(outformat) + '\n')
    f.write('filename             = ' + str(icfilename) + '\n')
    f.write('gadget_num_files     = ' + str(noutput) + '\n')
    f.write('gadget_spreadcoarse  = yes\n')
    f.write('\n')
    f.write('[poisson]' + '\n')
    f.write('fft_fine             = ' + str(fftfine) + '\n')
    f.write('accuracy             = ' + str(accuracy) + '\n')
    f.write('pre_smooth           = ' + str(presmooth) + '\n')
    f.write('post_smooth          = ' + str(postsmooth) + '\n')
    f.write('smoother             = ' + str(smoother)  + '\n')
    f.write('laplace_order        = ' + str(laplaceorder) + '\n')
    f.write('grad_order           = ' + str(gradorder) + '\n')
    f.close()


def make_music_file(master_music_cfg_dest,boxtype,seed,region_point_file,nrvir_label):
    haloid = float(master_music_cfg_dest.split("/")[-1].split("_")[0][1:])
    folder_single = master_music_cfg_dest.split("/")[-1].replace("cfg","")
    halo_label = folder_single.split("_")[0]
    lx_label = folder_single.split("LX")[1][1:2]

    f = open(master_music_cfg_dest,'w')
    f.write("[setup]\n")
    f.write("boxlength            = 100\n")
    f.write("zstart               = 127\n")
    f.write("levelmin             = 7\n")
    f.write("levelmin_TF          = 10\n")
    f.write("levelmax             = 1"+lx_label+"\n")
    f.write("padding              = 7\n")
    f.write("overlap              = 4\n")
    f.write("region               = " + boxtype + "\n")
    
    if boxtype == 'box':
        centx,centy,centz,dx,dy,dz = getcentext(region_point_file+".head")
        
        if "_BA_" in master_music_cfg_dest:
            extx=dx
            exty=dy
            extz=dz
        if "_BB_" in master_music_cfg_dest:
            extx=dx*1.2 
            exty=dy*1.2 
            extz=dz*1.2 
        if "_BC_" in master_music_cfg_dest:
            extx=dx*1.4 
            exty=dy*1.4 
            extz=dz*1.4 
        if "_BD_" in master_music_cfg_dest:
            extx=dx*1.6 
            exty=dy*1.6 
            extz=dz*1.6 

        f.write("ref_center           = " + str(centx) + "," + str(centy) + "," + str(centz)+"\n")
        f.write("ref_extent           = " + str(extx) + "," + str(exty) + "," + str(extz)+"\n")
    
    if boxtype == 'ellipsoid':
        #print nrvir_label
        if "_EA_" in master_music_cfg_dest:
            hipadding = 1.0
        if "_EB_" in master_music_cfg_dest:
            hipadding =	1.1
        if "_EC_" in master_music_cfg_dest:
            hipadding =	1.2      
        if "_ED_" in master_music_cfg_dest:
            hipadding = 1.3  

        if "_EX_" in master_music_cfg_dest:
            hipadding = 1.05    

        f.write("hipadding            = " + str(hipadding) + "\n")

    f.write("region_point_file    = " + region_point_file + "\n")
    f.write("align_top            = no\n")
    f.write("baryons              = no\n")
    f.write("use_2LPT             = no\n")
    f.write("use_2LLA             = no\n")
    f.write("periodic_TF          = yes\n")
    f.write("\n")
    f.write("[cosmology]\n")
    f.write("Omega_m              = 0.3175\n")
    f.write("Omega_L              = 0.6825\n")
    f.write("Omega_b              = 0.049\n")
    f.write("H0                   = 67.11\n")
    f.write("sigma_8              = 0.8344\n")
    f.write("nspec                = 0.9624\n")
    f.write("transfer             = eisenstein\n")
    f.write("\n")
    f.write("[random]\n")
    f.write("seed[10]              = 34567\n")
#    f.write("seed[11]              = " + str(seed) + "\n")

    for LX in xrange(1,int(lx_label)+1):
        seed = int(haloid*(float(LX % 10)))
        f.write("seed[1"+str(LX)+"]              = " + str(seed) + "\n")
        #print LX,seed

    f.write("\n")
    f.write("[output]\n")
    f.write("format               = gadget2_double\n")
    f.write("filename             = ./ics\n")
    f.write("gadget_num_files     = 8\n")
    f.write("gadget_spreadcoarse  = yes\n")
    f.write("\n")
    f.write("[poisson]\n")
    f.write("fft_fine             = yes\n")
    f.write("accuracy             = 1e-05\n")
    f.write("pre_smooth           = 3\n")
    f.write("post_smooth          = 3\n")
    f.write("smoother             = gs\n")
    f.write("laplace_order        = 6\n")
    f.write("grad_order           = 6\n")
    f.close()


def check_is_sorted(outpath,snap=0,hdf5=True):
    #TODO: option to check all snaps
    snap = str(snap).zfill(3)
    filename = outpath+'/outputs/snapdir_'+snap+'/snap_'+snap+'.0'
    if hdf5: filename += '.hdf5'
    h = rs.snapshot_header(filename)
    try:
        if h.sorted=='yes': return True
    except:
        return False
    
def get_size(start_path = '.'):
    total_size = 0
    for dirpath, dirnames, filenames in os.walk(start_path):
        for f in filenames:
            fp = os.path.join(dirpath, f)
            total_size += os.path.getsize(fp)
    return total_size
    
def checkmakedir(folder):
    if not os.path.isdir(folder):
        mkimagedir = 'mkdir -p ' + folder
        subprocess.call([mkimagedir],shell=True)
        
def getcurrentjobs():
    pipemyq = 'squeue -o "%i,%j,%t" -u bgriffen > currentqueue.out'
    subprocess.call(';'.join([pipemyq]),shell=True)
    lines = [line.strip() for line in open('currentqueue.out')]
    subprocess.call(["rm currentqueue.out"],shell=True)
    currentjobs = []
    jobstatus = []
    jobids = []
    for i in xrange(1,len(lines)):
        jobids.append(lines[i].split(",")[0])
        currentjobs.append(lines[i].split(",")[1])
        jobstatus.append(lines[i].split(",")[2])

    return currentjobs,jobids,jobstatus

def get_unrun_halos(base_halo_path):
    folder_paths_output = []
    for folder in base_halo_path:
        need_to_run = False
        for sub_dir in glob.glob(folder+"/H*BB*NV4"):
            if os.path.isdir(sub_dir+"/outputs/snapdir_255/"):
                need_to_run = True

        if not need_to_run:
            run_folder = glob.glob(folder+"/H*BB*LX11*NV4")
            #print run_folder
            folder_paths_output.append(run_folder[0])

    return folder_paths_output
    

def make_destination_folders_clean(base_path,suite_names,lx,nrvir):
    for folder in glob.glob(base_path+"H*"):
        haloid = folder.split("halos/")[-1].split("_")[0]
        new_folder_name = haloid + "_BB_Z127_P7_LN7_LX"+str(lx)+"_O4_NV"+str(nrvir)
        #folder_single = folder_path.split("/")[-1]
        for suite in suite_names:
            cmd_make_folder = "mkdir -p " + folder +  "/contamination_suite/" + new_folder_name + "_" + suite 
            subprocess.call([cmd_make_folder],shell=True)

def make_destination_folders(base_path,suite_names,lx,nrvir):
    special_attention = glob.glob(base_path+"special_attention/*")
    major_mergers = glob.glob(base_path+"massive_mergers/*")

    skip_list = []

    for halo in special_attention:
        skip_list.append(halo.split("/")[-1])

    for halo in major_mergers:
        skip_list.append(halo.split("/")[-1])

    for folder in glob.glob(base_path+"H*"):
        haloid = folder.split("halos/")[-1].split("_")[0]
        #folder_single = folder.split("/")[-1]
        for suite in suite_names:
            if haloid not in skip_list:
                new_folder_name = haloid + "_" + suite + "_Z127_P7_LN7_LX"+str(lx)+"_O4_NV"+str(nrvir)
                #old_folder_name = haloid + "_BB_Z127_P7_LN7_LX"+str(lx)+"_O4_NV"+str(nrvir)
                #print new_folder_name
                #cmd_rm_folder = "rm -rf " + folder +  "/contamination_suite/" + old_folder_name + "_" + suite 
                #subprocess.call([cmd_rm_folder],shell=True)
                cmd_make_folder = "mkdir -p " + folder + "/contamination_suite/" + new_folder_name
                subprocess.call([cmd_make_folder],shell=True)


def construct_music_cfg(lagr_file,master_music_cfg_orig,master_music_cfg_dest,typeic):

    xpos_arr = []
    ypos_arr = []
    zpos_arr = []

    with open(lagr_file) as infile:
        for line in infile:
            xpos_arr.append(float(line.split()[0]))
            ypos_arr.append(float(line.split()[1]))
            zpos_arr.append(float(line.split()[2]))

    lagrPos = np.array([xpos_arr,ypos_arr,zpos_arr]).T
 
    xmax = np.max(lagrPos[:,0]); xmin = np.min(lagrPos[:,0])
    dx = xmax-xmin
    centx = (xmax+xmin)/2.0
    ymax = np.max(lagrPos[:,1]); ymin = np.min(lagrPos[:,1])
    dy = ymax-ymin
    centy = (ymax+ymin)/2.0
    zmax = np.max(lagrPos[:,2]); zmin = np.min(lagrPos[:,2])
    dz = zmax-zmin
    centz = (zmax+zmin)/2.0
    
    need_extent = True

    if typeic == 'A':
       extx=dx
       exty=dy
       extz=dz
    if typeic == 'B':
        extx=dx*1.2 
        exty=dy*1.2 
        extz=dz*1.2 
    if typeic == 'C':
        extx=dx*1.4 
        exty=dy*1.4 
        extz=dz*1.4 
    if typeic == 'D':
        extx=dx*1.6 
        exty=dy*1.6 
        extz=dz*1.6 
    if typeic == 'CONVEX' or typeic == 'ELLIPSOID':
        need_extent = False

    fout = open(master_music_cfg_dest,'w')
    with open(master_music_cfg_orig) as infile:
        for line in infile:
            line = line.strip()
            if need_extent:
                if "ref_center" in line:
                    fout.write("ref_center           = "+str(centx)+","+str(centy)+","+str(centz)+"\n")
                if "ref_extent" in line:
                    fout.write("ref_extent           = "+str(extx)+","+str(exty)+","+str(extz)+"\n")
                    
            if "ref_" not in line and "region " not in line:
                fout.write(line+"\n")


    if "region " in line:
        if typeic == "CONVEX":
            fout.write("region           = convex_hull\n")
        if typeic == "ELLIPSOID":
            fout.write("region           = ellipsoid\n")
        if typeic != "ELLIPSOID" and typeic != "CONVEX":
            fout.write("region               = box\n")

    fout.close()
    
def makePBSicfile(cluster,runpath,ncores,haloid,nrvir,level,email=False):
    f1 = open(runpath + "smusic",'w')
    f1.write("#!/bin/sh\n")
    f1.write("#PBS -l nodes=1:ppn=" + str(ncores) + "\n")
    f1.write("#PBS -M brendan.f.griffen@gmail.com\n")
    f1.write("#PBS -N I" + str(haloid[0:3]) + "N" + str(nrvir) + "L" + str(level[1]) +"\n")
    f1.write("#PBS -m be\n")
    f1.write("\n")
    f1.write(". /opt/torque/etc/openmpi-setup.sh\n")
    f1.write("\n")
    f1.write("cd " + runpath + "\n")
    f1.write("\n")
    f1.write("./MUSIC ./" + runpath.split("/")[-2] + ".conf 1>OUTPUTmusic 2>ERRORmusic \n")
    f1.write("rm wnoise* temp*\n")
    f1.close()

    f = open(runpath + 'runfrombigbang.sh','w')
    f.write("#!/bin/bash\n")
    f.write("ssh antares << 'ENDSSH' \n")
    f.write("cd " + runpath + "\n")
    f.write("qsub " + runpath + "smusic" + " > submitlist\n")
    f.write("logout\n")
    f.close()

def makeSLURMicfile(cluster,runpath,ncores,haloid,nrvir,level,halotype,time=5000,memory=256,queue="general",email=False):
    f1 = open(runpath + "smusic",'w')
    f1.write("#!/bin/bash \n")
    f1.write("#SBATCH -o I" + str(haloid) + "B" + halotype + "N" + str(nrvir) + "L" + str(level[1]) + ".o%j \n")
    f1.write("#SBATCH -e I" + str(haloid) + "B" + halotype + "N" + str(nrvir) + "L" + str(level[1]) + ".e%j \n")
    #f1.write("#SBATCH -o I" + str(haloid[:5]) + "L" + str(level[1]) + ".o%j \n")
    #f1.write("#SBATCH -e I" + str(haloid[:5]) + "L" + str(level[1]) + ".e%j \n")

    f1.write("#SBATCH -N 1 -n 1\n")
    f1.write("#SBATCH --exclusive\n")
    f1.write("#SBATCH -p "+ queue + "\n")

    if "harvard" in cluster:
        if int(level) == 11:
            timein = 180
            mem = 64
        elif int(level) == 12:
            timein = 360
            mem = 128
        elif int(level) == 13:
            timein = 1440
            mem = 250
        elif int(level) == 14:
            timein = 1440
            mem = 505
        elif int(level) == 15:
            timein = 1440
            mem = 505

        f1.write("#SBATCH -t " + str(time) + "\n")
        f1.write("#SBATCH --mem="+str(memory)+"gb\n")

    f1.write("#SBATCH -J I" + str(haloid) + "B" + halotype + "N" + str(nrvir) + "L" + str(level[1]) + "\n")

    if email:
        f1.write("#SBATCH --mail-user=brendan.f.griffen@gmail.com \n")
        f1.write("#SBATCH --mail-type=FAIL\n")

    f1.write("\n")
    f1.write("export OMP_NUM_THREADS=" + str(ncores) + "\n")

    if "harvard" in cluster:
        f1.write("module purge \n")
        f1.write("module load -S centos6/binutils-2.23.2 \n")
        f1.write("module load -S centos6/gcc-4.8.0 \n")
        f1.write("module load -S centos6/gmp-5.1.1 \n")
        f1.write("module load -S centos6/openmpi-1.6.4_gcc-4.8.0 \n")
        f1.write("module load -S centos6/hdf5-1.8.11_gcc-4.8.0 \n")
        f1.write("module load -S centos6/fftw-3.3.2_openmpi-1.6.4_gcc-4.8.0 \n")
        f1.write("module load -S centos6/gsl-1.16_gcc-4.8.0 \n")

    f1.write("\n")
    f1.write("cd " + runpath + "\n")
    f1.write("\n")
    f1.write("./MUSIC ./" + runpath.split("/")[-2] + ".conf 1>OUTPUTmusic 2>ERRORmusic\n")
    f1.write("rm wnoise* temp*\n")
    f1.close()

def replacetextinfile(filein,findtext,replacewith):
    f = open(filein,'r')
    filedata = f.read()
    f.close()

    newdata = filedata.replace(findtext,replacewith)

    f = open(filein,'w')
    f.write(newdata)
    f.close()

def getcaterpillarcandidates():
    candidatelist = [   41792.,   265595.,  1940914.,  1452661.,   581230.,   265449., 1233057.,   329124.,   231858.,    93549.,   135990.,   388797., 1269360.,   796218.,  1940660.,   860386.,   197257.,  1848738., 1042812.,    94940.,  1798182.,  1879738.,  1848547.,    94768., 1940362.,   581380.,  1725363.,   919116.,  1268839.,  1327707., 795674.,  1725139.,  1764199.,   485763.,  1159679.,  1600215., 1475995.,   985988.,  1327150.,  1848216.,    94562.,   830980., 1292049.,  1354579.,    65895.,  1542739.,  1041422.,  1476079., 1507355.,  1422478.,  1232333.,    65777.,   706754.,  1129930., 1194662.,   830960.,  1764135.,  1042384.,   861128.,   952175., 649524.,   485321.,    94323.,    94093.,  1327014.,   356199., 1354250.,   682120.,  1042271.,   706597.,   795050.,   767828., 65263.,   616554.,  1939826.,  1326950.,   795187.,   299792., 327580.,  1421830.,  1194206.,   484859.,   264080.,   794721., 1665066.,   951984.,   889079.,   230667.,   985754.,  1506656., 674410.,  1194083.,  1763685.,   230604.,  1878536.,   889027., 1353966.,     4847.,   263605.,  1475260.,  1542569.,   134911., 1129405.,   917969.,   767570.,  1079498.,   579945.]
    return candidatelist

def plotxyzproj(ax1,ax2,ax3,pos,format='b-'):
    ax1.plot(pos[...,0],pos[...,1],format,markeredgewidth=0.0)
    ax2.plot(pos[...,0],pos[...,2],format,markeredgewidth=0.0)
    ax3.plot(pos[...,1],pos[...,2],format,markeredgewidth=0.0)

    ax1.set_xlabel('x-pos [Mpc/h]')
    ax1.set_ylabel('y-pos [Mpc/h]')
    ax2.set_xlabel('x-pos [Mpc/h]')
    ax2.set_ylabel('z-pos [Mpc/h]')
    ax3.set_xlabel('y-pos [Mpc/h]')
    ax3.set_ylabel('z-pos [Mpc/h]')

def plotxyzprojr(ax1,ax2,ax3,pos,radius,axislabels=True,**kwargs):
    xcirc,ycirc = drawcircle(pos[...,0],pos[...,1],radius)
    ax1.plot(xcirc,ycirc,**kwargs)
    xcirc,zcirc = drawcircle(pos[...,0],pos[...,2],radius)
    ax2.plot(xcirc,zcirc,**kwargs)
    ycirc,zcirc = drawcircle(pos[...,1],pos[...,2],radius)
    ax3.plot(ycirc,zcirc,**kwargs)
    
    if axislabels:
        ax1.set_xlabel('x-pos [Mpc/h]')
        ax1.set_ylabel('y-pos [Mpc/h]')
        ax2.set_xlabel('x-pos [Mpc/h]')
        ax2.set_ylabel('z-pos [Mpc/h]')
        ax3.set_xlabel('y-pos [Mpc/h]')
        ax3.set_ylabel('z-pos [Mpc/h]')

def getillustrismp(simtype):
    if "1" in simtype:
        mp = 4.4e6
    elif "2" in simtype:
        mp = 3.5e7
    elif "3" in simtype:
        mp = 2.8e8

    return mp

def getillustrisnsnaps(simtype):
    if "1" in simtype:
        nsnaps = 3975
    elif "2" in simtype:
        nsnaps = 2264
    elif "3" in simtype:
        nsnaps = 1425

    return nsnaps

def CorrectPos(pos, box):

    com=pos.sum()/len(pos)
    index=(abs(pos[:]-com)<box/2)

    if ~(np.alltrue(index)):
        index=(pos[:]<box/2)
        if (index.sum()>=len(index)/2):
            pos[~index]=pos[~index]-box
        else:
            pos[index]=pos[index]+box

def COM(posX,posY,posZ):

    tmpX=np.float64(posX)
    tmpY=np.float64(posY)
    tmpZ=np.float64(posZ)

    return tmpX.sum()/len(tmpX), tmpY.sum()/len(tmpY), tmpZ.sum()/len(tmpZ)
        
def plotquantl(ax,x,y,label):
    ax.plot(x,y,linestyle='-',linewidth=2,label=str(label))
    #ax.plot(x,y,marker='o',linewidth=10,label=str(label))

def plotquant(ax,x,y):
    ax.plot(x,y,linestyle='-',linewidth=2)

def plotlegend(ax,legendfontsize,location=1):
    handles, labels = ax.get_legend_handles_labels()
    leg = ax.legend(handles, labels,prop={'size':legendfontsize},ncol=1,loc=location)
    leg.draw_frame(False)


def imsave2(fname, arr, vmin=None, vmax=None, cmap=None, origin=None, z='dummy'):
    from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
    from matplotlib.figure import Figure
    fig = Figure(figsize=arr.shape[::-1], dpi=1, frameon=False)
    canvas = FigureCanvas(fig)
    fig.figimage(arr, cmap=cmap, vmin=vmin, vmax=vmax, origin=origin)
    fig.savefig(fname + '.png', format='PNG',dpi=1)
    
def tick_function(X):
    V = (1/X)-1
    return ["%.1f" % z for z in V]

def calcT(p):
    return 1 - (1./p)**2

def placetext(ax,xpos,ypos,teststr,fontweight,fontsize):
    xpos = float(xpos)
    ypos = float(ypos)
    ax.text(xpos, ypos,teststr,horizontalalignment='left',verticalalignment='center',transform = ax.transAxes,fontsize=fontsize,fontweight=fontsize)

def placetext_direct(ax,xpos,ypos,teststr,fontweight,fontsize):
    xpos = float(xpos)
    ypos = float(ypos)
    ax.text(xpos, ypos,teststr,
        horizontalalignment='center',
        verticalalignment='center',
        color='black',
        fontsize=fontsize,
        weight=fontweight)

#placetext_direct(ax3,0.,0.,mstring,'bold',12)
#placetext(ax1,0.05,0.9,"LX: "+LX,'bold',12)

def getxyzdeltamcut(resolution,icgeometry):
    if icgeometry == 'ellipsoid':
        if resolution == 'l11':
            x = 51.24
            y = 48.28
            z = 47.30
            m = 1.2E12
            mgroup = 1.4e12

        if resolution == 'l12':
            pass

    if icgeometry == 'box':
        if resolution == 'l11':
            x = 50.16
            y = 48.51
            z = 47.05
            m = 1.2E12
            mgroup = 1.4e12

        if resolution == 'l12':
            pass

    return x,y,z,m,mgroup

def getcentext(filename):
    with open(filename, 'r') as fp:
        lines = []
        for i in xrange(6):
            lines.append(fp.readline().strip('#\n'))

    #print lines
    lines = map(float, lines)
    return tuple(lines)

def getlagrxyz(filename):
    x, y, z = np.loadtxt(filename, comments='#', unpack=True)
    return x,y,z

def getcandidatelist(filename):
    listin = []
    for line in open(filename,'r'):
         li=line.strip()
         if not li.startswith("#"):
             line = line.partition('#')[0]
             listin.append(np.array(line.split(' ')[0:9]))
    
    listin = np.array(listin)
    listin = listin.astype(np.float)
    return listin
    
def load_matrix_from_file(f):
    """
    This function is to load an ascii format matrix (float numbers separated by
    whitespace characters and newlines) into a numpy matrix object.
 
    f is a file object or a file path.
    """
 
    import types
    import numpy
 
    if type(f) == types.StringType:
        fo = open(f, 'r')
        matrix = load_matrix_from_file(fo)
        fo.close()
        return matrix
    elif type(f) == types.FileType:
        file_content = f.read().strip()
        file_content = file_content.replace('\r\n', ';')
        file_content = file_content.replace('\n', ';')
        file_content = file_content.replace('\r', ';')
 
        return numpy.matrix(file_content)
 
    raise TypeError('f must be a file object or a file name.')
    
def linear_fit(xdata, ydata, ysigma=None):

    """
    Performs a linear fit to data.

    Parameters
    ----------
    xdata : An array of length N.
    ydata : An array of length N.
    sigma : None or an array of length N,
        If provided, it is the standard-deviation of ydata.
        This vector, if given, will be used as weights in the fit.

    Returns
    -------
    a, b   : Optimal parameter of linear fit (y = a*x + b)
    sa, sb : Uncertainties of the parameters
    """
    
    if ysigma is None:
        w = ones(len(ydata)) # Each point is equally weighted.
    else:
        w=1.0/(ysigma**2)

    sw = sum(w)
    wx = w*xdata # this product gets used to calculate swxy and swx2
    swx = sum(wx)
    swy = sum(w*ydata)
    swxy = sum(wx*ydata)
    swx2 = sum(wx*xdata)

    a = (sw*swxy - swx*swy)/(sw*swx2 - swx*swx)
    b = (swy*swx2 - swx*swxy)/(sw*swx2 - swx*swx)
    sa = sqrt(sw/(sw*swx2 - swx*swx))
    sb = sqrt(swx2/(sw*swx2 - swx*swx))

    if ysigma is None:
        chi2 = sum(((a*xdata + b)-ydata)**2)
    else:
        chi2 = sum((((a*xdata + b)-ydata)/ysigma)**2)
    dof = len(ydata) - 2
    rchi2 = chi2/dof
    
    #print 'results of linear_fit:'
    #print '   chi squared = ', chi2
    #print '   degrees of freedom = ', dof
    #print '   reduced chi squared = ', rchi2

    return a, b, sa, sb, rchi2, dof

def general_fit(f, xdata, ydata, p0=None, sigma=None, **kw):
    """
    Pass all arguments to curve_fit, which uses non-linear least squares
    to fit a function, f, to data.  Calculate the uncertaities in the
    fit parameters from the covariance matrix.
    """
    popt, pcov = curve_fit(f, xdata, ydata, p0, sigma, **kw)

    #calculate chi-squared
    if sigma is None or sigma.min() == 0.0:
        chi2 = sum(((f(xdata,*popt)-ydata))**2)
    else:
        chi2 = sum(((f(xdata,*popt)-ydata)/sigma)**2)
    #degrees of freedom
    dof = len(ydata) - len(popt)
    #reduced chi-squared
    rchi2 = chi2/dof
    # The uncertainties are the square roots of the diagonal elements
    punc = np.zeros(len(popt))
    for i in np.arange(0,len(popt)):
        punc[i] = np.sqrt(pcov[i,i])
    return popt, punc, rchi2, dof

def readvariable(filename):
    input_file = open(filename, 'r')
    index = array('d')
    index.fromstring(input_file.read())
    return index

def writevariable(filename,inputarray,vartype):
    output_file = open(filename, 'wb')
    writearray = array(vartype, inputarray)
    writearray.tofile(output_file)
    output_file.close()
    print "Writing to...",filename

def getbinnedxy(x,y,nbins=10):
    x = np.array(x)
    y = np.array(y)
    n, _ = np.histogram(x, bins=nbins)
    sy, _ = np.histogram(x, bins=nbins, weights=y)
    sy2, _ = np.histogram(x, bins=nbins, weights=y*y)
    mean = sy / n
    return (_[1:] + _[:-1])/2,mean

def plotbinnedxy(ax,x,y,labelin='data',nbins=10,color='r'):
    x = np.array(x)
    y = np.array(y)
    n, _ = np.histogram(x, bins=nbins)
    sy, _ = np.histogram(x, bins=nbins, weights=y)
    sy2, _ = np.histogram(x, bins=nbins, weights=y*y)

    mean = sy / n
    std = np.sqrt(sy2/n - mean*mean)
    ax.errorbar((_[1:] + _[:-1])/2, mean,linewidth=4, yerr=std, fmt=color,marker='.',capsize=0,linestyle='None',ecolor=color,markerfacecolor=color,label=labelin)
    #ax.errorbar((_[1:] + _[:-1])/2, mean, fmt=color,linestyle='-',linewidth=2,label=labelin)

def getystd(x,y,nbins):
    x = np.array(x)
    y = np.array(y)
    n, _ = np.histogram(x, bins=nbins)
    sy, _ = np.histogram(x, bins=nbins, weights=y)
    sy2, _ = np.histogram(x, bins=nbins, weights=y*y)

    mean = sy / n
    std = np.sqrt(sy2/n - mean*mean)
    return std

def residuals(p, y, x):
    err = y-pval(x,p)
    return err

def func(x, a, b):
    #print 10**x
    #return a + b*(x/10**14)
    #x = 10**x/10**14
    #return a + b*x 
    return a*x + b

def getbinned_meanstd(x,y,nbins):
    n, _ = np.histogram(x, bins=nbins)
    sy, _ = np.histogram(x, bins=nbins, weights=y)
    sy2, _ = np.histogram(x, bins=nbins, weights=y*y)
    color = 'r'
    mean = sy / n
    std = np.sqrt(sy2/n - mean*mean)
    return mean,std

def plotbinwithbestfit(ax,x,y,nbins,x0,color='r',linestyle='-',labelin='data',drawline=True,drawpoints=True):
    x = np.array(x)
    y = np.array(y)
    xbin,ybin = getbinnedxy(x,y,nbins)
    std = getystd(x,y,nbins)
    mask = np.isnan(std)
    std[mask] = 0.0
    if std.min() == 0.0:
        popt, punc, rchi2, dof = general_fit(func, xbin, ybin, x0)
    else:
        popt, punc, rchi2, dof = general_fit(func, xbin, ybin, x0,std)

    #print labelin,str('{:.2f}'.format(popt[0])), rchi2
    #print str('{:.2f}'.format(popt[0])) + ' & ' + str('{:.3f}'.format(punc[0])) + ' & ' + str('{:.2f}'.format(popt[1])) + ' & ' + str('{:.3f}'.format(punc[1])) + ' & ' + str('{:.2f}'.format(rchi2))

    #print " a:", '{:.2f}'.format(popt[0]),"+-",'{:.3f}'.format(punc[0]),", b:",'{:.2f}'.format(popt[1]),"+-",'{:.3f}'.format(punc[1])
    
    if drawpoints == True:
        ax.errorbar(xbin, ybin, yerr=std, fmt=color,marker='o',capsize=0,linestyle='None',ecolor=color,markerfacecolor=color,label=labelin)
        
    if drawline == True:
        ax.plot(xbin,func(xbin,popt[0],popt[1]),color=color,linestyle=linestyle,linewidth=3,label=labelin)
    
    return popt[0],punc[0],popt[1],punc[1]

def plotxy(x,y):
    fig = plt.figure(1, figsize=(10.0,10.0))
    plt.plot(x,y)
    plt.show()

def scatterxy(x,y,s=1,c='black',xmin=0.0, xmax=0.0, ymin=0.0, ymax=0.0):
    fig = plt.figure(1, figsize=(10.0,10.0))
    ax=fig.add_subplot(1,1,1)
    ax.scatter(x,y,s=s, c=c)
    if (xmin!=0.0) | (xmax!=0.0) | (ymin!=0.0) | (ymax!=0.0):
        ax.set_xlim([xmin,xmax])
        ax.set_ylim([ymin,ymax])
        plt.show()

def addsubtitle(ax,string,color='black',box=False,boxcolor='white'):
    if box:
        ax.text(0.95, 0.05,string,
            horizontalalignment='right',
            verticalalignment='bottom',
            color=color,
            weight='bold',
            transform = ax.transAxes,
            bbox={'facecolor':boxcolor, 'alpha':0.5, 'pad':10})
    else:
        ax.text(0.95, 0.05,string,
            horizontalalignment='right',
            verticalalignment='bottom',
            color=color,
            weight='bold',
            transform = ax.transAxes)


def create31fig(figx,figy,xmin,xmax,ymin,ymax,xlabel,ylabel,title=None):
    fig = plt.figure(figsize=(figx,figy))
    ax1 = fig.add_subplot(311)
    ax2 = fig.add_subplot(312)
    ax3 = fig.add_subplot(313)
    plt.subplots_adjust(hspace=0.08)
    plt.subplots_adjust(wspace=0.08)
    #ax1.set_xticklabels([])
    #ax2.set_xticklabels([])
    #xticklabels = ax1.get_xticklabels()+ ax2.get_xticklabels()
    #plt.setp(xticklabels, visible=False)
    ax1.set_title(title)
    ax2.set_ylabel(ylabel, fontsize = 20)
    ax3.set_xlabel(xlabel, fontsize = 20)
    ax1.set_xlim([xmin,xmax])
    ax2.set_xlim([xmin,xmax])
    ax3.set_xlim([xmin,xmax])
    ax1.set_ylim([ymin,ymax])
    ax2.set_ylim([ymin,ymax])
    ax3.set_ylim([ymin,ymax])
    return fig,ax1,ax2,ax3

def create13fig(size,xmin,xmax,ymin,ymax,xlabel1,ylabel1,xlabel2,ylabel2,xlabel3,ylabel3,title=None,fontsize=12):
    fig = plt.figure(figsize=(size+15,size))
    ax1 = fig.add_subplot(131)
    ax2 = fig.add_subplot(132)
    ax3 = fig.add_subplot(133)
    plt.subplots_adjust(hspace=0.08)
    plt.subplots_adjust(wspace=0.1)
    #ax1.set_xticklabels([])
    #ax2.set_xticklabels([])
    #xticklabels = ax1.get_xticklabels()+ ax2.get_xticklabels()
    #plt.setp(xticklabels,3 visible=False)
    ax2.set_title(title)
    ax1.set_xlabel(ylabel1, fontsize = fontsize)
    ax1.set_ylabel(ylabel1, fontsize = fontsize)
    ax2.set_xlabel(ylabel2, fontsize = fontsize)
    ax2.set_ylabel(ylabel2, fontsize = fontsize) 
    ax3.set_xlabel(ylabel3, fontsize = fontsize)
    ax3.set_ylabel(ylabel3, fontsize = fontsize) 
    ax1.set_xlim([xmin,xmax])
    ax2.set_xlim([xmin,xmax])
    ax3.set_xlim([xmin,xmax])
    ax1.set_ylim([ymin,ymax])
    ax2.set_ylim([ymin,ymax])
    ax3.set_ylim([ymin,ymax])
    return fig,ax1,ax2,ax3

def create3x3fig(size,xmin,xmax,ymin,ymax,xlabel,ylabel):
    fig = plt.figure(figsize=(size,size))
    ax1 = fig.add_subplot(331)
    ax2 = fig.add_subplot(332)
    ax3 = fig.add_subplot(333)

    ax4 = fig.add_subplot(334)
    ax5 = fig.add_subplot(335)
    ax6 = fig.add_subplot(336)

    ax7 = fig.add_subplot(337)
    ax9 = fig.add_subplot(339)

    plt.subplots_adjust(hspace=0.05)
    plt.subplots_adjust(wspace=0.05)

    ticklabels = ax1.get_xticklabels() + ax2.get_xticklabels() + ax3.get_xticklabels() + ax4.get_xticklabels() + \
                 ax5.get_xticklabels() + ax6.get_xticklabels() + ax9.get_xticklabels() + \
                 ax1.get_yticklabels() + ax2.get_yticklabels() + ax3.get_yticklabels() + ax4.get_yticklabels() + \
                 ax5.get_yticklabels() + ax6.get_yticklabels() + ax9.get_yticklabels()

    plt.setp(ticklabels, visible=False)

    ax7.set_ylabel(ylabel)
    ax7.set_xlabel(xlabel)

    ax1.set_xlim([xmin,xmax])
    ax3.set_xlim([xmin,xmax])
    ax4.set_xlim([xmin,xmax])
    ax5.set_xlim([xmin,xmax])
    ax6.set_xlim([xmin,xmax])
    ax7.set_xlim([xmin,xmax])
    ax9.set_xlim([xmin,xmax])

    ax1.set_ylim([ymin,ymax])
    ax3.set_ylim([ymin,ymax])
    ax4.set_ylim([ymin,ymax])
    ax5.set_ylim([ymin,ymax])
    ax6.set_ylim([ymin,ymax])
    ax7.set_ylim([ymin,ymax])
    ax9.set_ylim([ymin,ymax])


    return fig,ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax9
    
def drawcircle(x,y,r):
    try:
        phi = np.linspace(0.0,2*np.pi,100)
        na=np.newaxis
        x_line = x[na,:]+r[na,:]*np.sin(phi[:,na])
        y_line = y[na,:]+r[na,:]*np.cos(phi[:,na])
    except:
        x = float(x)
        y = float(y)
        r = float(r)
        phi = np.linspace(0.0,2*np.pi,100)
        na=np.newaxis
        x_line = x+r*np.sin(phi)
        y_line = y+r*np.cos(phi)

    return x_line,y_line

def getillustrispath():
    import platform
    node = platform.node()
    if "harvard" in node:
        basepath = '/n/home01/bgriffen/data/'
    if node == "bigbang.mit.edu":
        basepath = '/bigbang/data/bgriffen/'
    return basepath

    return basepath

def makelegend(ax,line=0,location=1):
    handles, labels = ax.get_legend_handles_labels()
    
    if len(labels) >3: 
        handles = [h for h in handles]
        labels = [h for h in labels]
        handles = handles[0:3]
        labels = labels[0:3]
    else:
        handles = [h for h in handles]
        labels = [h for h in labels]
    
    if len(labels) >= 2:
        handles[0], handles[1] = handles[1], handles[0]
        labels[0], labels[1] = labels[1], labels[0]

    ax.legend(handles, labels, loc=location,numpoints=1,prop={'size':16})
    # test1.py executed as script
    # do something
    
def forceAspect(ax,aspect=1):
    im = ax.get_images()
    extent =  im[0].get_extent()
    ax.set_aspect(abs((extent[1]-extent[0])/(extent[3]-extent[2]))/aspect)

def makecolormap():
    vals=5020
    rgbg=np.zeros([vals,4])
    rgbg[:,3]=1.0
    lambda2 = np.zeros(vals)
    lambda2 = np.arange(0.0,1.0,1.0/(vals))
    s=1.5
    gamma=0.90
    h=1.0
    r=1.5
    phi = 2*(3.14159)*(s/3.0 + r*lambda2)
    a = h*lambda2**gamma *( 1 - lambda2**gamma) / 2.0

    for color in range(0, vals):
        rgbg[color,0] = lambda2[color]**gamma - a[color] * 0.14871 * math.cos(phi[color]) + a[color] * 1.78277 * math.sin(phi[color])
        rgbg[color,1] = lambda2[color]**gamma - a[color] * 0.29227 * math.cos(phi[color]) - a[color] * 0.90649 * math.sin(phi[color])
        rgbg[color,2] = lambda2[color]**gamma + a[color] * 1.97249 * math.cos(phi[color])

    return col.LinearSegmentedColormap.from_list('newmap',rgbg,N=vals)

def cosmoconstant(cosmology):
#    if cosmology == 'WMAP1':
#        omegam = 
#        omegal = 
#        omegab = 
#        hubble = 
#        sigma8 = 
#        nspec = 
#    
#    if cosmology == 'WMAP3':
#        omegam = 
#        omegal = 
#        omegab = 
#        hubble = 
#        sigma8 = 
#        nspec = 
#  
#    if cosmology == 'WMAP5':
#        omegam = 
#        omegal = 
#        omegab = 
#        hubble = 
#        sigma8 = 
#        nspec = 
    
    if cosmology == 'WMAP7':
        omegam = 0.276
        omegal = 0.724
        omegab = 0.045
        hubble = 70.3
        sigma8 = 0.811
        nspec = 0.961
    
#    if cosmology == 'WMAP9':
#        omegam = 
#        omegal = 
#        omegab = 
#        hubble = 
#        sigma8 = 
#        nspec = 
    
    if cosmology == 'PLANCK':
        omegam = 0.3175
        omegal = 0.6825
        omegab = 0.0489991
        hubble = 67.11
        sigma8 = 0.8344
        nspec = 0.9624

    return omegam,omegal,omegab,hubble,sigma8,nspec


def create_mt_image(hpath,header,zoomid,parentid,image_path):
    
    fig1 = plt.figure(figsize=(22.0,12.0))
    ax1 = fig1.add_subplot(3,5,1)
    ax2 = fig1.add_subplot(3,5,2)
    ax3 = fig1.add_subplot(3,5,3)
    ax4 = fig1.add_subplot(3,5,4)
    ax5 = fig1.add_subplot(3,5,5)
    ax6 = fig1.add_subplot(3,5,6)
    ax7 = fig1.add_subplot(3,5,7)
    ax8 = fig1.add_subplot(3,5,8)
    ax9 = fig1.add_subplot(3,5,9)
    ax10 = fig1.add_subplot(3,5,10)
    ax11 = fig1.add_subplot(3,5,11)
    ax12 = fig1.add_subplot(3,5,12)
    ax13 = fig1.add_subplot(3,5,13)
    ax14 = fig1.add_subplot(3,5,14)
    ax15 = fig1.add_subplot(3,5,15)

    plt.subplots_adjust(hspace=0.1,wspace=0.3)

    zoomid = int(zoomid)

    if 'parent' in hpath:
        cat = MT.MTCatalogue(hpath + '/trees',indexbyrsid=False,haloids=[zoomid])
    else:
        cat = MT.MTCatalogue(hpath + '/trees',indexbyrsid=False,haloids=[zoomid],version=4)

    tree = cat[0]

    mainbranch = tree.getMainBranch()
    
    scale = mainbranch['scale']
    rvir = mainbranch['rvir']
    posX = mainbranch['posX']/header.hubble
    posY = mainbranch['posY']/header.hubble
    posZ = mainbranch['posZ']/header.hubble
    mvir = mainbranch['mvir']/header.hubble
    mall = mainbranch['m200c_all']/header.hubble
    m200b = mainbranch['m200b']/header.hubble
    vmax = mainbranch['vmax']
    vrms = mainbranch['vrms']
    rs = mainbranch['rs']
    xoff = mainbranch['xoff']
    voff = mainbranch['voff']
    virrat = mainbranch['T/|U|']
    spin = mainbranch['spin']
    spinb = mainbranch['spin_bullock']
    pecvx = mainbranch['pecVX']
    pecvy = mainbranch['pecVY']
    pecvz = mainbranch['pecVZ']
    jx = mainbranch['Jx']
    jy = mainbranch['Jy']
    jz = mainbranch['Jz']
    last_mmr = mainbranch['scale_of_last_MM']

    try:
        bta = mainbranch['b_to_a']
        cta = mainbranch['c_to_a']
    except:
        bta = mainbranch['b_to_a(500c)']
        cta = mainbranch['c_to_a(500c)']

    try:
        Ax = mainbranch['A[x]']
        Ay = mainbranch['A[y]']
        Az = mainbranch['A[z]']
    except:
        Ax = mainbranch['A[x](500c)']
        Ay = mainbranch['A[y](500c)']
        Az = mainbranch['A[z](500c)']

    ctb = bta*(1/cta)
    Tch = calcT(ctb)

    npart = mvir/(header.massarr[1]*10**10/header.hubble)

    normmvir = mvir/mvir[0]
    normvmax = vmax**2/vmax[0]**2

    #icand += 1

    #print "---------------------------------------"
    #print "         Candidate:",icand
    print "---------------------------------------"
    print "     rockstar id: %i" % (parentid)
    print "  merger tree id: %i" % (mainbranch['id'][0])
    print "           x-pos:",'{:.2f}'.format(posX[0]), "   \ [Mpc]"
    print "           y-pos:",'{:.2f}'.format(posY[0]), "   \ [Mpc]"
    print "           z-pos:",'{:.2f}'.format(posZ[0]), "   \ [Mpc]"
    print "            vmax:",'{:.2f}'.format(vmax[0]),"  \ [km/s]"
    print "     virial mass:",'{0:.2e}'.format(mvir[0]),"\ [Msol]"
    print "   virial radius:",'{:.2f}'.format(rvir[0]),"  \ [kpc]"
    print "---------------------------------------"

    plotquantl(ax1,scale,np.log10(mvir),'virial')
    #print scale
    #plt.show()

    plotquantl(ax1,scale,np.log10(mall),'+unbound')
    plotquantl(ax1,scale,np.log10(m200b),'m200')
    plotquant(ax2,scale,normmvir)
    plotquant(ax2,scale,last_mmr)

    plotquantl(ax3,scale,vmax,'max')
    plotquantl(ax3,scale,vrms,'rms')
    plotquant(ax4,scale,normvmax)
    plotquantl(ax5,scale,rvir,'virial')  
    plotquantl(ax5,scale,rs,'scale')
    plotquantl(ax6,scale,spin,'normal')
    plotquantl(ax6,scale,spinb,'bullock')
    plotquantl(ax7,scale,xoff,'position')
    plotquantl(ax7,scale,voff,'velocity')
    plotquant(ax8,scale,virrat)
    plotquantl(ax9,scale,np.log10(npart),parentid)
    plotquantl(ax10,scale,bta,'b/a')
    plotquantl(ax10,scale,cta,'c/a')
    plotquant(ax11,scale,Tch)
    plotquantl(ax12,scale,pecvx,'Vx')
    plotquantl(ax12,scale,pecvy,'Vy')
    plotquantl(ax12,scale,pecvz,'Vz')
    plotquantl(ax13,scale,jx,'Jx')
    plotquantl(ax13,scale,jy,'Jy')
    plotquantl(ax13,scale,jz,'Jz')
    plotquantl(ax14,scale,Ax,'Ax')
    plotquantl(ax14,scale,Ay,'Ay')
    plotquantl(ax14,scale,Az,'Az')

    if os.path.isfile("H"+str(parentid)+"_subhalo_mass_fraction.txt"):
        shmf = np.loadtxt("H"+str(parentid)+"_subhalo_mass_fraction.txt",delimiter=",")
        snap_shmf = shmf[:,0]
        scale_shmf = shmf[:,1]
        shmf_shmf = shmf[:,2]
        #print scale_shmf,shmf_shmf
        plotquant(ax15,scale_shmf,shmf_shmf)

    #if icand % 3 != 0:
    #    ncols += 1
    
    new_tick_locations = np.array([.2, .5, .9])
    
    ax11.text(0.5, 0.05,'oblate', horizontalalignment='center', verticalalignment='center', color='black', fontsize=11, transform = ax11.transAxes)
    ax11.text(0.5, 0.95,'prolate', horizontalalignment='center', verticalalignment='center', color='black', fontsize=11, transform = ax11.transAxes)
    
    #ax15.text(0.5, 0.5,'spare', horizontalalignment='center', verticalalignment='center', color='black', fontsize=11, transform = ax15.transAxes)
    
    ax1.set_xlim([0,1])
    ax2.set_xlim([0,1])
    ax3.set_xlim([0,1])
    ax4.set_xlim([0,1])
    ax5.set_xlim([0,1])
    ax6.set_xlim([0,1])
    ax7.set_xlim([0,1])
    ax8.set_xlim([0,1])
    ax9.set_xlim([0,1])
    ax10.set_xlim([0,1])
    ax11.set_xlim([0,1])
    ax12.set_xlim([0,1])
    ax13.set_xlim([0,1])
    ax14.set_xlim([0,1])
    ax15.set_xlim([0,1])
    #ax15.set_ylim([0,0.3])
    ax11.set_ylim([0,1])
    
    plotlegend(ax1,10,location=4)
    plotlegend(ax3,10,location=4)
    plotlegend(ax5,10,location=7)
    plotlegend(ax6,10)
    plotlegend(ax7,10)
    plotlegend(ax10,10,location=4)
    plotlegend(ax12,10)
    plotlegend(ax13,10)
    plotlegend(ax14,10)
    
    #plt.show()

    new_tick_locations = np.array([0.,.2, .4, 0.6, 0.8, 1.])
    
    ax1top = ax1.twiny()
    ax1top.set_xticklabels(tick_function(new_tick_locations))
    
    ax2top = ax2.twiny()
    ax2top.set_xticklabels(tick_function(new_tick_locations))
    
    ax3top = ax3.twiny()
    ax3top.set_xticklabels(tick_function(new_tick_locations))
    
    ax4top = ax4.twiny()
    ax4top.set_xticklabels(tick_function(new_tick_locations))
    
    ax5top = ax5.twiny()
    ax5top.set_xticklabels(tick_function(new_tick_locations))
    
    ax1top.set_xlabel(r'$\mathrm{redshift}$',size=14)
    ax2top.set_xlabel(r'$\mathrm{redshift}$',size=14)
    ax3top.set_xlabel(r'$\mathrm{redshift}$',size=14)
    ax4top.set_xlabel(r'$\mathrm{redshift}$',size=14)
    ax5top.set_xlabel(r'$\mathrm{redshift}$',size=14)
    
    ax1.set_xticks([])
    ax2.set_xticks([])
    ax3.set_xticks([])
    ax4.set_xticks([])
    ax5.set_xticks([])
    ax6.set_xticks([])
    ax7.set_xticks([])
    ax8.set_xticks([])
    ax9.set_xticks([])
    ax10.set_xticks([])
    
    ax1.set_ylabel(r'$\mathrm{log_{10}\ M(z)\ [M_\odot]}$',size=14)
    ax2.set_ylabel(r'$\mathrm{M_v(z)/M_v(z=0)}$',size=14)
    ax3.set_ylabel(r'$\mathrm{V(z)\ [km/s]}$',size=14)
    ax4.set_ylabel(r'$\mathrm{V_{max}(z)^2/V_{max}(z=0)^2}$',size=14)
    ax5.set_ylabel(r'$\mathrm{radius\ [kpc]}$',size=14)
    ax6.set_ylabel(r'$\mathrm{spin}$',size=14)
    ax7.set_ylabel(r'$\mathrm{offset}$',size=14)
    ax8.set_ylabel(r'$\mathrm{virial\ ratio}$',size=14)
    ax9.set_ylabel(r'$\mathrm{log_{10}\ number\ of\ particles}$',size=14)
    ax10.set_ylabel(r'$\mathrm{axis\ ratios}$',size=14)
    ax11.set_ylabel(r'$\mathrm{triaxiality\ parameter}$',size=14)
    ax12.set_ylabel(r'$\mathrm{peculiar\ V_x,\ V_y,\ V_z\ [km/s]}$',size=14)
    ax13.set_ylabel(r'$\mathrm{J_x, J_y, J_z}$',size=14)
    ax14.set_ylabel(r'$\mathrm{ellipticity\ axis}$',size=14)
    ax15.set_ylabel(r'$\mathrm{subhalo\ mass\ fraction}$',size=14)

    ax11.set_xlabel(r'$\mathrm{scale\ factor}$',size=14)
    ax12.set_xlabel(r'$\mathrm{scale\ factor}$',size=14)
    ax13.set_xlabel(r'$\mathrm{scale\ factor}$',size=14)
    ax14.set_xlabel(r'$\mathrm{scale\ factor}$',size=14)
    ax15.set_xlabel(r'$\mathrm{scale\ factor}$',size=14)

    fig1.savefig(image_path + 'H'+str(parentid)+'_all_mt_properties.png',bbox_inches='tight')
    plt.close(fig1)

    # -----------------------------------------------------------------

    fig2 = plt.figure(figsize=(15.0,6.0))
    ax1a = fig2.add_subplot(1,2,1)
    ax2a = fig2.add_subplot(1,2,2)

    plotquant(ax2a,scale,normvmax)
    plotquantl(ax1a,scale,normmvir,parentid)
    
    ax1a.set_xlim([0,1])
    ax2a.set_xlim([0,1])
    ax1a.set_ylim([0,1.1])
    ax2a.set_ylim([0,1.1])

    ax1a.set_ylabel(r'$\mathrm{M_v(z)/M_v(z=0)}$',size=14)
    ax2a.set_ylabel(r'$\mathrm{V_{max}(z)^2/V_{max}(z=0)^2}$',size=14)
    
    ax1atop = ax1a.twiny()
    ax1atop.set_xticklabels(tick_function(new_tick_locations))
    
    ax2atop = ax2a.twiny()
    ax2atop.set_xticklabels(tick_function(new_tick_locations))
    
    ax1atop.set_xlabel(r'$\mathrm{redshift}$',size=14)
    ax2atop.set_xlabel(r'$\mathrm{redshift}$',size=14)
    ax1a.set_xlabel(r'$\mathrm{scale\ factor}$',size=14)
    ax2a.set_xlabel(r'$\mathrm{scale\ factor}$',size=14)

    handles, labels = ax1a.get_legend_handles_labels()
    ax1a.legend(handles, labels, fontsize=9,loc=4,frameon=False)
    fig2.savefig(image_path + 'all_stacked_mvir_vmax.png',bbox_inches='tight')
    
def get_subhalo_mass_fraction(halodata,haloid):
    mhost = float(halodata.ix[haloid]['mgrav'])
    subhalos = halodata.get_all_subhalos_within_halo(haloid)
    if len(subhalos) > 0:
        submass = np.array(subhalos['mgrav']).sum()
        fsub = submass/mhost
    else:
        fsub = 0.0

    return fsub,len(subhalos)
    
def get_min_distance_to_contam(halopath,part_type=2):
    halox = htils.get_quant_zoom(halopath,'x')
    haloy = htils.get_quant_zoom(halopath,'y')
    haloz = htils.get_quant_zoom(halopath,'z')
    pos = rs.read_block(halopath + "/outputs/snapdir_255/snap_255","POS ",parttype=part_type)
    head = rs.snapshot_header(halopath + "/outputs/snapdir_255/snap_255.0.hdf5")
    dx = halox - pos[:,0]
    dy = haloy - pos[:,1]
    dz = haloz - pos[:,2]
    R = np.sqrt(dx**2+dy**2+dz**2)*1000./head.hubble
    return np.min(R)
    
def get_force_res(parameter_file):
    with open(parameter_file) as f:
        for line in f:
            if "SofteningHaloMaxPhys" in line:
                force_res = line.split()[-1]
    return force_res

def get_ncores(node_name):
    if node_name == "AMD":
        ncores = 64
    if node_name == "HyperNodes":
        ncores = 24
    if node_name == "RegNodes":
        ncores = 8
    return ncores

def write_SLURM_file(job_name,node_name,cfg_file,rsdir,halo_path,restart=False):
    f = open(halo_path + 'rockstar.sbatch','w')
    f.write('#!/bin/sh\n')
    f.write('#SBATCH -J ' + job_name + '\n')
    f.write('#SBATCH -o rockstar.o\n')
    f.write('#SBATCH -e rockstar.e\n')
    f.write('#SBATCH -p ' + node_name + '\n')
    f.write('#SBATCH -N 1\n')
    f.write('#SBATCH -t infinite\n')
    f.write('#SBATCH --exclusive\n')
    f.write('\n')

    f.write('rsdir=/home/bgriffen/data/lib/Rockstar-0.99.9-RC3/\n')
    f.write('exe=$rsdir/rockstar\n')
    sub.call(['mkdir -p ' + halo_path + 'halos/'],shell=True)
    f.write('outdir=' + halo_path + 'halos/\n')
    f.write('cd $rsdir\n')

    if restart:
        f.write('$exe -c $outdir/restart.cfg &\n')
    else:
        f.write('$exe -c $outdir/' + cfg_file + ' &\n')
    
    ncores = get_ncores(node_name)
    
    f.write('cd $outdir\n')
    #f.write('echo $outdir\n')
    f.write('''perl -e 'sleep 1 while (!(-e "auto-rockstar.cfg"))'\n''')
    f.write('srun -n ' + str(ncores) + ' $exe -c auto-rockstar.cfg')
    f.close()

def get_num_snaps(halo_path):
    dir_list = glob.glob(halo_path + 'outputs/snapdir_*')
    nsnaps = len(dir_list)
    return nsnaps

def get_nfiles(parameter_file):
    with open(parameter_file) as f:
        for line in f:
            if "NumFilesPerSnapshot" in line:
                nfiles = line.split()[-1]
    return nfiles

def write_rockstar_cfg(cfg_file,halo_path,num_blocks,num_writers,file_format,force_res):
    nsnaps = get_num_snaps(halo_path)

    f = open(halo_path + '/halos/rockstar.cfg','w')
    f.write('PARALLEL_IO=1\n')
    f.write('INBASE=' + halo_path + 'outputs' + '\n')
    f.write('OUTBASE=' + halo_path + 'halos\n')

    f.write('NUM_BLOCKS='+ str(num_blocks) + '\n')
    f.write('NUM_WRITERS=' + str(num_writers) + '\n')
    f.write('FORK_READERS_FROM_WRITERS=1\n')
    f.write('FILE_FORMAT="' + file_format + '"\n')
    #print file_format.lower()
    if file_format.lower() == "arepo":
        f.write('FILENAME=snapdir_<snap>/snap_<snap>.<block>.hdf5\n')
        f.write('AREPO_LENGTH_CONVERSION = 1\n')
        #f.write('AREPO_LENGTH_CONVERSION = 1e-3\n')
        f.write('AREPO_MASS_CONVERSION = 1e+10\n')

    if file_format.lower() == "gadget":
        f.write('FILENAME=snapdir_<snap>/snap_<snap>.<block>\n')
        f.write('GADGET_LENGTH_CONVERSION = 1\n')
        f.write('GADGET_MASS_CONVERSION = 1e+10\n')
        
    #f.write('NUM_SNAPS=' + str(nsnaps) + '\n')
    f.write('SNAPSHOT_NAMES= ' + halo_path + 'halos/snapshotlist.dat\n')
    f.write('FULL_PARTICLE_CHUNKS=0\n')
    f.write('FORCE_RES=' + str(force_res) + '\n')
    f.write('FULL_PARTICLE_BINARY=' + str(num_writers) + '\n')
    f.write('OUTPUT_FORMAT="BINARY"\n')
    f.write('MASS_DEFINITION="vir"\n')
    f.write('DELETE_BINARY_OUTPUT_AFTER_FINISHED=1\n')
    f.close()

def run_rockstar(halo_path,rsdir,ctrees,node_name="RegNodes",cfg_file="rockstar.cfg",file_format="AREPO",job_name="rockstar"):

    parameter_file = halo_path + "param.txt"
    
    print "Running:",halo_path
    nsnaps = get_num_snaps(halo_path)
    job_name_list = []
    if not os.path.isfile(halo_path+"halos/halos_"+str(nsnaps-1)+".0.fullbin") \
        and not os.path.isfile(halo_path + "halos/trees/tree_0_0_0.dat"):

        force_res = get_force_res(parameter_file)
        num_blocks = get_nfiles(parameter_file)
        num_writers = get_ncores(node_name)
        
        write_SLURM_file(job_name,node_name,cfg_file,rsdir,halo_path)
        write_rockstar_cfg(cfg_file,halo_path,num_blocks,num_writers,file_format,force_res)
        
        f = open(halo_path + "/halos/snapshotlist.dat","w")
        for snap in xrange(0,nsnaps):
            f.write(str(snap).zfill(3)+'\n')
        f.close()

        current_jobs,jobids,jobstatus = getcurrentjobs()
        run_rockstar = "cd " + halo_path + "; sbatch rockstar.sbatch"
        
        if job_name not in current_jobs and job_name not in job_name_list:
            sub.call([run_rockstar],shell=True)

        job_name_list.append(job_name)

def run_consistent_trees(halo_path,rsdir,ctrees):
    nsnaps = get_num_snaps(halo_path)
    if not os.path.isfile(halo_path + "halos/trees/tree_0_0_0.dat"):
        cmd1 = "perl " + rsdir + "/scripts/gen_merger_cfg.pl " + halo_path + "halos/rockstar.cfg"
        cmd2 = "cd " + ctrees
        cmd3 = "perl do_merger_tree.pl " + halo_path + "halos/outputs/merger_tree.cfg"
        sub.call(';'.join([cmd1,cmd2,cmd3]), shell=True)

    zfillhalos = 3
    
    if not os.path.isdir(halo_path + "halos/halos_"+str(nsnaps-1)):
        print "WILL BE RUNNING TO SNAPSHOT:",nsnaps
        for i in range(0,nsnaps):
            print "Moving snapshot...",i
            cmd_make_halos_dir = "mkdir -p " + halo_path + "halos/halos_" + str(i).zfill(zfillhalos)
            cmd_mv_rockstar_files = "mv " + halo_path + "halos/halos_" + str(i).zfill(zfillhalos) + ".* " + halo_path + "halos/halos_" + str(i).zfill(zfillhalos)
            cmd_mv_outlists = "mv " + halo_path + "halos/out_" + str(i) + ".list " + halo_path + "halos/halos_" + str(i).zfill(zfillhalos)
            sub.call(';'.join([cmd_make_halos_dir,cmd_mv_rockstar_files,cmd_mv_outlists]), shell=True)
    
def run_parents_list(halo_path,rsdir):
    nsnaps = get_num_snaps(halo_path)
    zfillhalos = 3
    if not os.path.isfile(halo_path + "halos/halos_"+str(nsnaps-1).zfill(zfillhalos)+"/parents.list"):
        print "CONSTRUCTING PARENTS LIST"
        for i in range(0,nsnaps):
            print "Constructing snapshot...",i
            cmd_find_parents = rsdir + "/util/find_parents " + halo_path + "halos/halos_" + str(i).zfill(zfillhalos) + "/out_" + str(i) + ".list > " + halo_path + "halos/halos_" + str(i).zfill(zfillhalos) + "/parents.list"
            sub.call([cmd_find_parents],shell=True)

def run_convert_mt(halo_path):
    if not os.path.isfile(halo_path + "halos/trees/treeindex.csv"):
        print "CREATING MERGER TREE INDEX FILE"
        MT.convertmt(halo_path + 'halos/trees/',version=4)
        print "MERGER TREE READY FOR ACTION!"

def mask_positions(x,y,z,center,delta,want_quad=False):
    if want_quad:
        cond1h = (np.array(x) >= 0) & (np.array(x) <= center[0] + delta)
        cond2h = (np.array(y) >= 0) & (np.array(y) <= center[1] + delta)
        cond3h = (np.array(z) >= 0) & (np.array(z) <= center[2] + delta)
    else:
        cond1h = (np.array(x) >= center[0] - delta) & (np.array(x) <= center[0] + delta)
        cond2h = (np.array(y) >= center[1] - delta) & (np.array(y) <= center[1] + delta)
        cond3h = (np.array(z) >= center[2] - delta) & (np.array(z) <= center[2] + delta)

    return cond1h & cond2h & cond3h

def adjust_render_map(map):
    #ma=map.max()/2.
    #mi=
    #ma=np.log10(ma)

    map=np.log10(map)
    ma=map.max()
    mi=map.min()
    print "Min heatmap:",mi
    print "Max heatmap:",ma
    #floor = ma*0.8
    #map[map>floor] = floor
    #map[map>ma]=ma
    return np.array(map)
    #map_proj = np.array(map).T

def disigmoidScaling(values, steepnessFactor=1, ref=None):
    ''' Sigmoid scaling in which values around a reference point are flattened
    arround a reference point

    Scaled value y is calculated as 
        y = sign(v - d)(1 - exp(-((x - d)/s)**2)))
    where v is the original value,  d is the referenc point and s is the 
    steepness factor
    '''
    if ref is None:
        mn = np.min(values)
        mx = np.max(values)
        ref = mn + (mx - mn) / 2.0

    sgn = np.sign(values - ref)
    term1 = ((values - ref)/steepnessFactor) ** 2
    term2 = np.exp(- term1) 
    term3 = 1.0 - term2 
    return sgn * term3

def render(ax,map,cmap,img_width):
    ax.imshow(map,origin='lower',cmap=cmap, aspect='auto',interpolation='gaussian',extent=[-img_width/2., img_width/2., -img_width/2., img_width/2.])

def make_map(filename,pos,hsmlvals,mass,img_res,img_width,px,py,cmap):

    center = np.array([0,0,0], dtype="float32")
    map = SphMap.CalcDensProjection(pos, hsmlvals, mass, mass, img_res, img_res, img_width, img_width, img_width, center[0], center[1], center[2], px, py, 3, 0)
    map = adjust_render_map(map).T
    np.savetxt(filename,map)

    #map = disigmoidScaling(map, steepnessFactor=4)

def get_projection_bool(projection):
    px = projection.find("x")+1
    py = projection.find("y")+1
    pz = projection.find("z")+1

    idx_all = ["x","y","z"]
    p_all = np.array([px,py,pz])
    pabs = np.array([0,1,2])
    p_mask = [p_all != 0]
    p_use = pabs[p_mask]
    return p_use[0],p_use[1]

def get_pos_mass_hsml(hpath,zoomid,snapshot,img_width_nrvir,img_width=0,nth_halo=0,want_just_img_width=False):

    header = htils.get_halo_header(hpath,snap=snapshot)

    halos = htils.load_rscat(hpath,snapshot)
    pos_host = np.array([halos.ix[zoomid]['posX'],halos.ix[zoomid]['posY'],halos.ix[zoomid]['posZ']])

    img_width = img_width_nrvir*2.*float(halos.ix[zoomid]['rvir'])/1000.

    if want_just_img_width:
        return img_width

    print "Number of particles: in snapshot",header.npart[1]

    t0 = time.clock()
    pos = htils.load_partblock(hpath,255,"POS ",parttype=1).astype('float32')
    mass = htils.load_partblock(hpath,255,"MASS",parttype=1).astype('float32')*1e10/header.hubble
    ids = htils.load_partblock(hpath,255,"ID  ",parttype=1)
    t1 = time.clock()

    print "[ > reading positions, masses, smoothing lengths: %3.2f minutes ]" % ((t1-t0)/60.)
    
    print "Getting image width from halo's rvir..."


    if nth_halo != 0:
        idx = np.argsort(np.array(halos['mvir']))
        zoomid = int(np.array(halos['id'])[idx[-nth_halo]])

    #print "Mvir: %3.2e Msun" % (float(halos.ix[zoomid]['mvir'])/header.hubble)
    #print "Rvir: %3.2f kpc" % (halos.ix[zoomid]['rvir'])
    #rint "Npart: %i" % (halos.ix[zoomid]['npart'])

    

    print "Getting image width from input..."

    #cat = htils.load_scat(hpath)
    #pos_host = cat.sub_pos[14]
    #img_width = 10./1000.

    t0 = time.clock()
    pos -= pos_host
    center = np.array([0,0,0], dtype="float32")
    mask_pos = mask_positions(pos[:,0],pos[:,1],pos[:,2],center,img_width)
    pos = pos[mask_pos]
    mass = mass[mask_pos]
    ids = ids[mask_pos]
    t1 = time.clock()
    #print "NSUM:",np.sum(mask_pos)
    print "[ > adjusting positions: %3.2f minutes ]" % ((t1-t0)/60.)

    t0 = time.clock()
    print "Do smoothing lengths exist?",os.path.isdir(hpath+"/outputs/hsmldir_255")
    if os.path.isdir(hpath+"/outputs/hsmldir_255"):
        hsml = readhsml.hsml_file(hpath+"/outputs/",255)
        hsmlvals = hsml.Hsml[ids]
    else:
        if not os.path.isdir(hpath+"/analysis/hsmls.dat"):
            print "Calculating smoothing lengths..."
            hsmlvals=ca.CalcHsml(pos, mass, 4, float(0), float(0.0), float(1.0), float(1.0), float(0.0))
            np.save(hpath+"/analysis/hsmls.dat",hsmlvals)
        else:
            print "Loading them from file: analysis/hsml.dat"
            hsmlvals = np.load(hpath+"/analysis/hsmls.dat")
        

    t1 = time.clock()


    
    #print "LEN(pos)",len(pos)
    #print "LEN(mass)",len(mass)
    #print "LEN(hsmlvals)",len(hsmlvals)

    print "[ > getting smoothing lengths: %3.2f minutes ]" % ((t1-t0)/60.)

    print "Image width: %3.2f [kpc]" % (img_width*1000.)
    print "Number of particles to be plotted:",len(pos)

    return pos,mass,hsmlvals,img_width

def overlay_rockstar_single(ax1,projection,hpath,zoomid,snapshot,render_type,img_width,nth_halo,sorted_by='vmax',want_label_on_halos=True):

    px,py = get_projection_bool(projection)

    halos = htils.load_rscat(hpath,snapshot)
    #rsid = htils.load_zoomid(hpath)

    x = np.array(halos['posX'])
    y = np.array(halos['posY'])
    z = np.array(halos['posZ'])

   # mask_pos = mask_positions(x,y,z,center,img_width)

    rvir = np.array(halos['rvir'])/1000.
    header = htils.get_halo_header(hpath,snap=snapshot)

    if sorted_by=='vmax': circles_cut = np.array(halos['vmax'])

    if sorted_by=='mvir': circles_cut = np.array(halos['mvir'])/header.hubble

    #if img_width == 0:

    halos = htils.load_rscat(hpath,snapshot)
    #img_width = img_width_nrvir*2.*float(halos.ix[rsid]['rvir'])/1000.

    if nth_halo == 0:
        pos_host = np.array([halos.ix[zoomid]['posX'],halos.ix[zoomid]['posY'],halos.ix[zoomid]['posZ']])
    else:
        idx = np.argsort(np.array(halos['mvir']))
        zoomid = int(np.array(halos['id'])[idx[-nth_halo]])
        pos_host = np.array([halos.ix[zoomid]['posX'],halos.ix[zoomid]['posY'],halos.ix[zoomid]['posZ']])

    #else:
    
    #cat = htils.load_scat(hpath)
    #pos_host = cat.sub_pos[14]
    #img_width = 10./1000.

    dx = x - pos_host[0]
    dy = y - pos_host[1]
    dz = z - pos_host[2]


    if sorted_by=='vmax': mask_halos = mask_positions(dx,dy,dz,np.array([0,0,0]),img_width) & (circles_cut >= 30)

    if sorted_by=='mvir': mask_halos = mask_positions(dx,dy,dz,np.array([0,0,0]),img_width) & (circles_cut >= 1e9)

    dx = dx[mask_halos]
    dy = dy[mask_halos]
    dz = dz[mask_halos]
    rvir = rvir[mask_halos]

    if px == 0: dx_use = dx
    if px == 1: dx_use = dy
    if px == 2: dx_use = dz

    if py == 0: dy_use = dx
    if py == 1: dy_use = dy
    if py == 2: dy_use = dz

    if want_label_on_halos:
        for item,xi,yi in zip(circles_cut[mask_halos],dx_use,dy_use):
            if xi > -img_width/2. and xi < img_width/2. and yi > -img_width/2. and yi < img_width/2.:
                mstring = "%3.0f" % (item)
                ax1.text(xi, yi,mstring, horizontalalignment='center',verticalalignment='center',weight='bold',fontsize=12,color='cyan')
        

    xcirc,ycirc = drawcircle(dx_use,dy_use,rvir)
    ax1.plot(xcirc,ycirc,color='cyan',linewidth=2)

    ax1.set_xlim([-img_width/2., img_width/2.])
    ax1.set_ylim([-img_width/2., img_width/2.])

    return ax1

def overlay_rockstar(ax1,ax2,ax3,hpath,zoomid,snapshot,render_type,img_width,nth_halo,sorted_by='vmax',want_label_on_halos=False):

    halos = htils.load_rscat(hpath,snapshot)
    #rsid = htils.load_zoomid(hpath)

    x = np.array(halos['posX'])
    y = np.array(halos['posY'])
    z = np.array(halos['posZ'])

   # mask_pos = mask_positions(x,y,z,center,img_width)

    rvir = np.array(halos['rvir'])/1000.
    header = htils.get_halo_header(hpath,snap=snapshot)


    #if img_width == 0:

    #halos = htils.load_rscat(hpath,snapshot)
    #img_width = img_width_nrvir*2.*float(halos.ix[rsid]['rvir'])/1000.

    if nth_halo == 0:
        pos_host = np.array([halos.ix[zoomid]['posX'],halos.ix[zoomid]['posY'],halos.ix[zoomid]['posZ']])
    else:
        idx = np.argsort(np.array(halos['mvir']))
        zoomid = int(np.array(halos['id'])[idx[-nth_halo]])
        pos_host = np.array([halos.ix[zoomid]['posX'],halos.ix[zoomid]['posY'],halos.ix[zoomid]['posZ']])

    #else:
    
    #cat = htils.load_scat(hpath)
    #pos_host = cat.sub_pos[14]
    #img_width = 10./1000.

    dx = x - pos_host[0]
    dy = y - pos_host[1]
    dz = z - pos_host[2]

    mask_halos_pos = mask_positions(dx,dy,dz,np.array([0,0,0]),img_width)

    subset = np.array(halos[sorted_by])
    if sorted_by=='vmax': 
        mask_halos = (mask_halos_pos) & (subset >= 30) #& (subset <= 25)

    if sorted_by=='mvir':
        mask_halos = (mask_halos_pos) & (subset >= 1e9)

    #mask_test = (subset >= 30)

    #print np.array(halos['vmax'])[mask_test]


    dx = dx[mask_halos]
    dy = dy[mask_halos]
    dz = dz[mask_halos]
    rvir = rvir[mask_halos]
    vmax = subset[mask_halos]

    #print vmax

    if want_label_on_halos:
        for item,xi,yi in zip(vmax,dx,dy):
            if xi > -img_width/2. and xi < img_width/2. and yi > -img_width/2. and yi < img_width/2.:
                mstring = "%3.0f" % (item)
                ax1.text(xi, yi, mstring, horizontalalignment='center',verticalalignment='center',weight='bold',fontsize=12,color='cyan')
        
        for item,xi,yi in zip(vmax,dx,dz):
            if xi > -img_width/2. and xi < img_width/2. and yi > -img_width/2. and yi < img_width/2.:
                mstring = "%3.0f" % (item)
                ax2.text(xi, yi, mstring, horizontalalignment='center',verticalalignment='center',weight='bold',fontsize=12,color='cyan')
        
        for item,xi,yi in zip(vmax,dy,dz):
            if xi > -img_width/2. and xi < img_width/2. and yi > -img_width/2. and yi < img_width/2.:
                mstring = "%3.0f" % (item)
                ax3.text(xi, yi, mstring, horizontalalignment='center',verticalalignment='center',weight='bold',fontsize=12,color='cyan')
    

    xcirc,ycirc = drawcircle(dx,dy,rvir)
    ax1.plot(xcirc,ycirc,color='cyan',linewidth=2)

    xcirc,zcirc = drawcircle(dx,dz,rvir)
    ax2.plot(xcirc,zcirc,color='cyan',linewidth=2)

    ycirc,zcirc = drawcircle(dy,dz,rvir)
    ax3.plot(ycirc,zcirc,color='cyan',linewidth=2)

    ax1.set_xlim([-img_width/2., img_width/2.])
    ax1.set_ylim([-img_width/2., img_width/2.])
    ax2.set_xlim([-img_width/2., img_width/2.])
    ax2.set_ylim([-img_width/2., img_width/2.])
    ax3.set_xlim([-img_width/2., img_width/2.])
    ax3.set_ylim([-img_width/2., img_width/2.])

    return ax1,ax2,ax3

def create_zoom_render(hpath,zoomid,lx=11,nth_halo=0,img_filename='density_projection.png',snapshot=255,cmap_name='hot',img_width_nrvir=1,img_width_force=0,img_res=512,render_type='single',projection='xy',resolution_list=[11,12,13,14],include_rockstar=False,include_subfind=False,want_label_on_halos=False):
    
    if "cat" in cmap_name:
        cmap = makecolormap()
    else:
        cmap = cmap_name

    t0_total =  time.clock()
    pid = int(hpath.split("/")[-2][1:])
    
    print
    print "==============================="
    print "Parent ID:",pid
    if nth_halo == 0: print "Zoom ID:",zoomid
    if nth_halo != 0: print "Doing nth:",nth_halo
    print "LX:",lx
    print "Snapshot:",snapshot
    print "X*Rvir:",img_width_nrvir
    print "Image Resolution:",img_res
    print "Render Type:",render_type
    print "Projection:",projection
    print "Forced Image Width [kpc]:",img_width_force
    print "Using halosdir: halos_bound/"
    print "Overlay Rockstar?",include_rockstar
    print "Overlay Rockstar Labels?",want_label_on_halos
    print "Overlay SUBFIND?",include_subfind
    print "Image Output:",img_filename
    print "==============================="

    header = htils.get_halo_header(hpath,snap=snapshot)

    if render_type == 'single':
        pos,mass,hsmlvals,img_width = get_pos_mass_hsml(hpath,zoomid,snapshot,img_width_nrvir,img_width_force,nth_halo)
        fig, ax1 = plt.subplots(1, 1, figsize=(12,12))
        px,py = get_projection_bool(projection)

        file_name = hpath+"/analysis/H"+str(pid)+"_LX"+str(lx)+"_DIM"+str(img_res)+"_MAP_KPC"+str(int(img_width*1000.))+"_P"+str(int(px))+str(int(py))+"_SN"+str(snapshot)+".dat"

        if not os.path.isfile(file_name): make_map(file_name,pos,hsmlvals,mass,img_res,img_width,px,py,cmap)
 
        map = np.loadtxt(file_name)
        render(ax1,map,cmap,img_width)

        if include_rockstar: overlay_rockstar_single(ax1,projection,hpath,zoomid,snapshot,render_type,img_width,nth_halo,sorted_by='vmax',want_label_on_halos=want_label_on_halos)
 
        ax1.set_xlabel(projection[0] + " [Mpc]")
        ax1.set_ylabel(projection[1] + " [Mpc]")
        ax1.text(0.05, 0.95,'LX:'+str(lx), horizontalalignment='left',verticalalignment='center',transform = ax1.transAxes,weight='bold',fontsize=14,color='w')
        
    elif render_type == 'multi':
 
        fig, axs = plt.subplots(1, 3, figsize=(15,5))
        fig.subplots_adjust(hspace=0,wspace=0)

        ax1 = axs[0]
        ax2 = axs[1]
        ax3 = axs[2]

        hpath = htils.get_hpath_lx(pid,lx)
        #zoomid = htils.load_zoomid(hpath)

        img_width = get_pos_mass_hsml(hpath,zoomid,snapshot,img_width_nrvir,img_width_force,nth_halo,want_just_img_width=True)
        file_name_xy = hpath+"/analysis/H"+str(pid)+"_LX"+str(lx)+"_DIM"+str(img_res)+"_MAP_KPC"+str(int(img_width*1000.))+"_P01_SN"+str(snapshot)+".dat"
        file_name_xz = hpath+"/analysis/H"+str(pid)+"_LX"+str(lx)+"_DIM"+str(img_res)+"_MAP_KPC"+str(int(img_width*1000.))+"_P02_SN"+str(snapshot)+".dat"
        file_name_yz = hpath+"/analysis/H"+str(pid)+"_LX"+str(lx)+"_DIM"+str(img_res)+"_MAP_KPC"+str(int(img_width*1000.))+"_P12_SN"+str(snapshot)+".dat"

        t0 = time.clock()
        
        #if not os.path.isfile(file_name_xy) or not os.path.isfile(file_name_xz) or not os.path.isfile(file_name_yz):
        pos,mass,hsmlvals,img_width = get_pos_mass_hsml(hpath,zoomid,snapshot,img_width_nrvir,img_width_force,nth_halo)

        t1 = time.clock()
        
        print "[ reading positions, masses, smoothing lengths: %3.2f minutes ]" % ((t1-t0)/60.)            
        
        t0 = time.clock()
        if not os.path.isfile(file_name_xy): make_map(file_name_xy,pos,hsmlvals,mass,img_res,img_width,0,1,cmap) # xy
        if not os.path.isfile(file_name_xz): make_map(file_name_xz,pos,hsmlvals,mass,img_res,img_width,0,2,cmap) # xz
        if not os.path.isfile(file_name_yz): make_map(file_name_yz,pos,hsmlvals,mass,img_res,img_width,1,2,cmap) # yz

        map_xy = np.loadtxt(file_name_xy)
        render(ax1,map_xy,cmap,img_width)
        map_xz = np.loadtxt(file_name_xz)
        render(ax2,map_xz,cmap,img_width)
        map_yz = np.loadtxt(file_name_yz)
        render(ax3,map_yz,cmap,img_width)
        
        t1 = time.clock()
        print "[ rendering to figure: %3.2f minutes ]" % ((t1-t0)/60.)
        
        ax1.set_ylabel('(y,z,z) [Mpc]')
        ax1.set_xlabel('x [Mpc]')
        ax2.set_xlabel('x [Mpc]')
        ax3.set_xlabel('y [Mpc]')

        ax2.set_yticklabels([])
        ax3.set_yticklabels([])
        ax2.set_yticks([])
        ax3.set_yticks([])

        if include_rockstar: overlay_rockstar(ax1,ax2,ax3,hpath,zoomid,snapshot,render_type,img_width,nth_halo,sorted_by='vmax',want_label_on_halos=want_label_on_halos)
        if include_subfind: overlay_rockstar(ax1,ax2,ax3,hpath,zoomid,snapshot,render_type,img_width,nth_halo,sorted_by='vmax')

        ax1.text(0.9, 0.05,'XY', horizontalalignment='left',verticalalignment='center',transform = ax1.transAxes,weight='bold',fontsize=14,color='w')
        ax2.text(0.9, 0.05,'XZ', horizontalalignment='left',verticalalignment='center',transform = ax2.transAxes,weight='bold',fontsize=14,color='w')
        ax3.text(0.9, 0.05,'YZ', horizontalalignment='left',verticalalignment='center',transform = ax3.transAxes,weight='bold',fontsize=14,color='w')
        ax1.text(0.05, 0.9,'LX:'+str(lx), horizontalalignment='left',verticalalignment='center',transform = ax1.transAxes,weight='bold',fontsize=14,color='w')
    
    elif render_type == 'convergence':
        fig, axs = plt.subplots(4, 3, figsize=(15,15))
        fig.subplots_adjust(hspace=0,wspace=0)

        for axi,lx in enumerate(resolution_list):
            print
            print "Running > %i" % (lx)

            ax1 = axs[axi,0]
            ax2 = axs[axi,1]
            ax3 = axs[axi,2]
            
            hpath = htils.get_hpath_lx(pid,lx)
            zoomid = htils.load_zoomid(hpath)

            img_width = get_pos_mass_hsml(hpath,zoomid,snapshot,img_width_nrvir,img_width_force,nth_halo,want_just_img_width=True)

            file_name_xy = hpath+"/analysis/H"+str(pid)+"_LX"+str(lx)+"_DIM"+str(img_res)+"_MAP_KPC"+str(int(img_width*1000.))+"_P01_SN"+str(snapshot)+".dat"
            file_name_xz = hpath+"/analysis/H"+str(pid)+"_LX"+str(lx)+"_DIM"+str(img_res)+"_MAP_KPC"+str(int(img_width*1000.))+"_P02_SN"+str(snapshot)+".dat"
            file_name_yz = hpath+"/analysis/H"+str(pid)+"_LX"+str(lx)+"_DIM"+str(img_res)+"_MAP_KPC"+str(int(img_width*1000.))+"_P12_SN"+str(snapshot)+".dat"

            t0 = time.clock()

            if not os.path.isfile(file_name_xy) or not os.path.isfile(file_name_xz) or not os.path.isfile(file_name_yz):
                pos,mass,hsmlvals,img_width = get_pos_mass_hsml(hpath,zoomid,snapshot,img_width_nrvir,img_width_force,nth_halo)

            t1 = time.clock()
            print "[ reading positions, masses, smoothing lengths: %3.2f minutes ]" % ((t1-t0)/60.)

            
            t0 = time.clock()

            if not os.path.isfile(file_name_xy): make_map(file_name_xy,pos,hsmlvals,mass,img_res,img_width,0,1,cmap) # xy
            if not os.path.isfile(file_name_xz): make_map(file_name_xz,pos,hsmlvals,mass,img_res,img_width,0,2,cmap) # xz
            if not os.path.isfile(file_name_yz): make_map(file_name_yz,pos,hsmlvals,mass,img_res,img_width,1,2,cmap) # yz

            map_xy = np.loadtxt(file_name_xy)
            render(ax1,map_xy,cmap,img_width)

            map_xz = np.loadtxt(file_name_xz)
            render(ax2,map_xz,cmap,img_width)

            map_yz = np.loadtxt(file_name_yz)
            render(ax3,map_yz,cmap,img_width)

            t1 = time.clock()
            print "[ rendering to figure: %3.2f minutes ]" % ((t1-t0)/60.)

            ax1.set_ylabel('(y,z,z) [Mpc]')

            if '14' not in hpath:
                ax1.set_xticklabels([])
                ax2.set_xticklabels([])
                ax3.set_xticklabels([])
            else:
                ax1.set_xlabel('x [Mpc]')
                ax2.set_xlabel('x [Mpc]')
                ax3.set_xlabel('y [Mpc]')

            ax2.set_yticklabels([])
            ax3.set_yticklabels([])
            
            t0 = time.clock()
            if include_rockstar: overlay_rockstar(ax1,ax2,ax3,hpath,zoomid,snapshot,render_type,img_width,nth_halo,sorted_by='vmax',want_label_on_halos=want_label_on_halos)
            if include_subfind: overlay_rockstar(ax1,ax2,ax3,hpath,zoomid,snapshot,render_type,img_width,nth_halo,sorted_by='vmax',want_label_on_halos=want_label_on_halos)
            t1 = time.clock()
            print "[ overlaying halos: %3.2f minutes ]" % ((t1-t0)/60.)

            ax1.text(0.9, 0.05,'XY', horizontalalignment='left',verticalalignment='center',transform = ax1.transAxes,weight='bold',fontsize=14,color='w')
            ax2.text(0.9, 0.05,'XZ', horizontalalignment='left',verticalalignment='center',transform = ax2.transAxes,weight='bold',fontsize=14,color='w')
            ax3.text(0.9, 0.05,'YZ', horizontalalignment='left',verticalalignment='center',transform = ax3.transAxes,weight='bold',fontsize=14,color='w')
            ax1.text(0.05, 0.9,'LX:'+str(lx), horizontalalignment='left',verticalalignment='center',transform = ax1.transAxes,weight='bold',fontsize=14,color='w')
    
    img_path = "/nfs/blank/h4231/bgriffen/Dropbox/work/projects/caterpillar/convergence/"
    img_filename = "H"+str(pid)+"_LX"+str(lx)+"_DIM"+str(img_res)+"_KPC"+str(int(img_width*1000.))+"_CM"+cmap_name.upper()+"_"+render_type.upper()+".png"
    print "writing out:",img_filename
    img_writeout = img_path + img_filename
    fig.savefig(img_writeout,bbox_inches='tight')
    #ax1.text(0.05, 0.95,'LX:'+str(lx), horizontalalignment='left',verticalalignment='center',transform = ax1.transAxes,weight='bold',fontsize=14,color='w')
    t1_total =  time.clock()

    print "[ TOTAL TIME FOR RENDER: %3.2f minutes ]" % ((t1_total-t0_total)/60.)

def create_fullbox_render(sim_path,snapshot=127,img_res=2048,lx=14,projection='all'):

    header = htils.get_halo_header(sim_path,snap=snapshot)

    pos = htils.load_partblock(sim_path,127,"POS ",parttype=1).astype('float32')
    mass = htils.load_partblock(sim_path,127,"MASS",parttype=1).astype('float32')*1e10/header.hubble
    ids = htils.load_partblock(sim_path,127,"ID  ",parttype=1)

    center = np.array([header.boxwidth/2.,header.boxwidth/2.,header.boxwidth/2.], dtype="float32")

    hsml = readhsml.hsml_file(sub_diri+"/outputs/",255)
    hsmlvals = hsml.Hsml[ids]

    #(ax1,file_name+"XY.dat",pos,mass,hsmlvals,nbins,box_width,0,1,mode,periodic_bool,center,want_quad=WANT_QUAD)
    map = SphMap.CalcDensProjection(pos, hsmlvals, mass, mass, img_res, img_res, header.boxwidth, header.boxwidth, img_width, center[0], center[1], center[2], dim0, dim1, 3, 0)

def get_subfind_particle_ids(subfind_path,groupid,snapnum):
    subfind_props = itils.loadSnapSubset(subfind_path,snapnum,groupid,0,1,["ParticleIDs"])
    return subfind_props["ParticleIDs"]

def is_match_by_ids(a,b,threshold=0.5):

    result_a = filter(set(a).__contains__, b)
    pc_a = len(result_a)/float(len(a))

    result_b = filter(set(b).__contains__, a)
    pc_b = len(result_b)/float(len(b))

    # display the number of subfind particles, number of rockstar particles and their respective match percentage

    print "SUB_N: %i ROCK_N: %i SUB_M: %3.3f ROCK_M: %3.3f" % (len(a),len(b),pc_a,pc_b)

    if pc_a >= threshold and pc_b >= threshold:
        return True
    else:
        return False
        
def transform_coordinates(x,y,z,vx,vy,vz,theta,phi):
    x_ = -z  * np.sin(theta)   + (x *  np.cos(phi) + y * np.sin(phi)) * np.cos(theta)
    y_ = -x  * np.sin(phi)     +  y  * np.cos(phi)
    z_ =  z  * np.cos(theta)   + (x *  np.cos(phi) + y * np.sin(phi)) * np.sin(theta)
    vx_ = -vz  * np.sin(theta) + (vx * np.cos(phi) + vy *np.sin(phi)) * np.cos(theta)
    vy_ = -vx  * np.sin(phi)   +  vy * np.cos(phi)
    vz_ =  vz  * np.cos(theta) + (vx * np.cos(phi) + vy *np.sin(phi)) * np.sin(theta)

    return x_,y_,z_,vx_,vy_,vz_
 
