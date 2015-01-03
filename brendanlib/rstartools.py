import subprocess as sub
import glob,sys,os
import readsnapshots.readsnapHDF5 as rhdf5
import mergertrees.MTCatalogue as MT
import numpy as np
import brendanlib.grifflib as glib

def get_subhalo_mass_fraction(halodata,haloid):
    mhost = float(halodata.data[halodata.data['id'] == haloid]['mvir'])
    subhalos = halodata.get_all_subhalos_within_halo(haloid)
    if len(subhalos) > 0:
        submass = np.array(subhalos['mvir']).sum()
        fsub = submass/mhost
    else:
        fsub = 0.0

    return fsub
    
def get_min_distance_to_contam(halopath,halox,haloy,haloz,part_type=2):
    pos = rhdf5.read_block(halopath + "/outputs/snapdir_255/snap_255","POS ",parttype=part_type)
    head = rhdf5.snapshot_header(halopath + "/outputs/snapdir_255/snap_255.0.hdf5")
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

        current_jobs,jobids,jobstatus = glib.getcurrentjobs()
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
