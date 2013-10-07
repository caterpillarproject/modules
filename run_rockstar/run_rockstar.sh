#!/bin/bash
if [ $# != 1 ]; then
    echo "ERROR: takes exactly one argument"
    echo "Usage: run_rockstar.sh my/rrs/file.rrs"
    exit 0
fi
echo "Using RRS file: $1"

## Read in rrs file parameters
while read line; do
    eval $line
done < $1

echo "Going into ${rrsdir}"
cd ${rrsdir}

##############################################################
## If you ran this script and rockstar completed but the 
## "finishing" step messed up (making folder structure/merger tree stuff),
## uncomment this and the SKIPTOFINISH below to ad-hoc comment out 
## the part of the script that runs rockstar.
#<<SKIPTOFINISH
##############################################################

##############################################################
## Create directory where temporary scripts will be stored
if [ ! -e rrs_scripts ]; then
    mkdir rrs_scripts
fi
## Create directory with the auto-generated rockstar .cfg file
if [ ! -e rrs_cfg ]; then
    mkdir rrs_cfg
fi
## Create the cfg file
nnodes=${#nodearr[@]}
let numwriters=nnodes*8
cat > rrs_cfg/${jobname}.cfg <<EOF
PARALLEL_IO = 1
INBASE = ${indir}
OUTBASE  = ${outdir}
FILENAME = snapdir_<snap>/snap_<snap>.<block>
NUM_BLOCKS = ${numblocks}
NUM_SNAPS = ${numsnaps}
NUM_WRITERS = ${numwriters}
FORK_READERS_FROM_WRITERS   = 1
FORK_PROCESSORS_PER_MACHINE = 8
FILE_FORMAT = "GADGET2"
FORCE_RES = ${forceres}
FULL_PARTICLE_CHUNKS = 1         #Print particles for 1 of the 8 tasks
STARTING_SNAP = ${startsnap}
MASS_DEFINITION = "200c"
BOX_SIZE = 100
EOF

config=${rrsdir}/rrs_cfg/${jobname}.cfg
echo "Created cfg file $config"

##############################################################
##### Write temporary scripts
##### Will overwrite any previously existing files!
##### Use the jobname variable in the config file to avoid this if necessary
## main.sh
cat > rrs_scripts/${jobname}_main.sh <<EOF
#!/bin/bash
ssh compute-0-${nodearr[0]} << 'ENDSSH'
echo "I AM THE MAIN ONE! DO YOU HEAR ME?"
$rsdir/rockstar -c $config
echo "I AM THE MAIN ONE! I ROAR AGAIN!"
ENDSSH
EOF
chmod +x rrs_scripts/${jobname}_main.sh

## worker_xx.sh
for node in "${nodearr[@]}"; do
    cat > rrs_scripts/${jobname}_worker_${node}.sh << EOF
#!/bin/bash
ssh compute-0-$node << 'ENDSSH'
echo "I AM COMPUTE NODE $node! HEAR ME ROAR!"
$rsdir/rockstar -c $outdir/auto-rockstar.cfg
echo "I AM COMPUTE NODE $node! I ROAR AGAIN!"
ENDSSH
EOF
    chmod +x rrs_scripts/${jobname}_worker_${node}.sh
done

##### Run Rockstar!
echo "Starting to run Rockstar..."
./rrs_scripts/${jobname}_main.sh &
sleep 10
let lastnode=$nnodes-1
for node in "${nodearr[@]}"; do
    if [ $node == ${nodearr[lastnode]} ]; then
	./rrs_scripts/${jobname}_worker_${node}.sh
	#echo "compute-0-$node"
    else
	#echo "compute-0-$node &"
	./rrs_scripts/${jobname}_worker_${node}.sh &
    fi
done
##############################################################
#SKIPTOFINISH
##############################################################

##############################################################
##### Finish Rockstar
##############################################################
echo "Making Merger Tree and generating parents.list"
perl $rsdir/scripts/gen_merger_cfg.pl $config
cd $ctrees
#make
perl do_merger_tree.pl $outdir/outputs/merger_tree.cfg

cd $outdir
## make folders
c=0
let lastsnap=$numsnaps-1
while [ $c -le $lastsnap ]
do
	echo "Making halos_0$c.."
	mkdir "halos_$c"
	mv halos_$c.* halos_$c
	mv out_$c.* halos_$c
	(( c++ ))
done


## generate parents.list
for ((i=0;i<=$lastsnap;++i))
do
    $rsdir/util/find_parents $outdir'/halos_'$i'/out_'$i'.list' 100 > $outdir'/halos_'$i'/parents.list'
done
