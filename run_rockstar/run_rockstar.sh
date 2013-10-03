#!/bin/bash
if [ $# != 1 ]; then
    echo "ERROR: takes exactly one argument"
    echo "Usage: run_rockstar.sh my/run/rockstar/file"
    exit 0
fi
echo "Using RRS file: $1"

## Read in rrs file parameters
while read line; do
    eval $line
done < $1

#<<COMMENT1
## Create directory where temporary scripts will be stored
if [ ! -e tmp_scripts ]; then
    mkdir tmp_scripts
fi

##### Write temporary scripts
##### Will overwrite any previously existing files!
##### Use the jobname variable in the config file to avoid this if necessary
## main.sh
cat > tmp_scripts/${jobname}_main.sh <<EOF
#!/bin/bash
ssh compute-0-${nodearr[0]} << 'ENDSSH'
$rsdir/rockstar -c $config
ENDSSH
EOF
chmod +x tmp_scripts/${jobname}_main.sh

## worker_xx.sh
for node in "${nodearr[@]}"; do
    cat > tmp_scripts/${jobname}_worker_${node}.sh << EOF
#!/bin/bash
ssh compute-0-$node << 'ENDSSH'
$rsdir/rockstar -c $outdir/auto-rockstar.cfg
ENDSSH
EOF
    chmod +x tmp_scripts/${jobname}_worker_${node}.sh
done

##### Run Rockstar!
./tmp_scripts/${jobname}_main.sh &
sleep 10
nnodes=${#nodearr[@]}
let lastnode=$nnodes-1
for node in "${nodearr[@]}"; do
    if [ $node == ${nodearr[lastnode]} ]; then
	./tmp_scripts/${jobname}_worker_${node}.sh
	#echo "compute-0-$node"
    else
	#echo "compute-0-$node &"
	./tmp_scripts/${jobname}_worker_${node}.sh &
    fi
done
#COMMENT1

##### Finish Rockstar
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
    $rsdir/util/find_parents $outdir'/halos_'$i'/out_'$i'.list' > $outdir'/halos_'$i'/parents.list'
done
