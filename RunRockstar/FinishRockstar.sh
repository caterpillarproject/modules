#!/bin/bash

echo "Making Merger Tree and generating parents.list"

## Load directories
while read line
do
    eval $line
done < directories

perl $rsdir/scripts/gen_merger_cfg.pl $rsdir$config
cd $ctrees
make
perl do_merger_tree.pl $outfile/outputs/merger_tree.cfg

cd $outfile
## make folders
c=0
while [ $c -le 63 ]
do
	echo "Making halos_0$c.."
	mkdir "halos_$c"
	mv halos_$c.* halos_$c
	mv out_$c.* halos_$c
	(( c++ ))
done


## generate parents.list
for i in {0..63}
do
    $rsdir/util/find_parents $outfile'/halos_'$i'/out_'$i'.list' > $outfile'/halos_'$i'/parents.list'
done
