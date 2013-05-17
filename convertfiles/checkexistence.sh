#!/bin/bash
basepath="/spacebase/data/AnnaGroup/caterpillarparent/planckparent/outputs/"
echo "READING FROM: "$basepath
for i in {0..63}
do
	if [ $i -lt '10' ]; 
	then
	    file='/snap_00'$i
        outfile=$basepath'snapdir_00'$i$file
    fi
    if [ $i -gt '9' ]; 
    then
        file='/snap_0'$i
        outfile=$basepath'snapdir_0'$i$file
    fi
    for j in {0..3}
    do
    tmpfile=$outfile'.'$j
        if [ -a $tmpfile ];
        then
        echo "Gadget file 'snapdir_0$i$file.$j' exists"
        else
        echo "Gadget file 'snapdir_0$i$file.$j' does not exist"
        fi
    done
done