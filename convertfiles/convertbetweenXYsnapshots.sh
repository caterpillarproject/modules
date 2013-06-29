# run bash  bash convertbetweenXYsnapshots.sh < 0 10 
# to convert snapshots 0 to 10

file='/spacebase/data/AnnaGroup/caterpillarparent/planckparent/outputs/'

read $sni
read $snf

echo 'Converting snapshots between '$sni' and '$snf
echo 'Reading from:'$file

for i in {sni..snf}
do
    if [ $i -lt '10' ]; 
    then
	    outfile=$file'snapdir_00'$i
    fi
    
    if [ $i -gt '9' ]; 
    then
        outfile=$file'snapdir_0'$i
    fi
    
    python HDF5toGadget2.py $file $i $outfile
    
    echo 'Created: '$outfile

done
