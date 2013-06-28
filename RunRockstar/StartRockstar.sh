while read line
do
    eval $line
done < directories
$rsdir/rockstar -c $rsdir$config