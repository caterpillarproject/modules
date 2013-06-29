## Load directories
while read line
do
    eval $line
done < directories

$rsdir/rockstar -c $outfile/auto-rockstar.cfg
