#!/bin/bash

#for FILE in mzero.cie m-00.cie m-05.cie m-10.cie m-15.cie m-20.cie m-30.cie m+05.cie
#do
#    wget "http://www.mso.anu.edu.au/~ralph/data/cool/"$FILE
#done
#
#for FILE in ff55 ff65 ff75 ff85 zf55 zf65 zf75 zf85
#do
#    wget "http://www.mso.anu.edu.au/~ralph/data/cool/pk6"$FILE".neq"
#done
#
#for FILE in 05 10 15 20 30
#do
#    wget "http://www.mso.anu.edu.au/~ralph/data/cool/pk6ff75m-"$FILE".neq"
#done

for FILE in 05 10 15 20 30
do
    wget "http://www.mso.anu.edu.au/~ralph/data/cool/pk6zf75m-"$FILE".neq"
done
