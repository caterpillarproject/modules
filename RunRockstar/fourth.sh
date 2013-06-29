#!/bin/bash

node4=compute-0-9

ssh $node4 <<'ENDSSH'
cd /spacebase/data/AnnaGroup/modules/RunRockstar
bash RockstarHelper.sh
ENDSSH
