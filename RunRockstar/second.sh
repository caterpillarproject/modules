#!/bin/bash

node2=compute-0-7

ssh $node2 <<'ENDSSH'
cd /spacebase/data/AnnaGroup/modules/RunRockstar
bash RockstarHelper.sh
ENDSSH
