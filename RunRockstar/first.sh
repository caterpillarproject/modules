#!/bin/bash

node1=compute-0-6

ssh $node1 <<'ENDSSH'
cd /spacebase/data/AnnaGroup/modules/RunRockstar
bash RockstarHelper.sh
ENDSSH
