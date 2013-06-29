#!/bin/bash

node3=compute-0-8

ssh $node3 <<'ENDSSH'
cd /spacebase/data/AnnaGroup/modules/RunRockstar
bash RockstarHelper.sh
ENDSSH
