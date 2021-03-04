#!/bin/bash

DIRECTORY=$1
PREFIX=$2

WDIR="/lustre/scratch/npaulat/ESAT"
cd $WDIR

IDS=()
while read ID; do IDS+=($ID); done < esat_batch_${PREFIX}_ids.txt

cd $DIRECTORY

for i in {0..299}; do cd ${IDS[$i]}; bash ../../ESAT_flag_stats.sh ${IDS[$i]}; cd ..; done
