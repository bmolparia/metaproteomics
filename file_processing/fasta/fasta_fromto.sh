#!/bin/bash

a=$1
b=$2
c=$((b+a))
d=$3
#echo $a
#echo $b
#echo $c
awk -v f=$a -v t=$b -v RS='>' -v ORS='>' 'BEGIN{printf ">"} NR>t{exit} NR>f && NR<=t' $d | head -n-1
