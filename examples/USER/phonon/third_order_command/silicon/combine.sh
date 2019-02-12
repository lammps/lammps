#!/bin/bash

#This script takes one argument
#The argument is the base name for the split up tensor
#The script then combines and sorts the tensor
#$1 file name

echo "$1"
[ -f $1 ] && rm $1

for i in $(ls ./$1*); do
    cat $i >> temp
    rm $i
done

sort temp | sort -s -n -k 3 | sort -s -n -k 1 > $1
rm temp
