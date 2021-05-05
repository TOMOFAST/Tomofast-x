#!/bin/sh -f

#FILENAME1=start_grav.txt
FILENAME1=true_model_grav.txt
FILENAME2=${FILENAME1}_formatted

awk '{\
if (NR > 1) \
    { print $1 " " $2 " " $3 " " $4 " " $5 " " $6 " " $7 " " 1 " " 1 " " NR-1 }  \
else \ 
    { print $0 }}' $FILENAME1 > $FILENAME2
