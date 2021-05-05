#!/bin/sh
# runtests.sh --
#    Bourne shell script to run unit tests.
#if test -f runtests.log ; then
#    rm runtests.log
#fi

echo ALL >ftnunit.run

chk=1

#until test ! -f ftnunit.lst -a $chk -eq 0 ; do
#until [ ! -f ftnunit.lst ] && [ $chk -eq 0 ] ; do
while [ -f ftnunit.lst ] || [ $chk -eq 1 ] ; do
    chk=0
#    $1 $2 $3 $4 $5 $6 $7 $8 $9 >>runtests.log 2>&1

# VO VO:
    ./tomofast3D | tee runtests.log
done

#rm ftnunit.run


