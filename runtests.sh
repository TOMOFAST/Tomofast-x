#!/bin/sh
# Bourne shell script to run the unit tests.

echo ALL >ftnunit.run

chk=1

while [ -f ftnunit.lst ] || [ $chk -eq 1 ] ; do
    chk=0
    ./tomofastx | tee runtests.log
done


