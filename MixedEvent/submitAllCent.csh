#!/bin/bash

for i in {1..9}
do
    echo "Submitting job for cent$i with i=$i"
    ./runAll.csh all cent${i}_VzBin10_EPBin12_Buff350_makeRealPair $i
done
