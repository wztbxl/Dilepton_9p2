#!/bin/bash

rm -rf CSFinal_iter*.pdf
rm -rf getMean_CS_iter*.pdf
rm -rf CS_real.root

for i in 10 1 2 3 4 5 6 7 8
do
    root -b -l -q "drawIter.C($i)"
done
