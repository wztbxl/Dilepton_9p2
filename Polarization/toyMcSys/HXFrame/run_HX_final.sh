#!/bin/bash
rm -rf getMean_HX_iter*.pdf
rm -rf HXFinal_iter*.pdf
rm -rf HX_real.root

for i in  1 2 3 4 5 6 7 8 
do
    root -b -l -q "drawIter.C($i)"
done
