#!/bin/bash

wdir=$PWD

# HX frame
parstheta=(-0.260 0.480 0.630 -0.047)
parsphi=(-0.034 0.008 -0.124 -0.300)

for((i=0; i<4; i++))
do
    ltheta=${parstheta[$i]}
    lphi=${parsphi[$i]}
    hadd Output/ToyMC_Data_Iter0_$ltheta\_$lphi\_0.0.root output/Iter0.Loop*.Data_$ltheta\_$lphi\_0.0.root
    rm -rf output/Iter0.Loop*.Data_$ltheta\_$lphi\_0.0.root
done

