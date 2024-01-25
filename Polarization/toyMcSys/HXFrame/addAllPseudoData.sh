#!/bin/bash

wdir=$PWD

# HX frame
parstheta=(-0.021 -0.331 -0.169 0.117)
parsphi=(-0.059 0.125 0.254 -0.012)

for((i=0; i<4; i++))
do
    ltheta=${parstheta[$i]}
    lphi=${parsphi[$i]}
    hadd Output/ToyMC_Data_Iter0_$ltheta\_$lphi\_0.0.root output/Iter0.Loop*.Data_$ltheta\_$lphi\_0.0.root
    rm -rf output/Iter0.Loop*.Data_$ltheta\_$lphi\_0.0.root
done

