#!/bin/bash

root -b -q -l "ana_data0.C(3, 15, 10, -2.0)"
root -b -q -l "ana_data1.C(3, 15, 10, -2.0)"
root -b -q -l "ana_data2.C(3, 15, 10, -2.0)"
root -b -q -l "ana_data3.C(3, 15, 10, -2.0)"
root -b -q -l "ana_data4.C(3, 15, 10, -2.0)"
root -b -q -l "ana_data5.C(3, 15, 10, -2.0)"

hadd Data_3_15_10_-2.0/All_3_15_10_-2.0.root Data_3_15_10_-2.0/Rebin*.root

for((i=0; i<12; i++))
do
    mkdir  File$i
done

root -b -q -l "getSigFile.C(3, 15, 10, -2.0)"


