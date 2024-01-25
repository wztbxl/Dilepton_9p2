#!/bin/bash

nLoop=$1
repeat=$2
nExpr=$3
ltheta=$4
lphi=$5
wdir=Loop$nLoop.$repeat.$ltheta.$lphi
mkdir $wdir
cp testData.C $wdir
cd $wdir
root -b -l -q "testData.C+($nExpr, $ltheta, $lphi, 0, $nLoop, $repeat)"
cd ..
rm -rf $wdir
