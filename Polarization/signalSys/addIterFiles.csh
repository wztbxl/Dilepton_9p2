#!/bin/bash
date

pwd=$PWD

for((f=0; f<12; f++))
do
fname=File$f
fiter=Iter_f$f
dir="$pwd/$fname/$fiter"
cd $dir
rm -rf *All.root
rm -rf final\_$f
hadd final\_$f\.root CSPt*Iter*.root HXPt*Iter*.root 
done

