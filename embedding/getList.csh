#!/bin/bash

ifile=100

echo "datalist" > data.list

for ((ifil0=100;ifile<155;ifile++))
do
    cd /star/u/wangzhen/run20/Dielectron/embedding/output/out_Electron_$ifile
    ls | sed "s:^:`pwd`/:" >> /star/u/wangzhen/run20/Dielectron/embedding/mc/data.list
    cd /star/u/wangzhen/run20/Dielectron/embedding/output/out_Positron_$ifile
    ls | sed "s:^:`pwd`/:" >> /star/u/wangzhen/run20/Dielectron/embedding/mc/data.list
done
