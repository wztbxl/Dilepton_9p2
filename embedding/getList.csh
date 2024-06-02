#!/bin/bash

ifile=100

echo "datalist" > data.list

for ((ifil0=100;ifile<119;ifile++))
do
    echo "open /star/u/wangzhen/run20/Dielectron/embedding/output/out_Electron_$ifile"
    cd /star/u/wangzhen/run20/Dielectron/embedding/output/out_Electron_$ifile
    rm ./data.list
    ls | sed "s:^:`pwd`/:" >> /star/u/wangzhen/run20/Dielectron/embedding/mc/data.list
    cd /star/u/wangzhen/run20/Dielectron/embedding/output/out_Positron_$ifile
    rm ./data.list
    ls | sed "s:^:`pwd`/:" >> /star/u/wangzhen/run20/Dielectron/embedding/mc/data.list
done
