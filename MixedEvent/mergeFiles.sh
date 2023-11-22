#!/bin/bash

if [ $# -ne 1 ]; then
     echo "Please input one arguement!"
     echo "arg1 for the file tag"
   exit 1
fi

#$1 should be the tag of output dir path
tag="$1"

cd /star/u/wangzhen/Hadd/hadd
echo $PWD

# loop centrality from 1-9
for i in {1..9}
do
    #build the path
    target_dir="/star/u/wangzhen/run20/Dielectron/MixedEvent/output_all/cent$i"

    #add the tag to the dir name
    if [ -n "$tag" ]; then
        target_dir="$target_dir"_"$tag"
    fi

    echo "Running command for cent$i:" $target_dir
    ./runhadd.csh "$target_dir" 50
done
