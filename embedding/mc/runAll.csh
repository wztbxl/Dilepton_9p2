#!/bin/bash
date

if [ $# -ne 2 ]; then
     echo "Please input two arguement!"
	 exit 1
fi

dir="/star/u/wangzhen/run20/Dielectron/embedding/mc"
echo $dir

if [ ! -d $dir/output/output_$1/$2 ]; then
     mkdir -p $dir/output_$1/$2
fi

if [ ! -d $dir/script/script_$1/$2 ]; then
     mkdir -p $dir/script_$1/$2
fi 

if [ ! -d  $dir/log/log_$1/$2 ]; then
     mkdir -p $dir/log_$1/$2
fi

if [ ! -d output/output_$1/$2 ]; then
     ln -s $dir/output_$1/$2 ./
fi

if [ ! -d script/script_$1/$2 ]; then
     ln -s $dir/script_$1/$2 ./
fi

if [ ! -d log/log_$1/$2 ]; then
     ln -s $dir/log_$1/$2 ./
fi

rm -rf job/*
  
ifile=0
for FILE in `cat /star/u/wangzhen/QA/wangzhen/embedding/myEmbedding/basicQA/mc/datalist_$1`
do
     echo $FILE
     cp run.con job/runAll_$1_$ifile.job
     cp ./run.csh script_$1/$1_$ifile.csh
 
     echo "./produceQA $FILE output_$1/$2/$ifile">>script_$1/$1_$ifile.csh

     echo "Executable       = script_$1/$1_$ifile.csh">>job/runAll_$1_$ifile.job
     echo "Output           = log_$1/$1_$ifile.out">>job/runAll_$1_$ifile.job
     echo "Error            = log_$1/$1_$ifile.err">>job/runAll_$1_$ifile.job
     echo "Log              = log_$1/$1_$ifile.olog">>job/runAll_$1_$ifile.job
     echo  "Queue" >>job/runAll_$1_$ifile.job
     echo  "     " >>job/runAll_$1_$ifile.job
     condor_submit job/runAll_$1_$ifile.job

     let "ifile+=1";
done


