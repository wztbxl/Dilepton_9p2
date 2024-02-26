#!/bin/bash
date

if [ $# -ne 1 ]; then
     echo "Please input one arguement!"
	 exit 1
fi

dir="/star/u/wangzhen/QA/wangzhen/embedding/myEmbedding/basicQA/data"
echo $dir

if [ ! -d $dir/output_$1 ]; then
     mkdir -p $dir/output_$1
fi

if [ ! -d $dir/script_$1 ]; then
     mkdir -p $dir/script_$1
fi 

if [ ! -d  $dir/log_$1 ]; then
     mkdir -p $dir/log_$1
fi

if [ ! -d output_$1 ]; then
     ln -s $dir/output_$1 ./
fi

if [ ! -d script_$1 ]; then
     ln -s $dir/script_$1 ./
fi

if [ ! -d log_$1 ]; then
     ln -s $dir/log_$1 ./
fi

rm -rf job/*
  
ifile=0
for FILE in `cat /star/u/wangzhen/QA/wangzhen/embedding/myEmbedding/basicQA/data/datalist_$1`
do
     echo $FILE
     cp run.con job/runAll_$1_$ifile.job
     cp ./run.csh script_$1/$1_$ifile.csh
 
     echo "./produceQA $FILE output_$1/$ifile">>script_$1/$1_$ifile.csh

     echo "Executable       = script_$1/$1_$ifile.csh">>job/runAll_$1_$ifile.job
     echo "Output           = log_$1/$1_$ifile.out">>job/runAll_$1_$ifile.job
     echo "Error            = log_$1/$1_$ifile.err">>job/runAll_$1_$ifile.job
     echo "Log              = log_$1/$1_$ifile.olog">>job/runAll_$1_$ifile.job
     echo  "Queue" >>job/runAll_$1_$ifile.job
     echo  "     " >>job/runAll_$1_$ifile.job
     condor_submit job/runAll_$1_$ifile.job

     let "ifile+=1";
done


