#!/bin/bash
date

if [ $# -ne 2 ]; then
     echo "Please input two arguement!"
     echo "arg1 for the datalist tag"
     echo "arg2 for the production tag"
	 exit 1
fi

dir=$(pwd)
echo $dir

if [ ! -d $dir/output_$1/$2 ]; then
     mkdir -p $dir/output_$1/$2
fi

if [ ! -d $dir/script_$1/$2 ]; then
     mkdir -p $dir/script_$1/$2
fi 

if [ ! -d  $dir/log_$1/$2 ]; then
     mkdir -p $dir/log_$1/$2
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

  
ifile=0
for FILE in `cat datalist_$1`
do
     echo $FILE
     cp ./run.csh script_$1/$2/$1_${ifile}.csh
    cp run.con jobs/runAll_${ifile}_${1}.job
 
     echo "./analysis $FILE output_$1/$2/${ifile} $3">>script_$1/$2/$1_${ifile}.csh

     echo "Executable       = script_$1/$2/$1_${ifile}.csh">>jobs/runAll_${ifile}_${1}.job
     echo "Output           = log_$1/$2/$1_${ifile}.out">>jobs/runAll_${ifile}_${1}.job
     echo "Error            = log_$1/$2/$1_$ifile.err">>jobs/runAll_${ifile}_${1}.job
     echo "Log              = log_$1/$2/$1_${ifile}.olog">>jobs/runAll_${ifile}_${1}.job
     echo  "Queue" >>jobs/runAll_${ifile}_${1}.job
     echo  "     " >>jobs/runAll_${ifile}_${1}.job
    
     echo "script_$1/$2/$1_${ifile}.csh"
    echo "jobs/runAll_${ifile}_${1}.job"  
    condor_submit jobs/runAll_${ifile}_${1}.job
     let "ifile+=1";
done
