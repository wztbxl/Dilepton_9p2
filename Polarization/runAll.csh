#!/bin/bash
date

dir="/star/u/liuzhen/lfsupc/WebPage/code/dimuon/Polarization/submit"
echo $dir

if [ -d $dir ]; then
    rm -rf $dir/*
fi

mkdir -p $dir/output
mkdir -p $dir/script
mkdir -p $dir/log

cp -v analysis run.csh $dir

cd $dir

ifile=0
for FILE in `cat /star/u/liuzhen/lfsupc/WebPage/code/dimuon/Polarization/datalist.lis`
do
     echo $FILE
     filename=$(basename $FILE)
     name=`echo $filename | sed 's/subDataList.//' | sed 's/.list//'`
     cp run.csh script/$name.csh

     echo "./analysis $FILE output/$name">>script/$name.csh
 
     job=submit_$name.job
     touch $job

     echo "Universe     = vanilla" >> $job
     echo "Notification = never" >> $job
     echo "Requirements = (CPU_Type != \"crs\") && (CPU_Experiment == \"star\")" >> $job
     echo "Initialdir   = $dir" >> $job
     echo "GetEnv       = True" >> $job
     echo "+Experiment  = \"star\"" >> $job
     echo "+Job_Type    = \"cas\"" >> $job
     echo  "     " >> $job

     echo "Executable       = script/$name.csh">>$job
     echo "Output           = log/$name.out">>$job
     echo "Error            = log/$name.err">>$job
     echo "Log              = log/$name.log">>$job
     echo  "Queue" >>$job
     echo  "     " >>$job
      
     condor_submit $job
     let "ifile+=1";
done

mv submit_*.job script
