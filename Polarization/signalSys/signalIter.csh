#!/bin/bash
date

for((f=0; f<12; f++))
do

fname=File$f
fiter=Iter_f$f
dir="/star/u/liuzhen/Run15/Run15pp200/submitjob/Polarization/signalSys/$fname/$fiter"
echo $dir
mkdir $dir

if [ -d $dir ]; then
    rm -rf $dir/*
    fi

    mkdir -p $dir/output
    mkdir -p $dir/script
    mkdir -p $dir/log

    cp -v /star/u/liuzhen/Run15/Run15pp200/submitjob/Polarization/signalSys/iter_code/analysis /star/u/liuzhen/Run15/Run15pp200/submitjob/Polarization/signalSys/iter_code/run.csh /star/u/liuzhen/Run15/Run15pp200/submitjob/Polarization/signalSys/iter_code/allFile.list /star/u/liuzhen/Run15/Run15pp200/submitjob/Polarization/signalSys/iter_code/ana_* $dir

    cd $dir
    cp ../systematic_Signal.root .
    echo "./analysis allFile.list 1 3 15 10 -2.0 0 0">> run.csh
    echo 'root -b -q -l "ana_ite.C(1)"' >> run.csh
    echo "for((i=2; i<=8; i++))" >> run.csh
    echo "do" >> run.csh
    echo "for((f=0; f<2; f++))" >> run.csh
    echo "do" >> run.csh
    echo " for((p=0; p<4; p++))" >> run.csh
    echo "do" >> run.csh
    echo "./analysis allFile.list \$i 3 15 10 -2.0 \$f \$p">> run.csh
    echo "done" >> run.csh
    echo "done" >> run.csh
    echo 'root -b -q -l "ana_ite.C($i)"' >> run.csh
    echo "done" >> run.csh

    job=submit_$f.job
    touch $job

    echo "Universe     = vanilla" >> $job
    echo "Notification = never" >> $job
    echo "Requirements = (CPU_Type != \"crs\") && (CPU_Experiment == \"star\")" >> $job
    echo "Initialdir   = $dir" >> $job
    echo "GetEnv       = True" >> $job
    echo "+Experiment  = \"star\"" >> $job
    echo "+Job_Type    = \"cas\"" >> $job
    echo  "     " >> $job

    echo "Executable       = run.csh">>$job
    echo "Output           = log/$f.out">>$job
    echo "Error            = log/$f.err">>$job
    echo "Log              = log/$f.log">>$job
    echo  "Queue" >>$job
    echo  "     " >>$job

    condor_submit $job
    mv submit_*.job script
    cd /star/u/liuzhen/Run15/Run15pp200/submitjob/Polarization/signalSys
    done
    
