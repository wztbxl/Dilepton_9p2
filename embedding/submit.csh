#! /bin/tcsh

if($#argv != 3) then
   echo "Please infile TWO arguments!"
   echo "argv0 -- the particle species. For me, it should be 'Electron' or 'Positron' "
   echo "argv1 -- the production number ID. For me, it should be '100', or '200' "
   echo "argv2 -- the request ID "
   exit 0
endif

set myPath = /star/u/wangzhen/run20/Dielectron/embedding
#set myPath = /star/u/wangzhen/QA/wangzhen/embedding/myEmbedding/basicQA
#set infileDir = /star/embed/embedding/AuAu54_production_2017/$1_$2_20183001/P18ic.SL18c/2017
#set infileDir = /star/embed/embedding/AuAu54_production_2017/$1_$2_20183001/P18ic.SL18c_embed/2017
#set infileDir = /star/embed/embedding/27GeV_production_2018/$1_$2_20193201/P19ib.SL19b/2018
#set infileDir = /star/embed/embedding/AuAu54_production_2017/$1_$2_20191701/P18ic.SL18c_embed/2017
#set infileDir = /star/embed/embedding/AuAu54_production_2017/$1_$2_20183001/P18ic.SL18c_embed/2017
#set infileDir = /star/embed/embedding/production_isobar_2018/$1_$2_20214217/P20ic.SL20c/2018
# set infileDir = /star/embed/embedding/production_isobar_2018/$1_$2_$3/P20ic.SL20c/2018
set infileDir = /star/data105/embedding/production_7p7GeV_2021//$1_$2_$3/P22ib.SL22b/2021
#set infileDir = /star/data18/embedding/AuAu54_production_2017/$1_$2_20183001/P18ic.SL18c_embed/2017
#/star/embed/embedding/AuAu54_production_2017/Electron_100_20183001/
set outfileDir = ./output/out_$1_$2
set logDir = ./log/log_$1_$2
set scriptDir = script/script_$1_$2

if(! -d $myPath/$outfileDir) then
      mkdir $myPath/$outfileDir
endif

if(! -d $myPath/$logDir) then
      mkdir $myPath/$logDir
endif

if(! -d $myPath/$scriptDir) then
      mkdir $myPath/$scriptDir
endif

#if(! -d $outfileDir) then
#      ln -s $myPath/$outfileDir ./
#endif

#if(! -d $logDir) then
#      ln -s $myPath/$logDir ./
#endif

#if(! -d $scriptDir) then
#      ln -s $myPath/$scriptDir ./
#endif

rm -rf $scriptDir/*
rm -rf $logDir/*
rm -rf $outfileDir/*.root
rm -rf job/*

#cp run.con job/runAll$1_$2.job
@ nfile=0



foreach file (`find $infileDir/*/* -name 'st_physics_adc*.MuDst.root'`)
#foreach file (`ls $infileDir/*/st_physics_adc*.event.root`)
  echo " *** $nfile ***"

 cp run.con job/runAll$1_$2_$nfile.job

  set baseName = `basename $file`
  set log = `basename $file`
  echo $baseName

  cp run.csh run_tmp.csh
  echo "root4star -b <<EOF">>run_tmp.csh
  echo ".O2">>run_tmp.csh
  echo -n '.x doEvent.C(1e9,"'>>run_tmp.csh
  echo -n $file>>run_tmp.csh
  echo -n '","'>>run_tmp.csh
  echo -n $outfileDir>>run_tmp.csh
  echo '")'>>run_tmp.csh
  echo ".q">>run_tmp.csh
  echo "EOF">>run_tmp.csh
  mv run_tmp.csh script/script_$1_$2/$baseName.csh

  echo "Executable       = script/script_$1_$2/$baseName.csh">>job/runAll$1_$2_$nfile.job
  echo "Output           = $logDir/$baseName.out">>job/runAll$1_$2_$nfile.job
  echo "Error            = $logDir/$baseName.err">>job/runAll$1_$2_$nfile.job
  echo "Log              = $logDir/$baseName.log">>job/runAll$1_$2_$nfile.job
  echo  "Queue" >>job/runAll$1_$2_$nfile.job
  echo  "     " >>job/runAll$1_$2_$nfile.job

  condor_submit job/runAll$1_$2_$nfile.job

  @ nfile++

end

  #condor_submit job/runAll$1_$2.job

