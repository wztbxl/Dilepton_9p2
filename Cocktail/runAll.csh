#!/bin/bash
date
#$1 = meson code $2 = decay model $3 = meson name
 dir="/star/u/wangzhen/QA/wangzhen/Cocktail"
 echo $dir
 echo $1
 echo $2
 echo $3

 if [ ! -d $dir/script ]; then
      mkdir -p $dir/script
 fi
 
 if [ ! -d  $dir/log ]; then
      mkdir -p $dir/log
 fi

 if [ ! -d $dir/jobs ]; then
     mkdir -p $dir/jobs
 fi
 
 if [ ! -d script ]; then
      ln -s $dir/script ./
 fi
 
 if [ ! -d log ]; then
      ln -s $dir/log ./
 fi


 ifile=0
 for ((ifile=0;ifile<300;ifile++)) 
 do
     cp ./run.csh script/runT$3_$2_$ifile.csh
     chmod +x script/runT$3_$2_$ifile.csh
     cp run.con jobs/runAll$3$ifile.job 
     echo "./mesondecay 080 $1 $2 $3 $ifile" >>script/runT$3_$2_$ifile.csh  
#     echo "./mesondecay 080 221 3 eta $ifile" >>script/runT$1_$ifile.csh
#     echo "./mesondecay 080 223 3 omega $ifile" >>script/runT$1_$ifile.csh
#submit less particle to reduce time
#     echo "./mesondecay 080 331 3 etaprim $ifile" >>script/runT$1_$ifile.csh
#     echo "./mesondecay 080 333 3 phi $ifile" >>script/runT$1_$ifile.csh
#     echo "./mesondecay 080 221 2 eta $ifile" >>script/runT$1_$ifile.csh
#     echo "./mesondecay 080 223 2 omega $ifile" >>script/runT$1_$ifile.csh
#     echo "./mesondecay 080 333 2 phi $ifile" >>script/runT$1_$ifile.csh
#     echo "./mesondecay 080 443 2 jpsi $ifile" >>script/runT$1_$ifile.csh
#     echo "./mesondecay 080 100443 2 psi $ifile" >>script/runT_$ifile.csh
#     echo "./mesondecay 113 2 rho $ifile" >>script/runT_$ifile.csh
#     echo "./mesondecay 1 2 photon $ifile" >>script/runT_$ifile.csh
 
     echo "Executable       = script/runT$3_$2_$ifile.csh">>jobs/runAll$3$ifile.job
     echo "Output           = log/runT$3_$2_$ifile.out">>jobs/runAll$3$ifile.job
     echo "Error            = log/runT$3_$2_$ifile.err">>jobs/runAll$3$ifile.job
     echo "Log              = log/runT$3_$2_$ifile.olog">>jobs/runAll$3$ifile.job
     echo  "Queue" >>jobs/runAll$3$ifile.job
     echo  "     " >>jobs/runAll$3$ifile.job
     echo  "ifile== $ifile" 
     condor_submit jobs/runAll$3$ifile.job

 done

