#!/bin/bash
date

#$1 is run number $2 is nJobs $3 is the error tag

echo submitting pre-crashed jobs for "$1" 
cp sched${1}.session.xml sched${1}.session.xml.old$3
tag=0
ifile=0
nbadjob=0
while [ "$ifile" -le $2 ];
do  
    grep Goodbye /star/u/wangzhen/run20/Dielectron/Efficiecny/log/$1_$ifile.out
    if [ $? -eq 0 ]; then
      echo "/star/u/wangzhen/run20/Dielectron/Efficiecny/output/$1_$ifile.root exit, skip ..."
    else
      echo "re-submit crashed job for miniTree/rootfiles_PicoDst/$1_$ifile.root"
      cp -f ./log/$1_$ifile.out ./log/$1_$ifile.out.old$3
      cp -f ./log/$1_$ifile.err ./log/$1_$ifile.err.old$3
      star-submit -r $ifile sched$1.session.xml
      let "nbadjob+=1"
      echo -n
     fi

    let "ifile+=1";
done

echo "total have $nbadjob need to be resubmit"
