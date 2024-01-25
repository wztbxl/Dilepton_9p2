#!/bin/bash
date

dca=(3 1)
nhitsfit=(15 25)
nhitsdedx=(10 15)
nsigmapi=(-1.5 -2.0 -2.5)

for((a=0; a<2; a++))
do
for((f=0; f<2; f++))
do
for((d=0; d<2; d++))
do
for((p=0; p<3; p++))
do

indca=${dca[$a]}
inhitsfit=${nhitsfit[$f]}
inhitsdedx=${nhitsdedx[$d]}
insigmapi=${nsigmapi[$p]}
fcuts=Data\_$indca\_$inhitsfit\_$inhitsdedx\_$insigmapi

mkdir $fcuts
cd $fcuts
dir=$PWD
echo $dir

if [ -d $dir ]; then
    rm -rf $dir/*
fi

mkdir -p $dir/output
mkdir -p $dir/script
mkdir -p $dir/log

cp -v ../../analysis ../../run.csh $dir

ifile=0
for FILE in `ls ../../DataList/*.list`
do
     echo $FILE
     filename=$(basename $FILE)
     name=`echo $filename | sed 's/subDataList.//' | sed 's/.list//'`
     cp run.csh script/$name.csh

     echo "./analysis $FILE output/$name $indca $inhitsfit $inhitsdedx $insigmapi">>script/$name.csh
 
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
cd ../
done
done
done
done

