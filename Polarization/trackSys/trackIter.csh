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
fIter=Iter\_$indca\_$inhitsfit\_$inhitsdedx\_$insigmapi
mkdir $fcuts/$fIter

dir=$PWD/$fcuts/$fIter
echo $dir

if [ -d $dir ]; then
    rm -rf $dir/*
    fi

    mkdir -p $dir/output
    mkdir -p $dir/script
    mkdir -p $dir/log

    cp -v $PWD/iter_code/analysis $PWD/iter_code/run.csh $PWD/iter_code/allFile.list $PWD/iter_code/ana_* $dir

    cd $dir
    cp ../Result*.root systematic_Signal.root
    echo "./analysis allFile.list 1 $indca $inhitsfit $inhitsdedx $insigmapi 0 0">> run.csh
    echo 'root -b -q -l "ana_ite.C(1)"' >> run.csh
    echo "for((i=2; i<=10; i++))" >> run.csh
    echo "do" >> run.csh
    echo "for((f=0; f<2; f++))" >> run.csh
    echo "do" >> run.csh
    echo " for((p=0; p<4; p++))" >> run.csh
    echo "do" >> run.csh
    echo "./analysis allFile.list \$i $indca $inhitsfit $inhitsdedx $insigmapi \$f \$p">> run.csh
    echo "done" >> run.csh
    echo "done" >> run.csh
    echo 'root -b -q -l "ana_ite.C($i)"' >> run.csh
    echo "done" >> run.csh

    job=submit_$fIter.job
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
    echo "Output           = log/$fIter.out">>$job
    echo "Error            = log/$fIter.err">>$job
    echo "Log              = log/$fIter.log">>$job
    echo  "Queue" >>$job
    echo  "     " >>$job

    condor_submit $job
    mv submit_*.job script
    cd ../../
    done
    done
    done
    done
