#!/bin/bash
date

dca=(3 1)
nhitsfit=(15 25)
nhitsdedx=(10 15)
nsigmapi=(-1.5 -2.0 -2.5)

pwd=$PWD

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

                dir="$pwd/$fcuts/$fIter"
                cd $dir
                rm -rf *All.root
                rm -rf final\_$indca\_$inhitsfit\_$inhitsdedx\_$insigmapi
                hadd final\_$indca\_$inhitsfit\_$inhitsdedx\_$insigmapi\.root CSPt*Iter*.root HXPt*Iter*.root 
            done
        done
    done
done
