#!/bin/bash

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

root -b -q -l "ana_data0.C($indca, $inhitsfit, $inhitsdedx, $insigmapi)"
root -b -q -l "ana_data1.C($indca, $inhitsfit, $inhitsdedx, $insigmapi)"
root -b -q -l "ana_data2.C($indca, $inhitsfit, $inhitsdedx, $insigmapi)"
cd Data\_$indca\_$inhitsfit\_$inhitsdedx\_$insigmapi
hadd All\_$indca\_$inhitsfit\_$inhitsdedx\_$insigmapi.root Rebin*.root
cd ../
root -b -q -l "average.C($indca, $inhitsfit, $inhitsdedx, $insigmapi)"

done
done
done
done


