#!/bin/bash
date

 if [ $# -ne 1 ]; then
    echo "Please input the centrality parameter, which should be \"080\", \"010\", \"1040\", \"4080\", \"4060\", \"6080\", \"6070\", or \"7080\""
    exit 0
 fi
 
 cen=$1
 
 if [ $cen != "080" -a $cen != "010" -a $cen != "1040" -a $cen != "4080" -a $cen != "4060" -a $cen != "6080" -a $cen != "6070" -a $cen != "7080" ]; then
    echo "Please check the centality string, which should be \"080\", \"010\", \"1040\", \"4080\", \"4060\", \"6080\", \"6070\", or \"7080\" !"
    exit 0
 fi

 origDir=$PWD
 echo ${origDir}

 #outDir="output/Cen${1}/"
 outDir="output/Cen${1}/"
 outDir1="output/TEMP/Cen${1}_R0"

#  for meson in "eta" 
 for meson in "eta" "etaprim" "omega" "phi" "pi0"
#
 do
   echo ${outDir}${meson}
   cd ${origDir}/${outDir}${meson}
   mkdir -p hadddalitz
   mkdir -p hadddalitz/raw
   mv ${meson}dalitz_*.root hadddalitz/raw
   cd hadddalitz
   ~/Scripts/Hadd/hadd.sh 0 raw round0 50
   

 done

for meson in "omega" "phi" "jpsi"

 do
   echo ${outDir}${meson}
   cd ${origDir}/${outDir}${meson}
   mkdir -p haddtwobody
   mkdir -p haddtwobody/raw
   mv ${meson}2ee_*.root haddtwobody/raw
   cd haddtwobody
   ~/Scripts/Hadd/hadd.sh 0 raw round0 50
   

 done

# cd ${origDir}/output/Cen${cen}/pi0
# mkdir -p hadddalitz
# mkdir -p hadddalitz/raw
# mv    pi0dalitz_*.root hadddalitz/raw
# cd    hadddalitz 
# ~/script/Hadd/hadd.sh 0 raw round0 50
#
# cd ${origDir}/output/Cen${cen}/eta
# mkdir -p hadddalitz
# mkdir -p hadddalitz/raw
# mv    *dalitz_*.root hadddalitz/raw
# cd    hadddalitz 
# ~/script/Hadd/hadd.sh 0 raw round0 50
#
# cd ${origDir}/output/Cen${cen}/etaprim
# mkdir -p hadddalitz
# mkdir -p hadddalitz/raw
# mv    *dalitz_*.root hadddalitz/raw
# cd    hadddalitz 
# ~/script/Hadd/hadd.sh 0 raw round0 50
#
# cd ${origDir}/output/Cen${cen}/omega
# mkdir -p haddtwobody hadddalitz
# mkdir -p haddtwobody/raw hadddalitz/raw
# mv    *2ee_*.root   haddtwobody/raw
# mv    *dalitz_*.root hadddalitz/raw
# cd    haddtwobody
# ~/script/Hadd/hadd.sh 0 raw round0 50
# cd    ../hadddalitz 
# ~/script/Hadd/hadd.sh 0 raw round0 50
#
# cd ${origDir}/output/Cen${cen}/phi
# mkdir -p haddtwobody hadddalitz
# mkdir -p haddtwobody/raw hadddalitz/raw
# mv    *2ee_*.root   haddtwobody/raw
# mv    *dalitz_*.root hadddalitz/raw
# cd    haddtwobody
# ~/script/Hadd/hadd.sh 0 raw round0 50
# cd    ../hadddalitz 
# ~/script/Hadd/hadd.sh 0 raw round0 50
# 
# cd ${origDir}/output/Cen${cen}/jpsi
# mkdir -p haddtwobody
# mkdir -p haddtwobody/raw
# mv    *2ee_*.root haddtwobody/raw
# cd    haddtwobody 
# ~/script/Hadd/hadd.sh 0 raw round0 50
#
# cd ${origDir}/output/Cen${cen}/psi
# mkdir -p haddtwobody
# mkdir -p haddtwobody/raw
# mv    *2ee_*.root haddtwobody/raw
# cd    haddtwobody 
# ~/script/Hadd/hadd.sh 0 raw round0 50

# cd ${origDir}/output/Cen${cen}/virtualphoton
# mkdir -p haddtwobody
# mkdir -p haddtwobody/raw
# mv    *2ee_*.root haddtwobody/raw
# cd    haddtwobody 
# ~/script/Hadd/hadd.sh 0 raw round0 50
cd ${origDir}/output/Cen${cen}/


 cd ${origDir}
