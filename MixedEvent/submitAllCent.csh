#!/bin/bash

echo "before submitting jobs, you Should check the file tag is want you want!!!!!!!!!!!!!!!"

read -p "Do you want to continue? (yes/no): " user_input
if [ "$user_input" = "yes" ]; then
    echo "Continuing with the next steps..."

    for i in {1..9}
    do
        echo "Submitting job for cent$i with i=$i"
        ./runAll.csh all cent${i}_VzBin10_EPBin12_Buff500_mixedOnly $i
    done
else
    echo "Aborted. Script will not continue."
    exit 1
fi
