#!/bin/bash

rm -rf PDF/*.pdf
for((i=1; i<8; i++))
do
    rm -rf *Final_iter$i.pdf
    rm -rf getMean_*_iter$i.pdf
done
