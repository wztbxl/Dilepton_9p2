#!/bin/bash

rm -f test/runtestLog
rm -f test/test.root
root4star -b -l <<EOF >& test/runtestLog
.O2
.x doEvent.C 
.q
EOF
~    
