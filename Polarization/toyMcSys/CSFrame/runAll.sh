#!/bin/bash

wdir=$PWD

# CS frame
parstheta=(-0.260 0.480 0.630 -0.047)
parsphi=(-0.034 0.008 -0.124 -0.300)

iteration=2

for((i=0; i<4; i++))
do
    ltheta=${parstheta[$i]}
    lphi=${parsphi[$i]}
    echo $ltheta $lphi
    for((nLoop=1; nLoop<=500; nLoop++))
    do
	script=submit.Iter$iteration.Loop$nLoop.theta$ltheta.phi$lphi.con
	echo 'Universe       = vanilla' >$script
	echo 'Notification   = never'>> $script
	echo 'Requirements   = (CPU_Type != "crs") && (CPU_Experiment == "star")'>>$script
	echo '+Experiment    = "star"'>>$script
	echo 'Priority       = +10'>>$script
	echo '+Job_Type      = "cas"'>>$script
	echo 'GetEnv         = true'>>$script
	echo "Executable     = $wdir/run.sh">>$script
	echo "Arguments      = $nLoop $iteration 1e7 $ltheta $lphi">>$script
	echo "Log            = $wdir/log/Iter$iteration.Loop$nLoop.log">>$script
	echo "Output         = $wdir/log/Iter$iteration.Loop$nLoop.out">>$script
	echo "Error          = $wdir/log/Iter$iteration.Loop$nLoop.err">>$script
	echo 'Queue'>>$script
    condor_submit $script
	echo "submit Iter$iteration Loop$nLoop complete"
    done
done

mv *.con SubmitFile
