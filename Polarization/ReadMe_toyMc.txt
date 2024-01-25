********************************************************************************
* I use toyMC to estimate the bias of the measured polarization parameters.    *
* The correction factors are then applied to get the final results.            *
* This part is described in analysis note Sect. 8.1.                           *
********************************************************************************

The central value of input lambda_theta and lambda_phi in toyMC study are the results of the 24 track cut variations. The lambda_theta and lambda_phi is changed by +/-0.1 and +/-0.2 separately.  
The results from toyMC study are in dimuon/Polarization/toyMcCorr：
  CS_Phi_iter.root CS_Theta_iter.root HX_Phi_iter.root HX_Theta_iter.root 

I have copied result_track_cuts.txt to this directory. (dimuon/Polarization/trackSys/result_track_cuts.txt)
Then:
   - root -b -q -l drawComAll.C 
   The useful out put files are: correct_factor_Cs.txt and correct_factor_Hx.txt


*******************************************************
* The following is to introduce the toyMC simulation. *
*******************************************************
I have set up a set of lambda_theta and lambda_phi as an example, according to the trackSys/result_track_cuts.txt. 

I use 500 events in my analysis. Here, we can use 5 events for a quick test.
In HXFrames:
 I. Generate 5 pseudo-data
   In testData.C: change line73 myDir; line77 SetSeed(0)
   In runAllPseudoData.sh: change line19 nLoop<=5
   - ./runAllPseudoData.sh 
   when job finished 
   - ./addAllPseudoData.sh 
 II. Generate pseudo-efficiency and do the iteration procedure
   In testData.C, change line77 SetSeed(1000)  
   (a) In runAll.sh, line9: iteration=1, line16: nLoop<=1
   (b) After job finished: 
        In defs.h, line3: change gIter=1; line4 nExpr=5
        Then do: root -b -q -l ana_JpsiPol.C
   (c) In runAll.sh, line 9: iteration=2, line 16: nLoop<=5
   (d) After job finished: 
        In defs.h, line 3: change gIter=2
        Then do: root -b -q -l ana_JpsiPol.C
       Repeat (c) and (d) change iteration and gIter from 2 to 8; -> For iterative procedure. 
       I have put the results of 500 events in ./Resultfiles/Run15_pp200.JpsiPolPar.root 
 III. Analysis final results
   Using root file ./Resultfiles/Run15_pp200.JpsiPolPar.root  
   * Notice: need to use the root file from CVS. * 
   - ./run_HX_final.sh
   The generated HX_compare_8.root will be used in the procedure of estimating the bias. 

The procedures are the same in directories of HXFrames (for helicity frame) and CSFrame (for Collins-soper frame). 


    


++++++++++++++++++++++++++ Parameters explanation ++++++++++++++++++++++++++++++
Not necessarily to read when doing code QA.

In testData.C, the meaning of arguments are:
- nExpr: number of events run in toyMC for each pseudo-experiment
- nLoop: index of pseudo-experiment. For example, if one would like to run 500 pseudo-experiments, this arguemnt needs to be set from 1 to 500 in runAll.sh
- iter: index of iteration: i) iter = 0: generate pseudo-data with realistic statistics (nExpr in runAllPseudoData.sh); ii) iter = 1: generate pseudo-embedding using polarization parametes equal to 0; iii) 2 <= iter < 10: generate pseudo-embedding using polarization parameters extracted from previous iteration
- LamThe, LamPhi, LamThePhi: input polarization parameters. They are not used if iter is set to [1,10)

In runAllPseudoData.sh and runAllCorr.sh:
- parstheta, parsphi: arrays of input polarization parameters
- iteration: index of iteration
- nLoop: index of pseudo-experiment
- nExpr is set in line "echo "Arguments = $nLoop $iteration 1e7 $ltheta $lphi">>$script"

+++++ Analyze pseudo-data +++++
This is done using ana_JpsiPol.C, and defs.h is also needed.

In defs.h:
- gIter: index of iteration under analysis
- nExpr: number of pseudo-experiments to be analysed. It can be set to a small number for testing purpose. 
- gUseMean: use mean of the distribution of the polarization parameters instead of fitting to get the averge polarization parameters. Set it to 0 unless needed, e.g. fitting does not work well.

In ana_JpsiPol.C, the following functions should be run in sequence:
- fit_test_EffCorr(): get the polarization parameters for a given iteration index set in defs.h for all pseudo-experiments


