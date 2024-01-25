Use the miniTree files to measure the J/psi polarization

++++ the main procedure +++++
Part I. Calculate the cos_theta and phi distribution 
   - Using the analysis.cxx. 
   In each sub-directories using the executable file analysis which is compiled using GNUmakefile

Part II. Fit to get the J/psi yields in 10 cos_theta bins and 15 phi bins for each J/psi pT interval.
   - Using the ana_data*.C 
   The fitting procedures are varied and using the mean of the raw counts and the raw cos_theta and phi distribution

Part III. Correct the cos_theta and phi distribution using iterative procedure
   - The code are in sub-directories iter_code
++++++++++++++++++++++++++++++

Notice: The final results are the mean of total 24 track cut variations

++++++++++ detailed procedure ++++++++++++++++++++
Systematic uncertainties:
++ For tracking and muon identification:++
  1) Get raw counts: 
     - cd trackSys
     - ./trackSys.csh
    after job finish
     - ./addFiles.csh
    To get the number of J/psi in different angular bins and in different dimuon pT bins
     - ./getResult.sh
        The codes ana_data0.C, ana_data1.C, ana_data2.C and average.C are used. The ROOT version I used is 5.34/36. 
        ana_data*.C is for fitting procedure which takes time. This step generate Rebin*.root 
        average.C is to get the mean value. The porpose is described in AnalysisNote chapter 5.2.
        ------------------------------------------------------------------------------------------------------------------
        + You can just test a set of track quality cut: eg: dca < 3; nHitsFit > 15; nHitsDedx > 10; -2 < nSigmaPi < 3    +
        + - root -b -q -l "ana_data0.C(3, 15, 10, -2)"                                                                   +
        + - root -b -q -l "ana_data1.C(3, 15, 10, -2)"                                                                   + 
        + - root -b -q -l "ana_data1.C(3, 15, 10, -2)"                                                                   +
        + - hadd Data_3_15_10_-2/All_3_15_10_-2.root Data_3_15_10_-2/Rebin*.root                                         +
        + - root -b -q -l "average.C(3, 15, 10, -2)"                                                                     +
        ------------------------------------------------------------------------------------------------------------------
        The plots can be found in: /star/u/liuzhen/Run15/Run15pp200/submitjob/Polarization/trackSys_fitting_plots.zip
        The output root files can be found in: /star/u/liuzhen/Run15/Run15pp200/submitjob/Polarization/trackSys/Data*
        The meaning of the variables are described in: /star/u/liuzhen/lfsupc/Run15/final_signal_extraction_plots/readme.txt 
        The purpose is described in AnalysisNote chapter 7.2.
  2) Iterative procedure:
     In directory: trackSys
      - ./trackIter.csh
     after job finish
      - ./addIterFiles.csh
      This procedure will generate eg: Data_3_15_10_-2/Iter_3_15_10_-2/final_3_15_10_-2.root
     after histogram added
     root -b -q -l draw_Final_Track.C
     The final root files can be found in: /star/u/liuzhen/Run15/Run15pp200/submitjob/Polarization/trackSys/Data*/Iter*
++++

++ For signal extraction:  ++
Use the track quality cuts: nHitsFit > 15; nHitsDedx > 10; dca < 3; -2 < nSigmaPi < 3
  1) Get raw counts:
      in the directory Polarization/signalSys:
      - mkdir Data_3_15_10_-2.0
      - cp ../trackSys/Data_3_15_10_-2.0/Data_3_15_10_-2.0.root Data_3_15_10_-2.0
      - ./runAll.sh
      The fitting results can be found in: /star/u/liuzhen/Run15/Run15pp200/submitjob/Polarization/signalSys/Data_3_15_10_-2.0/Rebin* 
      I have already put the output root file All_3_15_10_-2.0.root in Data_3_15_10_-2.0. 
    
  2) Iterative procedure:
      - ./signalIter.csh
      after job finish, change the directory in addIterFiles.csh
      The output root files can be find in: /star/u/liuzhen/Run15/Run15pp200/submitjob/Polarization/signalSys/File*/Iter*
      - ./addIterFiles.csh
      This procedure will generate eg: File*/Iter_f*/final_*.root 
      - root -b -q -l draw_Final_signal.C
++++

<Notice> 
the detailed root files is in:
    /star/u/liuzhen/Run15/Run15pp200/submitjob/Polarization/
    trackSys_rootfiles.zip    signalSys_rootfiles.zip

