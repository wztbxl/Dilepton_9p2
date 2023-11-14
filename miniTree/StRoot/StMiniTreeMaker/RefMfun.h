// temporary using the refmult 3 for the centrality defination for 7.7 GeV
// base on the analysis at https://drupal.star.bnl.gov/STAR/system/files/collab_meet_spring_2023.pdf
// Vz and more correction need to be added
// 9 for 080, 1 for 70-80%
Int_t GetCentralityRefMult3(Double_t refMultCorr){
  if(refMultCorr > 308){
    return 9;
  }else if (refMultCorr > 257){
    return 8;
  }else if (refMultCorr > 179){
    return 7;
  }else if (refMultCorr > 121){
    return 6;
  }else if (refMultCorr > 79){
    return 5;
  }else if (refMultCorr > 49){
    return 4;
  }else if (refMultCorr > 28){
    return 3;
  }else if (refMultCorr > 15){
    return 2;
  }else if (refMultCorr > 7){ 
    return 1;
  }else{
    return 0;
  }
} 

Int_t GetCentrality(Double_t refMultCorr){
  if(refMultCorr > 228){
    return 9;//0-5%
  }else if (refMultCorr > 191){
    return 8;//5-10%
  }else if (refMultCorr > 137){
    return 7;//10-20%
  }else if (refMultCorr > 97){
    return 6;//20-30
  }else if (refMultCorr > 68){
    return 5;//30-40%
  }else if (refMultCorr > 45){
    return 4;//40-50%
  }else if (refMultCorr > 29){
    return 3;//50-60%
  }else if (refMultCorr > 18){
    return 2;
  }else if (refMultCorr > 10){
    return 1;
  }else{
    return 0;
  }
} 