//refMult Function from qian
/*
 * in his definition centrality 9 is 0-5%
 * 1 is 70-80%
 * and 0 is event that refmult < 7
 */
Double_t GetRefMultCorr(const Int_t RefMult, const Double_t z);
Double_t GetWeight(Double_t refMultCorr);
Int_t GetCentrality(Double_t refMultCorr);

Double_t GetRefMultCorr(const Int_t RefMult, const Double_t z)
{
	const Double_t par_z[7] = {435.9, -0.02413, -0.003707, 0.0002204, 1.487e-5, -2.95e-7, -1.867e-8};
	const Double_t RefMult_ref = par_z[0];
	const Double_t  RefMult_z = par_z[0] + par_z[1]*z + par_z[2]*z*z + par_z[3]*z*z*z + par_z[4]*z*z*z*z + par_z[5]*z*z*z*z*z + par_z[6]*z*z*z*z*z*z; // Parametrization of mean RefMult vs. z_vertex position
	Double_t Hovno =1.0;
	if(RefMult_z > 0.0)
	{
		Hovno = (RefMult_ref + par_z[7])/RefMult_z;
	}
	Double_t RefMult_d = (Double_t)(RefMult) + gRandom->Rndm(); // random sampling over bin width -> avoid peak structures in corrected distribution
	return RefMult_d*Hovno;
}

Double_t GetWeight(Double_t refMultCorr)
{
	Double_t Weight = 1.0;

	Double_t par0 =3.905;
	Double_t par1 = -204.4;
	Double_t par2 = 1.851;
	Double_t par3 = 24.32;
	Double_t par4 = -0.01746;
	Double_t A    = 0;
	Double_t par6 = 6405;
	Double_t par7 = 3.743e-05;

	// Additional z-vertex dependent correction
	//const Double_t A = ((1.27/1.21))/(30.0*30.0); // Don't ask...
	//const Double_t A = (0.05/0.21)/(30.0*30.0); // Don't ask...
	Double_t mVz = 30.0; // this has to be modified...

	if(
			refMultCorr != -(par3/par2) // avoid denominator = 0
			&& refMultCorr<=70         // normalization stop at refMultCorr=300
	  )
	{
		Weight = par0 + par1/(par2*refMultCorr + par3) + par4*(par2*refMultCorr + par3) + par6/pow(par2*refMultCorr+par3 ,2) + par7*pow(par2*refMultCorr+par3 ,2); // Parametrization of MC/data RefMult ratio
		//Weight = par0 + par1/(par2*refMultCorr + par3) + par4*(par2*refMultCorr + par3) ; // Parametrization of MC/data RefMult ratio
		Weight = Weight + (Weight-1.0)*(A*mVz*mVz); // z-dependent weight correction
	}

	return Weight;
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
