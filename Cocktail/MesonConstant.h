#ifndef MESONCONSTANT_H 
#define MESONCONSTANT_H  

enum ParticleTypes {
	virtualphoton = 0,
	photon	= 1,	
	pi0		= 111, 
	eta		= 221,
	rho		= 113,
	omega	= 223,
	phi		= 333,
	etaprim	= 331,
	jpsi	= 443,
	psi		= 100443
};

enum DecayMode{
	twobody = 2,
	dalitz  = 3
};

//unit GeV/c^2
const Double_t Massphoton	= 0.0;
const Double_t Masselectron	= 0.000511;
const Double_t Masspi0 		= 0.1349766;
const Double_t Masseta		= 0.547853;
const Double_t Massrho		= 0.77549;
const Double_t Massomega  	= 0.78265;
const Double_t Massphi  	= 1.019455;
const Double_t Massetaprim 	= 0.95778;
const Double_t Massjpsi  	= 3.096916;
const Double_t Masspsi  	= 3.68609;

//full width unit GeV/c^2
const Double_t Widthphoton	= 0.;
const Double_t Widthpi0 	= 0.;
const Double_t Widtheta		= 1.3e-6;
const Double_t Widthrho		= 0.1491;
const Double_t Widthomega 	= 8.49e-3;
const Double_t Widthphi  	= 4.26e-3;
const Double_t Widthetaprim	= 1.94e-4;
const Double_t Widthjpsi  	= 9.29e-5;
const Double_t Widthpsi  	= 3.04e-4;

#endif
