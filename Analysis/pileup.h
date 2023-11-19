// TF1* f_upper;
// TF1* f_lower;
// f_upper = new TF1("f_upper","pol5",0,350);
// f_upper->SetParameters(6.32816,0.689232,-0.00185181,6.31563e-06,-8.29481e-09);
// f_lower = new TF1("f_lower","pol5",0,350);
// f_lower->SetParameters(-5.20165,0.144438,0.00186397,-1.28471e-05,4.28608e-08);

bool pileupRejection( double vz, int refmult, int NtofMatch )
{
    // std::cout << "start pileupRejection vz = " << vz << " refmult = " << refmult << "nTofMatch = " << NtofMatch  << " " << std::endl;
double b0=0,b1=0,b2=0,b3=0,b4=0;
double d0=0,d1=0,d2=0,d3=0,d4=0;

if(vz<-87.0)// -145,-87 cm
			{
				b0=25.6055790979197, b1=2.02528136596901, b2=-0.0058370984051939, b3=2.59602314466234e-05, b4=-5.3014743584261e-08;
				d0=-17.7059596791057, d1=0.614538168662738, d2=0.00534180935164814, d3=-1.79582873880806e-05, d4=1.01623054170579e-08;

			}
			else if(vz<-29.0)// -87,-29 cm
			{
				b0=23.0160060308621, b1=1.61885832757588, b2=-0.00275873189631398, b3=1.31262550392554e-05, b4=-2.94368020941846e-08;
				d0=-17.3591842617911, d1=0.796170989774258, d2=0.000670722514533827, d3=3.26258075150876e-06, d4=-1.60611460182112e-08;

			}
			else if(vz<29.0)// -29,29 cm
			{
				b0=16.4277056306649, b1=1.71652229539398, b2=-0.00406847684302521, b3=1.65203560938885e-05, b4=-2.96250329214512e-08;
				d0=-15.7887025834219, d1=0.789786364309292, d2=-0.000637115144252616, d3=1.00019972792727e-05, d4=-2.45208851616324e-08;

			}
			else if(vz<29.0)// 29,87 cm
			{
				b0=21.2024767158778, b1=1.70521848381614, b2=-0.00352260930859763, b3=1.60905730948817e-05, b4=-3.37443468806432e-08;
				d0=-17.1166088395929, d1=0.814739436616432, d2=0.000227197779215977, d3=6.55397838050604e-06, d4=-2.28812912596058e-08;

			}
			else// 87,145 cm
			{
				b0=26.0970905882739, b1=1.88889714311734, b2=-0.00195374948885512, b3=-6.14244087431038e-06, b4=1.99930095058841e-08;
				d0=-15.6624325989392, d1=0.52385751891358, d2=0.00794996911844969, d3=-4.09239155250494e-05, d4=6.40163739983216e-08;

			}

bool limit1 = false;
bool limit2 = false;

if(refmult<d0+
            d1*NtofMatch+
			d2*NtofMatch*NtofMatch+
			d3*NtofMatch*NtofMatch*NtofMatch+
			d4*NtofMatch*NtofMatch*NtofMatch*NtofMatch){return false;}
if(refmult>b0+
            b1*NtofMatch+
			b2*NtofMatch*NtofMatch+
			b3*NtofMatch*NtofMatch*NtofMatch+
			b4*NtofMatch*NtofMatch*NtofMatch*NtofMatch){return false;}

    // std::cout << " limit 1 " << limit1
    //           << " limit 2 " << limit2
    //           << (limit1 && limit2) << std::endl;
    // return (limit1 && limit2);
    return true;
}

// bool nPi_K_P_rejection(int refmult, int nPi_K_P )
// {
	
// 	if ( nPi_K_P >= f_lower->Eval(refmult) && nPi_K_P < f_upper->Eval(refmult))
// 	{
// 		return kTRUE; // pass the cut
// 	} else return kFALSE; // did not pass the cut
// }