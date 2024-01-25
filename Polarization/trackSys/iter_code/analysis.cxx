#include "headers.h"

Bool_t debug = false;
//----------//
int main(int argc, char** argv)
{
    if(argc!=9) {
        cout<<"input parameters: file.list nIterate dca nHitsFit nHitsDedx nsigmaPi frame pt"<<endl;
        return -1;
    }

    TString inFile;
    char outFile[1024];

    if(argc==9){
        inFile = argv[1];
        nIterate = atoi(argv[2]);
        cout<<"nIterate = "<<nIterate<<endl;
        mMaxDca = atof(argv[3]);
        mMinNHitsFit = atoi(argv[4]);
        mMinNHitsDedx = atoi(argv[5]);
        mMinNSigmaPiCut = atof(argv[6]);
        Frame = atoi(argv[7]);
        nPt = atof(argv[8]);
        if(nIterate==1)sprintf(outFile,"efficiency%d",nIterate);
        else{
            sprintf(outFile,"efficiencyFrame%dPt%dIter%d", Frame, nPt, nIterate);
        }
    }

    //---open files--and--add to chain---//
    TChain *chain = new TChain("EmbedTree");
    Int_t ifile=0;
    char filename[512];
    ifstream *inputStream = new ifstream;
    inputStream->open(inFile.Data());
    if (!(inputStream)) {
        printf("can not open list file\n");
        return 0;
    }
    for(;inputStream->good();){
        inputStream->getline(filename,512);
        if(inputStream->good()) {
            TFile *ftmp = new TFile(filename);
            if(!ftmp||!(ftmp->IsOpen())||!(ftmp->GetNkeys())) {
                cout<<"something wrong: "<< filename <<endl;
            } else {
                if(debug) cout<<"read in "<<ifile<<"th file: "<<filename<<endl;
                chain->Add(filename);
                ifile++;
            }
            delete ftmp;
        }
    }
    delete inputStream;

    myRandom = new TRandom3();
    //TDatime *clock = new TDatime();
    //myRandom->SetSeed(clock->GetTime());
    myRandom->SetSeed(1000);

    bookHistograms();

    //-- input parameters --//
    if(nIterate==1){
        lamThe=0;
        lamPhi=0;
    }
    else{
        TFile *fPar = TFile::Open(Form("%sPt%dIter%d.root",frameName[Frame], nPt, nIterate-1));
        TH1D *hLambdaTheta = (TH1D*)fPar->Get(Form("hIter%sCosPt%d", frameName[Frame], nPt));
        lamThe = hLambdaTheta->GetBinContent(nIterate-1);
        if(lamThe>1) lamThe = 1;
        else if(lamThe<-1) lamThe = -1;
        TH1D *hLambdaPhi = (TH1D*)fPar->Get(Form("hIter%sPhiPt%d", frameName[Frame], nPt));
        lamPhi = hLambdaPhi->GetBinContent(nIterate-1);
        if(lamPhi>1) lamPhi = 1;
        else if(lamPhi<-1) lamPhi = -1;
        //cout<<"****"<<Form("%sPt%dIter%d.root",frameName[Frame], nPt, nIterate-1)<<"****"<<hLambdaTheta<<"LambdaTheta:"<<lamThe<<endl;
    }
    cout<<"input lambda_theta and lambd_phi: "<<lamThe<<";  "<<lamPhi<<endl;

    TFile* fEff = TFile::Open("/star/u/liuzhen/Run15/Run15pp200/submitjob/Polarization/trackSys/iter_code/Run15pp200Efficiency.root");
    TH1D* hMuonTrigEff = (TH1D*) fEff->Get("hMuonTrigEff");
    TH2D* hTpcResp = (TH2D*) fEff->Get("hTpcResp");
    TH1D* hMtdRespEffEmbed = (TH1D*) fEff->Get("hMtdRespEffEmbed");

    //TF1* f = new TF1("f","(3.4948*TMath::Power(TMath::Exp((-0.395305)*x)+x/2.91793,(-8.46161)))*4*TMath::Pi()*x",0,20);
    TF1* f = new TF1("f","([0]*TMath::Power(TMath::Exp([1]*x)+x/[2],[3]))*4*TMath::Pi()*x",0,20);
    f->SetParameters(7.287, -0.235, 3.412, -10.533);
    f->SetNpx(200);
    TH1D* hMassShape = (TH1D*) f->GetHistogram();

    EmbedTree *event = new EmbedTree(chain);
    Int_t nEvts = chain->GetEntries();
    cout<<nEvts<<" events"<<endl;
    int perevt = nEvts/100;
    for(Int_t iEvt=0; iEvt<nEvts; iEvt++) {
        if(iEvt%perevt ==0) {
            cout<<"Looking at "<<((float)iEvt/nEvts)*100<<" % .."<<endl;
        }
        event->GetEntry(iEvt);
        hEvent->Fill(0.5);
        if(!passEvent(event)) continue;
        if(debug) cout<<"passEvent"<<endl;
        Int_t nMcTracks = event -> nMcTracks;
        if(debug) cout<<nMcTracks<<endl;

        Int_t runID = event->runID;
        int nEmbedJpsi = 0;
        for(Int_t i=0; i<nMcTracks; i++){
            //  muon1 mc
            Double_t mcpt1 = event -> mcpt[i];
            Int_t mcPkey1 = event->mcPkey[i];
            Double_t mcphi1 = event -> mcphi[i];
            Double_t mceta1 = event -> mceta[i];
            Int_t mcCharge1 = event->mccharge[i];
            if(mcpt1<=0) continue;
            if(event->mcGeantId[i]!=5 && event->mcGeantId[i]!=6) continue;

            // muon1 rc
            Int_t rcCharge1 = event->rcCharge[i];
            Int_t rcBackleg1 = event->rcBackleg[i];
            Int_t rcModule1 = event->rcModule[i];
            Int_t rcNHitsFit1 = event -> rcNHitsFit[i];
            Int_t rcNHitsDedx1 = event -> rcNHitsDedx[i];
            Int_t rcNHitsPoss1 = event -> rcNHitsPoss[i];
            Double_t rcNHitsFrac1 = (1.0*rcNHitsFit1)/(1.0*rcNHitsPoss1);
            Double_t rcDca1 = event -> rcDca[i];
            Double_t rceta1 = event -> rceta[i];
            Double_t rcphi1 = event -> rcphi[i];
            Double_t rcpt1 = event -> rcpt[i];
            Double_t rcNsigmaPi1 = event -> rcNSigmaPi[i];
            Double_t rcDy1 = event -> rcDy[i];
            Double_t rcDz1 = event -> rcDz[i];
            Double_t rcDtof1 = event -> rcDtof[i];

            for(Int_t j=i+1; j<nMcTracks; j++){
                // muon2 mc
                Double_t mcpt2 = event -> mcpt[j];
                Int_t mcPkey2 = event->mcPkey[j];
                Double_t mcphi2 = event -> mcphi[j];
                Double_t mceta2 = event -> mceta[j];
                Int_t mcCharge2 = event->mccharge[j];
                if(mcpt2<=0) continue;
                if(event->mcGeantId[j]!=5 && event->mcGeantId[j]!=6) continue;
                if( mcCharge1 * mcCharge2 > 0 ) continue;
                if( mcPkey1 != mcPkey2 ) continue;

                mcMuon1.SetPtEtaPhiM(mcpt1,mceta1,mcphi1,muMass);
                mcMuon2.SetPtEtaPhiM(mcpt2,mceta2,mcphi2,muMass);
                mcJpsi = mcMuon1 + mcMuon2;
                if(TMath::Abs(mcJpsi.Rapidity())>mPairYCut) continue;

                hJpsiSpec[0][0]->Fill(mcJpsi.Pt());
                hJpsiSpec[0][1]->Fill(mcJpsi.Pt(), weightFunc(mcJpsi.Pt(), hMassShape));

                nEmbedJpsi++;

                if(mcCharge1>0) {
                    calPolarization(mcMuon1, mcJpsi, mhMcMPtCosPhi, mhMcMPtCosPhiCS,  weightFunc(mcJpsi.Pt(), hMassShape), lamThe, lamPhi);
                }
                else if(mcCharge2>0) {
                    calPolarization(mcMuon2, mcJpsi, mhMcMPtCosPhi, mhMcMPtCosPhiCS,  weightFunc(mcJpsi.Pt(), hMassShape), lamThe, lamPhi);
                }

                //  muon2 rc
                Int_t rcBackleg2 = event->rcBackleg[j];
                Int_t rcModule2 = event->rcModule[j];
                Int_t rcNHitsFit2 = event -> rcNHitsFit[j];
                Int_t rcNHitsDedx2 = event -> rcNHitsDedx[j];
                Int_t rcNHitsPoss2 = event -> rcNHitsPoss[j];
                Double_t rcNHitsFrac2 = (1.0*rcNHitsFit2)/(1.0*rcNHitsPoss2);
                Double_t rcDca2 = event -> rcDca[j];
                Double_t rceta2 = event -> rceta[j];
                Double_t rcphi2 = event -> rcphi[j];
                Double_t rcpt2 = event -> rcpt[j];
                Double_t rcNsigmaPi2 = event -> rcNSigmaPi[j];
                Double_t rcDy2 = event -> rcDy[j];
                Double_t rcDz2 = event -> rcDz[j];
                Double_t rcDtof2 = event -> rcDtof[j];
                Int_t rcCharge2 = event->rcCharge[j];

                rcMuon1.SetPtEtaPhiM(rcpt1,rceta1,rcphi1,muMass);
                rcMuon2.SetPtEtaPhiM(rcpt2,rceta2,rcphi2,muMass);
                rcJpsi = rcMuon1 + rcMuon2;
                if(TMath::Abs(rcJpsi.Rapidity())>mPairYCut) continue;

                if(rcCharge1*rcCharge2>0) continue;
                if(mcPkey1!=mcPkey2) continue;

                if(!passTrack(rcpt1,rcNHitsDedx1,rcNHitsFit1,rcNHitsFrac1,rcDca1,rceta1,rcphi1)) continue;
                if(!passTrack(rcpt2,rcNHitsDedx2,rcNHitsFit2,rcNHitsFrac2,rcDca2,rceta2,rcphi2)) continue;
                if(debug) cout<<"pass track"<<endl;

                //  smear
                rcpt1 = smearPt(rcpt1); rcpt2 = smearPt(rcpt2);

                // tpc response
                if(!calTpcRespEff(rcpt1, rceta1, rcphi1, hTpcResp)) continue;
                if(!calTpcRespEff(rcpt2, rceta2, rcphi2, hTpcResp)) continue;
                if(debug) cout<<"tpc response"<<endl;

                hJpsiSpec[1][0]->Fill(rcJpsi.Pt());
                hJpsiSpec[1][1]->Fill(rcJpsi.Pt(),weightFunc(mcJpsi.Pt(), hMassShape));

                //  have MTD hit
                if(rcModule1==0 || rcBackleg1==0) continue;
                if(rcModule2==0 || rcBackleg2==0) continue;//already have eta cuts
                if(!checkMuCandidate( rcpt1, rcNsigmaPi1, rcDy1, rcDz1, rcDtof1, 1 )) continue;
                if(!checkMuCandidate ( rcpt2, rcNsigmaPi2, rcDy2, rcDz2, rcDtof2, 1 )) continue;

                if(debug) cout<<"have MTD hit"<<endl;

                if(!isActiveModule(runID, rcBackleg1, rcModule1)) continue;
                if(!isActiveModule(runID, rcBackleg2, rcModule2)) continue;

                hHitMap->Fill(rcModule1, rcBackleg1);

                TH1D* hMtdRespEff1 = (TH1D*)fEff->Get(Form("MtdRespEffCosmic_BL%d_Mod%d",rcBackleg1,rcModule1));
                TH1D* hMtdRespEff2 = (TH1D*)fEff->Get(Form("MtdRespEffCosmic_BL%d_Mod%d",rcBackleg2,rcModule2));

                if(!calMtdResponseEff(rcpt1, hMtdRespEff1, hMtdRespEffEmbed)) continue;
                if(!calMtdResponseEff(rcpt2, hMtdRespEff2, hMtdRespEffEmbed)) continue;

                hJpsiSpec[2][0]->Fill(rcJpsi.Pt());
                hJpsiSpec[2][1]->Fill(rcJpsi.Pt(),weightFunc(mcJpsi.Pt(), hMassShape));

                //  efficiency correction
                if(!calMtdTriggerEff(rcpt1,hMuonTrigEff)) continue;
                if(!calMtdTriggerEff(rcpt2,hMuonTrigEff)) continue;

                hJpsiSpec[3][0]->Fill(rcJpsi.Pt());
                hJpsiSpec[3][1]->Fill(rcJpsi.Pt(),weightFunc(mcJpsi.Pt(), hMassShape));

                rcMuon1.SetPtEtaPhiM(rcpt1,rceta1,rcphi1,muMass);
                rcMuon2.SetPtEtaPhiM(rcpt2,rceta2,rcphi2,muMass);
                rcJpsi = rcMuon1 + rcMuon2;

                hHitEtaVsPhi->Fill(rceta1, rcphi1);

                if(rcCharge1>0) {
                    calPolarization(rcMuon1, rcJpsi, mhRcMPtCosPhi, mhRcMPtCosPhiCS,  weightFunc(mcJpsi.Pt(), hMassShape), lamThe, lamPhi);
                }
                else if(rcCharge2>0) {
                    calPolarization(rcMuon2, rcJpsi, mhRcMPtCosPhi, mhRcMPtCosPhiCS,  weightFunc(mcJpsi.Pt(), hMassShape), lamThe, lamPhi);
                }
            }//muon2
        }//muon1
        hNEmbedJpsi->Fill(nEmbedJpsi);

    }//event loop

    writeHistograms(outFile);
    fEff->Close();
    delete chain;
    cout<<"end of program"<<endl;
    return 0;
}

//-----------------------------------------------
void bookHistograms()
{
    hEvent = new TH1D("hEvent","# of event",5,0,5);
    hEvent->GetXaxis()->SetBinLabel(1,"all event");
    hEvent ->GetXaxis()->SetBinLabel(2,"TPC Zero Vertex Cut");
    hEvent ->GetXaxis()->SetBinLabel(3,"TPC Vr Cut");
    hEvent ->GetXaxis()->SetBinLabel(4,"TPC Vz Cut");
    hEvent ->GetXaxis()->SetBinLabel(5,"Vz diff Cut");

    hEtaVsPhi = new TH2D("hEtaVsPhi","#varphi VS #eta of reconstructed tracks in embedding;#eta;#varphi",40,-1,1,48,-PI,PI);
    hHitMap = new TH2D("hHitMap","reconstructed MC track hit map in embedding; module; backleg",7,0,7,32,0,32);
    hHitEtaVsPhi = new TH2D("hHitEtaVsPhi","#varphi VS #eta of reconstructed tracks in embedding hit with MTD;#eta;#varphi",40,-1,1,48,-PI,PI);

    const Int_t dim = 4;
    //const Int_t nBins[dim]={200,400,40,96};
    const Int_t nBins[dim]={200,400,40,60};
    const Double_t lowBins[dim]={2,0,-1,0.};
    const Double_t highBins[dim]={4,20,1,TMath::Pi()};

    mhMcMPtCosPhi = new THnSparseF("mhMcMPtCosPhi","Mc HF UnlikeSignPair;m_{#mu#mu} (GeV/c^{2});p_{T} (GeV/c);cos#theta;#phi",dim,nBins,lowBins,highBins);
    mhMcMPtCosPhi->Sumw2();
    mhMcMPtCosPhiCS = new THnSparseF("mhMcMPtCosPhiCS","Mc CS UnlikeSignPair;m_{#mu#mu} (GeV/c^{2});p_{T} (GeV/c);cos#theta;#phi",dim,nBins,lowBins,highBins);
    mhMcMPtCosPhiCS->Sumw2();
    mhRcMPtCosPhi = new THnSparseF("mhRcMPtCosPhi","Rc HF UnlikeSignPair;m_{#mu#mu} (GeV/c^{2});p_{T} (GeV/c);cos#theta;#phi",dim,nBins,lowBins,highBins);
    mhRcMPtCosPhi->Sumw2();
    mhRcMPtCosPhiCS = new THnSparseF("mhRcMPtCosPhiCS","Rc CS UnlikeSignPair;m_{#mu#mu} (GeV/c^{2});p_{T} (GeV/c);cos#theta;#phi",dim,nBins,lowBins,highBins);
    mhRcMPtCosPhiCS->Sumw2();

    //  efficiency spectrum
    const char* nSbs[4] = {"MC", "Tpc", "MtdMth", "MtdTrig"};
    for(Int_t i=0; i<4; i++){
        hJpsiSpec[i][0] = new TH1D(Form("hJpsiSpec_%s", nSbs[i]), Form("hJpsiSpec_%s", nSbs[i]), 20, 0, 10);
        hJpsiSpec[i][0]->Sumw2();
        hJpsiSpec[i][1] = new TH1D(Form("hJpsiSpec_%s_w", nSbs[i]), Form("hJpsiSpec_%s_w", nSbs[i]), 20, 0, 10);
        hJpsiSpec[i][1]->Sumw2();
    }

    hNEmbedJpsi = new TH1D("hNEmbedJpsi", "hNEmbedJpsi; # of embedded j/psi per event; N", 50, 0, 50);

}

//-----------------------------------------------
Bool_t passEvent(EmbedTree *evt)
{
    Float_t tpcVx = evt -> tpcVx;
    Float_t tpcVy = evt -> tpcVy;
    Float_t tpcVz = evt -> tpcVz;
    Float_t tpcVr = sqrt(tpcVx*tpcVx+tpcVy*tpcVy);
    if (!(fabs(tpcVx)>0 || fabs(tpcVy)>0 || fabs(tpcVz)>0)) return false;
    hEvent->Fill(1.5);
    if (tpcVr > mMaxVr) return false;
    hEvent->Fill(2.5);
    if (fabs(tpcVz)>mMaxVz) return false;
    hEvent->Fill(3.5);
    return true;
}

//-----------------------------------------------
Bool_t passTrack(Double_t pt, Int_t nHitsDedx, Int_t nHitsFit, Double_t nHitsFrac, Double_t dca, Double_t eta, Double_t phi)
{
    if(pt < mMinTrkPt   || pt > mMaxTrkPt) return false;
    if (mMaxDca<1e4 && dca>mMaxDca) return false;
    if(nHitsFit<mMinNHitsFit) return false;
    if (nHitsDedx < mMinNHitsDedx) return false;
    if(1.*nHitsFrac<mMinFitHitsFaction) return false;
    if(TMath::Abs(eta)>mMaxTrkEta) return false;
    return true;

    return true;
}

//-----------------------------------------------
Double_t smearPt(double pT)
{
    pT = pT * myRandom->Gaus(1., mSmearPtRes * sqrt(pT));
    return pT;
}

//-----------------------------------------------
Double_t weightFunc(Double_t pT, TH1D* hShape)
{
    Double_t w = hShape->GetBinContent(hShape->FindBin(pT));
    return w;
}

//-------------------------------------------------------
Double_t weightPola( double ptheta, double pphi, double costheta, double phi )
{
    double w = (1 + ptheta*costheta*costheta + pphi*(1-costheta*costheta)*TMath::Cos(2*phi))/((3+ptheta)/3.);
    //double w = 1 + ptheta*costheta*costheta + pphi*(1-costheta*costheta)*TMath::Cos(2*phi);
    return w;
}

//-------------------------------------------------------
Bool_t checkMuCandidate ( Double_t pt, Double_t nSigmaPion, Double_t dy, Double_t dz, Double_t dtof, Int_t matchFlag = 1 )
{
    if(debug) cout<<"in check muon candidate"<<endl;
    if(
            pt > mMuPtCut &&
            mMinNSigmaPiCut < nSigmaPion && nSigmaPion < mMaxNSigmaPiCut &&
            checkMuDeltaY(pt, dy) &&
            checkMuDeltaZ(pt, dz) &&
            checkMuDeltaTof(dtof) &&
            matchFlag > 0
      )return true;
    else return false;
}
//-------------------------------------------------------
Bool_t checkMuDeltaY(Double_t pt, Double_t dy)
{
    // Pt dependent cut
    double sigy   =  -17.6867 + 18.4528*exp(0.637142/pt);
    double marg   = 0.;
    if(pt>3.) marg = 0.5;
    if( fabs(dy)  <= (mMuDYSigCut+marg)*sigy ) return true;
    else return false;
}
//-------------------------------------------------------
Bool_t checkMuDeltaZ(Double_t pt, Double_t dz)
{
    // Pt dependent cut
    double sigz   =  -32.6793 + 32.6034*exp(0.444217/pt);
    double marg   = 0.;
    if(pt>3.) marg = 0.5;
    if( fabs(dz)  <= (mMuDZSigCut+marg)*sigz ) return true;
    else return false;
}
//-------------------------------------------------------
Bool_t checkMuDeltaTof(Double_t dtof)
{
    return true;
    //if( mMuMindTofCut <= dtof && dtof  <= mMuMaxdTofCut ) return true;
    //else return false;
}

//-----------------------------------------------
Bool_t  calTpcRespEff(Double_t pt, Double_t eta, Double_t phi, TH2D* hEff)
{
    Double_t newPhi = phi;
    if(phi<0) newPhi = phi + 2*PI;
    if(eta>0.2) return true;
    if(pt<1) pt = 1;
    if(pt>=10) pt = 10;
    Double_t eff = hEff->GetBinContent(hEff->GetXaxis()->FindBin(pt), hEff->GetYaxis()->FindBin(newPhi));
    Double_t p = myRandom->Uniform(1);
    if(p<eff) return true;
    else return false;
}

//-----------------------------------------------
Bool_t  calMtdTriggerEff(Double_t pt, TH1D* hEff)
{
    Double_t eff = hEff->GetBinContent(hEff->FindBin(pt));
    Double_t p = myRandom->Uniform(1);
    if(p<eff) return true;
    else return false;
}

//-----------------------------------------------
Bool_t calMtdResponseEff(Double_t pt, TH1D* hCosmic, TH1D* hEmb)
{
    Double_t eff = hCosmic->GetBinContent(hCosmic->FindBin(pt))/hEmb->GetBinContent(hEmb->FindBin(pt));
    Double_t p = myRandom->Uniform(1);
    if(p<eff) return true;
    else return false;
}

//-----------------------------------------------
Bool_t isActiveModule(Int_t mRunId, Int_t backleg, Int_t module)
{
    Bool_t isActive = kTRUE;
    if( 16139054 <= mRunId && mRunId <= 16149025 ){
        if(13 <= backleg && backleg <= 22) isActive = kFALSE;
        if(13 <= backleg && backleg <= 17 && module == 3) isActive = kTRUE;
    }
    else if( 16154003 <= mRunId && mRunId <= 16156034 ){
        if(backleg==15 || backleg==16 || backleg==17) isActive = kFALSE;
    }
    return isActive;
}

//-----------------------------------------------
void calPolarization(TLorentzVector iVec,TLorentzVector vec, THnSparse * hn, THnSparse *hnCS, double w, double pTheta, double pPhi )
{
    //iVec:positron TLorentzVector  vec:JPSI LorentzVector
    TLorentzVector mu(iVec);//positron
    TLorentzVector Proton1(0.,0.,100.,100.),Proton2(0.,0.,-100.,100.);//pp200GeV

    TVector3 XXHX,YYHX,ZZHX;//hilicity frame XYZaxis
    ZZHX = vec.Vect();
    YYHX = vec.Vect().Cross(Proton1.Vect());
    XXHX = YYHX.Cross(ZZHX);

    TVector3 jpsi = vec.BoostVector();

    mu.Boost(-1*jpsi);

    Float_t theta = mu.Angle(ZZHX);
    Float_t phi = TMath::ATan2((mu.Vect().Dot(YYHX.Unit())),(mu.Vect().Dot(XXHX.Unit())));
    if(phi<0) phi = -1.*phi;
    Float_t cosTheta = TMath::Cos(theta);
    if(isnan(cosTheta)) return;
    if(isnan(phi)) return;
    Double_t fill[]={vec.M(),vec.Pt(),cosTheta,phi};
    Double_t polar = weightPola(pTheta, pPhi, cosTheta, phi);
    hn->Fill(fill,w*polar);

    Proton1.Boost(-1*jpsi);
    Proton2.Boost(-1*jpsi);
    TVector3 XX,YY,ZZ;//Collins-Soper frame
    ZZ = Proton1.Vect()*(1/Proton1.Vect().Mag())-Proton2.Vect()*(1/Proton2.Vect().Mag());
    YY = Proton1.Vect().Cross(Proton2.Vect());
    XX = YY.Cross(ZZ);
    Float_t thetaCS = mu.Angle(ZZ);
    Float_t phiCS = TMath::ATan2((mu.Vect().Dot(YY.Unit())),(mu.Vect().Dot(XX.Unit())));
    if(phiCS<0) phiCS = -1.*phiCS;
    Float_t cosThetaCS = TMath::Cos(thetaCS);
    if(isnan(cosThetaCS)) return;
    if(isnan(phi)) return;
    Double_t fillCS[]={vec.M(),vec.Pt(),cosThetaCS,phiCS};
    Double_t polarCS = weightPola(pTheta, pPhi, cosThetaCS, phiCS);
    hnCS->Fill(fillCS,w*polarCS);
}

//-----------------------------------------------
void writeHistograms(char* outFile)
{
    char buf[1024];
    sprintf(buf,"%s.histo.root",outFile);
    cout<<"Writing histograms into "<<buf<<endl;
    TFile *mFile = new TFile(buf,"recreate");

    mFile->cd();
    hEvent->Write();
    hEtaVsPhi->Write();
    hHitMap->Write();
    hHitEtaVsPhi->Write();

    mhMcMPtCosPhi->Write();
    mhMcMPtCosPhiCS->Write();
    mhRcMPtCosPhi->Write();
    mhRcMPtCosPhiCS->Write();

    for(Int_t i=0; i<4; i++){
        hJpsiSpec[i][0]->Write();
        hJpsiSpec[i][1]->Write();
    }

    hNEmbedJpsi->Write();

}
