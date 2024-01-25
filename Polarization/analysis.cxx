#include "headers.h"
#include "MINIEVENT.h"
#include "constant.h"

//-- declare function --//
void bookHistograms();
Bool_t passEvent(MINIEVENT* event);
Bool_t isTPCTrackOK(Double_t pt, Double_t eta, Int_t nHitsFit, Int_t nHitsMax, Int_t nHitsDedx, Double_t dca);
Bool_t checkMuCandidate(Double_t pt, Double_t nSigmaPion, Double_t dy, Double_t dz, Double_t dtof, Int_t matchFlag);
Bool_t checkMuDeltaY(Double_t pt, Double_t dy);
Bool_t checkMuDeltaZ(Double_t pt, Double_t dz);
Bool_t checkMuDeltaTof(Double_t dtof);
Bool_t isActiveModule(Int_t mRunId, Int_t backleg, Int_t module);
void calPolarization(TLorentzVector iVec,TLorentzVector vec, THnSparse * hn, THnSparse *hnCS);
void writeHistograms(char* outFile);

const Bool_t debug = false;

int main(int argc, char** argv)
{
    if(argc!=1&&argc!=3&&argc!=7) return -1;
    TString inFile="test.list";
    char outFile[1024];
    sprintf(outFile,"test");
    if(argc==3){
        inFile = argv[1];
        sprintf(outFile,"%s",argv[2]);
    }
    if(argc==7){
        inFile = argv[1];
        sprintf(outFile,"%s",argv[2]);
        mDcaCut = atof(argv[3]);
        mNHitsFitCut = atoi(argv[4]);
        mNHitsDedxCut = atoi(argv[5]);
        mMinNSigmaPiCut = atof(argv[6]);
    }

    //+---------------------------------+
    //| open files and add to the chain |
    //+---------------------------------+
    TChain *chain = new TChain("qaTree");

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
                cout<<"something wrong"<<endl;
            } else {
                cout<<"read in "<<ifile<<"th file: "<<filename<<endl;
                chain->Add(filename);
                ifile++;
            }
            delete ftmp;
        }
    }
    delete inputStream;

    bookHistograms();

    //+-------------+
    //| loop events |
    //+-------------+
    MINIEVENT *event = new MINIEVENT(chain);
    Int_t nEvts = chain->GetEntries();
    cout<<nEvts<<" events"<<endl;
    int perevt = nEvts/100;
    for(int ievt=0;ievt<nEvts;ievt++) {
        if(ievt%perevt ==0) cout<<"Looking at "<<((float)ievt/nEvts)*100<<" % .."<<endl;
        event->GetEntry(ievt);
        hEvent->Fill(0.5);
        if(!passEvent(event)) continue;
        Int_t runId = event->runId;
        Int_t npTracks = event->npTracks;
        for(int i=0;i<npTracks;i++){
            Int_t nHitsFit1 = event -> nHitsFit[i];
            Int_t nHitsPoss1 = event -> nHitsPoss[i];
            Int_t nHitsDedx1 = event -> nHitsDedx[i];
            Double_t nSigmaPion1 = event -> nSigmaPion[i];
            Double_t pt1 = event -> pt[i];
            Double_t eta1 = event -> eta[i];
            Double_t phi1 = event -> phi[i];
            Double_t dca1 = event -> dca[i];
            Int_t charge1 = event -> charge[i];
            Int_t matchFlag1 = event -> matchFlag[i];
            Int_t module1 = event -> module[i];
            Int_t backleg1 = event -> backleg[i];
            Double_t dy1 = event -> dy[i];
            Double_t dz1 = event -> dz[i];
            Double_t dtof1 = event -> dtof[i];
            if(!isTPCTrackOK(pt1,eta1,nHitsFit1,nHitsPoss1,nHitsDedx1,dca1)) continue;
            if(!isActiveModule(runId, backleg1, module1)) continue;
            if(!checkMuCandidate(pt1, nSigmaPion1, dy1, dz1, dtof1, matchFlag1)) continue;

            for(int j=i+1; j<npTracks;j++){
                Int_t nHitsFit2 = event -> nHitsFit[j];
                Int_t nHitsPoss2 = event -> nHitsPoss[j];
                Int_t nHitsDedx2 = event -> nHitsDedx[j];
                Double_t nSigmaPion2 = event -> nSigmaPion[j];
                Double_t pt2 = event -> pt[j];
                Double_t eta2 = event -> eta[j];
                Double_t phi2 = event -> phi[j];
                Double_t dca2 = event -> dca[j];
                Int_t charge2 = event -> charge[j];
                Int_t matchFlag2 = event -> matchFlag[j];
                Int_t module2 = event -> module[j];
                Int_t backleg2 = event -> backleg[j];
                Double_t dy2 = event -> dy[j];
                Double_t dz2 = event -> dz[j];
                Double_t dtof2 = event -> dtof[j];
                if(!isTPCTrackOK(pt2,eta2,nHitsFit2,nHitsPoss2,nHitsDedx2,dca2)) continue;
                if(!isActiveModule(runId, backleg2, module2)) continue;
                if(!checkMuCandidate(pt2, nSigmaPion2, dy2, dz2, dtof2, matchFlag2)) continue;

                muon1.SetPtEtaPhiM(pt1, eta1, phi1, Mmuon);
                muon2.SetPtEtaPhiM(pt2, eta2, phi2, Mmuon);
                Jpsi = muon1 + muon2;
                if(TMath::Abs(Jpsi.Rapidity())>mPairYCut) continue;

                if(charge1*charge2<0){
                    if(charge1>0) {
                        calPolarization(muon1,Jpsi,mhULMPtCosPhi, mhULMPtCosPhiCS);
                        hMuonPtEtaPhi->Fill(muon1.Pt(), muon1.Eta(), muon1.Phi());
                    }
                    else {
                        calPolarization(muon2,Jpsi,mhULMPtCosPhi, mhULMPtCosPhiCS);
                    }
                }
                else if(charge1*charge2>0) {
                    calPolarization(muon1,Jpsi,mhLSMPtCosPhi,mhLSMPtCosPhiCS);
                    hMuonPtEtaPhi->Fill(muon2.Pt(), muon2.Eta(), muon2.Phi());
                }
            }
        }
    }//event loop

    writeHistograms(outFile);
    delete chain;

    cout<<"end of program"<<endl;
    return 0;

}//main


//------------------------------------------------------
Bool_t passEvent(MINIEVENT* event)
{
    Double_t vx = event->tpcVx;
    Double_t vy = event->tpcVy;
    Double_t vz = event->tpcVz;
    Double_t vr = sqrt( vx*vx + vy*vy );

    if (!(fabs(vx)>0 || fabs(vy)>0 || fabs(vz)>0)) return false;
    hEvent->Fill(1.5);
    if( fabs(vr) > mVxyCut ) return false;
    hEvent->Fill(2.5);
    if( fabs(vz) > mVzCut ) return false;
    hEvent->Fill(3.5);
    return true;
}
//-------------------------------------------------------
Bool_t isTPCTrackOK(Double_t pt, Double_t eta, Int_t nHitsFit, Int_t nHitsMax, Int_t nHitsDedx, Double_t dca)
{
    if(
            pt >= mTpcMinPtCut  &&
            pt <= mTpcMaxPtCut  &&
            (float) nHitsFit/nHitsMax  > mNHitsFitMaxCut  &&
            nHitsFit  >= mNHitsFitCut  &&
            nHitsDedx >= mNHitsDedxCut &&
            fabs(dca) <= mDcaCut  &&
            TMath::Abs(eta)<mMuEtaCut
      )return true;
    else return false;
}
//-------------------------------------------------------
Bool_t checkMuCandidate ( Double_t pt, Double_t nSigmaPion, Double_t dy, Double_t dz, Double_t dtof, Int_t matchFlag = 0 )
{
    if(
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
    if( mMuMindTofCut <= dtof && dtof  <= mMuMaxdTofCut ) return true;
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
void calPolarization(TLorentzVector iVec,TLorentzVector vec, THnSparse * hn, THnSparse *hnCS)
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
    Double_t fill[]={vec.M(),vec.Pt(),cosTheta,phi};
    hn->Fill(fill);

    Proton1.Boost(-1*jpsi);
    Proton2.Boost(-1*jpsi);
    TVector3 XX,YY,ZZ;//Collins-Soper frame
    ZZ = Proton1.Vect()*(1/Proton1.Vect().Mag())-Proton2.Vect()*(1/Proton2.Vect().Mag());
    YY = Proton1.Vect().Cross(Proton2.Vect());
    XX = YY.Cross(ZZ);
    Float_t thetaCS = mu.Angle(ZZ);
    Float_t phiCS = TMath::ATan2((mu.Vect().Dot(YY.Unit())),(mu.Vect().Dot(XX.Unit())));
    if(phiCS<0) phiCS  = -1.* phiCS;
    Float_t cosThetaCS = TMath::Cos(thetaCS);
    Double_t fillCS[]={vec.M(),vec.Pt(),cosThetaCS,phiCS};
    hnCS->Fill(fillCS);
}

//-----------------------------------------------------
void bookHistograms()
{
    hEvent = new TH1D("hEvent","# of event",5,0,5);
    hEvent->GetXaxis()->SetBinLabel(1,"all event");
    hEvent ->GetXaxis()->SetBinLabel(2,"TPC Zero Vertex Cut");
    hEvent ->GetXaxis()->SetBinLabel(3,"TPC Vr Cut");
    hEvent ->GetXaxis()->SetBinLabel(4,"TPC Vz Cut");

    hEtaVsPhi = new TH2D("hEtaVsPhi","#varphi VS #eta ;#eta;#varphi",40,-1,1,48,-PI,PI);
    hHitMap = new TH2D("hHitMap","; module; backleg",7,0,7,32,0,32);
    hHitEtaVsPhi = new TH2D("hHitEtaVsPhi","#varphi VS #eta of hit with MTD;#eta;#varphi",40,-1,1,48,-PI,PI);

    const Int_t dim = 4;
    const Int_t nBins[dim]={220,400,10,15};
    const Double_t lowBins[dim]={2,0,-1, 0};
    const Double_t highBins[dim]={4.2,20,1,TMath::Pi()};
    mhULMPtCosPhi = new THnSparseF("mhULMPtCosPhi","HF UnlikeSignPair;m_{#mu#mu} (GeV/c^{2});p_{T} (GeV/c);cos#theta;#phi",dim,nBins,lowBins,highBins);
    mhULMPtCosPhi->Sumw2();
    mhLSMPtCosPhi = new THnSparseF("mhLSMPtCosPhi","HF UnlikeSignPair;m_{#mu#mu} (GeV/c^{2});p_{T} (GeV/c);cos#theta;#phi",dim,nBins,lowBins,highBins);
    mhLSMPtCosPhi->Sumw2();
    mhULMPtCosPhiCS = new THnSparseF("mhULMPtCosPhiCS","CS UnlikeSignPair;m_{#mu#mu} (GeV/c^{2});p_{T} (GeV/c);cos#theta;#phi",dim,nBins,lowBins,highBins);
    mhULMPtCosPhiCS->Sumw2();
    mhLSMPtCosPhiCS = new THnSparseF("mhLSMPtCosPhiCS","CS UnlikeSignPair;m_{#mu#mu} (GeV/c^{2});p_{T} (GeV/c);cos#theta;#phi",dim,nBins,lowBins,highBins);
    mhLSMPtCosPhiCS->Sumw2();

    hMuonPtEtaPhi = new TH3D("hMuonPtEtaPhi", "#mu^{+} track info ;p_{T} (GeV/c); #eta; #phi ", 400, 0, 20, 40, -1, 1, 96, -2*PI, 2*PI);
    hMuonPtEtaPhi->Sumw2();

}

//-------------------------------------------------------
void writeHistograms(char* outFile)
{
    TFile *mFile = new TFile(Form("%s.histo.root",outFile),"recreate");
    mFile->cd();

    hEvent->Write();
    hEtaVsPhi->Write();
    hHitMap->Write();
    hHitEtaVsPhi->Write();

    mhULMPtCosPhi->Write();
    mhLSMPtCosPhi->Write();
    mhULMPtCosPhiCS->Write();
    mhLSMPtCosPhiCS->Write();

    hMuonPtEtaPhi->Write();

}
