//--------------------------------------------------------
void pdfAction(TCanvas *c, TPDF *ps)
{
    ps->On();
    c->Update();
    c->cd();
    ps->NewPage();
    ps->Off();
}

//--------------------------------------------------------
TLatex* drawLatex(Double_t x, Double_t y, TString text, Int_t textFont, Double_t textSize, Int_t colorIndex)
{
    TLatex *latex = new TLatex(x,y,text.Data());
    latex->SetNDC();
    latex->SetTextFont(textFont);
    latex->SetTextSize(textSize);
    latex->SetTextColor(colorIndex);
    latex->Draw("same");
    return latex;
}

//--------------------------------------------------------
TLegend* drawLegend(TLegend* leg, Int_t textFont=132, Double_t textSize=0.04)
{
    leg->SetTextFont(textFont);
    //leg->SetTextAlign(11);
    //align = 10*HorizontalAlign + VerticalAlign (1=left adjusted, 2=centered, 3=right adjusted)
    leg->SetTextSize(textSize);
    leg->SetBorderSize(0);
    leg->SetFillColor(0);
    return leg;
}


//--------------------------------------------------------
TLine* drawLine(Double_t xlow,Double_t ylow, Double_t xup, Double_t yup, Double_t lineWidth, Int_t lineStyle,Int_t lineColor)
{
    TLine *l1 = new TLine(xlow,ylow,xup,yup);
    l1->SetLineWidth(lineWidth);
    l1->SetLineColor(lineColor);
    l1->SetLineStyle(lineStyle);
    l1->Draw("same");
    return l1;
}

//--------------------------------------------------------
const Double_t maxDiff = 1.e-10;

TH1D *rebHisto(TH1D *oldHisto, TString name, Int_t nBinsX, Double_t *BinX, TString NorX = "X")
{
    TH1D *newHisto = new TH1D(name.Data(),name.Data(),nBinsX,BinX);
    
    Int_t oldNBins = oldHisto->GetNbinsX();
    Int_t newNBins = newHisto->GetNbinsX();
    if(newNBins>oldNBins){
        cout<<" Failed to rebin \""<<oldHisto->GetName()<<"\"! --- NBinsX are more than original histogram!"<<endl;
        return NULL;
    }
    
    newHisto->GetXaxis()->SetTitle(oldHisto->GetXaxis()->GetTitle());
    
    Bool_t NorFlag = kTRUE;
    Double_t oldBinWidth = 1.e-3;
    for(Int_t i=0; i<oldNBins; i++){
        if(oldBinWidth>oldHisto->GetBinWidth(i+1))
            oldBinWidth = oldHisto->GetBinWidth(i+1);
    }
    for(Int_t i=0;i<newNBins;i++){
        Double_t newBinCenter = newHisto->GetBinCenter(i+1);
        Double_t newBinWidth = newHisto->GetBinWidth(i+1);
        Double_t newBinLow = newBinCenter-newBinWidth/2.+oldBinWidth/2.;
        Double_t newBinHi  = newBinCenter+newBinWidth/2.+oldBinWidth/2.;
        
        Int_t oldBinLow = oldHisto->FindBin(newBinLow);
        Int_t oldBinHi  = oldHisto->FindBin(newBinHi)-1;
        Int_t oldBinDiff = oldBinHi-oldBinLow+1;
        
        if(
           TMath::Abs(newBinCenter-newBinWidth/2.-oldHisto->GetXaxis()->GetBinLowEdge(oldBinLow)) > maxDiff
           || TMath::Abs(newBinCenter+newBinWidth/2.-oldHisto->GetXaxis()->GetBinUpEdge(oldBinHi)) > maxDiff
           ){
            cout<<" Failed to rebin \""<<oldHisto->GetName()<<"\"! --- X-axis bin boundary is different from original histogram!"<<endl;
            return NULL;
        }
        
        Double_t newBinContent  = oldHisto->Integral(oldBinLow,oldBinHi);
        Double_t newBinErr = 0.;
        for(Int_t j=oldBinLow;j<=oldBinHi;j++){
            newBinErr += pow(oldHisto->GetBinError(j),2);
        }
        newBinErr = sqrt(newBinErr);
        if(NorX.CompareTo("X")==0){
            newHisto->SetBinContent(i+1,newBinContent/oldBinDiff);
            newHisto->SetBinError(i+1,newBinErr/oldBinDiff);
        }
        else{
            newHisto->SetBinContent(i+1,newBinContent);
            newHisto->SetBinError(i+1,newBinErr);
            NorFlag = kFALSE;
        }
    }	
    
    if(!NorFlag){
        cout<<"The \""<<newHisto->GetName()<<"\" is not normalized!"<<endl;
        cout<<"If you want to normalize the \""<<newHisto->GetName()<<"\", the NorXY argument should be \"X\"or\"\" !"<<endl;
    }
    newHisto->SetMinimum(0);
    return newHisto;
}

//__________________________________________________
TH2D *rebHisto(TH2D *oldHisto, TString name, Int_t nBinsX, Double_t *BinX, Int_t nBinsY, Double_t *BinY, TString RebXY="XY", TString NorXY="N")
{
    TH2D *newHisto;
    if(RebXY.CompareTo("X")==0){//rebinX
        nBinsY = oldHisto->GetNbinsY();
        Double_t LowY = oldHisto->GetYaxis()->GetBinLowEdge(1);
        Double_t UpY = oldHisto->GetYaxis()->GetBinUpEdge(nBinsY);
        newHisto = new TH2D(name.Data(),name.Data(),nBinsX,BinX,nBinsY,LowY,UpY);
    }
    else if(RebXY.CompareTo("Y")==0 || RebXY.CompareTo("")==0){//rebinY
        nBinsX = oldHisto->GetNbinsX();
        Double_t LowX = oldHisto->GetXaxis()->GetBinLowEdge(1);
        Double_t UpX = oldHisto->GetXaxis()->GetBinUpEdge(nBinsX);
        newHisto = new TH2D(name.Data(),name.Data(),nBinsX,LowX,UpX,nBinsY,BinY);
    }
    else if(RebXY.CompareTo("XY")==0){//rebinXY
        newHisto = new TH2D(name.Data(),name.Data(),nBinsX,BinX,nBinsY,BinY);
    }
    else{
        cout<<"Failed to rebin \""<<oldHisto->GetName()<<"\" --- The RebXY parameter should be \"X\", \"Y\"or\"\", \"XY\"!"<<endl;
        return NULL;
    }
    
    newHisto->GetXaxis()->SetTitle(oldHisto->GetXaxis()->GetTitle());
    newHisto->GetYaxis()->SetTitle(oldHisto->GetYaxis()->GetTitle());
    
    Int_t oldNBinsX = oldHisto->GetNbinsX();
    Int_t oldNBinsY = oldHisto->GetNbinsY();
    Int_t newNBinsX = newHisto->GetNbinsX();
    Int_t newNBinsY = newHisto->GetNbinsY();
    if(newNBinsX>oldNBinsX){
        cout<<" Failed to rebin \""<<oldHisto->GetName()<<"\"! --- NBinsX are more than original histogram!"<<endl;
        return NULL;
    }
    if(newNBinsY>oldNBinsY){
        cout<<" Failed to rebin \""<<oldHisto->GetName()<<"\"! --- NBinsY are more than original histogram!"<<endl;
        return NULL;
    }
    
    Bool_t NorFlag = kTRUE;
    Double_t oldBinWidthX = 1.e-3;
    for(Int_t i=0; i<oldNBinsX; i++){
        if(oldBinWidthX>oldHisto->GetXaxis()->GetBinWidth(i+1))
            oldBinWidthX = oldHisto->GetXaxis()->GetBinWidth(i+1);
    }
    Double_t oldBinWidthY = 1.e-3;
    for(Int_t i=0; i<oldNBinsY; i++){
        if(oldBinWidthY>oldHisto->GetYaxis()->GetBinWidth(i+1))
            oldBinWidthY = oldHisto->GetYaxis()->GetBinWidth(i+1);
    }
    for(Int_t i=0;i<newNBinsX;i++){
        for(Int_t j=0;j<newNBinsY;j++){
            Double_t newBinCenterX = newHisto->GetXaxis()->GetBinCenter(i+1);
            Double_t newBinWidthX = newHisto->GetXaxis()->GetBinWidth(i+1);
            Double_t newBinLowX = newBinCenterX-newBinWidthX/2.+oldBinWidthX/2.;
            Double_t newBinHiX  = newBinCenterX+newBinWidthX/2.+oldBinWidthX/2.;
            Double_t newBinCenterY = newHisto->GetYaxis()->GetBinCenter(j+1);
            Double_t newBinWidthY = newHisto->GetYaxis()->GetBinWidth(j+1);
            Double_t newBinLowY = newBinCenterY-newBinWidthY/2.+oldBinWidthY/2.;
            Double_t newBinHiY  = newBinCenterY+newBinWidthY/2.+oldBinWidthY/2.;
            
            Int_t oldBinLowX = oldHisto->GetXaxis()->FindBin(newBinLowX);
            Int_t oldBinHiX  = oldHisto->GetXaxis()->FindBin(newBinHiX)-1;
            Int_t oldBinDiffX = oldBinHiX-oldBinLowX+1;
            Int_t oldBinLowY = oldHisto->GetYaxis()->FindBin(newBinLowY);
            Int_t oldBinHiY  = oldHisto->GetYaxis()->FindBin(newBinHiY)-1;
            Int_t oldBinDiffY = oldBinHiY-oldBinLowY+1;
            
            if(
               TMath::Abs(newBinCenterX-newBinWidthX/2.-oldHisto->GetXaxis()->GetBinLowEdge(oldBinLowX)) > maxDiff
               || TMath::Abs(newBinCenterX+newBinWidthX/2.-oldHisto->GetXaxis()->GetBinUpEdge(oldBinHiX)) > maxDiff
               ){
                cout<<" Failed to rebin \""<<oldHisto->GetName()<<"\"! --- X-axis bin boundary is different from original histogram!"<<endl;
                return NULL;
            }
            
            if(
               TMath::Abs(newBinCenterY-newBinWidthY/2.-oldHisto->GetYaxis()->GetBinLowEdge(oldBinLowY)) > maxDiff
               || TMath::Abs(newBinCenterY+newBinWidthY/2.-oldHisto->GetYaxis()->GetBinUpEdge(oldBinHiY)) > maxDiff
               ){
                cout<<" Failed to rebin \""<<oldHisto->GetName()<<"\"! --- Y-axis bin boundary is different from original histogram!"<<endl;
                return NULL;
            }
            
            Double_t newBinContent  = oldHisto->Integral(oldBinLowX,oldBinHiX,oldBinLowY,oldBinHiY);
            Double_t newBinErr = 0.;
            for(Int_t binx=oldBinLowX;binx<=oldBinHiX;binx++){
                for(Int_t biny=oldBinLowY;biny<=oldBinHiY;biny++){
                    newBinErr += pow(oldHisto->GetBinError(binx,biny),2.);
                }
            }
            newBinErr = sqrt(newBinErr);
            
            if(NorXY.CompareTo("X")==0){
                newHisto->SetBinContent(i+1,j+1,newBinContent/oldBinDiffX);
                newHisto->SetBinError(i+1,j+1,newBinErr/oldBinDiffX);
            }
            else if(NorXY.CompareTo("Y")==0 || NorXY.CompareTo("")==0){
                newHisto->SetBinContent(i+1,j+1,newBinContent/oldBinDiffY);
                newHisto->SetBinError(i+1,j+1,newBinErr/oldBinDiffY);
            }
            else if(NorXY.CompareTo("XY")==0){
                newHisto->SetBinContent(i+1,j+1,newBinContent/oldBinDiffX/oldBinDiffY);
                newHisto->SetBinError(i+1,j+1,newBinErr/oldBinDiffX/oldBinDiffY);
            }
            else{
                newHisto->SetBinContent(i+1,j+1,newBinContent);
                newHisto->SetBinError(i+1,j+1,newBinErr);
                NorFlag = kFALSE;
            }
        }
    }	
    
    if(!NorFlag){
        cout<<"The \""<<newHisto->GetName()<<"\" is not normalized!"<<endl;
        cout<<"If you want to normalize the \""<<newHisto->GetName()<<"\", the NorXY argument should be \"X\", \"Y\"or\"\", \"XY\"!"<<endl;
    }
    return newHisto;
}

//--------------------------------------------------------
TH1D* divideHisto(TH1D* h1, TH1D* h2, double a=1., double b=1., const char* p="B")
{
    TH1D* h = new TH1D(*h1);
    h->Divide(h1, h2, a, b, p);
    return h;
}

//--------------------------------------------------------
TH1D* reflectHisto(TH1D* oldHisto, TString name, TString title)
{
    Int_t oldNBins = oldHisto->GetNbinsX();
    if(oldNBins%2!=0){
        cout<<" Failed to rebin \""<<oldHisto->GetName()<<"\"! --- NBinsX are odd!"<<endl;
        return NULL;
    }
    double oldBinWidth = oldHisto->GetBinWidth(1);
    double newLowEdge = oldHisto->GetBinCenter(oldNBins/2) + oldBinWidth/2.;
    double newHighEdge = oldHisto->GetBinCenter(oldNBins) + oldBinWidth/2.;
    
    TH1D* newHisto = new TH1D(name.Data(), title.Data(), oldNBins/2, newLowEdge, newHighEdge);
    
    for(int i=1; i<=oldNBins/2; i++){
        double oldContent1 = oldHisto->GetBinContent(i);
        double oldError1 = oldHisto->GetBinError(i);
        double oldContent2 = oldHisto->GetBinContent(oldNBins+1-i);
        double oldError2 = oldHisto->GetBinError(oldNBins+1-i);
        
        double newContent = oldContent1 + oldContent2;
        double newError = sqrt( pow(oldError1,2) + pow(oldError2,2) );
        
        newHisto->SetBinContent(oldNBins/2 + 1 - i, newContent);
        newHisto->SetBinError(oldNBins/2 + 1 -i, newError);
    }
    
    newHisto->SetMinimum(0);
    return newHisto;
}

//--------------------------------------------------------
TH1D* normHisto(TH1D* h)
{
    int nBins = h->GetNbinsX();
    for(int i=1; i<=nBins; i++){
        double cont = h->GetBinContent(i);
        double contErr = h->GetBinError(i);
        double width = h->GetBinWidth(i);
        //cout<<"cont: "<<cont<<" contErr: "<<contErr<<" width: "<< width <<endl;
        h->SetBinContent(i, cont/width);
        h->SetBinError(i, contErr/width);
    }
    return h;
}

