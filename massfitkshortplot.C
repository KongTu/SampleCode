#include <iostream>

#include "TF1.h"
#include "TH1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TAxis.h"
//#include "RooFit.h"
//#include "RooRealVar.h"
//#include "RooPlot.h"
//#include "RooDataHist.h"
//#include "RooGaussian.h"
//#include "RooPolynomial.h"
//#include "RooAddPdf.h"
#include "TString.h"

#include <vector>

void massfitkshortplot()
{
    gStyle->SetMarkerSize(0.8);
    TFile* f1 = TFile::Open("HighMultCorrelation220_260_ref_v14.root");
    using namespace RooFit;
    gStyle->SetOptTitle(kFALSE);

    TCanvas* cc1 = new TCanvas("cc1","cc1",1600,900);
    cc1->Divide(4,3);
    
    TH1D* massks[11];
    
    TLatex* tex = new TLatex();
    tex->SetNDC();
    
    char label_pT[11][200]={"0.2<p_{T}<0.4 GeV","0.4<p_{T}<0.6 GeV","0.6<p_{T}<0.8 GeV","0.8<p_{T}<1.0 GeV","1.0<p_{T}<1.4 GeV","1.4<p_{T}<1.8 GeV","1.8<p_{T}<2.2 GeV","2.2<p_{T}<2.8 GeV","2.8<p_{T}<3.6 GeV","3.6<p_{T}<4.6 GeV","4.6<p_{T}<6.0 GeV"};
    char label_energy[200]={"CMS pPb"};
    char label_n[200]={"220 #leq N^{offline}_{trk} < 260"};
    
    tex->SetTextSize(tex->GetTextSize()*1.1);
    
    for(int i=0;i<11;i++)
    {
        massks[i] = (TH1D*)f1->Get(Form("demo/masskshort_pt%d",i));
        RooRealVar x("x","mass",0.43,0.565);
        RooDataHist data("data","dataset",x,massks[i]);
        RooPlot* xframe = x.frame(270);
        xframe->GetXaxis()->SetTitle("invariant mass (GeV)");
        xframe->GetYaxis()->SetTitle("Candidates / 0.0005 GeV");
        xframe->GetXaxis()->CenterTitle(1);
        xframe->GetYaxis()->CenterTitle(1);
        xframe->GetXaxis()->SetTitleSize(xframe->GetXaxis()->GetTitleSize()*1.4);
        xframe->GetYaxis()->SetTitleSize(xframe->GetYaxis()->GetTitleSize()*1.4);
        data->plotOn(xframe,Name("data"));
        RooRealVar mean("mean","mean",0.50,0.49,0.51);
        RooRealVar sigma1("sigma1","sigma1",0.01,0.001,0.04);
        RooRealVar sigma2("sigma2","sigma2",0.01,0.001,0.04);
        RooRealVar sig1("sig1","signal1",10,0,10000000);
        RooRealVar sig2("sig2","signal2",10,0,10000000);
        RooRealVar a("a","a",0,-100000,100000);
        RooRealVar b("b","b",0,-100000,100000);
        RooRealVar cp("cp","cp",0,-100000,100000);
        RooRealVar d("d","d",0,-100000,100000);
        RooGaussian gaus1("gaus1","gaus1",x,mean,sigma1);
        RooGaussian gaus2("gaus2","gaus2",x,mean,sigma2);
        RooPolynomial poly("poly","poly",x,RooArgList(a,b,cp,d));
        RooRealVar polysig("polysig","polysig",10,0,10000000);
        RooAddPdf sum("sum","sum",RooArgList(gaus1,gaus2,poly),RooArgList(sig1,sig2,polysig));
        
        x->setRange("cut",0.44,0.56);

        sum->fitTo(data,Range("cut"));
        sum->fitTo(data,Range("cut"));
        sum->fitTo(data,Range("cut"));
        sum->fitTo(data,Range("cut"));
        sum->fitTo(data,Range("cut"));
        sum->fitTo(data,Range("cut"));
        
        
        sum->plotOn(xframe,Name("sum"),NormRange("cut"),LineWidth(0.5),LineColor(kRed));
        sum->plotOn(xframe,Components(poly),NormRange("cut"),LineStyle(kDashed),LineWidth(0.5),LineColor(kRed));
        cc1->cd(i+1);
        xframe->Draw();
        tex->DrawLatex(0.60,0.80,label_pT[i]);
    }
    
    cc1->cd(1);
    tex->DrawLatex(0.20,0.85,label_energy);
    tex->SetTextSize(tex->GetTextSize()*1.0);
    tex->DrawLatex(0.20,0.75,label_n);

    cc1->Print("masspeak_ks.pdf");
}
