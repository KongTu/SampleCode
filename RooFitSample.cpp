TH1D* massks  = f1.Get("demo/masskshortnew0.2_0.4");
    RooRealVar x("x","mass",0.44,0.56);
    RooDataHist data("data","dataset",x,massks );
    RooPlot* xframe = x.frame(240);
    data.plotOn(xframe,Name("data"));
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
    sum->fitTo(data);
    sum->fitTo(data);
    sum->fitTo(data);
    sum->fitTo(data);
    sum->fitTo(data);
    sum->fitTo(data);
    
    sum->plotOn(xframe,Name("sum"));
    sum->plotOn(xframe,Components(poly),LineStyle(kDashed));
    cc1->cd(1);
    xframe->Draw();
    double chi2  = xframe->chiSquare("sum","data");
    double meanf  = mean->getVal();
    double meanfe  = mean->getError();
    double sigmaf1  = sigma1->getVal();
    double sigmaf2  = sigma2->getVal();
    double bkgf  = polysig->getVal();
    double sigf1  = sig1->getVal();
    double sigf2  = sig2->getVal();
    double sigwf1  = sigf1 /(sigf1 +sigf2 );
    double sigwf2  = sigf2 /(sigf1 +sigf2 );
    double c1 = a->getVal();
    double c2 = b->getVal();
    
    double sigmaf  = sqrt(sigmaf1 **2*sigwf1  + sigmaf2 **2*sigwf2 );
    double massmin  = meanf  - 2*sigmaf ;
    double massmax  = meanf  + 2*sigmaf ;
    
    TLine* l1 = new TLine(massmin ,0,massmin ,500);
    TLine* l2 = new TLine(massmax ,0,massmax ,500);
    l1.Draw("same");
    l2.Draw("same");
    
    int nmin  = massks ->GetXaxis()->FindBin(massmin );
    int nmax  = massks ->GetXaxis()->FindBin(massmax );
    int anmin  = massks ->GetXaxis()->FindBin(0.44);
    int anmax  = massks ->GetXaxis()->FindBin(0.56);
    
    double awyh1  = massks ->Integral(anmin ,nmin );
    double awyh2  = massks ->Integral(nmax ,anmax );
    double awyh  = awyh1  + awyh2 ;
    double totyh  = massks ->Integral(nmin ,nmax );
    
    x->setRange("cut",massmin ,massmax );
    RooAbsReal* ibkg  = poly->createIntegral(x,NormSet(x),Range("cut"));
    RooAbsReal* isig1  = gaus1->createIntegral(x,NormSet(x),Range("cut"));
    RooAbsReal* isig2  = gaus2->createIntegral(x,NormSet(x),Range("cut"));
    double ibkgf  = ibkg ->getVal();
    double bkgfe  = polysig->getError();
    double isig1f  = isig1 ->getVal();
    double isig2f  = isig2 ->getVal();
    
    double bkgy  = ibkgf *bkgf ;
    double bkgye  = ibkgf *bkgfe ;
    double sigy1  = isig1f *sigf1 ;
    double sigy2  = isig2f *sigf2 ;
    double sigy  = sigy1  + sigy2 ;
    double toty  = bkgy  + sigy ;
    
    double abkgy  = (1-ibkgf )*bkgf ;
    double asigy1  = (1-isig1f )*sigf1 ;
    double asigy2  = (1-isig2f )*sigf2 ;
    double asigy  = asigy1  + asigy2 ;
    double awy  = abkgy  + asigy ;
    
    double sigfrac  = sigy /toty ;
    double bkgfrac  = bkgy /toty ;
    
    double sigyh  = totyh  - bkgy ;
    double sigfrach  = sigyh /totyh ;
    double bkgfrach  = bkgy /totyh ;
    
    double signif  = sigyh /sqrt(bkgy );
