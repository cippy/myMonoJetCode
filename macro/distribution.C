//ROOT header files
#include <TROOT.h>
#include <TAttFill.h>
#include <TAxis.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TF1.h>
#include <TFile.h>
#include <TFitResult.h>
#include <TGraphErrors.h>
#include <THStack.h>
#include <TH1.h>
#include <TH1D.h>
#include <TKey.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TMatrixDSym.h>
#include <TMultiGraph.h>
#include <TPad.h>
#include <TPaveStats.h>
#include <TPaveText.h>
#include <TStyle.h>
#include <TVirtualFitter.h>

//C or C++ header files
#include <stdio.h>
#include <stdlib.h>
#include <cstdlib> //as stdlib.h
#include <cstdio>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <iomanip> 
#include <algorithm>  // to use the "reverse" function to reverse the order in the array

#include <Rtypes.h> // to use kColor

using namespace std;

void myAddOverflowInLastBin(TH1D *h) {  // this is the function in functionsForAnalysis.cc, but I copied here to simplify compilation

  // to avoid problems regarding memory leak for not deleting htemp in previous function, I sum directly the content of overflow bin in last bin

  Int_t lastBinNumber = h->GetNbinsX();
  Int_t overflowBinNumber = 1 + lastBinNumber;
  Double_t lastBinContent = h->GetBinContent(lastBinNumber);
  Double_t overflowBinContent = h->GetBinContent(overflowBinNumber);
  Double_t lastBinError = h->GetBinError(lastBinNumber);
  Double_t overflowBinError = h->GetBinError(overflowBinNumber);

  // add content of overflow bin in last bin and set error as square root of sum of error squares (with the assumption that they are uncorrelated)
  h->SetBinContent(lastBinNumber, lastBinContent + overflowBinContent);
  h->SetBinError(lastBinNumber, sqrt(lastBinError * lastBinError + overflowBinError * overflowBinError));

}


void myRebinHisto(TH1D *h, const Int_t rebinFactor = 1) {

  if (rebinFactor != 1) {
    h->Rebin(rebinFactor);
    if ( (h->GetNbinsX() % rebinFactor) != 0) myAddOverflowInLastBin(h);
  }

}

void setSampleName(const Int_t signalRegion0_controlRegion1, vector<string> &sampleName, vector<string> &MC_TexLabel, const Int_t z0_w1 = 0, const Int_t mumu0_ee1 = 0) {

  // for Control region, mumu0_ee1 says if we use muon (0) or electron (1)

  if (signalRegion0_controlRegion1 == 0) {

    //sampleName.push_back("GJets");
    sampleName.push_back("DYJetsToLL");
    sampleName.push_back("QCD");
    sampleName.push_back("Diboson");
    sampleName.push_back("Top");
    sampleName.push_back("WJetsToLNu");
    sampleName.push_back("ZJetsToNuNu");
       
    //MC_TexLabel.push_back("#gamma + jets");
    MC_TexLabel.push_back("Z(ll)+jets");
    MC_TexLabel.push_back("QCD");
    MC_TexLabel.push_back("Diboson");
    MC_TexLabel.push_back("t#bar{t},single t");
    MC_TexLabel.push_back("W(l#nu)+jets");
    MC_TexLabel.push_back("Z(#nu#nu)+jets");
    
  } else {

    if(z0_w1 == 0) {

      sampleName.push_back("ZJetsToNuNu");
      //sampleName.push_back("GJets");
      sampleName.push_back("QCD");
      sampleName.push_back("WJetsToLNu");
      sampleName.push_back("Diboson");
      sampleName.push_back("Top");
      sampleName.push_back("DYJetsToLL");

      MC_TexLabel.push_back("Z(#nu#nu)+jets");
      //MC_TexLabel.push_back("#gamma + jets");
      MC_TexLabel.push_back("QCD");
      MC_TexLabel.push_back("W(l#nu)+jets");
      MC_TexLabel.push_back("Diboson");
      MC_TexLabel.push_back("t#bar{t},single t");

      if (mumu0_ee1 == 0) MC_TexLabel.push_back("Z(#mu#mu)+jets");
      else MC_TexLabel.push_back("Z(ee)+jets");

    } else if(z0_w1 == 1) {

      sampleName.push_back("ZJetsToNuNu");
      sampleName.push_back("DYJetsToLL");
      //sampleName.push_back("GJets");
      sampleName.push_back("QCD");      
      sampleName.push_back("Diboson");
      sampleName.push_back("Top");   
      sampleName.push_back("WJetsToLNu");

      MC_TexLabel.push_back("Z(#nu#nu)+jets");
      MC_TexLabel.push_back("Z(ll)+jets");
      //MC_TexLabel.push_back("#gamma + jets");
      MC_TexLabel.push_back("QCD");     
      MC_TexLabel.push_back("Diboson");
      MC_TexLabel.push_back("t#bar{t},single t");

      if (mumu0_ee1 == 0) MC_TexLabel.push_back("W(#mu#nu)+jets");
      else MC_TexLabel.push_back("W(e#nu)+jets");

    }

  }

}



void setHistColor(vector<Int_t> &histColor, const Int_t nSamples) {

  Int_t colorList[] = {kCyan, kViolet, kBlue, kRed, /*kOrange+1, */kYellow, kGreen};  // the first color is for the main object. This array may contain more values than nSamples

  for (Int_t i = 0; i < nSamples; i++) {   // now color are assigned in reverse order (the main contribution is the last object in the sample array)

    histColor.push_back(colorList[i]);

  }

  std::reverse(histColor.begin(), histColor.end());   // reverse order in the array

}



void setDistribution(const Int_t mumu0_ee1, const string var, string &hvarName, string &xAxisName) {

  if ( !(strcmp("metBin",var.c_str())) ) { 

    hvarName = "HYieldsMetBin";
    xAxisName = "recoil [GeV]"; //"#slash{E}_{T} [GeV]";

  } else if ( !(strcmp("met",var.c_str())) ) {

    hvarName = "HmetNoLepDistribution"; 
    //xAxisName = "#slash{E}_{T} [GeV]";
    xAxisName = "recoil [GeV]";

  } else if ( !(strcmp("ht",var.c_str())) ) {

    hvarName = "HhtDistribution"; 
    xAxisName = "#slash{H}_{T} [GeV]";

  } else if ( !(strcmp("nvtx",var.c_str())) ) {

    hvarName = "HvtxDistribution"; 
    xAxisName = "N(vertices)";

  } else if ( !(strcmp("njets",var.c_str())) ) {

    hvarName = "HnjetsDistribution"; 
    xAxisName = "N(jets)";

  } else if ( !(strcmp("j1pt",var.c_str())) ) {

    hvarName = "Hjet1ptDistribution"; 
    xAxisName = "leading jet p_{T} [GeV]";

  } else if ( !(strcmp("j2pt",var.c_str())) ) {

    hvarName = "Hjet2ptDistribution"; 
    xAxisName = "second jet p_{T} [GeV]";

  } else if ( !(strcmp("j1eta",var.c_str())) ) {
    
    hvarName = "Hjet1etaDistribution"; 
    xAxisName = "leading jet #eta";

  } else if ( !(strcmp("prunedMass",var.c_str())) ) {
    
    hvarName = "HprunedMassDistribution"; 
    xAxisName = "pruned mass [GeV]";

  } else if ( !(strcmp("tau2OverTau1",var.c_str())) ) {
    
    hvarName = "Htau2OverTau1Distribution"; 
    xAxisName = "#tau_{2} / #tau_{1}";

  } else if ( !(strcmp("invMass",var.c_str())) ) {

    hvarName = "HinvMass";
    if (mumu0_ee1 == 0) xAxisName = "m(#mu#mu) [GeV]";
    else if (mumu0_ee1 == 1) xAxisName = "m(ee) [GeV]";

  } else if ( !(strcmp("zpt",var.c_str())) ) {

    hvarName = "HzptDistribution";
    xAxisName = "p_{T}(Z) [GeV]";

  } else if ( !(strcmp("gammapt",var.c_str())) ) {

    hvarName = "HgammaptDistribution";
    xAxisName = "p_{T}(#gamma) [GeV]";

  }else if ( !(strcmp("Mt",var.c_str())) ) {

    hvarName = "HtransverseMass";
    xAxisName = "m_{T} [GeV]";

  } else if ( !(strcmp("lep1pt",var.c_str())) ) {

    hvarName = "Hlep1ptDistribution"; 
    if (mumu0_ee1 == 0) xAxisName = "leading muon p_{T} [GeV]";
    else if (mumu0_ee1 == 1) xAxisName = "leading electron p_{T} [GeV]";

  } else if ( !(strcmp("lep2pt",var.c_str())) ) {

    hvarName = "Hlep2ptDistribution"; 
    if (mumu0_ee1 == 0) xAxisName = "second muon p_{T} [GeV]";
    else if (mumu0_ee1 == 1) xAxisName = "second electron p_{T} [GeV]";

  } else {
    
    cout << " Variable " << var << " not available. End of programme." << endl;
    exit(EXIT_FAILURE);

  }

}

// ================================================




void distribution(const string folderNameWithRootFiles = "", 
		  const Int_t signalRegion0_controlRegion1 = 0, 
		  const Int_t z0_w1 = 0,
		  const Int_t mu0_e1 = 0,
		  const Int_t data0_noData1 = 0, 
		  const Int_t monoJ0_monoV1 = 0,
		  const string var = "met", 
		  const Int_t yAxisLog_flag = 0, 
		  const Int_t MCpoissonUncertainty_flag = 0, 
		  const Double_t xAxisMin = 0, 
		  const Double_t xAxisMax = -1, 
		  const Double_t yAxisMin = 0, 
		  const Double_t yAxisMax = -1, 
		  const Int_t binDensity_flag = 0, 
		  const Int_t MCnormalizedToData_flag = 0,
		  const Int_t rebinFactor = 1,
		  const Double_t minRatioY = 0.0,
		  const Double_t maxRatioY = 2.0) 
{

  // if signalRegion0_controlRegion1 == 0 (default), will do met distribution in the monojet signal region, else it will do the control region

  // z0_w1 is for CR to choose between Z or W to leptons

  // mu0_e1 is for lepton flavour in CS (not used in SR)

  // data0_noData1 is to use or not a data file to compared with MC

  // yAxisLog_flag is to choose whether or not to use log scale in Y axis (default is 0, that is, normal scale)
  // xAxisMin and xAxisMax are the ranges for x Axis. Default values are used if xAxisMin > xAxisMax (otherwise user values are used)

  TH1::SetDefaultSumw2();            //all the following histograms will automatically call TH1::Sumw2() 

  gROOT->SetStyle("Plain");  // to have white legend (on my pc it's already white, but in tier2 it appears grey)
  gStyle->SetHistTopMargin(0.);
  gStyle->SetFillColor(10);

  string filenameExtension = ".root";
  // string fileDirectoryPath = "spring15_25ns_rootfiles/";
  string fileDirectoryPath = "/cmshome/ciprianim/CMSSW721/output/monojet/" + folderNameWithRootFiles + "/";

  string plotDirectoryPath = fileDirectoryPath;
  // string plotDirectoryPath = "/cmshome/ciprianim/CMSSW721/pdfsFromAnalysis/plots/monojet/met_distribution/";
  //string plotDirectoryPath = "./distributions/";
  string plotFileExtension = ".pdf"; //could be png
  string suffix;

  if (mu0_e1 == 0) {
    if (z0_w1 == 0) suffix = "mumu";
    if (z0_w1 == 1) suffix = "munu";
  } else if (mu0_e1 == 1) {
    if (z0_w1 == 0) suffix = "ee";
    if (z0_w1 == 1) suffix = "enu";
  } else {
    cout << "Error: mu0_e1 must be 0 or 1. End of programme." << endl;
    exit(EXIT_FAILURE);
  }

  TH1D* hvar = NULL;   // to get histogram from file
  string hvarName;          // name of histogram to take from file
  string xAxisName;        // name of X axis when plotting distribution. It is a tex string (with ROOT standard), e.g. "#slash{E}_{T} [GeV]" for MET

  // ===== TO BE MODIFIED =====

  // hvarName = "HmetNoLepDistribution";
  // xAxisName = "#slash{E}_{T} [GeV]";

  setDistribution(mu0_e1, var, hvarName, xAxisName);
  if (monoJ0_monoV1 == 1) hvarName += "_monoV";
    
  // =====================

  vector<TH1D*> hMC;
  TH1D* hdata = NULL;

  vector<string> sampleName;
  vector<string> MC_TexLabel;
  if (data0_noData1 == 1) setSampleName(signalRegion0_controlRegion1, sampleName, MC_TexLabel, z0_w1, mu0_e1);
  else setSampleName(signalRegion0_controlRegion1, sampleName, MC_TexLabel, z0_w1, mu0_e1);

  string data_TexLabel = "data";

  Int_t nFiles = (Int_t) sampleName.size();

  vector<Int_t> histColor;
  setHistColor(histColor, nFiles);

  string filenameBase;
  string canvasName;
  
  if (signalRegion0_controlRegion1 == 0) {
    
    canvasName = var + "_monojetSR";
    filenameBase = "monojet_SR_";    
    
  } else {
    
    if(z0_w1 == 0) {
      
      canvasName = var + "_z" + suffix + "jetsCR";
      if (mu0_e1 == 0) filenameBase = "zmumujets_CR_";
      else if (mu0_e1 == 1) filenameBase = "zeejets_CR_";
       
    } else if(z0_w1 == 1) {

      canvasName = var + "_w" + suffix + "jetsCR";
      if (mu0_e1 == 0) filenameBase = "wmunujets_CR_";
      else if (mu0_e1 == 1) filenameBase = "wenujets_CR_";
      
    }

  }

  if (monoJ0_monoV1 == 1) canvasName += "_monoV";

  vector<string> MCfileName;
  for (Int_t i = 0; i < nFiles; i++) {
    MCfileName.push_back(fileDirectoryPath + filenameBase + sampleName[i] + filenameExtension);
  }

  for(Int_t i = 0; i < nFiles; i++) {

    //cout<<"fileName : "<<MCfileName[i]<<endl;

    TFile* f = TFile::Open(MCfileName[i].c_str(),"READ");
    if (!f || !f->IsOpen()) {
      cout<<"*******************************"<<endl;
      cout<<"Error opening file \""<<MCfileName[i]<<"\".\nApplication will be terminated."<<endl;
      cout<<"*******************************"<<endl;
      exit(EXIT_FAILURE);
    }

    //cout << "check 1 " << endl;    

    hvar = (TH1D*)f->Get(hvarName.c_str());
    if (!hvar) {
      cout << "Error: histogram not found in file ' " << MCfileName[i] << "'. End of programme." << endl;
      exit(EXIT_FAILURE);
    }
    hMC.push_back( (TH1D*)hvar->Clone() );

  } 

  // ==== FILE NAME WILL HAVE TO BE MODIFIED, NOW IT IS JUST A TEST =====

  string datafileName = fileDirectoryPath;

  if (data0_noData1 == 0) {

    if (signalRegion0_controlRegion1 == 0) {

      datafileName += "monojet_SR_DATA.root";

    } else {

      if(z0_w1 == 0) {

       if (mu0_e1 == 0) datafileName += "zmumujets_CR_DATA.root";
       else if (mu0_e1 == 1) datafileName += "zeejets_CR_DATA.root";

     } else if(z0_w1 == 1) {

       if (mu0_e1 == 0) datafileName += "wmunujets_CR_DATA.root";
       else if (mu0_e1 == 1) datafileName += "wenujets_CR_DATA.root";

     }

      

    }

  }

  // ==== opening data file ======

  if (data0_noData1 == 0) {

    //cout<<"fileName : "<<datafileName<<endl;

    TFile* f = TFile::Open(datafileName.c_str(),"READ");
    if (!f || !f->IsOpen()) {
      cout<<"*******************************"<<endl;
      cout<<"Error opening file \""<<datafileName<<"\".\nApplication will be terminated."<<endl;
      cout<<"*******************************"<<endl;
      exit(EXIT_FAILURE);
    }

    hvar = (TH1D*)f->Get(hvarName.c_str());

    if (!hvar) {
      cout << "Error: histogram not found in file ' " << datafileName << "'. End of programme." << endl;
      exit(EXIT_FAILURE);
    }
    hdata = (TH1D*)hvar->Clone();

  }

  // ===============================

  THStack* hMCstack = new THStack("hMCstack","");
  Double_t stackNorm = 0.0;

  for (Int_t j = 0; j < nFiles; j++) {

    for (Int_t i = 1; i <= hMC[j]->GetNbinsX(); i++) {

      if (MCpoissonUncertainty_flag == 1) {

	hMC[j]->SetBinError(i,sqrt(hMC[j]->GetBinContent(i)));
	
      }

    }

    hMC[j]->SetFillColor(histColor[j]);
    stackNorm += hMC[j]->Integral();

  }

  // if (data0_noData1 == 0) {

  //   cout << "==================================================" << endl;
  //   cout << "Before any scaling" << endl;
  //   cout << "Events in data: " << hdata->Integral() << endl;
  //   cout << "Events in MC   : " << stackNorm << endl;
  //   cout << "==================================================" << endl;
  
  // }

  // loop again on MC histograms to scale them and then fill the thstack

  for (Int_t j = 0; j < nFiles; j++) {

    if (data0_noData1 == 0 && MCnormalizedToData_flag != 0) {

      Double_t dataNorm = hdata->Integral();

      if (binDensity_flag != 0) hMC[j]->Scale(dataNorm/stackNorm,"width");
      else hMC[j]->Scale(dataNorm/stackNorm);

    } else if (binDensity_flag != 0) hMC[j]->Scale(1.0,"width");  // option width divides by bin width and manages the correct error setting

    myRebinHisto(hMC[j], rebinFactor);
    hMCstack->Add(hMC[j]);

  }


  if (data0_noData1 == 0) {

    myRebinHisto(hdata, rebinFactor);

    if (binDensity_flag != 0) hdata->Scale(1.0,"width");

    // if (MCnormalizedToData_flag != 0 || binDensity_flag != 0) {
    //   cout << "==================================================" << endl;
    //   cout << "After scaling" << endl;
    //   if (binDensity_flag != 0) {
    // 	cout << "Events in data: " << hdata->Integral("width") << endl;
    // 	cout << "Events in MC   : " << ((TH1D*) hMCstack->GetStack()->Last())->Integral("width") << endl;
    //   } else {
    // 	cout << "Events in data: " << hdata->Integral() << endl;
    // 	cout << "Events in MC   : " << ((TH1D*) hMCstack->GetStack()->Last())->Integral() << endl;
    //   }
    //   cout << "==================================================" << endl;
    // }

  }


  // now here we go with the canvas

  TH1D * ratioplot = NULL; // will use it for the ratio plots
  //TH1D *h1overMCstack = NULL; //will be used to create ratioplot

  TPad *subpad_1 = NULL;  // will use it to access specific subpad in canvas
  TPad *subpad_2 = NULL; 

  TCanvas *c;
  if (data0_noData1 == 0) c = new TCanvas(canvasName.c_str(), (var + " distribution").c_str(), 700, 700);
  else c = new TCanvas(canvasName.c_str(), (var + " distribution").c_str());
  TLegend *leg = new TLegend(0.7,0.6,0.99,0.94);  

  // if there are data, split canvas to draw the dta/MC ratio plot

  if (data0_noData1 == 0) {

    subpad_1 = new TPad("pad_1","",0.0,0.28,1.0,1.0);
    if (yAxisLog_flag) subpad_1->SetLogy();
    //subpad_1->SetBottomMargin(0);
    subpad_2 = new TPad("pad_2","",0.0,0.0,1.0,0.32);
    subpad_2->SetGridy();
    //subpad_2->SetTopMargin(0);
    subpad_2->SetBottomMargin(0.3);
    subpad_1->Draw();
    subpad_2->Draw();

    subpad_1->cd();

  } else if (yAxisLog_flag) c->SetLogy();

  
  hMCstack->Draw("HIST");
  //if (yAxisMin > 0) hMCstack->SetMinimum(yAxisMin);
  TH1D* stackCopy = (TH1D*)hMCstack->GetStack()->Last();
  // TH1D* stackCopy = (TH1D*)(((TH1D*)hMCstack->GetStack()->Last())->DrawCopy("E2 SAME"));
  // stackCopy->SetFillColor(kBlack);
  // stackCopy->SetFillStyle(3144);

  if (yAxisMin < yAxisMax) {  // default option [0,-1] implies default axis
 
    hMCstack->GetYaxis()->SetRangeUser(yAxisMin,yAxisMax);
    if (data0_noData1 == 0) subpad_1->Update();  // to be done after Draw() to access pad parameters such as default axis range
    else c->Update(); 

  } else if (yAxisMin > yAxisMax && yAxisMax == -1) {  //if only lower bound is given (and the upper is left as -1) use default upper bound

    if (data0_noData1 == 0) {
      hMCstack->GetYaxis()->SetRangeUser(yAxisMin,subpad_1->GetY2());
      subpad_1->Update();  // to be done after Draw() to access pad parameters such as default axis range      
    } else {
      hMCstack->GetYaxis()->SetRangeUser(yAxisMin,c->GetY2());
      c->Update();  
    }

  }

  // if (data0_noData1 == 0) {

  //   if (yAxisMin > 0) hMCstack->GetYaxis()->SetRangeUser(yAxisMin, subpad_1->GetY2());
  //   else hMCstack->GetYaxis()->SetRangeUser(yAxisMin, c->GetY2());

  // }

  //hMCstack->SetMaximum(4000.0);

  if (data0_noData1 == 1) {    //  when using data ( == 0) the x axis will not have labels (they will only be below in the ratio plot)
    hMCstack->GetXaxis()->SetTitle(xAxisName.c_str());
    hMCstack->GetXaxis()->SetTitleSize(0.06);
    hMCstack->GetXaxis()->SetTitleOffset(0.6);
  }

  Double_t stackXmax = 0.0; // will be used to save the maximum value of x Axis for MC stack

  if (xAxisMin < xAxisMax) {
    if (data0_noData1 == 0) subpad_1->Update();  // to be done after Draw() to access pad parameters such as default axis range
    else c->Update();  
    hMCstack->GetXaxis()->SetRangeUser(xAxisMin,xAxisMax);
  } else if (xAxisMin > xAxisMax && xAxisMax == -1) {  //if only lower bound is given (and the upper is left as -1) use default upper bound

    stackXmax = hMCstack->GetXaxis()->GetXmax();
    //cout << "stackXmax: " << stackXmax << endl;

    if (data0_noData1 == 0) {
      subpad_1->Update();  // to be done after Draw() to access pad parameters such as default axis range  
      //hMCstack->GetXaxis()->SetRangeUser(xAxisMin,subpad_1->GetX2());
      //cout << "subpad_1->GetX2(): " << subpad_1->GetX2() << endl;
    } else {
      c->Update();  
      //hMCstack->GetXaxis()->SetRangeUser(xAxisMin,c->GetX2());
    }

    hMCstack->GetXaxis()->SetRangeUser(xAxisMin,stackXmax);   // N.B. SetRangeUser cannot set axis ranges outside the original coordinates: if they are 5-10, I can set them 6-9, 5-8 ecc... but if I do 5-12, the range will be set to 5-10

  }

  if (binDensity_flag == 0) hMCstack->GetYaxis()->SetTitle("events");
  else hMCstack->GetYaxis()->SetTitle("events / GeV");
  hMCstack->GetYaxis()->SetTitleSize(0.06);
  hMCstack->GetYaxis()->SetTitleOffset(0.8);
  hMCstack->GetYaxis()->CenterTitle();
  for (Int_t j = (nFiles-1); j >= 0; j--) {
    leg->AddEntry(hMC[j],Form("%s",MC_TexLabel[j].c_str()),"lf");
  }

  if (data0_noData1 == 0) {    
    hdata->SetMarkerStyle(8); // large dot
    hdata->Draw("EX0 SAME"); //X0 doesn't draw x error
    leg->AddEntry(hdata,Form("%s",data_TexLabel.c_str()),"p");
  }

  gStyle->SetStatStyle(0);
  leg->Draw(); 
  leg->SetMargin(0.3); 
  leg->SetBorderSize(0);
  //leg->SetFillStyle(0);  // transparent legend

  if (data0_noData1 == 0) { // if using data, draw the ratio plot

    subpad_2->cd();
    ratioplot = new TH1D(*stackCopy);  // to have ratioplot with the same x axis range as stack, I copy it from stackcopy created above, then I substitute bin content with hdata/stackcopy
    ratioplot->Divide(hdata,stackCopy);
    ratioplot->SetStats(0);
    ratioplot->SetTitle("");
    ratioplot->GetXaxis()->SetLabelSize(0.10);
    ratioplot->GetXaxis()->SetTitle(xAxisName.c_str());
    ratioplot->GetXaxis()->SetTitleSize(0.14);
    ratioplot->GetXaxis()->SetTitleOffset(0.8);
    ratioplot->GetYaxis()->SetLabelSize(0.10);
    ratioplot->GetYaxis()->SetTitle("data / MC");
    ratioplot->GetYaxis()->SetTitleSize(0.15);
    ratioplot->GetYaxis()->SetTitleOffset(0.3);
    ratioplot->GetYaxis()->CenterTitle();
    ratioplot->GetYaxis()->SetRangeUser(minRatioY,maxRatioY);
    ratioplot->GetYaxis()->SetNdivisions(011);
    ratioplot->SetMarkerStyle(8);  //medium dot
    if (xAxisMin < xAxisMax) ratioplot->GetXaxis()->SetRangeUser(xAxisMin,xAxisMax); 
    else if (xAxisMin > xAxisMax && xAxisMax == -1) {
      Double_t ratioplotXmax = ratioplot->GetXaxis()->GetXmax();
      //cout << "ratioplotXmax: " << ratioplotXmax << endl;
      subpad_2->Update();
      ratioplot->GetXaxis()->SetRangeUser(xAxisMin,ratioplotXmax); 
      //ratioplot->GetXaxis()->SetRangeUser(xAxisMin,stackXmax); //ok, it works
    }
    ratioplot->DrawCopy("E");

  }

  //c->SaveAs( (plotDirectoryPath + c->GetName() + plotFileExtension).c_str() );
  c->SaveAs( (plotDirectoryPath + c->GetName() + ".pdf").c_str() );
  c->SaveAs( (plotDirectoryPath + c->GetName() + ".png").c_str() );

}

//==========================================================================
// useless, doing this in analyzer
// void fillSystematic(TFile *f, vector<Double_t> &systVec, const string& norm, const string& up, const string& down) {

//   TH1D *hnorm = NULL;
//   TH1D *hup = NULL;
//   TH1D *hdown = NULL;

//   hnorm = (TH1D*)f->Get(norm.c_str());
//   //cout << h->GetName() << endl;
//   if (!hnorm) {
//     cout << "Error in computeSystematicError(): histogram not found in file. End of programme." << endl;
//     exit(EXIT_FAILURE);
//   }
//   hup = (TH1D*)f->Get(up.c_str());
//   //cout << h->GetName() << endl;
//   if (!hup) {
//     cout << "Error in computeSystematicError(): histogram not found in file. End of programme." << endl;
//     exit(EXIT_FAILURE);
//   }
//   hdown = (TH1D*)f->Get(down.c_str());
//   //cout << h->GetName() << endl;
//   if (!hdown) {
//     cout << "Error in computeSystematicError(): histogram not found in file. End of programme." << endl;
//     exit(EXIT_FAILURE);
//   }
  
//   if (TH1::CheckConsistency(hnorm,hup) && TH1::CheckConsistency(hnorm,hdown)) {

//     Int_t nbins = hnorm->GetNbinsX();
//     Double_t maxdiff = 0.0;
//     Double_t updiff = 0.0;
//     Double_t downdiff = 0.0;

//     for (Int_t i = 1; i <= nbins; i++) {

//       updiff = fabs(hnorm->GetBinContent(i) - hup->GetBinContent(i));
//       downdiff = fabs(hnorm->GetBinContent(i) - hdown->GetBinContent(i)); 
//       maxdiff = (updiff > downdiff) ? updiff : downdiff; 
//       systVec.push_back(maxdiff);

//     }

//   } else {
//     cout << "Error in fillSystematic(): TH1::CheckConsistency() returned false. End of programme." << endl;
//     exit(EXIT_FAILURE);
//   }  

// }

//==========================================================================
// useless, doing this in analyzer
// void computeSystematicError(TFile *f, vector<Double_t> &qcdRenScale, vector<Double_t> &qcdFacScale, vector<Double_t> &qcdPdf, vector<Double_t> &ewk) {

//   // here the file is read inside fillSystematic and histograms are used to compute systematics bin by bin for each scale factor
    
//   fillSystematic(f, qcdRenScale, "HYieldsMetBin", "HYieldsMetBin_qcdRenScaleUp", "HYieldsMetBin_qcdRenScaleDown");
//   fillSystematic(f, qcdFacScale, "HYieldsMetBin", "HYieldsMetBin_qcdFacScaleUp", "HYieldsMetBin_qcdFacScaleDown");
//   fillSystematic(f, qcdPdf, "HYieldsMetBin", "HYieldsMetBin_qcdPdfUp", "HYieldsMetBin_qcdPdfDown");
//   fillSystematic(f, ewk, "HYieldsMetBin", "HYieldsMetBin_ewkUp", "HYieldsMetBin_ewkDown");

// }

void makeTransferFactor(const string folderNameWithRootFilesSR = "",
			const string folderNameWithRootFilesCR = "",
			const Int_t bkgToEstimate_z0_w1 = 0,
			const Int_t monoJ0_monoV1 = 0,
			//const Int_t signalRegion0_controlRegion1 = 0, 
			//const Int_t z0_w1 = 0,
			//const Int_t mu0_e1 = 0,
			//const Int_t data0_noData1 = 0, 
			//const string var = "met", 
			//const Int_t yAxisLog_flag = 0, 
			//const Int_t MCpoissonUncertainty_flag = 0, 
			const Double_t xAxisMin = 0, 
			const Double_t xAxisMax = -1, 
			const Double_t yAxisMin = 0, 
			const Double_t yAxisMax = -1, 
			//const Int_t binDensity_flag = 0, 
			//const Int_t MCnormalizedToData_flag = 0,
			const Int_t numDenCorrelation = 0,
			const Int_t rebinFactor = 1)
{

  // numDenCorrelation is referred to systematic uncertainties between numerator and denominator in the ratio and for now it is an integer with 0 or 1 (0 for no correlation, 1 for 100% correlation). Roughly speaking, no correlation means propagating the whole systematic uncertainties in ratio (actually the total uncertainties, systematic and statistical added in quadrature), while total correlation would mean that in the ratio the uncertainties cancel out (only the systematics, because statistical uncertainties are supposed to be (totally) uncorrelated).

  // It's not trivial how to propagate the correlation on uncertainties, while I get to it, this variable will becaome a double in [0,1] 

  // monoJ0_monoV1 is used to decide to make transfer factors with monojet (default) or monoV events. To use monojet, set flag to 1;

  // the following flag is used to manage Z(vv)/W(lv) scale factor computation where both sample are taken from signal region (folderNameWithRootFilesCR == folderNameWithRootFilesSR)
  Int_t zvv_wlv_flag = 0;
  if (folderNameWithRootFilesCR == folderNameWithRootFilesSR) zvv_wlv_flag = 1;

  TH1::SetDefaultSumw2();            //all the following histograms will automatically call TH1::Sumw2() 

  gROOT->SetStyle("Plain");  // to have white legend (on my pc it's already white, but in tier2 it appears grey)
  gStyle->SetFillColor(10);

  string filenameExtension = ".root";
  // string fileDirectoryPath = "spring15_25ns_rootfiles/";
  string fileDirectoryPathSR = "/cmshome/ciprianim/CMSSW721/output/monojet/" + folderNameWithRootFilesSR + "/";
  string fileDirectoryPathCR = "/cmshome/ciprianim/CMSSW721/output/monojet/" + folderNameWithRootFilesCR + "/";
  string plotDirName;

  vector<string> suffixSR;   // to build plot name (will be e.g. znunu_zmumu_trasnferFactor.pdf)
  vector<string> texSuffixSR;  // to build Y axis name of plot with tex format
  suffixSR.push_back("znunu"); texSuffixSR.push_back("Z(#nu#nu)");
  suffixSR.push_back("wlnu");   texSuffixSR.push_back("W(l#nu)");

  vector<string> suffixCR;   // to build plot name (will be e.g. znunu_zmumu_trasnferFactor.pdf)
  vector<string> texSuffixCR;  // to build Y axis name of plot with tex format
  suffixCR.push_back("zmumu"); texSuffixCR.push_back("Z(#mu#mu)");
  suffixCR.push_back("zee");   texSuffixCR.push_back("Z(ee)");
  suffixCR.push_back("wmunu"); texSuffixCR.push_back("W(#mu#nu)");
  suffixCR.push_back("wenu");  texSuffixCR.push_back("W(#e#nu)");

  Int_t indexRightSuffixSR = -1;
  Int_t indexRightSuffixCR = -1;

  string canvasName = "";
  string texYTitle = "";

  string bkgSampleSR = "";

  if (bkgToEstimate_z0_w1 == 0) {
    indexRightSuffixSR = 0;
    bkgSampleSR = "ZJetsToNuNu";
  } else if (bkgToEstimate_z0_w1 == 1) {
    indexRightSuffixSR = 1;
    bkgSampleSR = "WJetsToLNu";
  }

  for (UInt_t i = 0; i < suffixCR.size(); i++) {
    size_t foundSuffix;
    foundSuffix = fileDirectoryPathCR.find(suffixCR.at(i));
    if (foundSuffix != string::npos) indexRightSuffixCR = i;
  }

  if (zvv_wlv_flag) {
    plotDirName = "wlnu";
    texYTitle = "Z(#nu#nu) / W(l#nu)";
    canvasName = "znunu_wlnu_transferFactor";
  } else {
    plotDirName = suffixCR.at(indexRightSuffixCR); 
    texYTitle = texSuffixSR.at(indexRightSuffixSR) + " / " + texSuffixCR.at(indexRightSuffixCR);
    canvasName =  suffixSR.at(indexRightSuffixSR) + "_" + suffixCR.at(indexRightSuffixCR) + "_transferFactor";
  }

  if (monoJ0_monoV1 == 1) canvasName += "_monoV"; 

  string calibEleSuffix = "_CalibEle";  // used to select proper directory where to save plot, depending on input file name
  size_t useCalibEle = fileDirectoryPathCR.find(calibEleSuffix);
  if (useCalibEle != string::npos) plotDirName += calibEleSuffix;
  
  string plotDirectoryPath = fileDirectoryPathSR;
  //  if (monoJ0_monoV1 == 1) plotDirectoryPath = "/cmshome/ciprianim/CMSSW721/output/monojet/SR_CR_ratio/monoV/" + plotDirName + "/";
  //else plotDirectoryPath = "/cmshome/ciprianim/CMSSW721/output/monojet/SR_CR_ratio/monojet/" + plotDirName + "/";
  string plotFileExtension = ".pdf";

  string filenameSR = fileDirectoryPathSR + "monojet_SR_" + bkgSampleSR + filenameExtension;

  string filenameCR = fileDirectoryPathCR;
  if (zvv_wlv_flag) {
    filenameCR += "monojet_SR_WJetsToLNu" + filenameExtension;
  } else {
    filenameCR += suffixCR.at(indexRightSuffixCR) + "jets_CR_"; 
    if (plotDirName.find("z") != string::npos) filenameCR += ("DYJetsToLL" + filenameExtension);
    else if (plotDirName.find("w") != string::npos) filenameCR += ("WJetsToLNu" + filenameExtension); 
  }

  TH1D* hvar = NULL;   // to get histogram from file
  TH1D* hsyst = NULL;
  //string hvarName;          // name of histogram to take from file
  //string xAxisName;        // name of X axis when plotting distribution. It is a tex string (with ROOT standard), e.g. "#slash{E}_{T} [GeV]" for MET
  string hvarName = "HYieldsMetBin";
  string hsystName = "HSyst_total";
  if (monoJ0_monoV1 == 1) {
    hvarName += "_monoV";
    hsystName += "_monoV";
  }
  string xAxisName = "recoil [GeV]";

  TH1D* hSR = NULL;
  TH1D* hCR = NULL;
  TH1D* hSRsyst = NULL;
  TH1D* hCRsyst = NULL;

  TFile* f = TFile::Open(filenameSR.c_str(),"READ");
  if (!f || !f->IsOpen()) {
    cout<<"*******************************"<<endl;
    cout<<"Error opening file \""<<filenameSR<<"\".\nApplication will be terminated."<<endl;
    cout<<"*******************************"<<endl;
    exit(EXIT_FAILURE);
  } else {

    cout << "Opening -> " << filenameSR << endl;
    hvar = (TH1D*)f->Get(hvarName.c_str());
    hsyst = (TH1D*)f->Get(hsystName.c_str());
    //cout << hvar->GetName() << endl;
    if (!hvar || !hsyst) {
      cout << "Error: histogram not found in file ' " << filenameSR << "'. End of programme." << endl;
      exit(EXIT_FAILURE);
    }
    hSR = (TH1D*)hvar->Clone();
    hSRsyst = (TH1D*)hsyst->Clone();

    //f->Close();  //if I close the file I lose the histogram and the clone (but cannot understand why the clone should be gone)

  }

  TFile* f2 = TFile::Open(filenameCR.c_str(),"READ");
  if (!f2 || !f2->IsOpen()) {
    cout<<"*******************************"<<endl;
    cout<<"Error opening file \""<<filenameCR<<"\".\nApplication will be terminated."<<endl;
    cout<<"*******************************"<<endl;
    exit(EXIT_FAILURE);
  } else {

    cout << "Opening -> " << filenameCR << endl;
    hvar = (TH1D*)f2->Get(hvarName.c_str());
    hsyst = (TH1D*)f2->Get(hsystName.c_str());
    //cout << hvar->GetName() << endl;
    if (!hvar || !hsyst) {
      cout << "Error: histogram not found in file ' " << filenameCR << "'. End of programme." << endl;
      exit(EXIT_FAILURE);
    }
    hCR = (TH1D*)hvar->Clone();
    hCRsyst = (TH1D*)hsyst->Clone();

    //f2->Close();  //if I close the file I lose the histogram and the clone (but cannot understand why the clone should be gone)
    
  }

  TH1D* hTransferFactor = NULL;
  TH1D* hTransferFactor_totUnc = NULL;
 
  hTransferFactor = (TH1D*)hSR->Clone("hTransferFactor"); 
  hTransferFactor->Divide(hCR);  // now the uncertainty is only statistical

  // now that hSR and hCR have been used with inly statistical uncertainty, compute their total uncertainty and use it in ratio

  for (Int_t i = 1; i <= (hSR->GetNbinsX()+1); i++) {

    //computing total uncertainty summing stat and syst in quadrature (histogram with systematic has uncertainty as bin content)
    Double_t totUnc = 0.0;
    // ---> handling correlation, to be revised  <----
    if (numDenCorrelation == 0) {
      totUnc = hSR->GetBinError(i) * hSR->GetBinError(i) + hSRsyst->GetBinContent(i) * hSRsyst->GetBinContent(i); 
      hSR->SetBinError(i, sqrt(totUnc));
      totUnc = hCR->GetBinError(i) * hCR->GetBinError(i) + hCRsyst->GetBinContent(i) * hCRsyst->GetBinContent(i);
      hCR->SetBinError(i, sqrt(totUnc));
    } // else if (numDenCorrelation == 1) {   // in this case just coonsider the statistical uncertainty (default)
    //   // if total correlation holds, systematic will cancel out in ratio
    //   totUnc = hSR->GetBinError(i) * hSR->GetBinError(i); 
    //   hSR->SetBinError(i, sqrt(totUnc));
    //   totUnc = hCR->GetBinError(i) * hCR->GetBinError(i);
    //   hCR->SetBinError(i, sqrt(totUnc));
    // }

  }

  // now there is a systematic uncertainty related to the (charged) lepton ID. Tis is 2% for each lepton, thus 2% for W ad 4% for Z in control regions only.
  // this uncertainty can be summed in quadrature to the others
  // this is applied to the denominator only if it comes from a control region (for Zvv/Wlv both numerator and denominator are taken from the signal region)
  // this uncertainty is set after implementing all the other uncertainties deriving from change in QCD or EWK scale (done just above)
  if (!zvv_wlv_flag) {

    Double_t lepIDuncertainty = 0.0; // 2% for 1 lep CR, 4% for 2 lep CR, which is 0.02 * binContent or 0.04 * binContent
    if (plotDirName.find("z") != string::npos) lepIDuncertainty = 0.04 * 0.04;       // computing square of the factor to use it directly below
    else if (plotDirName.find("w") != string::npos) lepIDuncertainty = 0.02 * 0.02;  // computing square of the factor to use it directly below

    for (Int_t i = 1; i <= (hCR->GetNbinsX()+1); i++) {

      //computing total uncertainty summing stat and syst in quadrature (histogram with systematic has uncertainty as bin content)
      Double_t totUnc = hCR->GetBinError(i) * hCR->GetBinError(i) + lepIDuncertainty * hCR->GetBinContent(i) * hCR->GetBinContent(i); 
      // note that here the bin error used in totUnc is the sum in quadrature of statistical and other systematic uncertainties   
      hCR->SetBinError(i, sqrt(totUnc));

    }

  }

  hTransferFactor_totUnc = (TH1D*)hSR->Clone("hTransferFactor_totUnc"); 
  hTransferFactor_totUnc->Divide(hCR);  // now the uncertainty is the total one

  // not used for metBin histograms
  // myRebinHisto(hTransferFactor_totUnc, rebinFactor);
  // myRebinHisto(hTransferFactor, rebinFactor);

  //cout << "CHECK" << endl;

  TCanvas *c = new TCanvas(canvasName.c_str(),"");
  TLegend *leg = new TLegend(0.7,0.7,0.9,0.89); 
  //cout << "CHECK YYY" << endl;
  hTransferFactor_totUnc->Draw("HIST");
  TH1D* tfCopy = (TH1D*)(hTransferFactor_totUnc->DrawCopy("E2 SAME"));
  tfCopy->SetFillColor(kCyan);
  
  hTransferFactor->SetMarkerStyle(8); // large dot
  hTransferFactor->Draw("E1 SAME"); //EX0 doesn't draw x error
  leg->AddEntry(tfCopy,"tot. unc.","lf");
  leg->AddEntry(hTransferFactor,"stat. unc.","p");
  
  //cout << "CHECK ZZZ" << endl;

  if (yAxisMin < yAxisMax) {  // default option [0,-1] implies default axis

    c->Update();  
    hTransferFactor_totUnc->GetYaxis()->SetRangeUser(yAxisMin,yAxisMax);

  } // else if (yAxisMin > yAxisMax && yAxisMax == -1) {  //if only lower bound is given (and the upper is left as -1) use default upper bound

  //   cout << "CHECK CCC" << endl;
  //   c->Update();  
  //   hTransferFactor_totUnc->GetYaxis()->SetRangeUser(yAxisMin,c->GetY2());
    
  //   cout << "CHECK BBB" << endl;

  // }

  //cout << "CHECK AAA" << endl;

  hTransferFactor_totUnc->SetTitle("");  // no title on the plot
  hTransferFactor_totUnc->SetStats(0); //no stat box
  hTransferFactor_totUnc->GetYaxis()->SetTitle((texYTitle + " transfer factor").c_str());
  hTransferFactor_totUnc->GetYaxis()->SetTitleSize(0.06);
  hTransferFactor_totUnc->GetYaxis()->SetLabelSize(0.045);
  hTransferFactor_totUnc->GetYaxis()->SetTitleOffset(0.8);
  hTransferFactor_totUnc->GetYaxis()->CenterTitle();

  hTransferFactor_totUnc->GetXaxis()->SetTitle(xAxisName.c_str());
  hTransferFactor_totUnc->GetXaxis()->SetTitleSize(0.06);
  hTransferFactor_totUnc->GetXaxis()->SetLabelSize(0.045);
  hTransferFactor_totUnc->GetXaxis()->SetTitleOffset(0.8);
  
  //cout << "CHECK IV" << endl;

  if (xAxisMin < xAxisMax) {

    c->Update();  
    hTransferFactor_totUnc->GetXaxis()->SetRangeUser(xAxisMin,xAxisMax);

  } else if (xAxisMin > xAxisMax && xAxisMax == -1) {  //if only lower bound is given (and the upper is left as -1) use default upper bound
    
    c->Update();   
    hTransferFactor_totUnc->GetXaxis()->SetRangeUser(xAxisMin,hTransferFactor_totUnc->GetXaxis()->GetXmax());   // N.B. SetRangeUser cannot set axis ranges outside the original coordinates: if they are 5-10, I can set them 6-9, 5-8 ecc... but if I do 5-12, the range will be set to 5-10

  }

  gStyle->SetStatStyle(0);
  leg->Draw(); 
  leg->SetMargin(0.3); 
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);

  //c->SaveAs( (plotDirectoryPath + c->GetName() + plotFileExtension).c_str() );
  c->SaveAs( (plotDirectoryPath + c->GetName() + ".pdf").c_str() );
  c->SaveAs( (plotDirectoryPath + c->GetName() + ".png").c_str() );
		

}



Int_t main(int argc, char* argv[]) {

  // distribution() function template
  /*distribution(const string folderNameWithRootFiles = "", 
	       const Int_t signalRegion0_controlRegion1 = 0, 
	       const Int_t z0_w1 = 0,
	       const Int_t mu0_e1 = 0,
	       const Int_t data0_noData1 = 0, 
	       const Int_t monoJ0_monoV1 = 0,
	       const string var = "met", 
	       const Int_t yAxisLog_flag = 0, 
	       const Int_t MCpoissonUncertainty_flag = 0, 
	       const Double_t xAxisMin = 0, 
	       const Double_t xAxisMax = -1, 
	       const Double_t yAxisMin = 0, 
	       const Double_t yAxisMax = -1, 
	       const Int_t binDensity_flag = 0, 
	       const Int_t MCnormalizedToData_flag = 0,
	       const Int_t rebinFactor = 1,
	       const Double_t minRatioY = 0.0,
	       const Double_t maxRatioY = 2.0)*/

  Int_t j0v1 = 0;  // passing from outside whether to make plots with monojet events (0) or monoV (1). A value different from 1 is considered as monojet

  if (argc > 1) {
    j0v1 = atoi(argv[1]);
  }


  // signal region
  distribution("SignalRegion_spring15_25ns_2p32fb_newSF",0,0,0,0,j0v1,"metBin",1,0,200,-1,0,-1,1,0);
  distribution("SignalRegion_spring15_25ns_2p32fb_newSF",0,0,0,0,j0v1,"met",1,0,200,-1,0,-1,0,0,2);
  distribution("SignalRegion_spring15_25ns_2p32fb_newSF",0,0,0,0,j0v1,"ht",1,0,0,-1,0,-1,0,0,2);
  distribution("SignalRegion_spring15_25ns_2p32fb_newSF",0,0,0,0,j0v1,"j1pt",1,0,100,-1,0,-1,0,0,2);
  distribution("SignalRegion_spring15_25ns_2p32fb_newSF",0,0,0,0,j0v1,"j1eta",0,0,0,-1,0,-1,0,0);
  distribution("SignalRegion_spring15_25ns_2p32fb_newSF",0,0,0,0,j0v1,"nvtx",0,0,0,-1,0,-1,0,0);
  distribution("SignalRegion_spring15_25ns_2p32fb_newSF",0,0,0,0,j0v1,"njets",1,0,0,-1,0,-1,0,0);
  if (j0v1 != 1) distribution("SignalRegion_spring15_25ns_2p32fb_newSF",0,0,0,0,j0v1,"j2pt",1,0,0,-1,0,-1,0,0,2);
  else {
    distribution("SignalRegion_spring15_25ns_2p32fb_newSF",0,0,0,0,j0v1,"prunedMass",0,0,0,-1,0,-1,0,0);
    distribution("SignalRegion_spring15_25ns_2p32fb_newSF",0,0,0,0,j0v1,"tau2OverTau1",0,0,0,-1,0,-1,0,0);
  }

  // Z->mumu control region
  distribution("ControlRegion_zmumu_spring15_25ns_2p32fb_newSF",1,0,0,0,j0v1,"metBin",1,0,200,-1,0,-1,1,0);
  distribution("ControlRegion_zmumu_spring15_25ns_2p32fb_newSF",1,0,0,0,j0v1,"met",1,0,200,-1,0,-1,0,0,2);
  distribution("ControlRegion_zmumu_spring15_25ns_2p32fb_newSF",1,0,0,0,j0v1,"ht",1,0,0,-1,0,-1,0,0,2);
  distribution("ControlRegion_zmumu_spring15_25ns_2p32fb_newSF",1,0,0,0,j0v1,"j1pt",1,0,100,-1,0,-1,0,0,2);
  distribution("ControlRegion_zmumu_spring15_25ns_2p32fb_newSF",1,0,0,0,j0v1,"j1eta",0,0,0,-1,0,-1,0,0);
  distribution("ControlRegion_zmumu_spring15_25ns_2p32fb_newSF",1,0,0,0,j0v1,"nvtx",0,0,0,-1,0,-1,0,0);
  distribution("ControlRegion_zmumu_spring15_25ns_2p32fb_newSF",1,0,0,0,j0v1,"njets",1,0,0,-1,0,-1,0,0);
  distribution("ControlRegion_zmumu_spring15_25ns_2p32fb_newSF",1,0,0,0,j0v1,"zpt",1,0,0,-1,0,-1,0,0,2);
  distribution("ControlRegion_zmumu_spring15_25ns_2p32fb_newSF",1,0,0,0,j0v1,"invMass",0,0,0,-1,0,-1,0,0);
  distribution("ControlRegion_zmumu_spring15_25ns_2p32fb_newSF",1,0,0,0,j0v1,"lep1pt",1,0,0,-1,0,-1,0,0,4);
  distribution("ControlRegion_zmumu_spring15_25ns_2p32fb_newSF",1,0,0,0,j0v1,"lep2pt",1,0,0,-1,0,-1,0,0,4);
  if (j0v1 != 1) distribution("ControlRegion_zmumu_spring15_25ns_2p32fb_newSF",1,0,0,0,j0v1,"j2pt",1,0,0,-1,0,-1,0,0,2);
  else {
    distribution("ControlRegion_zmumu_spring15_25ns_2p32fb_newSF",1,0,0,0,j0v1,"prunedMass",0,0,0,-1,0,-1,0,0);
    distribution("ControlRegion_zmumu_spring15_25ns_2p32fb_newSF",1,0,0,0,j0v1,"tau2OverTau1",0,0,0,-1,0,-1,0,0);
  }

  // // Z->ee control region
  // distribution("ControlRegion_zee_spring15_25ns_2p32fb_newSF",1,0,1,0,j0v1,"metBin",1,0,200,-1,0,-1,1,0);
  // distribution("ControlRegion_zee_spring15_25ns_2p32fb_newSF",1,0,1,0,j0v1,"met",1,0,200,-1,0,-1,0,0,2);
  // distribution("ControlRegion_zee_spring15_25ns_2p32fb_newSF",1,0,1,0,j0v1,"ht",1,0,0,-1,0,-1,0,0,2);
  // distribution("ControlRegion_zee_spring15_25ns_2p32fb_newSF",1,0,1,0,j0v1,"j1pt",1,0,100,-1,0,-1,0,0,2);
  // distribution("ControlRegion_zee_spring15_25ns_2p32fb_newSF",1,0,1,0,j0v1,"j1eta",0,0,0,-1,0,-1,0,0);
  // distribution("ControlRegion_zee_spring15_25ns_2p32fb_newSF",1,0,1,0,j0v1,"nvtx",0,0,0,-1,0,-1,0,0);
  // distribution("ControlRegion_zee_spring15_25ns_2p32fb_newSF",1,0,1,0,j0v1,"njets",1,0,0,-1,0,-1,0,0);
  // distribution("ControlRegion_zee_spring15_25ns_2p32fb_newSF",1,0,1,0,j0v1,"zpt",1,0,0,-1,0,-1,0,0,2);
  // distribution("ControlRegion_zee_spring15_25ns_2p32fb_newSF",1,0,1,0,j0v1,"invMass",0,0,0,-1,0,-1,0,0);
  // distribution("ControlRegion_zee_spring15_25ns_2p32fb_newSF",1,0,1,0,j0v1,"lep1pt",1,0,0,-1,0,-1,0,0,4);
  // distribution("ControlRegion_zee_spring15_25ns_2p32fb_newSF",1,0,1,0,j0v1,"lep2pt",1,0,0,-1,0,-1,0,0,4);
  // if (j0v1 != 1) distribution("ControlRegion_zee_spring15_25ns_2p32fb_newSF",1,0,1,0,j0v1,"j2pt",1,0,0,-1,0,-1,0,0,2);
  // else {
  //   distribution("ControlRegion_zee_spring15_25ns_2p32fb_newSF",1,0,1,0,j0v1,"prunedMass",0,0,0,-1,0,-1,0,0);
  //   distribution("ControlRegion_zee_spring15_25ns_2p32fb_newSF",1,0,1,0,j0v1,"tau2OverTau1",0,0,0,-1,0,-1,0,0);
  // }

  // // Z->ee control region (calibEle)
  // distribution("ControlRegion_zee_spring15_25ns_2p32fb_newSF_CalibEle",1,0,1,0,j0v1,"metBin",1,0,200,-1,0,-1,1,0);
  // distribution("ControlRegion_zee_spring15_25ns_2p32fb_newSF_CalibEle",1,0,1,0,j0v1,"met",1,0,200,-1,0,-1,0,0,2);
  // distribution("ControlRegion_zee_spring15_25ns_2p32fb_newSF_CalibEle",1,0,1,0,j0v1,"j1pt",1,0,100,-1,0,-1,0,0,2);
  // distribution("ControlRegion_zee_spring15_25ns_2p32fb_newSF_CalibEle",1,0,1,0,j0v1,"j1eta",0,0,0,-1,0,-1,0,0);
  // distribution("ControlRegion_zee_spring15_25ns_2p32fb_newSF_CalibEle",1,0,1,0,j0v1,"nvtx",0,0,0,-1,0,-1,0,0);
  // distribution("ControlRegion_zee_spring15_25ns_2p32fb_newSF_CalibEle",1,0,1,0,j0v1,"njets",1,0,0,-1,0,-1,0,0);
  // distribution("ControlRegion_zee_spring15_25ns_2p32fb_newSF_CalibEle",1,0,1,0,j0v1,"zpt",1,0,0,-1,0,-1,0,0,2);
  // distribution("ControlRegion_zee_spring15_25ns_2p32fb_newSF_CalibEle",1,0,1,0,j0v1,"invMass",0,0,0,-1,0,-1,0,0);
  // distribution("ControlRegion_zee_spring15_25ns_2p32fb_newSF_CalibEle",1,0,1,0,j0v1,"lep1pt",1,0,0,-1,0,-1,0,0,4);
  // distribution("ControlRegion_zee_spring15_25ns_2p32fb_newSF_CalibEle",1,0,1,0,j0v1,"lep2pt",1,0,0,-1,0,-1,0,0,4);
  // if (j0v1 != 1) distribution("ControlRegion_zee_spring15_25ns_2p32fb_newSF_CalibEle",1,0,1,0,j0v1,"j2pt",1,0,0,-1,0,-1,0,0,2);
  // else {
  //   distribution("ControlRegion_zee_spring15_25ns_2p32fb_newSF_CalibEle",1,0,1,0,j0v1,"prunedMass",0,0,0,-1,0,-1,0,0);
  //   distribution("ControlRegion_zee_spring15_25ns_2p32fb_newSF_CalibEle",1,0,1,0,j0v1,"tau2OverTau1",0,0,0,-1,0,-1,0,0);
  // }

  // W->munu control region
  distribution("ControlRegion_wmunu_spring15_25ns_2p32fb_newSF",1,1,0,0,j0v1,"metBin",1,0,200,-1,0,-1,1,0);
  distribution("ControlRegion_wmunu_spring15_25ns_2p32fb_newSF",1,1,0,0,j0v1,"met",1,0,200,-1,0,-1,0,0,2);
  distribution("ControlRegion_wmunu_spring15_25ns_2p32fb_newSF",1,1,0,0,j0v1,"ht",1,0,0,-1,0,-1,0,0,2);
  distribution("ControlRegion_wmunu_spring15_25ns_2p32fb_newSF",1,1,0,0,j0v1,"j1pt",1,0,100,-1,0,-1,0,0,2);
  distribution("ControlRegion_wmunu_spring15_25ns_2p32fb_newSF",1,1,0,0,j0v1,"j1eta",0,0,0,-1,0,-1,0,0);
  distribution("ControlRegion_wmunu_spring15_25ns_2p32fb_newSF",1,1,0,0,j0v1,"nvtx",0,0,0,-1,0,-1,0,0);
  distribution("ControlRegion_wmunu_spring15_25ns_2p32fb_newSF",1,1,0,0,j0v1,"njets",1,0,0,-1,0,-1,0,0);
  distribution("ControlRegion_wmunu_spring15_25ns_2p32fb_newSF",1,1,0,0,j0v1,"Mt",1,0,0,-1,0,-1,0,0);
  distribution("ControlRegion_wmunu_spring15_25ns_2p32fb_newSF",1,1,0,0,j0v1,"lep1pt",1,0,0,-1,0,-1,0,0,4);
  if (j0v1 != 1) distribution("ControlRegion_wmunu_spring15_25ns_2p32fb_newSF",1,1,0,0,j0v1,"j2pt",1,0,0,-1,0,-1,0,0,2);
  else {
    distribution("ControlRegion_wmunu_spring15_25ns_2p32fb_newSF",1,1,0,0,j0v1,"prunedMass",0,0,0,-1,0,-1,0,0);
    distribution("ControlRegion_wmunu_spring15_25ns_2p32fb_newSF",1,1,0,0,j0v1,"tau2OverTau1",0,0,0,-1,0,-1,0,0);
  }

  // // W->enu control region
  // distribution("ControlRegion_wenu_spring15_25ns_2p32fb_newSF",1,1,1,0,j0v1,"metBin",1,0,200,-1,0,-1,1,0);
  // distribution("ControlRegion_wenu_spring15_25ns_2p32fb_newSF",1,1,1,0,j0v1,"met",1,0,200,-1,0,-1,0,0,2);
  // distribution("ControlRegion_wenu_spring15_25ns_2p32fb_newSF",1,1,1,0,j0v1,"ht",1,0,0,-1,0,-1,0,0,2);
  // distribution("ControlRegion_wenu_spring15_25ns_2p32fb_newSF",1,1,1,0,j0v1,"j1pt",1,0,100,-1,0,-1,0,0,2);
  // distribution("ControlRegion_wenu_spring15_25ns_2p32fb_newSF",1,1,1,0,j0v1,"j1eta",0,0,0,-1,0,-1,0,0);
  // distribution("ControlRegion_wenu_spring15_25ns_2p32fb_newSF",1,1,1,0,j0v1,"nvtx",0,0,0,-1,0,-1,0,0);
  // distribution("ControlRegion_wenu_spring15_25ns_2p32fb_newSF",1,1,1,0,j0v1,"njets",1,0,0,-1,0,-1,0,0);
  // distribution("ControlRegion_wenu_spring15_25ns_2p32fb_newSF",1,1,1,0,j0v1,"Mt",1,0,0,-1,0,-1,0,0);
  // distribution("ControlRegion_wenu_spring15_25ns_2p32fb_newSF",1,1,1,0,j0v1,"lep1pt",1,0,0,-1,0,-1,0,0,4);
  // if (j0v1 != 1) distribution("ControlRegion_wenu_spring15_25ns_2p32fb_newSF",1,1,1,0,j0v1,"j2pt",1,0,0,-1,0,-1,0,0,2);
  // else {
  //   distribution("ControlRegion_wenu_spring15_25ns_2p32fb_newSF",1,1,1,0,j0v1,"prunedMass",0,0,0,-1,0,-1,0,0);
  //   distribution("ControlRegion_wenu_spring15_25ns_2p32fb_newSF",1,1,1,0,j0v1,"tau2OverTau1",0,0,0,-1,0,-1,0,0);
  // }

  // // W->enu control region (calibEle)
  // distribution("ControlRegion_wenu_spring15_25ns_2p32fb_newSF_CalibEle",1,1,1,0,j0v1,"metBin",1,0,200,-1,0,-1,1,0);
  // distribution("ControlRegion_wenu_spring15_25ns_2p32fb_newSF_CalibEle",1,1,1,0,j0v1,"met",1,0,200,-1,0,-1,0,0,2);
  // distribution("ControlRegion_wenu_spring15_25ns_2p32fb_newSF_CalibEle",1,1,1,0,j0v1,"ht",1,0,0,-1,0,-1,0,0,2);
  // distribution("ControlRegion_wenu_spring15_25ns_2p32fb_newSF_CalibEle",1,1,1,0,j0v1,"j1pt",1,0,100,-1,0,-1,0,0,2);
  // distribution("ControlRegion_wenu_spring15_25ns_2p32fb_newSF_CalibEle",1,1,1,0,j0v1,"j1eta",0,0,0,-1,0,-1,0,0);
  // distribution("ControlRegion_wenu_spring15_25ns_2p32fb_newSF_CalibEle",1,1,1,0,j0v1,"nvtx",0,0,0,-1,0,-1,0,0);
  // distribution("ControlRegion_wenu_spring15_25ns_2p32fb_newSF_CalibEle",1,1,1,0,j0v1,"njets",1,0,0,-1,0,-1,0,0);
  // distribution("ControlRegion_wenu_spring15_25ns_2p32fb_newSF_CalibEle",1,1,1,0,j0v1,"Mt",1,0,0,-1,0,-1,0,0);
  // distribution("ControlRegion_wenu_spring15_25ns_2p32fb_newSF_CalibEle",1,1,1,0,j0v1,"lep1pt",1,0,0,-1,0,-1,0,0,4);
  // if (j0v1 != 1) distribution("ControlRegion_wenu_spring15_25ns_2p32fb_newSF_CalibEle",1,1,1,0,j0v1,"j2pt",1,0,0,-1,0,-1,0,0,2);
  // else {
  //   distribution("ControlRegion_wenu_spring15_25ns_2p32fb_newSF_CalibEle",1,1,1,0,j0v1,"prunedMass",0,0,0,-1,0,-1,0,0);
  //   distribution("ControlRegion_wenu_spring15_25ns_2p32fb_newSF_CalibEle",1,1,1,0,j0v1,"tau2OverTau1",0,0,0,-1,0,-1,0,0);
  // }


  // void makeTransferFactor(const string folderNameWithRootFilesSR = "",
  // 			const string folderNameWithRootFilesCR = "",
  // 			const Int_t bkgToEstimate_z0_w1 = 0,
  // 			const Int_t monoJ0_monoV1 = 0,
  // 			//const Int_t signalRegion0_controlRegion1 = 0, 
  // 			//const Int_t z0_w1 = 0,
  // 			//const Int_t mu0_e1 = 0,
  // 			//const Int_t data0_noData1 = 0, 
  // 			//const string var = "met", 
  // 			//const Int_t yAxisLog_flag = 0, 
  // 			//const Int_t MCpoissonUncertainty_flag = 0, 
  // 			const Double_t xAxisMin = 0, 
  // 			const Double_t xAxisMax = -1, 
  // 			const Double_t yAxisMin = 0, 
  // 			const Double_t yAxisMax = -1, 
  // 			//const Int_t binDensity_flag = 0, 
  // 			//const Int_t MCnormalizedToData_flag = 0,
  // 			const Int_t numDenCorrelation = 0,
  // 			const Int_t rebinFactor = 1)
  // {

  if (j0v1 == 1) {
    makeTransferFactor("SignalRegion_spring15_25ns_2p32fb_newSF","ControlRegion_zmumu_spring15_25ns_2p32fb_newSF",0,j0v1,0,-1,0,10,1,1);
    makeTransferFactor("SignalRegion_spring15_25ns_2p32fb_newSF","ControlRegion_wmunu_spring15_25ns_2p32fb_newSF",1,j0v1,0,-1,0,0.8,1,1);
    makeTransferFactor("SignalRegion_spring15_25ns_2p32fb_newSF","SignalRegion_spring15_25ns_2p32fb_newSF",0,j0v1,0,-1,0,10,1,1);
  } else {
    makeTransferFactor("SignalRegion_spring15_25ns_2p32fb_newSF","ControlRegion_zmumu_spring15_25ns_2p32fb_newSF",0,j0v1,0,-1,4,15,1,1);
    makeTransferFactor("SignalRegion_spring15_25ns_2p32fb_newSF","ControlRegion_wmunu_spring15_25ns_2p32fb_newSF",1,j0v1,0,-1,0,-1,1,1);
    makeTransferFactor("SignalRegion_spring15_25ns_2p32fb_newSF","SignalRegion_spring15_25ns_2p32fb_newSF",0,j0v1,0,-1,0,10,1,1);
  }

  return 0;

}


