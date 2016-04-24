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

void myGetEnvVariable(const string envVar, string& input) {

  // here we get environmental variable that we need to use the code, such as CMSSW_BASE                                                                               
  char* pPath;
  pPath = getenv (envVar.c_str());
  if (pPath!=NULL) {
    input = string(pPath);  // assign char* to string. Can also do --> string someString(char*);                  
    //cout << "Environment variable is : "<< input << endl;                                                                                                    
  }      


}

Int_t mySetMaxBinWidth(const TH1D *stackCopy, Double_t width) {

  Double_t maxBinWidth = 0.0;

  for (Int_t i = 1; i <= stackCopy->GetNbinsX(); i++) {
    if (stackCopy->GetBinWidth(i) > maxBinWidth) maxBinWidth = stackCopy->GetBinWidth(i);
  }
  width = maxBinWidth;
  if (fabs(width) < 0.00001) {
    cout << "WARNING: maximum bin width is 0. Returning 0 and setting width to 1." << endl;
    width = 1.0;
    return 0;
  } else return 1;

}

void setSampleName(const Int_t signalRegion0_controlRegion1, vector<string> &sampleName, vector<string> &MC_TexLabel, const Int_t z0_w1_g2 = 0, const Int_t mumu0_ee1 = 0) {

  // for Control region, mumu0_ee1 says if we use muon (0) or electron (1)
  // for trees with skim 1 lep Tight we don't have GJets samples in 76X

  if (signalRegion0_controlRegion1 == 0) {

    sampleName.push_back("GJets");
    sampleName.push_back("DYJetsToLL");
    sampleName.push_back("QCD");
    sampleName.push_back("Diboson");
    sampleName.push_back("Top");
    sampleName.push_back("WJetsToLNu");
    sampleName.push_back("ZJetsToNuNu");
       
    MC_TexLabel.push_back("#gamma + jets");
    MC_TexLabel.push_back("Z(ll)+jets");
    MC_TexLabel.push_back("QCD");
    MC_TexLabel.push_back("Diboson");
    MC_TexLabel.push_back("t#bar{t},single t");
    MC_TexLabel.push_back("W(l#nu)+jets");
    MC_TexLabel.push_back("Z(#nu#nu)+jets");
    
  } else {

    if(z0_w1_g2 == 0) {

      sampleName.push_back("ZJetsToNuNu");
      if (mumu0_ee1 == 0) sampleName.push_back("GJets");
      sampleName.push_back("QCD");
      sampleName.push_back("WJetsToLNu");
      sampleName.push_back("Diboson");
      sampleName.push_back("Top");
      sampleName.push_back("DYJetsToLL");

      MC_TexLabel.push_back("Z(#nu#nu)+jets");
      if (mumu0_ee1 == 0) MC_TexLabel.push_back("#gamma + jets");
      MC_TexLabel.push_back("QCD");
      MC_TexLabel.push_back("W(l#nu)+jets");
      MC_TexLabel.push_back("Diboson");
      MC_TexLabel.push_back("t#bar{t},single t");

      if (mumu0_ee1 == 0) MC_TexLabel.push_back("Z(#mu#mu)+jets");
      else MC_TexLabel.push_back("Z(ee)+jets");

    } else if(z0_w1_g2 == 1) {

      sampleName.push_back("ZJetsToNuNu");
      if (mumu0_ee1 == 0) sampleName.push_back("GJets");
      sampleName.push_back("DYJetsToLL");
      sampleName.push_back("QCD");      
      sampleName.push_back("Diboson");
      sampleName.push_back("Top");   
      sampleName.push_back("WJetsToLNu");

      MC_TexLabel.push_back("Z(#nu#nu)+jets");
      if (mumu0_ee1 == 0) MC_TexLabel.push_back("#gamma + jets");
      MC_TexLabel.push_back("Z(ll)+jets");
      MC_TexLabel.push_back("QCD");     
      MC_TexLabel.push_back("Diboson");
      MC_TexLabel.push_back("t#bar{t},single t");

      if (mumu0_ee1 == 0) MC_TexLabel.push_back("W(#mu#nu)+jets");
      else MC_TexLabel.push_back("W(e#nu)+jets");

    } else if(z0_w1_g2 == 2) {

      //sampleName.push_back("ZJetsToNuNu");
      sampleName.push_back("DYJetsToLL");
      //sampleName.push_back("Diboson");
      //sampleName.push_back("Top");   
      sampleName.push_back("WJetsToLNu");
      sampleName.push_back("QCD");      
      sampleName.push_back("GJets");

      //MC_TexLabel.push_back("Z(#nu#nu)+jets");
      MC_TexLabel.push_back("Z(ll)+jets");
      //MC_TexLabel.push_back("Diboson");
      //MC_TexLabel.push_back("t#bar{t},single t");
      MC_TexLabel.push_back("W(l#nu)+jets");
      MC_TexLabel.push_back("QCD");     
      MC_TexLabel.push_back("#gamma + jets");

    }

  }

}



void setHistColor(vector<Int_t> &histColor, const Int_t nSamples) {

  Int_t colorList[] = {kCyan, kViolet, kBlue, kRed, kYellow, kGreen, kOrange+1};  // the first color is for the main object. This array may contain more values than nSamples

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

    hvarName = "HrecoilDistribution"; 
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

  } else if ( !(strcmp("ph1pt",var.c_str())) ) {

    hvarName = "Hphoton1ptDistribution";
    xAxisName = "leading photon p_{T} [GeV]";

  } else if ( !(strcmp("ph1eta",var.c_str())) ) {

    hvarName = "Hphoton1etaDistribution";
    xAxisName = "leading photon #eta [GeV]";

  } else if ( !(strcmp("Mt",var.c_str())) ) {

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
		  const Int_t z0_w1_g2 = 0,
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

  // z0_w1_g2 is for CR to choose among Z or W to leptons or photons

  // mu0_e1 is for lepton flavour in CR (not used in SR)

  // data0_noData1 is to use or not a data file to compared with MC

  // yAxisLog_flag is to choose whether or not to use log scale in Y axis (default is 0, that is, normal scale)
  // xAxisMin and xAxisMax are the ranges for x Axis. Default values are used if xAxisMin > xAxisMax (otherwise user values are used)

  TH1::SetDefaultSumw2();            //all the following histograms will automatically call TH1::Sumw2() 

  gROOT->SetStyle("Plain");  // to have white legend (on my pc it's already white, but in tier2 it appears grey)
  gStyle->SetHistTopMargin(0.);
  gStyle->SetFillColor(10);

  string filenameExtension = ".root";
  // string fileDirectoryPath = "spring15_25ns_rootfiles/";
  //string fileDirectoryPath = "/cmshome/ciprianim/CMSSW721/output/monojet/" + folderNameWithRootFiles + "/";

  // here we get environmental variable that we need to use the code, such as CMSSW_BASE         
  string working_cmssw_path = "";
  myGetEnvVariable("CMSSW_BASE",working_cmssw_path);
  working_cmssw_path += "/src"; // now this string is $CMSSW_BASE/src      
  string fileDirectoryPath = working_cmssw_path + "/myMonoJetCode/output/monojet/76X/Fall15_25ns/lumi_2p32fb/" + folderNameWithRootFiles + "/"; 

  string plotDirectoryPath = fileDirectoryPath;
  // string plotDirectoryPath = "/cmshome/ciprianim/CMSSW721/pdfsFromAnalysis/plots/monojet/met_distribution/";
  string suffix;

  if (z0_w1_g2 == 0) {
    if (mu0_e1 == 0) suffix = "zmumu";
    else if (mu0_e1 == 1) suffix = "zee";
    else {
      cout << "Error: mu0_e1 must be 0 or 1. End of programme." << endl;
      exit(EXIT_FAILURE);
    }
  } else if (z0_w1_g2 == 1) {
    if (mu0_e1 == 0) suffix = "wmunu";
    else if (mu0_e1 == 1) suffix = "wenu";
    else {
      cout << "Error: mu0_e1 must be 0 or 1. End of programme." << endl;
      exit(EXIT_FAILURE);
    }
  } else if (z0_w1_g2 == 2) {
    suffix = "gamma";
  } else {
    cout << "Error: z0_w1_g2 must be 0, 1 or 2. End of programme." << endl;
    exit(EXIT_FAILURE);
  }

  TH1D* hvar = NULL;   // to get histogram from file
  string hvarName;          // name of histogram to take from file
  string xAxisName;        // name of X axis when plotting distribution. It is a tex string (with ROOT standard), e.g. "#slash{E}_{T} [GeV]" for MET

  setDistribution(mu0_e1, var, hvarName, xAxisName);
  if (monoJ0_monoV1 == 1) hvarName += "_monoV";

  vector<TH1D*> hMC;
  TH1D* hdata = NULL;

  vector<string> sampleName;
  vector<string> MC_TexLabel;
  if (data0_noData1 == 1) setSampleName(signalRegion0_controlRegion1, sampleName, MC_TexLabel, z0_w1_g2, mu0_e1);
  else setSampleName(signalRegion0_controlRegion1, sampleName, MC_TexLabel, z0_w1_g2, mu0_e1);

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
    
    canvasName = var + "_" + suffix + "jetsCR";
    filenameBase = suffix + "jets_CR_";

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

  string datafileName = fileDirectoryPath;

  if (data0_noData1 == 0) {

    if (signalRegion0_controlRegion1 == 0) datafileName += "monojet_SR_DATA.root";    
    else datafileName += suffix + "jets_CR_DATA.root";

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

  }


  // now here we go with the canvas

  TH1D * ratioplot = NULL; // will use it for the ratio plots
  TPad *subpad_1 = NULL;  // will use it to access specific subpad in canvas
  TPad *subpad_2 = NULL; 
  TLegend *leg = new TLegend(0.7,0.6,0.99,0.94);  
  TCanvas *c = NULL;
  if (data0_noData1 == 0) c = new TCanvas(canvasName.c_str(), (var + " distribution").c_str(), 700, 700);
  else c = new TCanvas(canvasName.c_str(), (var + " distribution").c_str());

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
    if (yAxisLog_flag) subpad_1->SetLogy();

    subpad_1->cd();

  } else if (yAxisLog_flag) c->SetLogy();

  
  hMCstack->Draw("HIST");
  //if (yAxisMin > 0) hMCstack->SetMinimum(yAxisMin);
  TH1D* stackCopy = (TH1D*)hMCstack->GetStack()->Last();   // used to make ratioplot without affecting the plot
  // TH1D* stackCopy = (TH1D*)(((TH1D*)hMCstack->GetStack()->Last())->DrawCopy("E2 SAME"));
  // stackCopy->SetFillColor(kBlack);
  // stackCopy->SetFillStyle(3144);

  // will use the following value to set the maximum Y axis value
  Double_t maxYvalue =  hMCstack->GetMaximum();
  if (data0_noData1 == 0 && (hMCstack->GetMaximum() < hdata->GetMaximum())) maxYvalue = hdata->GetMaximum();
  Double_t minYvalue =  hMCstack->GetMinimum();
  if (data0_noData1 == 0 && (hMCstack->GetMinimum() > hdata->GetMinimum())) minYvalue = hdata->GetMinimum();
  if (fabs(minYvalue) < 0.001) {
    minYvalue = 0.05;  // 0 can occur with data. Since data can be 1 or more (if not 0), then set minimum scale to 0.01 (arbitrary choice, might be not ideal if you have only MC and want to show lower values)
    if (binDensity_flag != 0) {
      // Double_t maxWidth = 1.0;
      // if (mySetMaxBinWidth(stackCopy, maxWidth)) minYvalue = 0.8 / maxWidth;  // not working as expected
      minYvalue = 0.001; 
      // in case bin density is shown, I saw that setting minYvalue to 0.05 is not ideal
      // This is because if I have empty bins (likely this happens when you have data and not just MC, although empty bins in MC might be present as well), then the minimum value is set to 0.05, but then I might not see bins which are not empty because their content (for data it would be >= 1) is also divided by bin width. In that case I set the minimum to be 20% lower than 1/MaxBinWidth, that is 0.8 * MaxBinWidth.
      // Note for future users/developers: I am aware this can all be made better and more generic, but I am forced to do things as fast as I can, so consider this to be a working temporary solution ;)
    }
  }

  if (!yAxisLog_flag) minYvalue = 0.0; // if using linear y axis, very likely the distribution will start from 0 (counting events)

  // ===> WARNING: seting the minimum to 0.01 could be bad if one has histograms normalized to 1

  //  TH1* histForAxis = hMCstack->GetHistogram();

  if (yAxisMin < yAxisMax) {  // default option [0,-1] implies default axis (this case is contemplated by the else condition below)
 
    if (data0_noData1 == 0) subpad_1->Update();  // to be done after Draw() to access pad parameters such as default axis range                                        
    else c->Update();
    //    hMCstack->GetYaxis()->SetRangeUser(yAxisMin,yAxisMax);
    // getting histogram used to draw axis (retrieved using THStack::GetHistogram() )                                                                                 
    // histForAxis->SetMaximum(yAxisMax);
    // if (fabs(yAxisMin) > 0.00001) histForAxis->SetMinimum(yAxisMin);  // if lower limit is set by user for Y axis, use it to set the minimum         
    // else histForAxis->SetMinimum(minYvalue);  // in lower limit left as default, use lowest possible value (for log scale, use 0.01 if it would be 0) 

    hMCstack->SetMaximum(yAxisMax);
    if (fabs(yAxisMin) > 0.00001) hMCstack->SetMinimum(yAxisMin);  // if lower limit is set by user for Y axis, use it to set the minimum                           
    else hMCstack->SetMinimum(minYvalue);  // if lower limit is left as default, use lowest possible value (for log scale, use 0.01 if it would be 0) 

    if (data0_noData1 == 0) subpad_1->Update();  // to be done after Draw() to access pad parameters such as default axis range
    else c->Update(); 

  } else if (yAxisMin > yAxisMax && (fabs(yAxisMax + 1.0) < 0.0001)) {  //if only lower bound is given (and the upper is left as -1) use default upper bound

    //Double_t stackYmax = hMCstack->GetYaxis()->GetXmax();

    if (data0_noData1 == 0) subpad_1->Update();  // to be done after Draw() to access pad parameters such as default axis range                                        
    else c->Update();

    // getting histogram used to draw axis (retrieved using THStack::GetHistogram() )                                                                               
    // if (yAxisLog_flag == 0) histForAxis->SetMaximum(1.2 * maxYvalue);  // for linear scale set upper Yaxis bound as 20% bigger than maximum histogram               
    // else histForAxis->SetMaximum(20 * maxYvalue);  // for log scale use 20 times bigger value for max Y                                                            
    // if (fabs(yAxisMin) > 0.00001) histForAxis->SetMinimum(yAxisMin);  // if only lower limit is set by user for Y axis, use it to set the minimum        
    // else histForAxis->SetMinimum(minYvalue);  // in any other case, use lowest possible value (for log scale, use 0.01 if it would be 0) 

    if (yAxisLog_flag == 0) hMCstack->SetMaximum(1.2 * maxYvalue);  // for linear scale set upper Yaxis bound as 20% bigger than maximum histogram                  
    else hMCstack->SetMaximum(10 * maxYvalue);  // for log scale use 10 times bigger value for max Y                                                                
    if (fabs(yAxisMin) > 0.00001) hMCstack->SetMinimum(yAxisMin);  // if only lower limit is set by user for Y axis, use it to set the minimum                      
    else hMCstack->SetMinimum(minYvalue);  // in any other case, use lowest possible value

    //hMCstack->GetYaxis()->SetRangeUser(yAxisMin,stackYmax);
    if (data0_noData1 == 0) subpad_1->Update();  // to be done after Draw() to access pad parameters such as default axis range
    else c->Update();  
    //hMCstack->GetYaxis()->SetRangeUser(yAxisMin,stackYmax);

  }

  // getting histogram used to draw axis (retrieved using THStack::GetHistogram() ) 
  // TH1* histForAxis = hMCstack->GetHistogram();
  // if (yAxisLog_flag == 0) histForAxis->SetMaximum(1.2 * maxYvalue);  // for linear scale set upper Yaxis bound as 20% bigger than maximum histogram
  // else histForAxis->SetMaximum(20 * maxYvalue);  // for log scale use 20 times bigger value for max Y
  // if (yAxisMin > yAxisMax && yAxisMax == -1) histForAxis->SetMinimum(yAxisMin);  // if only lower limit is set by user for Y axis, use it to set the minimum
  // else histForAxis->SetMinimum(minYvalue);  // in any other case, use lowest possible value (for log scale, use 0.01 if it would be 0)

  // //update pad just in case
  // if (data0_noData1 == 0) subpad_1->Update();  // to be done after Draw() to access pad parameters such as default axis range
  // else c->Update(); 
  // //hMCstack->SetMaximum(4000.0);

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

  } else if (xAxisMin > xAxisMax && (fabs(xAxisMax + 1.0) < 0.0001)) {  //if only lower bound is given (and the upper is left as -1) use default upper bound

    stackXmax = hMCstack->GetXaxis()->GetXmax();
    //cout << "stackXmax: " << stackXmax << endl;

    if (data0_noData1 == 0) subpad_1->Update();  // to be done after Draw() to access pad parameters such as default axis range 
    else c->Update();       //hMCstack->GetXaxis()->SetRangeUser(xAxisMin,c->GetX2());

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
    else if (xAxisMin > xAxisMax && (fabs(xAxisMax + 1.0) < 0.0001) ) {
      Double_t ratioplotXmax = ratioplot->GetXaxis()->GetXmax();
      //cout << "ratioplotXmax: " << ratioplotXmax << endl;
      subpad_2->Update();
      ratioplot->GetXaxis()->SetRangeUser(xAxisMin,ratioplotXmax); 
      //ratioplot->GetXaxis()->SetRangeUser(xAxisMin,stackXmax); //ok, it works
    }
    ratioplot->DrawCopy("E");

  }

  c->SaveAs( (plotDirectoryPath + c->GetName() + ".pdf").c_str() );
  c->SaveAs( (plotDirectoryPath + c->GetName() + ".png").c_str() );

}


// =========================================================


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

  // here we get environmental variable that we need to use the code, such as CMSSW_BASE                              
  string working_cmssw_path = "";
  myGetEnvVariable("CMSSW_BASE",working_cmssw_path);
  working_cmssw_path += "/src"; // now this string is $CMSSW_BASE/src      

  string filenameExtension = ".root";
  // string fileDirectoryPath = "spring15_25ns_rootfiles/";
  string fileDirectoryPathSR = working_cmssw_path + "/myMonoJetCode/output/monojet/76X/Fall15_25ns/lumi_2p32fb/" + folderNameWithRootFilesSR + "/";
  string fileDirectoryPathCR = working_cmssw_path + "/myMonoJetCode/output/monojet/76X/Fall15_25ns/lumi_2p32fb/" + folderNameWithRootFilesCR + "/";
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
  suffixCR.push_back("wenu");  texSuffixCR.push_back("W(e#nu)");
  suffixCR.push_back("gamma");  texSuffixCR.push_back("#gamma");

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

  // decide things about canvas name, title (in tex format) depending on zvv_wlv_flag: if 1 (set) then we do Z(vv)/W(lv), otherwise we do usual scale factor SR/CR
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

  // won't use CalibEle variables, but keep that for now just in case
  // string calibEleSuffix = "_CalibEle";  // used to select proper directory where to save plot, depending on input file name
  // size_t useCalibEle = fileDirectoryPathCR.find(calibEleSuffix);
  // if (useCalibEle != string::npos) plotDirName += calibEleSuffix;
  
  // scale factors plot are saved in the same directory as SR plots. It's just a choice of convenience since they are applied to SR
  string plotDirectoryPath = fileDirectoryPathSR;
  string plotFileExtension = ".pdf";

  // create name for SR, e.g: <path_to_plots>/SignalRegion_25ns_2p32fb/monojet_SR_ZJetsToNuNu.root
  // this is the SR file to be opened: similar thing done for CR
  string filenameSR = fileDirectoryPathSR + "monojet_SR_" + bkgSampleSR + filenameExtension;

  // start CR root file name with the directory where it is (note, it might be the same as SR, in this case fileDirectoryPathCR == fileDirectoryPathSR)
  string filenameCR = fileDirectoryPathCR;
  if (zvv_wlv_flag) {
    filenameCR += "monojet_SR_WJetsToLNu" + filenameExtension;
  } else {
    filenameCR += suffixCR.at(indexRightSuffixCR) + "jets_CR_"; 
    if (plotDirName.find("z") != string::npos) filenameCR += ("DYJetsToLL" + filenameExtension);
    else if (plotDirName.find("w") != string::npos) filenameCR += ("WJetsToLNu" + filenameExtension); 
    else if (plotDirName.find("gamma") != string::npos) filenameCR += ("GJets" + filenameExtension); 
  }

  // now open files and get histograms
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

  // now that hSR and hCR have been used with only statistical uncertainty, compute their total uncertainty and use it in ratio

  for (Int_t i = 1; i <= (hSR->GetNbinsX()+1); i++) {

    //computing total uncertainty summing stat and syst in quadrature (histogram with systematic has uncertainty as bin content)
    // to stat. uncertainty (GetBinError) sum syst. one (given by the content of hsyst, and not the error of hsyst, since hsyst was devised as the syst unc of HYieldsMetBin histogram)
    Double_t totUnc = 0.0;
    // ---> handling correlation, to be revised  <----
    if (numDenCorrelation == 0) {
      // add uncert. in quadrature for numerator
      totUnc = hSR->GetBinError(i) * hSR->GetBinError(i) + hSRsyst->GetBinContent(i) * hSRsyst->GetBinContent(i); 
      // assign uncert. as square root of quadratic sum above
      hSR->SetBinError(i, sqrt(totUnc));
      // and here we go again for the denominator
      totUnc = hCR->GetBinError(i) * hCR->GetBinError(i) + hCRsyst->GetBinContent(i) * hCRsyst->GetBinContent(i);
      hCR->SetBinError(i, sqrt(totUnc));
    } // else if (numDenCorrelation == 1) {   // in this case just consider the statistical uncertainty (default)
    //   // if total correlation holds, systematic will cancel out in ratio
    //   totUnc = hSR->GetBinError(i) * hSR->GetBinError(i); 
    //   hSR->SetBinError(i, sqrt(totUnc));
    //   totUnc = hCR->GetBinError(i) * hCR->GetBinError(i);
    //   hCR->SetBinError(i, sqrt(totUnc));
    // }

  }


  // now there is a systematic uncertainty related to the (charged) lepton ID. This is 2% for each lepton, thus 2% for W ad 4% for Z in control regions only.
  // this uncertainty can be summed in quadrature to the others
  // this is applied to the denominator only if it comes from a control region (for Zvv/Wlv both numerator and denominator are taken from the signal region)
  // this uncertainty is set after implementing all the other uncertainties deriving from change in QCD or EWK scale (done just above)
  // ============================================================= 
  // NOTE: commented, since now it is implemented as a scale factor
  // =============================================================         
  // if (!zvv_wlv_flag) {

  //   Double_t lepIDuncertainty = 0.0; // 2% for 1 lep CR, 4% for 2 lep CR, which is 0.02 * binContent or 0.04 * binContent
  //   if (plotDirName.find("z") != string::npos) lepIDuncertainty = 0.04 * 0.04;       // computing square of the factor to use it directly below
  //   else if (plotDirName.find("w") != string::npos) lepIDuncertainty = 0.02 * 0.02;  // computing square of the factor to use it directly below

  //   for (Int_t i = 1; i <= (hCR->GetNbinsX()+1); i++) {

  //     //computing total uncertainty summing stat and syst in quadrature (histogram with systematic has uncertainty as bin content)
  //     Double_t totUnc = hCR->GetBinError(i) * hCR->GetBinError(i) + lepIDuncertainty * hCR->GetBinContent(i) * hCR->GetBinContent(i); 
  //     // note that here the bin error used in totUnc is the sum in quadrature of statistical and other systematic uncertainties   
  //     hCR->SetBinError(i, sqrt(totUnc));

  //   }

  // }

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

  // get maximum and minimum value of histogram including uncertainty, so to choose y axis range automatically 


  //Double_t maxYvalue =  hTransferFactor_totUnc->GetMaximum();
  //Double_t minYvalue =  hTransferFactor_totUnc->GetMinimum();
  Double_t maxYvalue =  0.0;  
  Double_t minYvalue =  10000.0;

  for (Int_t i = 0; i <= (hTransferFactor_totUnc->GetNbinsX()+1); i++) {
    Double_t tmp = hTransferFactor_totUnc->GetBinContent(i) + hTransferFactor_totUnc->GetBinError(i);
    if ( tmp > maxYvalue) maxYvalue = tmp;
    tmp = hTransferFactor_totUnc->GetBinContent(i) - hTransferFactor_totUnc->GetBinError(i);
    if (tmp < minYvalue) minYvalue = (tmp < 0.0) ? 0.0 : tmp;

  }

  // now that we have max and min for scale factors including the uncertainty band, set them to a slightly wider range to have a nice plot
  maxYvalue += 0.2 * maxYvalue;
  minYvalue -= 0.2 * maxYvalue;
  if (minYvalue < 0.0) minYvalue = 0.0;

  if (yAxisMin < yAxisMax) {  // default option [0,-1] implies default axis

    c->Update();  
    // hTransferFactor_totUnc->GetYaxis()->SetRangeUser(yAxisMin,yAxisMax);
    hTransferFactor_totUnc->SetMaximum(yAxisMax);
    if (fabs(yAxisMin) > 0.00001) hTransferFactor_totUnc->SetMinimum(yAxisMin);  // if lower limit is set by user for Y axis, use it to set the minimum   
    //else hTransferFactor_totUnc->SetMinimum(0.5 * minYvalue);  // if lower limit is left as default, use half of the lowest possible value
    c->Update();  

  } else if ( yAxisMin > yAxisMax && (fabs(yAxisMax + 1.0) < 0.0001) ) {  //if only lower bound is given (and the upper is left as -1) use default upper bound

    //   cout << "CHECK CCC" << endl;
    c->Update();  
    //hTransferFactor_totUnc->GetYaxis()->SetRangeUser(yAxisMin,c->GetY2());
    hTransferFactor_totUnc->SetMaximum(maxYvalue);  
    if (fabs(yAxisMin) > 0.00001) hTransferFactor_totUnc->SetMinimum(yAxisMin);  // if only lower limit is set by user for Y axis, use it to set the minimum   
    //else hTransferFactor_totUnc->SetMinimum(minYvalue);
    //   cout << "CHECK BBB" << endl;
    c->Update();  

  }

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
    c->Update();  

  } else if ( xAxisMin > xAxisMax && (fabs(xAxisMax + 1.0) < 0.0001) ) {  //if only lower bound is given (and the upper is left as -1) use default upper bound
    
    c->Update();   
    hTransferFactor_totUnc->GetXaxis()->SetRangeUser(xAxisMin,hTransferFactor_totUnc->GetXaxis()->GetXmax());   // N.B. SetRangeUser cannot set axis ranges outside the original coordinates: if they are 5-10, I can set them 6-9, 5-8 ecc... but if I do 5-12, the range will be set to 5-10
    c->Update();  

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
	       const Int_t z0_w1_g2 = 0,
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
  distribution("SignalRegion",0,0,0,0,j0v1,"metBin",1,0,200,-1,0,-1,1,0);
  distribution("SignalRegion",0,0,0,0,j0v1,"met",1,0,200,-1,0,-1,0,0,2);
  distribution("SignalRegion",0,0,0,0,j0v1,"ht",1,0,0,-1,0,-1,0,0,4);
  distribution("SignalRegion",0,0,0,0,j0v1,"j1pt",1,0,100,-1,0,-1,0,0,2);
  distribution("SignalRegion",0,0,0,0,j0v1,"j1eta",0,0,0,-1,0,-1,0,0);
  distribution("SignalRegion",0,0,0,0,j0v1,"nvtx",0,0,0,-1,0,-1,0,0);
  distribution("SignalRegion",0,0,0,0,j0v1,"njets",1,0,0,-1,0,-1,0,0);
  if (j0v1 != 1) distribution("SignalRegion",0,0,0,0,j0v1,"j2pt",1,0,0,-1,0,-1,0,0,2);
  else {
    distribution("SignalRegion",0,0,0,0,j0v1,"prunedMass",0,0,0,-1,0,-1,0,0);
    distribution("SignalRegion",0,0,0,0,j0v1,"tau2OverTau1",0,0,0,-1,0,-1,0,0);
  }

  // Z->mumu control region
  distribution("ControlRegion_zmumu",1,0,0,0,j0v1,"metBin",1,0,200,-1,0,-1,1,0);
  distribution("ControlRegion_zmumu",1,0,0,0,j0v1,"met",1,0,200,-1,0,-1,0,0,2);
  distribution("ControlRegion_zmumu",1,0,0,0,j0v1,"ht",1,0,0,-1,0,-1,0,0,4);
  distribution("ControlRegion_zmumu",1,0,0,0,j0v1,"j1pt",1,0,100,-1,0,-1,0,0,2);
  distribution("ControlRegion_zmumu",1,0,0,0,j0v1,"j1eta",0,0,0,-1,0,-1,0,0);
  distribution("ControlRegion_zmumu",1,0,0,0,j0v1,"nvtx",0,0,0,-1,0,-1,0,0);
  distribution("ControlRegion_zmumu",1,0,0,0,j0v1,"njets",1,0,0,-1,0,-1,0,0);
  distribution("ControlRegion_zmumu",1,0,0,0,j0v1,"zpt",1,0,0,-1,0,-1,0,0,4);
  distribution("ControlRegion_zmumu",1,0,0,0,j0v1,"invMass",0,0,0,-1,0,-1,0,0);
  distribution("ControlRegion_zmumu",1,0,0,0,j0v1,"lep1pt",1,0,0,-1,0,-1,0,0,4);
  distribution("ControlRegion_zmumu",1,0,0,0,j0v1,"lep2pt",1,0,0,-1,0,-1,0,0,4);
  if (j0v1 != 1) distribution("ControlRegion_zmumu",1,0,0,0,j0v1,"j2pt",1,0,0,-1,0,-1,0,0,2);
  else {
    distribution("ControlRegion_zmumu",1,0,0,0,j0v1,"prunedMass",0,0,0,-1,0,-1,0,0);
    distribution("ControlRegion_zmumu",1,0,0,0,j0v1,"tau2OverTau1",0,0,0,-1,0,-1,0,0);
  }

  // Z->ee control region
  distribution("ControlRegion_zee",1,0,1,0,j0v1,"metBin",1,0,200,-1,0,-1,1,0);
  distribution("ControlRegion_zee",1,0,1,0,j0v1,"met",1,0,200,-1,0,-1,0,0,2);
  distribution("ControlRegion_zee",1,0,1,0,j0v1,"ht",1,0,0,-1,0,-1,0,0,2);
  distribution("ControlRegion_zee",1,0,1,0,j0v1,"j1pt",1,0,100,-1,0,-1,0,0,2);
  distribution("ControlRegion_zee",1,0,1,0,j0v1,"j1eta",0,0,0,-1,0,-1,0,0);
  distribution("ControlRegion_zee",1,0,1,0,j0v1,"nvtx",0,0,0,-1,0,-1,0,0);
  distribution("ControlRegion_zee",1,0,1,0,j0v1,"njets",1,0,0,-1,0,-1,0,0);
  distribution("ControlRegion_zee",1,0,1,0,j0v1,"zpt",1,0,0,-1,0,-1,0,0,4);
  distribution("ControlRegion_zee",1,0,1,0,j0v1,"invMass",0,0,0,-1,0,-1,0,0);
  distribution("ControlRegion_zee",1,0,1,0,j0v1,"lep1pt",1,0,0,-1,0,-1,0,0,4);
  distribution("ControlRegion_zee",1,0,1,0,j0v1,"lep2pt",1,0,0,-1,0,-1,0,0,4);
  if (j0v1 != 1) distribution("ControlRegion_zee",1,0,1,0,j0v1,"j2pt",1,0,0,-1,0,-1,0,0,2);
  else {
    distribution("ControlRegion_zee",1,0,1,0,j0v1,"prunedMass",0,0,0,-1,0,-1,0,0);
    distribution("ControlRegion_zee",1,0,1,0,j0v1,"tau2OverTau1",0,0,0,-1,0,-1,0,0);
  }

  // W->munu control region
  distribution("ControlRegion_wmunu",1,1,0,0,j0v1,"metBin",1,0,200,-1,0,-1,1,0);
  distribution("ControlRegion_wmunu",1,1,0,0,j0v1,"met",1,0,200,-1,0,-1,0,0,2);
  distribution("ControlRegion_wmunu",1,1,0,0,j0v1,"ht",1,0,0,-1,0,-1,0,0,4);
  distribution("ControlRegion_wmunu",1,1,0,0,j0v1,"j1pt",1,0,100,-1,0,-1,0,0,2);
  distribution("ControlRegion_wmunu",1,1,0,0,j0v1,"j1eta",0,0,0,-1,0,-1,0,0);
  distribution("ControlRegion_wmunu",1,1,0,0,j0v1,"nvtx",0,0,0,-1,0,-1,0,0);
  distribution("ControlRegion_wmunu",1,1,0,0,j0v1,"njets",1,0,0,-1,0,-1,0,0);
  distribution("ControlRegion_wmunu",1,1,0,0,j0v1,"Mt",0,0,0,-1,0,-1,0,0);
  distribution("ControlRegion_wmunu",1,1,0,0,j0v1,"lep1pt",1,0,0,-1,0,-1,0,0,4);
  if (j0v1 != 1) distribution("ControlRegion_wmunu",1,1,0,0,j0v1,"j2pt",1,0,0,-1,0,-1,0,0,2);
  else {
    distribution("ControlRegion_wmunu",1,1,0,0,j0v1,"prunedMass",0,0,0,-1,0,-1,0,0);
    distribution("ControlRegion_wmunu",1,1,0,0,j0v1,"tau2OverTau1",0,0,0,-1,0,-1,0,0);
  }

  // W->enu control region
  distribution("ControlRegion_wenu",1,1,1,0,j0v1,"metBin",1,0,200,-1,0,-1,1,0);
  distribution("ControlRegion_wenu",1,1,1,0,j0v1,"met",1,0,200,-1,0,-1,0,0,2);
  distribution("ControlRegion_wenu",1,1,1,0,j0v1,"ht",1,0,0,-1,0,-1,0,0,4);
  distribution("ControlRegion_wenu",1,1,1,0,j0v1,"j1pt",1,0,100,-1,0,-1,0,0,2);
  distribution("ControlRegion_wenu",1,1,1,0,j0v1,"j1eta",0,0,0,-1,0,-1,0,0);
  distribution("ControlRegion_wenu",1,1,1,0,j0v1,"nvtx",0,0,0,-1,0,-1,0,0);
  distribution("ControlRegion_wenu",1,1,1,0,j0v1,"njets",1,0,0,-1,0,-1,0,0);
  distribution("ControlRegion_wenu",1,1,1,0,j0v1,"Mt",0,0,0,-1,0,-1,0,0);
  distribution("ControlRegion_wenu",1,1,1,0,j0v1,"lep1pt",1,0,0,-1,0,-1,0,0,4);
  if (j0v1 != 1) distribution("ControlRegion_wenu",1,1,1,0,j0v1,"j2pt",1,0,0,-1,0,-1,0,0,2);
  else {
    distribution("ControlRegion_wenu",1,1,1,0,j0v1,"prunedMass",0,0,0,-1,0,-1,0,0);
    distribution("ControlRegion_wenu",1,1,1,0,j0v1,"tau2OverTau1",0,0,0,-1,0,-1,0,0);
  }

  // Gamma control region
  distribution("ControlRegion_gamma",1,2,0,0,j0v1,"metBin",1,0,200,-1,0,-1,1,0);
  distribution("ControlRegion_gamma",1,2,0,0,j0v1,"met",1,0,200,-1,0,-1,0,0,2);
  distribution("ControlRegion_gamma",1,2,0,0,j0v1,"ht",1,0,0,-1,0,-1,0,0,4);
  distribution("ControlRegion_gamma",1,2,0,0,j0v1,"j1pt",1,0,100,-1,0,-1,0,0,2);
  distribution("ControlRegion_gamma",1,2,0,0,j0v1,"j1eta",0,0,0,-1,0,-1,0,0);
  distribution("ControlRegion_gamma",1,2,0,0,j0v1,"nvtx",0,0,0,-1,0,-1,0,0);
  distribution("ControlRegion_gamma",1,2,0,0,j0v1,"njets",1,0,0,-1,0,-1,0,0);
  distribution("ControlRegion_gamma",1,2,0,0,j0v1,"ph1pt",1,0,150,-1,0,-1,0,0,2);
  distribution("ControlRegion_gamma",1,2,1,0,j0v1,"ph1eta",0,0,0,-1,0,-1,0,0);
  if (j0v1 != 1) distribution("ControlRegion_gamma",1,2,0,0,j0v1,"j2pt",1,0,0,-1,0,-1,0,0,2);
  else {
    distribution("ControlRegion_gamma",1,2,0,0,j0v1,"prunedMass",0,0,0,-1,0,-1,0,0);
    distribution("ControlRegion_gamma",1,2,0,0,j0v1,"tau2OverTau1",0,0,0,-1,0,-1,0,0);
  }


  // void makeTransferFactor(const string folderNameWithRootFilesSR = "",
  // 			const string folderNameWithRootFilesCR = "",
  // 			const Int_t bkgToEstimate_z0_w1 = 0,
  // 			const Int_t monoJ0_monoV1 = 0,
  // 			const Double_t xAxisMin = 0, 
  // 			const Double_t xAxisMax = -1, 
  // 			const Double_t yAxisMin = 0, 
  // 			const Double_t yAxisMax = -1, 
  // 			const Int_t numDenCorrelation = 0,
  // 			const Int_t rebinFactor = 1)
  // {

  if (j0v1 == 1) {
    makeTransferFactor("SignalRegion","ControlRegion_zmumu",0,j0v1,0,-1,0,-1,1,1);
    makeTransferFactor("SignalRegion","ControlRegion_zee",0,j0v1,0,-1,0,-1,1,1);
    makeTransferFactor("SignalRegion","ControlRegion_wmunu",1,j0v1,0,-1,0,-1,1,1);
    makeTransferFactor("SignalRegion","ControlRegion_wenu",1,j0v1,0,-1,0,0.-1,1,1);
    makeTransferFactor("SignalRegion","ControlRegion_gamma",0,j0v1,0,-1,0,0.-1,0,1);
    makeTransferFactor("SignalRegion","SignalRegion",0,j0v1,0,-1,0,-1,1,1);
  } else {
    makeTransferFactor("SignalRegion","ControlRegion_zmumu",0,j0v1,0,-1,0,-1,1,1);
    makeTransferFactor("SignalRegion","ControlRegion_zee",0,j0v1,0,-1,0,-1,1,1);
    makeTransferFactor("SignalRegion","ControlRegion_wmunu",1,j0v1,0,-1,0,-1,1,1);
    makeTransferFactor("SignalRegion","ControlRegion_wenu",1,j0v1,0,-1,0,-1,1,1);
    makeTransferFactor("SignalRegion","ControlRegion_gamma",0,j0v1,0,-1,0,-1,0,1);
    makeTransferFactor("SignalRegion","SignalRegion",0,j0v1,0,-1,0,-1,1,1);
  }

  return 0;

}


