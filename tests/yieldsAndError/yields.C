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
#include <TH1F.h>
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

Int_t checkFile(const TFile* f, const string fname) {

  if (!f || !f->IsOpen()) {
    cout<<"*******************************"<<endl;
    cout<<"Error opening file \""<<fname<<"\".\nSkipping."<<endl;
    cout<<"*******************************"<<endl;
    //exit(EXIT_FAILURE);
    return 1;
  }
  return 0;

}

//=====================================

void yields() {
  
  TFile *f = NULL;

  TH1F *htmp = NULL;
  vector<TH1F*> histo;

  vector<string> sample;
  sample.push_back("ZJetsToNuNu");
  sample.push_back("WJetsToLNu");
  sample.push_back("EWKZToNuNu");
  sample.push_back("EWKW");
  sample.push_back("Diboson");
  sample.push_back("Top");
  sample.push_back("GJets");
  sample.push_back("QCD");
  sample.push_back("DYJetsToLL");
  sample.push_back("EWKDYJetsLL");
  sample.push_back("ggH125");
  sample.push_back("vbfH125");

  vector<string> fname;
  string path = "/afs/cern.ch/work/m/mciprian/myCode/CMSSW_8_0_19/src/myMonoJetCode/output/vbfHiggsToInv/80X/lumi_2p3fb/SignalRegion_2015selection/vbfHiggsToInv_SR_";
  for (UInt_t i = 0; i < sample.size(); i++) {
    fname.push_back(path + sample[i] + ".root");
  }

  for (UInt_t i = 0; i < fname.size(); i++) {

    f = new TFile(fname[i].c_str(),"READ");
    if (checkFile(f,fname[i])) {
      histo.push_back(NULL);
    } else {

      htmp = (TH1F*)f->Get("HvbfTaggedJets_deltaEta");
      if (!htmp || htmp == NULL) {
	cout << "Error: histogram not found in file ' " << fname[i] << "'. End of programme." << endl;
	exit(EXIT_FAILURE);
      }
      histo.push_back((TH1F*) htmp->Clone());

    }

  }

  Double_t error = -1.0;
  Double_t integral = -1.0;

  for (UInt_t i = 0; i < histo.size(); i++) {
    if (histo[i] != NULL) integral = histo[i]->IntegralAndError(0, histo[i]->GetNbinsX()+1,error);
    cout << sample[i] << "\t" << integral << "\t" << error << endl;
  }

}


//=====================================

void SoverSqrtB(string srdir) {

  TFile *f = NULL;

  TH1F *htmp = NULL;
  vector<TH1F*> histoB;
  vector<TH1F*> histoS;

  vector<string> sampleB;
  sampleB.push_back("ZJetsToNuNu");
  sampleB.push_back("WJetsToLNu");
  sampleB.push_back("EWKZToNuNu");
  sampleB.push_back("EWKW");
  // samBple.push_back("Diboson");
  // sampleB.push_back("Top");
  // sampleB.push_back("GJets");
  // sampleB.push_back("QCD");
  // sampleB.push_back("DYJetsToLL");
  // sampleB.push_back("EWKDYJetsLL");

  vector<string> sampleS;
  sampleS.push_back("vbfH125");
  sampleS.push_back("ggH125");

  vector<string> fnameB;
  vector<string> fnameS;
  string path = "/afs/cern.ch/work/m/mciprian/myCode/CMSSW_8_0_19/src/myMonoJetCode/output/vbfHiggsToInv/80X/" + srdir + "/vbfHiggsToInv_SR_";
  for (UInt_t i = 0; i < sampleB.size(); i++) {
    fnameB.push_back(path + sampleB[i] + ".root");
  }
  for (UInt_t i = 0; i < sampleS.size(); i++) {
    fnameS.push_back(path + sampleS[i] + ".root");
  }

  for (UInt_t i = 0; i < fnameB.size(); i++) {

    f = new TFile(fnameB[i].c_str(),"READ");
    if (checkFile(f,fnameB[i])) {
      histoB.push_back(NULL);
    } else {

      htmp = (TH1F*)f->Get("HvbfTaggedJets_deltaEta");
      if (!htmp || htmp == NULL) {
	cout << "Error: histogram not found in file ' " << fnameB[i] << "'. End of programme." << endl;
	exit(EXIT_FAILURE);
      }
      histoB.push_back((TH1F*) htmp->Clone());

    }

  }

  for (UInt_t i = 0; i < fnameS.size(); i++) {

    f = new TFile(fnameS[i].c_str(),"READ");
    if (checkFile(f,fnameS[i])) {
      histoS.push_back(NULL);
    } else {

      htmp = (TH1F*)f->Get("HvbfTaggedJets_deltaEta");
      if (!htmp || htmp == NULL) {
	cout << "Error: histogram not found in file ' " << fnameS[i] << "'. End of programme." << endl;
	exit(EXIT_FAILURE);
      }
      histoS.push_back((TH1F*) htmp->Clone());

    }

  }

  Double_t background = 0.0;
  Double_t signal = 0.0;
  Double_t initialSevents = 10903. + 21890.;  //VBF + ggH
  //initialSevents.push_back();
  
  for (UInt_t i = 0; i < fnameS.size(); i++) {
    signal += histoS[i]->Integral();
  }

  for (UInt_t i = 0; i < fnameB.size(); i++) {
    background += histoB[i]->Integral();
  }

  Double_t effS = signal/initialSevents;
  cout << srdir << " : \t S/sqrt(B) = " << signal/sqrt(background) << "\t S/sqrt(S+B) = " << signal/sqrt(signal+background) 
       << "\tSignal efficiency = " << effS << endl;

}


//========================

void callSoverSqrtB() {

  SoverSqrtB("lumi_24p5fb/SignalRegion_looseSel_met130");
  SoverSqrtB("lumi_24p5fb/SignalRegion_looseSel_met150");
  SoverSqrtB("lumi_24p5fb/SignalRegion_looseSel_met175");
  SoverSqrtB("lumi_24p5fb/SignalRegion_looseSel_met200");

}
