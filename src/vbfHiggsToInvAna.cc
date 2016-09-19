#define vbfHiggsToInvAna_cxx
#include "EmanTreeAnalysis.h"
//#include "functionsForAnalysis.h"
//#include "myClasses.h"
//C or C++ header files
#include <stdio.h>
#include <stdlib.h>
#include <cstdlib> //as stdlib.h
#include <cstdio>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>      // std::istringstream ; to read array of numbers from a line in a file
#include <string>
#include <vector>
#include <iomanip> //for input/output manipulators
//ROOT header files
#include <TAxis.h>
#include <TCanvas.h>
#include <TEfficiency.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TGraphAsymmErrors.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TMatrixDSym.h>
#include <TPad.h>
#include <TPaveStats.h>
#include <TPaveText.h>
#include <TStyle.h>
#include <TVector2.h>
#include <TVector3.h>
#include <TVirtualFitter.h>

using namespace std;
using namespace myAnalyzerTEman;

#ifdef vbfHiggsToInvAna_cxx

//===============================================

vbfHiggsToInvAna::vbfHiggsToInvAna(TTree *tree) : AnalysisDarkMatter(tree) {
  //cout <<"check in constructor "<<endl;

  // initialize some variables with sensible values. They will be set later depending on config file
  // following variables, if any, are specific to vbf H->inv analysis  

  DELTAETAMIN_VBFJETS = 3.6;
  JMET_DPHI_MIN = 2.3;
  INVMASS_VBFJETS = 600.0;
  JET1PT_VBFJETS = 80.0;
  JET2PT_VBFJETS = 70.0;

}

//===============================================


void vbfHiggsToInvAna::setNumberParameterValue(const string parameterName, const Double_t value) {

  AnalysisDarkMatter::setNumberParameterValue(parameterName, value);

  if (parameterName == "DELTAETAMIN_VBFJETS") DELTAETAMIN_VBFJETS = value;
  else if (parameterName == "INVMASS_VBFJETS") INVMASS_VBFJETS = value;
  else if (parameterName == "JET1PT_VBFJETS") JET1PT_VBFJETS = value;
  else if (parameterName == "JET2PT_VBFJETS") JET2PT_VBFJETS = value;

}

//===============================================

void vbfHiggsToInvAna::setVarFromConfigFile() {

  AnalysisDarkMatter::setVarFromConfigFile();

}

//===============================================

void vbfHiggsToInvAna::setSelections() {

  AnalysisDarkMatter::setSelections();

  vbfTaggedJets_deltaEtaC.set("|Deta(j1,j2)|",Form("|Deta(j1,j2)| > %2.1f",DELTAETAMIN_VBFJETS));
  vbfTaggedJets_inVMassC.set("M(j1,j2)",Form("Inv.Mass(j1,j2) > %3.0f",INVMASS_VBFJETS));
  vbfTaggedJets_jetsPtC.set("VBF jets pT",Form("pT j1(j2) > %2.0f(%2.0f)",JET1PT_VBFJETS,JET2PT_VBFJETS));
  jetMetDphiMinC.set("Dphi(j,MET)",Form("Dphi(jets,MET) > %1.1f",JMET_DPHI_MIN),"Dphi between any jet (pT>30, |eta|<4.7) and MET");

  selection::checkMaskLength();

}

//===============================================

void vbfHiggsToInvAna::setHistograms() {

  AnalysisDarkMatter::setHistograms();

  HvbfTaggedJets_deltaEta = new TH1D("HvbfTaggedJets_deltaEta","",100,0.0,10.0);
  HvbfTaggedJets_invMass = new TH1D("HvbfTaggedJets_invMass","",100,600.0,1600.0);
  HvbfTaggedJets_jet1pt = new TH1D("HvbfTaggedJets_jet1pt","",120,0.0,1200.0);
  HvbfTaggedJets_jet2pt = new TH1D("HvbfTaggedJets_jet2pt","",100,0.0,100.0);
  HvbfTaggedJets_jet1eta = new TH1D("HvbfTaggedJets_jet1eta","",100,-5.0,5.0);
  HvbfTaggedJets_jet2eta = new TH1D("HvbfTaggedJets_jet2eta","",100,-5.0,5.0);

}

//===============================================

void vbfHiggsToInvAna::setHistogramLastBinAsOverFlow(const Int_t hasScaledHistograms = 0) {

  AnalysisDarkMatter::setHistogramLastBinAsOverFlow(hasScaledHistograms);

  myAddOverflowInLastBin(HvbfTaggedJets_invMass);
  myAddOverflowInLastBin(HvbfTaggedJets_jet1pt);
  myAddOverflowInLastBin(HvbfTaggedJets_jet2pt);

}

//===============================================

void vbfHiggsToInvAna::createSystematicsHistogram() {

  AnalysisDarkMatter::createSystematicsHistogram();  

}

//===============================================


void vbfHiggsToInvAna::fillEventMask(UInt_t & eventMask) {

  AnalysisDarkMatter::fillEventMask(eventMask);  
  
  eventMask += vbfTaggedJets_deltaEtaC.addToMask(vbfTaggedJet_deltaEta > DELTAETAMIN_VBFJETS);
  eventMask += vbfTaggedJets_inVMassC.addToMask(vbfTaggedJet_invMass > INVMASS_VBFJETS);
  eventMask += vbfTaggedJets_jetsPtC.addToMask(vbfTaggedJet_leadJetPt > JET1PT_VBFJETS && vbfTaggedJet_trailJetPt > JET2PT_VBFJETS);

}

//===============================================


#endif
