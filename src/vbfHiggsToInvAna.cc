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

  //  cout << "CHECK VBFANA yyy"<< endl;
  AnalysisDarkMatter::setSelections();
  //cout << "CHECK VBFANA xxx"<< endl;

  vbfTaggedJets_deltaEtaC.set("|Deta(j1,j2)|",Form("|Deta(j1,j2)| > %1.1f",DELTAETAMIN_VBFJETS));
  //cout << "CHECK VBFANA aaa"<< endl;
  vbfTaggedJets_inVMassC.set("M(j1,j2)",Form("Inv.Mass(j1,j2) > %3.0f",INVMASS_VBFJETS));
  //cout << "CHECK VBFANA bbb"<< endl;
  vbfTaggedJets_jetsPtC.set("VBF jets pT",Form("pT j1(j2) > %2.0f(%2.0f)",JET1PT_VBFJETS,JET2PT_VBFJETS));
  //cout << "CHECK VBFANA 3"<< endl;
  //  jetMetDphiMinC.set("Dphi(j,MET)",Form("Dphi(jets,MET) > %1.1f",JMET_DPHI_MIN),"Dphi between any jet (pT>30, |eta|<4.7) and MET");
  jetNoiseCleaningC.set("jet1 noise clean","noise cleaning","energy fractions (only for jet1 in |eta| < 2.5): CHF > 0.1; NHF < 0.8");
  jetMetDphiMinC.set("Dphi(j,MET)",Form("Dphi(jets,MET) > %1.1f",JMET_DPHI_MIN),"Dphi between first 4 jets (pT>30, |eta|<4.7) and MET");
  if (MET_FILTERS_FLAG != 0) metFiltersC.set("met filters","met filters","EcalDeadCellTriggerPrimitiveFilter, HBHENoiseFilter, HBHENoiseIsoFilter, goodVertices, eeBadScFilter, globalTightHalo2016Filter");
  
  selection::checkMaskLength();

}

//===============================================

void vbfHiggsToInvAna::setHistograms() {

  AnalysisDarkMatter::setHistograms();

  HjetMetDphiMinAllJets = new TH1D("HjetMetDphiMinAllJets","",32,0.0,3.2);
  HvbfTaggedJets_mT = new TH1D("HvbfTaggedJets_mT","",150,0.0,3000.0);
  HvbfTaggedJets_deltaEta = new TH1D("HvbfTaggedJets_deltaEta","",100,0.0,10.0);
  HvbfTaggedJets_invMass = new TH1D("HvbfTaggedJets_invMass","",450,500.0,5000.0);

}

//===============================================

void vbfHiggsToInvAna::setHistogramLastBinAsOverFlow(const Int_t hasScaledHistograms = 0) {

  AnalysisDarkMatter::setHistogramLastBinAsOverFlow(hasScaledHistograms);

  myAddOverflowInLastBin(HvbfTaggedJets_invMass);
  myAddOverflowInLastBin(HvbfTaggedJets_mT);

}

//===============================================

void vbfHiggsToInvAna::createSystematicsHistogram() {

  AnalysisDarkMatter::createSystematicsHistogram();  

}

//===============================================


void vbfHiggsToInvAna::fillEventMask(ULong64_t & eventMask) {

  AnalysisDarkMatter::fillEventMask(eventMask);  
  
  eventMask += vbfTaggedJets_deltaEtaC.addToMask((JetClean_eta[0] * JetClean_eta[1]) < 0.0 && fabs(JetClean_eta[0] - JetClean_eta[1]) > DELTAETAMIN_VBFJETS);
  eventMask += vbfTaggedJets_inVMassC.addToMask(vbfJetsInvMass() > INVMASS_VBFJETS);
  eventMask += vbfTaggedJets_jetsPtC.addToMask( JetClean_pt[0] > JET1PT_VBFJETS && JetClean_pt[1] > JET2PT_VBFJETS);
  if (fabs(JetClean_eta[0]) < 2.5) eventMask += jetNoiseCleaningC.addToMask(JetClean_leadClean[0] > 0.5);
  else eventMask += jetNoiseCleaningC.addToMask(1);  // consider cut passed if jet is not central
  if (MET_FILTERS_FLAG != 0) eventMask += metFiltersC.addToMask(Flag_EcalDeadCellTriggerPrimitiveFilter == 1 && Flag_HBHENoiseFilter == 1 && Flag_HBHENoiseIsoFilter == 1 && Flag_goodVertices == 1 && Flag_eeBadScFilter == 1 && Flag_globalTightHalo2016Filter == 1);
  //  eventMask += jetMetDphiMinC.addToMask(fabs(dphijmAllJets) > JMET_DPHI_MIN);
  eventMask += jetMetDphiMinC.addToMask(fabs(dphijm) > JMET_DPHI_MIN);

}

//===============================================

Double_t vbfHiggsToInvAna::vbfJetsMT() {

  return TMath::Sqrt( 2. * JetClean_pt[0] * JetClean_pt[1] * (1. - TMath::Cos(dphijj)) );

}

//===============================================

Double_t vbfHiggsToInvAna::vbfJetsInvMass() {

  TLVjet1.SetPtEtaPhiM(JetClean_pt[0],JetClean_eta[0],JetClean_phi[0],0.0);
  TLVjet2.SetPtEtaPhiM(JetClean_pt[1],JetClean_eta[1],JetClean_phi[1],0.0);
  return (TLVjet1 + TLVjet2).Mag();

}

//===============================================


#endif
