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
  HvbfTaggedJets_invMass = new TH1D("HvbfTaggedJets_invMass","",500,0.0,5000.0);

}

//===============================================

void vbfHiggsToInvAna::setHistogramLastBinAsOverFlow(const Int_t hasScaledHistograms = 0) {

  AnalysisDarkMatter::setHistogramLastBinAsOverFlow(hasScaledHistograms);

  myAddOverflowInLastBin(HvbfTaggedJets_invMass);
  myAddOverflowInLastBin(HvbfTaggedJets_mT);

}

//===============================================

void vbfHiggsToInvAna::setBranchStatusForAnalysis() {
  
  if (fChain == 0) return;
  
  fChain->SetBranchStatus("*",0);  
  // warning: in Emanuele's trees non integer values are float

  fChain->SetBranchStatus("nMu10V",1);  // # of muons passing loose selection
  fChain->SetBranchStatus("nEle10V",1);  // # of electrons passing loose selection for electron veto
  fChain->SetBranchStatus("nGamma15V",1);  // # of photons passing loose selection for photon veto
  fChain->SetBranchStatus("nTauClean18V",1);

  fChain->SetBranchStatus("nJetClean",1);    // # of jet with pt > 30 & eta < 4.7 and cleaning for against muons misidentified as PFjets   
  fChain->SetBranchStatus("JetClean_pt",1);  
  fChain->SetBranchStatus("JetClean_eta",1);  
  fChain->SetBranchStatus("JetClean_phi",1);  
  fChain->SetBranchStatus("JetClean_mass",1);  

  fChain->SetBranchStatus("Jet_chHEF",1);  // < 0.99 to reduce spikes until we have caloMET  

  fChain->SetBranchStatus("metNoMu_pt",1);
  //fChain->SetBranchStatus("metNoMu_eta",1);
  fChain->SetBranchStatus("metNoMu_phi",1);
  fChain->SetBranchStatus("htJet25",1);
  fChain->SetBranchStatus("met_pt",1);
  //fChain->SetBranchStatus("met_eta",1);                                                                                                                              
  fChain->SetBranchStatus("met_phi",1);

  fChain->SetBranchStatus("nVert",1);  // number of good vertices 

  fChain->SetBranchStatus("HLT_MonoJetMetNoMuMHT90",1);
  fChain->SetBranchStatus("HLT_MonoJetMetNoMuMHT120",1);
  fChain->SetBranchStatus("HLT_Met170",1);

  // met filters to be used (the config file has a parameter saying whether they should be used or not)
  fChain->SetBranchStatus("Flag_EcalDeadCellTriggerPrimitiveFilter",1);
  fChain->SetBranchStatus("Flag_HBHENoiseFilter",1);
  fChain->SetBranchStatus("Flag_HBHENoiseIsoFilter",1);
  fChain->SetBranchStatus("Flag_goodVertices",1);
  fChain->SetBranchStatus("Flag_eeBadScFilter",1);
  fChain->SetBranchStatus("Flag_globalTightHalo2016Filter",1);
  
  fChain->SetBranchStatus("dphijj",1);          // dphi between 1st and 2nd jet, 999 if second jet doesn't exist
  fChain->SetBranchStatus("dphijm",1);          // flag for dphi minimum between met and any of the jets in the event (using only the first four jets
  fChain->SetBranchStatus("nBTag15",1);  // for b-jet veto
  fChain->SetBranchStatus("JetClean_leadClean",1); // has new cleaning on energy fractions (added on 17 November 2015) 

  fChain->SetBranchStatus("dphijmAllJets",1);   
  fChain->SetBranchStatus("vbfTaggedJet_deltaEta",1);   
  fChain->SetBranchStatus("vbfTaggedJet_invMass",1);   
  fChain->SetBranchStatus("vbfTaggedJet_leadJetPt",1);   
  fChain->SetBranchStatus("vbfTaggedJet_trailJetPt",1);   
  fChain->SetBranchStatus("vbfTaggedJet_leadJetEta",1);   
  fChain->SetBranchStatus("vbfTaggedJet_trailJetEta",1);   

  fChain->SetBranchStatus("puw",1); // added on 21 Sept 2016, substituting vtxWeight
  fChain->SetBranchStatus("weight",1);   // modified since 17 November 2015: now it includes the whol weight, e.g. 1000*xsec*genWeight ...

  if (!ISDATA_FLAG) {
 
    fChain->SetBranchStatus("SF_BTag",1);
    //variables in sfFriend tree
    fChain->SetBranchStatus("SF_trig1lep",1);
    fChain->SetBranchStatus("SF_trigmetnomu",1);
    fChain->SetBranchStatus("SF_LepTightLoose",1);
    fChain->SetBranchStatus("SF_LepTightLooseUp",1);
    fChain->SetBranchStatus("SF_LepTightLooseDown",1);
    fChain->SetBranchStatus("SF_LepTight",1);
    fChain->SetBranchStatus("SF_LepTightUp",1);
    fChain->SetBranchStatus("SF_LepTightDown",1);
    fChain->SetBranchStatus("SF_NLO_QCD",1);
    fChain->SetBranchStatus("SF_NLO_QCD_renScaleUp",1);
    fChain->SetBranchStatus("SF_NLO_QCD_renScaleDown",1);
    fChain->SetBranchStatus("SF_NLO_QCD_facScaleUp",1);
    fChain->SetBranchStatus("SF_NLO_QCD_facScaleDown",1);
    fChain->SetBranchStatus("SF_NLO_QCD_pdfUp",1);
    fChain->SetBranchStatus("SF_NLO_QCD_pdfDown",1);
    fChain->SetBranchStatus("SF_NLO_EWK",1);
    fChain->SetBranchStatus("SF_NLO_EWK_up",1);
    fChain->SetBranchStatus("SF_NLO_EWK_down",1);

  }


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
  eventMask += vbfTaggedJets_jetsPtC.addToMask( nJetClean > 1 && JetClean_pt[0] > JET1PT_VBFJETS && JetClean_pt[1] > JET2PT_VBFJETS);
  if (fabs(JetClean_eta[0]) < 2.5) eventMask += jetNoiseCleaningC.addToMask(JetClean_leadClean[0] > 0.5 && Jet_chHEF[0] < 0.99);
  else eventMask += jetNoiseCleaningC.addToMask(1);  // consider cut passed if jet is not central
  if (MET_FILTERS_FLAG != 0) eventMask += metFiltersC.addToMask(Flag_EcalDeadCellTriggerPrimitiveFilter == 1 && Flag_HBHENoiseFilter == 1 
								&& Flag_HBHENoiseIsoFilter == 1 && Flag_goodVertices == 1 && Flag_eeBadScFilter == 1 
								&& Flag_globalTightHalo2016Filter == 1);
  //  eventMask += jetMetDphiMinC.addToMask(fabs(dphijmAllJets) > JMET_DPHI_MIN);
  eventMask += jetMetDphiMinC.addToMask(fabs(dphijm) > JMET_DPHI_MIN);

}

//===============================================

Double_t vbfHiggsToInvAna::vbfJetsMT() {

  return TMath::Sqrt( 2. * JetClean_pt[0] * JetClean_pt[1] * (1. - TMath::Cos(dphijj)) );

}

//===============================================

Double_t vbfHiggsToInvAna::vbfJetsInvMass() {

  TLVjet1.SetPtEtaPhiM(JetClean_pt[0],JetClean_eta[0],JetClean_phi[0],JetClean_mass[0]);
  TLVjet2.SetPtEtaPhiM(JetClean_pt[1],JetClean_eta[1],JetClean_phi[1],JetClean_mass[1]);
  return (TLVjet1 + TLVjet2).Mag();

}

//===============================================


#endif
