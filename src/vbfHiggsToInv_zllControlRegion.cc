#define vbfHiggsToInv_zllControlRegion_cxx
#include "EmanTreeAnalysis.h"
//#include "monojetAna.h" // already included in EmanTreeAnalysis.h
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
//my headers
#include "functionsForAnalysis.h"
#include "myClasses.h"

using namespace std;
using namespace myAnalyzerTEman;

#ifdef vbfHiggsToInv_zllControlRegion_cxx

vbfHiggsToInv_zllControlRegion::vbfHiggsToInv_zllControlRegion(TTree *tree) : vbfHiggsToInv_LeptonControlRegion(tree) {
  //cout <<"check in constructor "<<endl;
  //edimarcoTree_v3::Init(tree);
  //vbfHiggsToInv_LeptonControlRegion::Init(tree);  // could also be just Init(tree)

}

#endif

//===============================================

void vbfHiggsToInv_zllControlRegion::setSelections() {

  vbfHiggsToInv_LeptonControlRegion::setSelections();

  gammaLooseVetoC.set("photon veto","photons veto");
  oppChargeLeptonsC.set("OS/SF lep","OS/SF leptons");

  if (!ISDATA_FLAG && (GENLEP_TAG != 0)) {
    if (using_zlljets_MCsample_flag ) {
      genLepC.set("genLep",Form("%s generated",FLAVOUR));     
      recoGenLepMatchC.set("reco-gen match","reco-gen match (DR = 0.1)","only for zlljets: looks for matching of reco and gen particles");      
    }
    if (using_ztautaujets_MCsample_flag) genTauC.set("genTauC","taus generated"); 
  }

  if (fabs(LEP_PDG_ID) == 13) {  // if we have Z -> mumu do stuff...

    invMassC.set("M_mumu",Form("mass in [%3.0lf,%3.0lf]",DILEPMASS_LOW,DILEPMASS_UP));
    twoLepLooseC.set("2 loose mu",Form("2 loose %s",FLAVOUR));
    tightLepC.set(">0 tight mu",Form(">0 tight %s",FLAVOUR));

  } else if (fabs(LEP_PDG_ID) == 11) {   // if we have Z -> ee do different stuff...

    invMassC.set("M_ee",Form("mass in [%3.0lf,%3.0lf]",DILEPMASS_LOW,DILEPMASS_UP));
    twoLepLooseC.set("2 loose ele",Form("2 loose %s",FLAVOUR));
    tightLepC.set(">0 tight ele",Form(">0 tight %s",FLAVOUR));

  }

  selection::checkMaskLength();
  selection::printActiveSelections(cout); 

}

//===============================================

void vbfHiggsToInv_zllControlRegion::setMask() {

  analysisMask.setName(Form("%s control sample (inclusive)",CONTROL_SAMPLE));

  if (!ISDATA_FLAG && (GENLEP_TAG != 0)) {     
    if (using_zlljets_MCsample_flag ) {
      analysisMask.append(genLepC.get2ToId());
      analysisMask.append(recoGenLepMatchC.get2ToId());
    }
    else if (using_ztautaujets_MCsample_flag) analysisMask.append(genTauC.get2ToId());    
  }
  if ( HLT_FLAG != 0 ) analysisMask.append(HLTC.get2ToId());
  if (MET_FILTERS_FLAG != 0) analysisMask.append(metFiltersC.get2ToId());
  analysisMask.append(twoLepLooseC.get2ToId());
  analysisMask.append(tightLepC.get2ToId());
  analysisMask.append(oppChargeLeptonsC.get2ToId());     
  analysisMask.append(invMassC.get2ToId());
  analysisMask.append(lepLooseVetoC.get2ToId());
  analysisMask.append(gammaLooseVetoC.get2ToId());
  if (TAU_VETO_FLAG) analysisMask.append(tauLooseVetoC.get2ToId());
  if (B_VETO_FLAG) analysisMask.append(bjetVetoC.get2ToId());
  analysisMask.append(vbfTaggedJets_jetsPtC.get2ToId());
  analysisMask.append(vbfTaggedJets_inVMassC.get2ToId());
  analysisMask.append(vbfTaggedJets_deltaEtaC.get2ToId());
  analysisMask.append(jetMetDphiMinC.get2ToId());
  if (METNOLEP_START != 0) analysisMask.append(recoilC.get2ToId());
   
  analysisSelectionManager.SetMaskPointer(&analysisMask);

  if ( HLT_FLAG != 0 ) analysisSelectionManager.append(&HLTC);
  if (MET_FILTERS_FLAG != 0) analysisSelectionManager.append(&metFiltersC);
  analysisSelectionManager.append(&twoLepLooseC);
  analysisSelectionManager.append(&tightLepC);
  analysisSelectionManager.append(&oppChargeLeptonsC);  
  analysisSelectionManager.append(&invMassC);
  analysisSelectionManager.append(&lepLooseVetoC);
  analysisSelectionManager.append(&gammaLooseVetoC);
  if (TAU_VETO_FLAG) analysisSelectionManager.append(&tauLooseVetoC);
  if (B_VETO_FLAG) analysisSelectionManager.append(&bjetVetoC);
  analysisSelectionManager.append(&vbfTaggedJets_jetsPtC);
  analysisSelectionManager.append(&vbfTaggedJets_inVMassC);
  analysisSelectionManager.append(&vbfTaggedJets_deltaEtaC);
  analysisSelectionManager.append(&jetMetDphiMinC);
  if (METNOLEP_START != 0) analysisSelectionManager.append(&recoilC);
   
  // creating collection of pointers to mask used in the analysis                                                                                                      
  anaMasksPtrCollection.push_back(&analysisMask);


}

//===============================================

void vbfHiggsToInv_zllControlRegion::setHistograms() {

  vbfHiggsToInv_LeptonControlRegion::setHistograms();
  
  HinvMass = new TH1D("HinvMass","",NinvMassBins,DILEPMASS_LOW,DILEPMASS_UP);    // for MC it's done on Z->mumu or Z->ee at gen level
  HzptDistribution = new TH1D("HzptDistribution","",240,0.0,1200.0); 
  Hlep2ptDistribution = new TH1D("Hlep2ptDistribution","",200,0.0,1000.0);
  Hlep2etaDistribution = new TH1D("Hlep2etaDistribution","",100,-5.0,5.0);

  if (suffix == "DYJetsToLL") {
    hasScaledHistograms_flag = 1;
    setScaleFactorHistograms();
  }

}

//===============================================

void vbfHiggsToInv_zllControlRegion::setScaleFactorHistograms() {

  vbfHiggsToInv_LeptonControlRegion::setScaleFactorHistograms();

  // vbfHiggsToInv
  HYieldsMetBin_LepTightLooseUp = new TH1D("HYieldsMetBin_LepTightLooseUp","yields in bins of met; #slash{E}_{T};# of events",nMetBins,metBinEdgesVector.data());
  HYieldsMetBin_LepTightLooseDown = new TH1D("HYieldsMetBin_LepTightLooseDown","yields in bins of met; #slash{E}_{T};# of events",nMetBins,metBinEdgesVector.data());

  HSyst_LepTightLoose = new TH1D("HSyst_LepTightLoose","systematic uncertainty for lepton tight and loose id",nMetBins,metBinEdgesVector.data());

}


//===============================================

void vbfHiggsToInv_zllControlRegion::setHistogramLastBinAsOverFlow(const Int_t hasScaledHistograms = 0) {

  vbfHiggsToInv_LeptonControlRegion::setHistogramLastBinAsOverFlow(hasScaledHistograms);

  myAddOverflowInLastBin(Hlep2ptDistribution);
  myAddOverflowInLastBin(HzptDistribution);

  if (hasScaledHistograms) {
    myAddOverflowInLastBin(HYieldsMetBin_LepTightLooseUp);
    myAddOverflowInLastBin(HYieldsMetBin_LepTightLooseDown);
  }

}

//===============================================

void vbfHiggsToInv_zllControlRegion::setNumberParameterValue(const std::string parameterName, const Double_t value) {

  vbfHiggsToInv_LeptonControlRegion::setNumberParameterValue(parameterName, value);

  if (parameterName == "LEP2PT") LEP2PT = value;
  else if (parameterName == "LEP2ETA") LEP2ETA = value;
  else if (parameterName == "DILEPMASS_LOW") DILEPMASS_LOW = value;
  else if (parameterName == "DILEPMASS_UP") DILEPMASS_UP = value;
  else if (parameterName == "HLT_LEP2PT") HLT_LEP2PT = value;
  else if (parameterName == "HLT_LEP2ETA") HLT_LEP2ETA = value;

  if (!ISDATA_FLAG) {

    if (parameterName == "GENLEP2PT") GENLEP2PT = value;
    else if (parameterName == "GENLEP2ETA") GENLEP2ETA = value;
    else if (parameterName == "GEN_ZMASS_LOW") GEN_ZMASS_LOW = value;
    else if (parameterName == "GEN_ZMASS_UP") GEN_ZMASS_UP = value;

  }

}

//===============================================

void vbfHiggsToInv_zllControlRegion::setControlSampleSpecificParameter() {

  vbfHiggsToInv_LeptonControlRegion::setControlSampleSpecificParameter(); //check if defined before uncommenting

  invMassBinWidth = 0.5;  // invariant mass histogram's bin width in GeV
  NinvMassBins = (DILEPMASS_UP - DILEPMASS_LOW) / invMassBinWidth;
  // the following flag is needed to enable search for Z->ll at generator level. For MC samples different from DYJetsToLL I must not require 2 gen leptons from Z
  // unless it is Z->tautau, in which case I start from generated taus and apply selection (tau can produce muon or electron)
  if ( !ISDATA_FLAG && ( suffix == "DYJetsToLL" )  )  using_zlljets_MCsample_flag = 1; 
  else using_zlljets_MCsample_flag = 0;
 
  if ( !ISDATA_FLAG && ( suffix == "ZJetsToTauTau" )) using_ztautaujets_MCsample_flag = 1; 
  else using_ztautaujets_MCsample_flag = 0;  
  
  if (fabs(LEP_PDG_ID) == 13) {  // if we have Z -> mumu do stuff...

    strcpy(FLAVOUR,"muons");
    strcpy(LL_FLAVOUR,"mumu");
    strcpy(CONTROL_SAMPLE,"Z-->mumu");

  } else if (fabs(LEP_PDG_ID) == 11) {   // if we have Z -> ee do different stuff...

    strcpy(FLAVOUR,"electrons");
    strcpy(LL_FLAVOUR,"ee");
    strcpy(CONTROL_SAMPLE,"Z-->ee");

  }

}

//===============================================

void vbfHiggsToInv_zllControlRegion::setVarFromConfigFile() {

  vbfHiggsToInv_LeptonControlRegion::setVarFromConfigFile();
  setControlSampleSpecificParameter();

}

//===============================================

Double_t vbfHiggsToInv_zllControlRegion::computeEventWeight() {

  if (ISDATA_FLAG || unweighted_event_flag) return 1.0;
  else {

    //    Double_t tmp = LUMI * weight * puw * SF_BTag; //SF_BTag is in evVarFriend, not sfFriend
    Double_t tmp = LUMI * weight * puw; // no SF_BTag for now

    if (suffix == "ZJetsToNuNu" || suffix == "DYJetsToLL") tmp /= 1.23;
    else if (suffix == "WJetsToLNu") tmp/= 1.21; 
   
    if (hasSFfriend_flag != 0) { 
      if (fabs(LEP_PDG_ID) == 13) tmp *= SF_trigmetnomu;
      else if (fabs(LEP_PDG_ID) == 11) tmp *= SF_trig1lep;
      return tmp * SF_NLO_QCD * SF_NLO_EWK * SF_LepTightLoose;

    } else return tmp; 

  }

}

//===============================================

void vbfHiggsToInv_zllControlRegion::createSystematicsHistogram() {

  // first step is calling the vbfHiggsToInvAna method through vbfHiggsToInv_LeptonControlRegion
  vbfHiggsToInv_LeptonControlRegion::createSystematicsHistogram();
    
  //for vbfHiggsToInv

  // computing systematic uncertainties and saving them as histograms.
  myBuildSystematicsHistogram(HSyst_LepTightLoose, HYieldsMetBin, HYieldsMetBin_LepTightLooseUp, HYieldsMetBin_LepTightLooseDown);

  // define an empty histogram to sum uncertainties in a clean way
  TH1D *Htmp = new TH1D("Htmp","",nMetBins,metBinEdgesVector.data());
  vector<TH1D*> hptr;
  hptr.push_back(HSyst_LepTightLoose);

  HSyst_total->Multiply(HSyst_total);  // before adding anything we must get the square of this histogram, because we will add histograms in quadrature

  for (UInt_t i = 0; i < hptr.size(); i++) {
    Htmp->Multiply(hptr[i],hptr[i]); // square of bin content for each single systematic histogram
    HSyst_total->Add(Htmp);             // adding the squares
  }

  for (Int_t i = 0; i <= (HSyst_total->GetNbinsX() + 1); i++) {  // computing square root of each bin's content (from underflow to overflow bin, but they should be empty)
    HSyst_total->SetBinContent(i, sqrt(HSyst_total->GetBinContent(i)));

  }

  delete Htmp;
  

}

//===============================================                                                                                                                       

void vbfHiggsToInv_zllControlRegion::fillEventMask(UInt_t & eventMask) {

  vbfHiggsToInv_LeptonControlRegion::fillEventMask(eventMask);

  eventMask += invMassC.addToMask((mZ1 > DILEPMASS_LOW) && (mZ1 < DILEPMASS_UP));  
  eventMask += oppChargeLeptonsC.addToMask(LepGood_pdgId[0] == - LepGood_pdgId[1]); 
  eventMask += twoLepLooseC.addToMask(nLepLoose > 1.5 && nLepLoose < 2.5);
  // Warning: tightLep cut is slightly different between W and Z (1 or 2 leptons), so we keep this here        
  if (fabs(LEP_PDG_ID) == 11) {
    eventMask += tightLepC.addToMask(nLepTight > 0.5 && ptr_lepton_pt[0] > LEP1PT && fabs(LepGood_pdgId[0]) == 11);
    if (HLT_FLAG != 0) eventMask += HLTC.addToMask(HLT_SingleEl == 1);
  } else {
    eventMask += tightLepC.addToMask(nLepTight > 0.5 && fabs(LepGood_pdgId[0]) == 13);
    if (HLT_FLAG != 0) eventMask += HLTC.addToMask(HLT_MonoJetMetNoMuMHT90 > 0.5 || HLT_MonoJetMetNoMuMHT120 > 0.5 || HLT_Met170 > 0.5);
  }

}


//===============================================

void vbfHiggsToInv_zllControlRegion::loop(vector< Double_t > &yRow, vector< Double_t > &eRow, vector< Double_t > &uncRow, 
				vector< Double_t > &yRow_monoJ, vector< Double_t > &eRow_monoJ, vector< Double_t > &uncRow_monoJ,
				vector< Double_t > &yRow_monoV, vector< Double_t > &eRow_monoV, vector< Double_t > &uncRow_monoV)
{

   if (fChain == 0) return;

   fChain->SetBranchStatus("*",0);  
   // warning: in Emanuele's trees non integer values are float

   //fChain->SetBranchStatus("LHEorigWeight",1); // contains negative values: the weight in the event is weight*LHEorigWeight

   //fChain->SetBranchStatus("genWeight",1); 

   fChain->SetBranchStatus("nMu10V",1);  // # of muons passing loose selection
   fChain->SetBranchStatus("nEle10V",1);  // # of electrons passing loose selection for electron veto
   fChain->SetBranchStatus("nGamma15V",1);  // # of photons passing loose selection for photon veto
   fChain->SetBranchStatus("nMu20T",1);  // # of muons passing tight selection (isolation included)
   //added on 23/01/2016
   fChain->SetBranchStatus("nEle40T",1);   

   fChain->SetBranchStatus("nTauClean18V",1);

   fChain->SetBranchStatus("dphijj",1);          // dphi between 1st and 2nd jet, 999 if second jet doesn't exist
   fChain->SetBranchStatus("nJetClean30",1);    // # of jet with pt > 30 & eta < 2.5 and cleaning for against muons misidentified as PFjets   
   fChain->SetBranchStatus("JetClean_pt",1);  
   fChain->SetBranchStatus("JetClean_eta",1); 

   fChain->SetBranchStatus("dphijmAllJets",1);
   fChain->SetBranchStatus("vbfTaggedJet_deltaEta",1);
   fChain->SetBranchStatus("vbfTaggedJet_invMass",1);
   fChain->SetBranchStatus("vbfTaggedJet_leadJetPt",1);
   fChain->SetBranchStatus("vbfTaggedJet_trailJetPt",1);
   fChain->SetBranchStatus("vbfTaggedJet_leadJetEta",1);
   fChain->SetBranchStatus("vbfTaggedJet_trailJetEta",1);
 
   fChain->SetBranchStatus("nLepGood",1);
   fChain->SetBranchStatus("LepGood_pdgId",1);  // must be 13 for muons ( -13 for mu+), 11 for electrons and 15 for taus
   fChain->SetBranchStatus("LepGood_pt",1);
   fChain->SetBranchStatus("LepGood_eta",1);
   fChain->SetBranchStatus("LepGood_phi",1);   
   fChain->SetBranchStatus("LepGood_mass",1);
   //fChain->SetBranchStatus("LepGood_charge",1);
   //   fChain->SetBranchStatus("LepGood_tightId",1);
   //fChain->SetBranchStatus("LepGood_relIso04",1);
   // fChain->SetBranchStatus("ngenLep",1);         // not needed, using GenPart to study generator level quantities
   // fChain->SetBranchStatus("genLep_pdgId",1);
   // fChain->SetBranchStatus("genLep_pt",1);
   // fChain->SetBranchStatus("genLep_eta",1);
   // fChain->SetBranchStatus("genLep_phi",1);
   fChain->SetBranchStatus("mZ1",1);  // best m(ll) SF/OS

   fChain->SetBranchStatus("met_pt",1);
   //fChain->SetBranchStatus("met_eta",1);
   fChain->SetBranchStatus("met_phi",1);

   fChain->SetBranchStatus("metNoMu_pt",1);
   //fChain->SetBranchStatus("metNoMu_eta",1);
   fChain->SetBranchStatus("metNoMu_phi",1);
   fChain->SetBranchStatus("htJet25",1);

   fChain->SetBranchStatus("nVert",1);  // number of good vertices 

   fChain->SetBranchStatus("HLT_MonoJetMetNoMuMHT90",1);
   fChain->SetBranchStatus("HLT_MonoJetMetNoMuMHT120",1);
   fChain->SetBranchStatus("HLT_Met170",1); 
   fChain->SetBranchStatus("HLT_SingleEl",1);
   //fChain->SetBranchStatus("HLT_SinglePho",1); 

   // met filters to be used (the config file has a parameter saying whether they should be used or not)
   // fChain->SetBranchStatus("cscfilter",1);
   // fChain->SetBranchStatus("ecalfilter",1);
   // fChain->SetBranchStatus("hbheFilterNew25ns",1);
   // fChain->SetBranchStatus("hbheFilterIso",1);
   // fChain->SetBranchStatus("Flag_eeBadScFilter",1);

   //added on November 2015. These are new variables (except for weight, which has just changed in the definition)
   fChain->SetBranchStatus("nBTag15",1);  // for b-jet veto
   fChain->SetBranchStatus("dphijm",1);          // flag for dphi minimum between met and any of the jets in the event (using only the first four jets
   fChain->SetBranchStatus("weight",1);   // modified since 17 November 2015: now it includes the whol weight, e.g. 1000*xsec*genWeight ...
   fChain->SetBranchStatus("events_ntot",1);    // equivalent to SUMWEIGHTS for samples before 17 November 2015
   fChain->SetBranchStatus("JetClean_leadClean",1); // has new cleaning on energy fractions (added on 17 November 2015) 

   // // For mon-V categhory the following variables are needed.  -- > WARNING: this collection was made with |eta| < 2.4, not 2.5
   // fChain->SetBranchStatus("nFatJetClean",1);             // at least one for mono-V
   // fChain->SetBranchStatus("FatJetClean_pt",1);           // leading jet is required to be > 250
   // fChain->SetBranchStatus("FatJetClean_eta",1);          // just for the histogram
   // fChain->SetBranchStatus("FatJetClean_prunedMass",1);   // in 65-105 for leading jet in V-tag
   // fChain->SetBranchStatus("FatJetClean_tau1",1);         // tau2/tau1 < 0.6 (I guess for the leading jet)
   // fChain->SetBranchStatus("FatJetClean_tau2",1);

   if (!ISDATA_FLAG) {
     fChain->SetBranchStatus("nGenPart",1);
     fChain->SetBranchStatus("GenPart_pdgId",1);
     fChain->SetBranchStatus("GenPart_motherId",1);
     fChain->SetBranchStatus("GenPart_pt",1);
     fChain->SetBranchStatus("GenPart_eta",1);
     fChain->SetBranchStatus("GenPart_phi",1);
     fChain->SetBranchStatus("GenPart_mass",1);
     fChain->SetBranchStatus("GenPart_motherIndex",1);

     fChain->SetBranchStatus("puw",1); // added on 21 Sept 2016, substituting vtxWeight

     //added on 23/01/2016
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
     //fChain->SetBranchStatus("SF_NLO",1);
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

   //cout << "CHECK 1 " << endl;
   setVarFromConfigFile();
   //cout << "CHECK 2 " << endl;
   setSelections();
   //cout << "CHECK 3 " << endl;
   setMask();
   //cout << "CHECK 4 " << endl;

   //  TVector2 metNoLepTV, ele;
 
   TLorentzVector l1gen, l2gen, Zgen;     // gen level Z and l1,l2  (Z->(l1 l2)
   TLorentzVector l1reco, l2reco, Zreco;

   // following indices refer to the leading pair of OS/SF in the list of LepGood. They are initialized with 0 and 1 by default, but can be set with function
   // myGetPairIndexInArray (see functionsForAnalysis.cc for reference). 
   // When myGetPairIndexInArray() is called, the index of "correct" particles will be used. If they are not found (e.g. a pair of OS/SF is mismeasured as 2 mu+), 
   // indices are set as 0 and 1 (and following selection asking lep[0] and lep[1] to be OS or whatever will fail).
   // Int_t firstIndex = 0;
   // Int_t secondIndex = 1;

   Int_t firstIndexGen = 0;
   Int_t secondIndexGen = 1;
   Int_t recoLepFound_flag = 0;
   Int_t genLepFound_flag = 0;
   //Int_t recoGenMatchDR_flag = 0;
   Int_t genTauFound_flag = 0;
   Int_t Z_index = 0; 

   //Int_t HLT_passed_flag = 1;          // some computations (for e) require a trigger preselection, while others don't. The former will be done if the flag is set to 1
                                                       // it's set to 1 because if the trigger selection is not applied every event must be considered to be a "good" event having passed all preselections
                                                       // Actually in this code the trigger is not always necessary, but I keep it like this nonetheless.

   if (fabs(LEP_PDG_ID) == 13) {  // if we have Z -> mumu do stuff...
  
     ptr_nLepLoose = &nMu10V;                      // ask 2 muons
     ptr_nLep10V = &nEle10V;                         // veto on electrons
     ptr_nLepTight = &nMu20T;                                            
     ptr_metNoLepPt = &metNoMu_pt;               // for muons  get this variable from the tree 
     //ptr_metNoLepEta = &metNoMu_eta;               // for muons  get this variable from the tree 
     ptr_metNoLepPhi = &metNoMu_phi;         // for muons  get this variable from the tree

   } else if (fabs(LEP_PDG_ID) == 11) {   // if we have Z -> ee do different stuff...

     ptr_nLepLoose = &nEle10V;                      // ask 2 electrons
     ptr_nLep10V = &nMu10V;                         // veto on muons  
     //if (fChain->GetBranch("nEle40T")) ptr_nLepTight = &nEle40T;
     //else if (fChain->GetBranch("nEle20T")) ptr_nLepTight = &nEle20T;  // the most recent version for this variable si nEle40T, but older versions use nEle20T
     // or alternatively use the following methods
     //if (fChain->GetListOfBranches()->FindObject("nEle40T"))
     ptr_nLepTight = &nEle40T;
   
   }

   ptr_nRecoLepton = &nLepGood;
   ptr_lepton_pt = LepGood_pt;
   ptr_lepton_eta = LepGood_eta;
   ptr_lepton_phi = LepGood_phi;
   ptr_lepton_mass = LepGood_mass;


   
   cout << "Opening file " <<ROOT_FNAME<< " in folder " << outputFolder << endl;

   TFile *rootFile = new TFile((outputFolder + ROOT_FNAME).c_str(),"RECREATE");
   if (!rootFile || !rootFile->IsOpen()) {
     cout<<"Error: file \""<<outputFolder + ROOT_FNAME<<"\" was not opened."<<endl;
     exit(EXIT_FAILURE);
   }
 

   TH1::SetDefaultSumw2();            //all the following histograms will automatically call TH1::Sumw2() 
   //TH1::StatOverflows();                 //enable use of underflows and overflows for statistics computation 
   TVirtualFitter::SetDefaultFitter("Minuit");

   setHistograms();

   TH1D *HzlljetsInvMassMetBin[nMetBins];

   for (Int_t i = 0; i < nMetBins; i++) {

     HzlljetsInvMassMetBin[i] = new TH1D(Form("HzlljetsInvMassMetBin_met%2.0lfTo%2.0lf",metBinEdgesVector[i],metBinEdgesVector[i+1]),"",NinvMassBins,DILEPMASS_LOW,DILEPMASS_UP);

   } 

   Long64_t nentries = fChain->GetEntriesFast();
   cout<<"vbfHiggsToInv_zllControlRegion::loop()"<<endl;
   cout<<"nentries = "<<nentries<<endl;   

   Long64_t nbytes = 0, nb = 0;

   for (Int_t jentry=0; jentry<nentries; jentry++) {

     Long64_t ientry = LoadTree(jentry);
     if (ientry < 0) break;
     nb = fChain->GetEntry(jentry);   nbytes += nb;
     // if (Cut(ientry) < 0) continue;   

     if (jentry%500000 == 0) cout << jentry << endl;

     UInt_t eventMask = 0; 

     //cout << "CHECK -2 in zll"<< endl;

     newwgt = computeEventWeight();
     nTotalWeightedEvents += newwgt;  // counting events with weights
     
     //cout << "CHECK -1 in zll"<< endl;

     nLepLoose = *ptr_nLepLoose;
     nLep10V = *ptr_nLep10V;
     nLepTight = *ptr_nLepTight;
     nRecoLepton = *ptr_nRecoLepton;

     //cout << "CHECK in zll"<< endl;

     //Double_t ZgenMass;        // not used for now
     //Double_t ZtoLLGenPt = 0;    // filled below (only if running on MC DYJetsToLL)
     Double_t ZtoLLRecoPt = 0;   // filled below

     // genLepFound_flag is used when analysing DYJetsToLL in MC fo Z->mumu or Z->ee. For other MC samples it's not used.

     if (!ISDATA_FLAG && (GENLEP_TAG != 0)) {

       if (using_zlljets_MCsample_flag) {

	 //genLepFound_flag = 1;
	 genLepFound_flag = myPartGenAlgo(nGenPart, GenPart_pdgId, GenPart_motherId, LEP_PDG_ID, 23, firstIndexGen, secondIndexGen, Z_index, GenPart_motherIndex); 
     
	 if (genLepFound_flag != 0) {
	   eventMask += genLepC.addToMask(1);
	   l1gen.SetPtEtaPhiM(GenPart_pt[firstIndexGen],GenPart_eta[firstIndexGen],GenPart_phi[firstIndexGen],GenPart_mass[firstIndexGen]);
	   l2gen.SetPtEtaPhiM(GenPart_pt[secondIndexGen],GenPart_eta[secondIndexGen],GenPart_phi[secondIndexGen],GenPart_mass[secondIndexGen]);
	   Zgen = l1gen + l2gen;
	   //ZtoLLGenPt = Zgen.Pt();                             
	 }

       } else if (using_ztautaujets_MCsample_flag) {

	 genTauFound_flag = myPartGenAlgo(nGenPart, GenPart_pdgId, GenPart_motherId, 15, 23);
	 eventMask += genTauC.addToMask( genTauFound_flag );

       }

     }

     // look if the first two good leptons are OS/SF (flavour depending on config file)
     // the check for 2 OS/SF leptons in LepGood list is just useful to speed up things, but would not be necessary if other variables or computations that require this condition are only used at the end of selection (which already includes 2 OS/SF condition).
     // thus, we evaluate this right now and use this piece of information for later use

     //cout << "CHECK 2 in zll"<< endl;

     if ( (nRecoLepton >= 2) && (fabs(LepGood_pdgId[0]) ==  LEP_PDG_ID) && (LepGood_pdgId[0] == (-1) * LepGood_pdgId[1]) ) recoLepFound_flag = 1;
     else recoLepFound_flag = 0;

     if (recoLepFound_flag) {
       l1reco.SetPtEtaPhiM(ptr_lepton_pt[0],ptr_lepton_eta[0],ptr_lepton_phi[0],ptr_lepton_mass[0]);
       l2reco.SetPtEtaPhiM(ptr_lepton_pt[1],ptr_lepton_eta[1],ptr_lepton_phi[1],ptr_lepton_mass[1]);
       Zreco = l1reco + l2reco;
       ZtoLLRecoPt = Zreco.Pt();
     }

     if (fabs(LEP_PDG_ID) == 13) { 

       // if ( HLT_FLAG != 0) {

       // 	 if ( HLT_MonoJetMetNoMuMHT90 > 0.5 || HLT_MonoJetMetNoMuMHT120 > 0.5 || HLT_Met170 > 0.5) HLT_passed_flag = 1; 	 
       // 	 else HLT_passed_flag = 0; //continue;

       // }  // end of   if ( HLT_FLAG )

       metNoLepPt = *ptr_metNoLepPt;       
       metNoLepPhi = *ptr_metNoLepPhi; 

     } else if (fabs(LEP_PDG_ID) == 11) { 

       // if ( HLT_FLAG != 0 ) {

       // 	 if (HLT_SingleEl == 1) HLT_passed_flag = 1; 
       // 	 else HLT_passed_flag = 0;  //continue;

       // }  // end of   if ( HLT_FLAG )

       metNoLepTV.SetMagPhi(met_pt,met_phi);
       // summing just electrons from Z if found
       if (recoLepFound_flag) {
	 ele.SetMagPhi(ptr_lepton_pt[0],ptr_lepton_phi[0]);
	 metNoLepTV += ele;
	 ele.SetMagPhi(ptr_lepton_pt[1],ptr_lepton_phi[1]);
	 metNoLepTV += ele;
       }

       metNoLepPt = metNoLepTV.Mod();

     }


     //cout << "CHECK 3 in zll"<< endl;

     // beginning of eventMask building
     
     // genLepC added to mask above if ISDATA_FLAG == false (in order not to repeat here the check) 
     
     fillEventMask(eventMask); 
 
     // end of eventMask building

     // test matching of reco and gen lep for DY MC 
     if (using_zlljets_MCsample_flag && (GENLEP_TAG != 0)) {
    
       // now we require a DeltaR cut between gen-level and reco-level leptons (mu or e) from Z in Drell-Yan events to assess that lreco comes from lgen
       
       // since genLepFound_flag && recoLepFound_flag require 2 OS/SF leptons (with correct flavour) to be found in their lists, we have two possible cases:
       // 1) the first gen lepton and the first reco lepton have same charge
       // 2) the first gen lepton and the first reco lepton have opposite charge
       // if the first case is not satisfied, then the second will automatically be

       if (genLepFound_flag && recoLepFound_flag) {       

	 Double_t DeltaR_lreco_lgen_pair1 = 100.0;
	 Double_t DeltaR_lreco_lgen_pair2 = 100.0;
       
	 if(LepGood_pdgId[0] == GenPart_pdgId[firstIndexGen] && LepGood_pdgId[1] == GenPart_pdgId[secondIndexGen]) {
	 
	   DeltaR_lreco_lgen_pair1 = l1reco.DeltaR(l1gen);
	   DeltaR_lreco_lgen_pair2 = l2reco.DeltaR(l2gen);

	 } else {
	 
	   DeltaR_lreco_lgen_pair1 = l1reco.DeltaR(l2gen);
	   DeltaR_lreco_lgen_pair2 = l2reco.DeltaR(l1gen);
	 
	 }
       
	 if (DeltaR_lreco_lgen_pair1 < 0.1 && DeltaR_lreco_lgen_pair2 < 0.1) eventMask += recoGenLepMatchC.addToMask(1);
	 else eventMask += recoGenLepMatchC.addToMask(0);

       }

     }

     analysisMask.countEvents(eventMask,newwgt);

     // cout << "check 3 in zll" << endl;

     if ( ((eventMask & analysisMask.globalMask.back()) == analysisMask.globalMask.back()) ) {
       
       // this histogram holds the final yields in bins of MET
       HYieldsMetBin->Fill(metNoLepPt,newwgt);
	 
       HhtDistribution->Fill(htJet25,newwgt);
       HinvMass->Fill(mZ1,newwgt);
       HrecoilDistribution->Fill(metNoLepPt,newwgt);
       HzptDistribution->Fill(ZtoLLRecoPt,newwgt);
       HvtxDistribution->Fill(nVert,newwgt);
       HnjetsDistribution->Fill(nJetClean30,newwgt);
       HjetMetDphiMinDistribution->Fill(dphijmAllJets,newwgt);
       Hlep1ptDistribution->Fill(ptr_lepton_pt[0],newwgt);
       Hlep1etaDistribution->Fill(ptr_lepton_eta[0],newwgt);
       Hlep2ptDistribution->Fill(ptr_lepton_pt[1],newwgt);
       Hlep2etaDistribution->Fill(ptr_lepton_eta[1],newwgt);
       Hjet1etaDistribution->Fill(vbfTaggedJet_leadJetEta,newwgt);
       Hjet1ptDistribution->Fill(vbfTaggedJet_leadJetPt,newwgt);
       Hjet2etaDistribution->Fill(vbfTaggedJet_trailJetEta,newwgt);
       Hjet2ptDistribution->Fill(vbfTaggedJet_trailJetPt,newwgt);
       HvbfTaggedJets_deltaEta->Fill(vbfTaggedJet_deltaEta);
       HvbfTaggedJets_invMass->Fill(vbfTaggedJet_invMass);
       if (nJetClean30 > 1) Hj1j2dphiDistribution->Fill(dphijj,newwgt);
       	 
       if (hasScaledHistograms_flag) {

	 HYieldsMetBin_qcdRenScaleUp->Fill(metNoLepPt,(newwgt * SF_NLO_QCD_renScaleUp/ SF_NLO_QCD));
	 HYieldsMetBin_qcdRenScaleDown->Fill(metNoLepPt,(newwgt * SF_NLO_QCD_renScaleDown/ SF_NLO_QCD));
	 HYieldsMetBin_qcdFacScaleUp->Fill(metNoLepPt,(newwgt * SF_NLO_QCD_facScaleUp/ SF_NLO_QCD));
	 HYieldsMetBin_qcdFacScaleDown->Fill(metNoLepPt,(newwgt * SF_NLO_QCD_facScaleDown/ SF_NLO_QCD));
	 HYieldsMetBin_qcdPdfUp->Fill(metNoLepPt,(newwgt * SF_NLO_QCD_pdfUp/ SF_NLO_QCD));
	 HYieldsMetBin_qcdPdfDown->Fill(metNoLepPt,(newwgt * SF_NLO_QCD_pdfDown/ SF_NLO_QCD));
	 HYieldsMetBin_ewkUp->Fill(metNoLepPt,(newwgt * SF_NLO_EWK_up/ SF_NLO_EWK));
	 HYieldsMetBin_ewkDown->Fill(metNoLepPt,(newwgt * SF_NLO_EWK_down/ SF_NLO_EWK));
	 HYieldsMetBin_LepTightLooseUp->Fill(metNoLepPt,(newwgt * SF_LepTightLooseUp / SF_LepTightLoose));
	 HYieldsMetBin_LepTightLooseDown->Fill(metNoLepPt,(newwgt * SF_LepTightLooseDown / SF_LepTightLoose));

       }
	
     }

     // now entering analysis in bins of met

     if ((metNoLepPt > metBinEdgesVector.front()) && (metNoLepPt < metBinEdgesVector.back())) {

       Int_t bin = myGetBin(metNoLepPt,metBinEdgesVector.data(),nMetBins);

       if ( ((eventMask & analysisMask.globalMask.back()) == analysisMask.globalMask.back()) ) { 
 
	 HzlljetsInvMassMetBin[bin]->Fill(mZ1,newwgt); 

       }
	 
     }                      // end of    if ((metNoLepPt > metBinEdges[0]) && (metNoLepPt < metBinEdges[nMetBins])) 
       
   }                        // end of loop on entries

   mySpaces(cout,2);
   selection::printSelectionFlowAndYields(cout, LUMI, nTotalWeightedEvents, &analysisMask);

   mySpaces(cout,2);
   myPrintYieldsMetBinInStream(cout, HYieldsMetBin, metBinEdgesVector.data(), nMetBins);
   mySpaces(cout,2);
 
   cout<<"creating file '"<<TXT_FNAME<<"' in folder "<< outputFolder <<" ..."<<endl;

   ofstream myfile((outputFolder + TXT_FNAME).c_str(),ios::out);

   if ( !myfile.is_open() ) {

     cout<<"Error: unable to open file "<<TXT_FNAME<<" !"<<endl;
     exit(EXIT_FAILURE);
     
   }


   if (!ISDATA_FLAG && unweighted_event_flag) myfile << "======   Using unweighted events (w = 1)   ======" << endl;
   mySpaces(myfile,3);
   selection::printSelectionFlowAndYields(myfile, LUMI, nTotalWeightedEvents, &analysisMask);
   mySpaces(myfile,3);
   myPrintYieldsMetBinInStream(myfile, HYieldsMetBin, metBinEdgesVector.data(), nMetBins);

   myfile.close();


   if ((GENLEP_TAG != 0) && (using_zlljets_MCsample_flag == 1 || using_ztautaujets_MCsample_flag == 1)) {
     fillRowVector(nTotalWeightedEvents, analysisSelectionManager, analysisMask, yRow, eRow, uncRow,1);
     // fillRowVector(nTotalWeightedEvents, analysisSelectionManager_monoJ, analysisMask_monoJ, yRow_monoJ, eRow_monoJ, uncRow_monoJ,1);
     // fillRowVector(nTotalWeightedEvents, analysisSelectionManager_monoV, analysisMask_monoV, yRow_monoV, eRow_monoV, uncRow_monoV,1);
   } else {
     fillRowVector(nTotalWeightedEvents, analysisSelectionManager, analysisMask, yRow, eRow, uncRow,0);
     // fillRowVector(nTotalWeightedEvents, analysisSelectionManager_monoJ, analysisMask_monoJ, yRow_monoJ, eRow_monoJ, uncRow_monoJ,0);
     // fillRowVector(nTotalWeightedEvents, analysisSelectionManager_monoV, analysisMask_monoV, yRow_monoV, eRow_monoV, uncRow_monoV,0);
   }

   // fill last bin with overflow 
   setHistogramLastBinAsOverFlow(hasScaledHistograms_flag);
   if (hasScaledHistograms_flag) createSystematicsHistogram();

   rootFile->Write();

   rootFile->Close();
   delete rootFile;

   //creating a .tex file to build tables with data
   //myCreateTexTable(TEX_FNAME, outputFolder, LUMI,nTotalWeightedEvents, &analysisMask);
   myCreateTexTable(TEX_FNAME, outputFolder, LUMI,nTotalWeightedEvents, anaMasksPtrCollection);

   // end of tex file

   mySpaces(cout,2);

}




