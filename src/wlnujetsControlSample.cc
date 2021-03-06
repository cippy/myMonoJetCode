#define wlnujetsControlSample_cxx
#include "EmanTreeAnalysis.h"
//#include "AnalysisDarkMatter.h" // already included in EmanTreeAnalysis.h
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

#ifdef wlnujetsControlSample_cxx

wlnujetsControlSample::wlnujetsControlSample(TTree *tree) : monojet_LeptonControlRegion(tree) {
  //cout <<"check in constructor "<<endl;
  //edimarcoTree_v3::Init(tree);
  // suffix = "";
  // uncertainty = "";
  // configFileName = NULL;
  // ISDATA_FLAG = 0;
  // unweighted_event_flag = 0;
  // hasSFfriend_flag = 0;
  //AnalysisDarkMatter::Init(tree);  // could also be just Init(tree)

}

#endif

//===============================================

void wlnujetsControlSample::setSelections() {

  monojet_LeptonControlRegion::setSelections();

  gammaLooseVetoC.set("photon veto","photons veto");

  if (!ISDATA_FLAG && (GENLEP_TAG != 0)) {
    if (using_wlnujets_MCsample_flag) {
      genLepC.set("genLep",Form("%s generated",FLAVOUR));     
      recoGenLepMatchC.set("reco-gen match","reco-gen match (DR = 0.1)","only for wlnujets: looks for matching of reco and gen particles");      
    }

    if (using_wtaunujets_MCsample_flag) genTauC.set("genTau","tau generated"); 
  }

  if (fabs(LEP_PDG_ID) == 13) {  // if we have Z -> mumu do stuff...

    oneLepLooseC.set("1 loose mu",Form("1 loose %s",FLAVOUR));
    tightLepC.set("1 tight mu",Form("1 tight %s",FLAVOUR));

  } else if (fabs(LEP_PDG_ID) == 11) {   // if we have Z -> ee do different stuff...

    oneLepLooseC.set("1 loose ele",Form("1 loose %s",FLAVOUR));
    tightLepC.set("1 tight ele",Form("1 tight %s",FLAVOUR));
    metC.set("met","met > 50");
  
  }

  selection::checkMaskLength();
  selection::printActiveSelections(cout); 

}

//===============================================

void wlnujetsControlSample::setMask() {

  analysisMask.setName(Form("%s control sample (inclusive)",CONTROL_SAMPLE));

   if (!ISDATA_FLAG && (GENLEP_TAG != 0)) {
     if (using_wlnujets_MCsample_flag) {
       analysisMask.append(genLepC.get2ToId());
       analysisMask.append(recoGenLepMatchC.get2ToId());
     }
     if (using_wtaunujets_MCsample_flag) analysisMask.append(genTauC.get2ToId());
   }
   if ( HLT_FLAG != 0 ) analysisMask.append(HLTC.get2ToId());
   if (MET_FILTERS_FLAG != 0) analysisMask.append(metFiltersC.get2ToId());
   analysisMask.append(oneLepLooseC.get2ToId());
   analysisMask.append(tightLepC.get2ToId());
   analysisMask.append(lepLooseVetoC.get2ToId());
   if (TAU_VETO_FLAG) analysisMask.append(tauLooseVetoC.get2ToId());
   analysisMask.append(gammaLooseVetoC.get2ToId());
   analysisMask.append(bjetVetoC.get2ToId());
   if (fabs(LEP_PDG_ID) == 11) analysisMask.append(metC.get2ToId());
   if (METNOLEP_START != 0) analysisMask.append(recoilC.get2ToId());
   analysisMask.append(jet1C.get2ToId());
   analysisMask.append(jetNoiseCleaningC.get2ToId());
   analysisMask.append(jetMetDphiMinC.get2ToId());

   analysisSelectionManager.SetMaskPointer(&analysisMask);

   if ( HLT_FLAG != 0 ) analysisSelectionManager.append(&HLTC);
   if (MET_FILTERS_FLAG != 0) analysisSelectionManager.append(&metFiltersC);
   analysisSelectionManager.append(&oneLepLooseC);
   analysisSelectionManager.append(&tightLepC);
   analysisSelectionManager.append(&lepLooseVetoC);
   if (TAU_VETO_FLAG) analysisSelectionManager.append(&tauLooseVetoC);
   analysisSelectionManager.append(&gammaLooseVetoC);
   analysisSelectionManager.append(&bjetVetoC);
   if (fabs(LEP_PDG_ID) == 11) analysisSelectionManager.append(&metC);
   if (METNOLEP_START != 0) analysisSelectionManager.append(&recoilC);
   analysisSelectionManager.append(&jet1C);
   analysisSelectionManager.append(&jetNoiseCleaningC);
   analysisSelectionManager.append(&jetMetDphiMinC);

   // ========== Mono-J ==============

   analysisMask_monoJ.setName(Form("%s control sample",CONTROL_SAMPLE));
   
   analysisMask_monoJ.append(analysisMask.globalMask.back()); // all the common selections
   analysisMask_monoJ.append(noVtagC.get2ToId());
   
   analysisSelectionManager_monoJ.SetMaskPointer(&analysisMask_monoJ);
   
   analysisSelectionManager_monoJ.append("all cuts");
   analysisSelectionManager_monoJ.append(&noVtagC);

   // ========== Mono-V ==============

  analysisMask_monoV.setName("monoV signal selection");
  
  analysisMask_monoV.append(analysisMask.globalMask.back()); // all the common selections
  analysisMask_monoV.append(VtagC.get2ToId());

  analysisSelectionManager_monoV.SetMaskPointer(&analysisMask_monoV);

  analysisSelectionManager_monoV.append("all cuts");
  analysisSelectionManager_monoV.append(&VtagC);
  
  // creating collection of pointers to mask used in the analysis
  anaMasksPtrCollection.push_back(&analysisMask);
  anaMasksPtrCollection.push_back(&analysisMask_monoJ);
  anaMasksPtrCollection.push_back(&analysisMask_monoV);

}

//===============================================

void wlnujetsControlSample::setHistograms() {

  monojet_LeptonControlRegion::setHistograms();
    
  HtransverseMass = new TH1D("HtransverseMass","",40,0.0,200.0);
  HtransverseMass_monoV = new TH1D("HtransverseMass_monoV","",40,0.0,200.0);
  
  if (suffix == "WJetsToLNu") {
    hasScaledHistograms_flag = 1;
    setScaleFactorHistograms();
  }

}

//===============================================

void wlnujetsControlSample::setHistogramLastBinAsOverFlow(const Int_t hasScaledHistograms = 0) {

  monojet_LeptonControlRegion::setHistogramLastBinAsOverFlow(hasScaledHistograms);

  myAddOverflowInLastBin(HtransverseMass);

  myAddOverflowInLastBin(HtransverseMass_monoV);

}

//===============================================

void wlnujetsControlSample::setNumberParameterValue(const std::string parameterName, const Double_t value) {

  monojet_LeptonControlRegion::setNumberParameterValue(parameterName, value);

  // N.B.: if there are more checks, the first must be an "if", the successive can be "else if" or "if" (but "else if" is better because when parameter is ound the query ends

  //if (parameterName == "<PAR_NAME>") PAR_NAME = value;

  // if (!ISDATA_FLAG) {

  //   if (parameterName == "<PAR_NAME>") PAR_NAME = value;

  // }

}

//===============================================

void wlnujetsControlSample::setControlSampleSpecificParameter() {

  monojet_LeptonControlRegion::setControlSampleSpecificParameter(); //check if defined before uncommenting

  // the following flag is needed to enable search for W->lnu at generator level. For MC samples different from WJetsToLNu I must not require 2 gen leptons from Z
  // unless it is Z->tautau, in which case I start from generated taus and apply selection (tau can produce muon or electron)
  if ( !ISDATA_FLAG && ( suffix == "WJetsToLNu" )  )  using_wlnujets_MCsample_flag = 1; 
  else using_wlnujets_MCsample_flag = 0;
 
  if ( !ISDATA_FLAG && ( suffix == "WJetsToTauNu" )) using_wtaunujets_MCsample_flag = 1; 
  else using_wtaunujets_MCsample_flag = 0;  

  if (fabs(LEP_PDG_ID) == 13) {  // if we have Z -> mumu do stuff...

    strcpy(FLAVOUR,"muon");
    strcpy(LL_FLAVOUR,"munu");
    strcpy(CONTROL_SAMPLE,"W-->munu");

  } else if (fabs(LEP_PDG_ID) == 11) {   // if we have Z -> ee do different stuff...

    strcpy(FLAVOUR,"electron");
    strcpy(LL_FLAVOUR,"enu");
    strcpy(CONTROL_SAMPLE,"W-->enu");

  }
  
}

//===============================================

void wlnujetsControlSample::setVarFromConfigFile() {

  monojet_LeptonControlRegion::setVarFromConfigFile();
  setControlSampleSpecificParameter();

}

//===============================================

Double_t wlnujetsControlSample::computeEventWeight() {

  if (ISDATA_FLAG || unweighted_event_flag) return 1.0;
  else {

    Double_t tmp = LUMI * weight * vtxWeight * SF_BTag; //SF_BTag is in evVarFriend, not sfFriend

    if (hasSFfriend_flag != 0) { 

      // sf_nlo_weight = (*ptr_sf_nlo_QCD) * (*ptr_sf_nlo_EWK);
      // if (fabs(LEP_PDG_ID) == 13) tmp *= SF_trigmetnomu;
      // else if (fabs(LEP_PDG_ID) == 11) tmp *=  SF_trig1lep;
      // return tmp * sf_nlo_weight * SF_LepTight;
      if (fabs(LEP_PDG_ID) == 13) tmp *= SF_trigmetnomu;
      else if (fabs(LEP_PDG_ID) == 11) tmp *=  SF_trig1lep;
      if (suffix == "ZJetsToNuNu" || suffix == "DYJetsToLL") tmp /= 1.23;
      else if (suffix == "WJetsToLNu") tmp/= 1.21; 
      return tmp * SF_NLO_QCD * SF_NLO_EWK * SF_LepTight;


    } else return tmp; 

  }

}

//===============================================

void wlnujetsControlSample::loop(vector< Double_t > &yRow, vector< Double_t > &eRow, vector< Double_t > &uncRow, 
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
   fChain->SetBranchStatus("nEle20T",1);  // # of electrons passing tight selection (isolation included)
   //fChain->SetBranchStatus("nTau18V",1);
   fChain->SetBranchStatus("nTauClean18V",1);

   fChain->SetBranchStatus("dphijj",1);          // dphi between 1st and 2nd jet, 999 if second jet doesn't exist
   fChain->SetBranchStatus("nJetClean30",1);    // # of jet with pt > 30 & eta < 2.5 and cleaning for against muons misidentified as PFjets   
   fChain->SetBranchStatus("JetClean_pt",1);  
   fChain->SetBranchStatus("JetClean_eta",1); 
 
   fChain->SetBranchStatus("nLepGood",1);
   fChain->SetBranchStatus("LepGood_pdgId",1);  // must be 13 for muons ( -13 for mu+), 11 for electrons and 15 for taus
   fChain->SetBranchStatus("LepGood_pt",1);
   fChain->SetBranchStatus("LepGood_eta",1);
   fChain->SetBranchStatus("LepGood_phi",1);   
   fChain->SetBranchStatus("LepGood_mass",1);
   //fChain->SetBranchStatus("LepGood_charge",1);
   fChain->SetBranchStatus("LepGood_tightId",1);
   fChain->SetBranchStatus("LepGood_relIso04",1);
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
   fChain->SetBranchStatus("HLT_SingleEl",1);
 
   // met filters to be used (the config file has a parameter saying whether they should be used or not)
   fChain->SetBranchStatus("cscfilter",1);
   fChain->SetBranchStatus("ecalfilter",1);
   fChain->SetBranchStatus("hbheFilterNew25ns",1);
   fChain->SetBranchStatus("hbheFilterIso",1);
   fChain->SetBranchStatus("Flag_eeBadScFilter",1);

   //added on November 2015. These are new variables (except for weight, which has just changed in the definition)
   fChain->SetBranchStatus("nBTag15",1);  // for b-jet veto
   fChain->SetBranchStatus("dphijm",1);          // flag for dphi minimum between met and any of the jets in the event (using only the first four jets
   fChain->SetBranchStatus("weight",1);   // modified since 17 November 2015: now it includes the whol weight, e.g. 1000*xsec*genWeight ...
   fChain->SetBranchStatus("events_ntot",1);    // equivalent to SUMWEIGHTS for samples before 17 November 2015
   fChain->SetBranchStatus("JetClean_leadClean",1); // has new cleaning on energy fractions (added on 17 November 2015) 

   // For mon-V categhory the following variables are needed.  -- > WARNING: this collection was made with |eta| < 2.4, not 2.5
   fChain->SetBranchStatus("nFatJet",1);             // at least one for mono-V
   fChain->SetBranchStatus("FatJet_pt",1);           // leading jet is required to be > 250
   fChain->SetBranchStatus("FatJet_eta",1);          // just for the histogram
   fChain->SetBranchStatus("FatJet_mass",1);
   fChain->SetBranchStatus("FatJet_prunedMass",1);   // in 65-105 for leading jet in V-tag
   fChain->SetBranchStatus("FatJet_tau1",1);         // tau2/tau1 < 0.6 (I guess for the leading jet)
   fChain->SetBranchStatus("FatJet_tau2",1);

   //added on 23/01/2016
   fChain->SetBranchStatus("nEle40T",1);   
   fChain->SetBranchStatus("nCalibEle",1);
   fChain->SetBranchStatus("CalibEle_pt",1);
   fChain->SetBranchStatus("CalibEle_energy",1);
   fChain->SetBranchStatus("CalibEle_eta",1);
   fChain->SetBranchStatus("CalibEle_phi",1);
   fChain->SetBranchStatus("CalibEle_mass",1);

   if (!ISDATA_FLAG) {
     fChain->SetBranchStatus("nGenPart",1);
     fChain->SetBranchStatus("GenPart_pdgId",1);
     fChain->SetBranchStatus("GenPart_motherId",1);
     fChain->SetBranchStatus("GenPart_pt",1);
     fChain->SetBranchStatus("GenPart_eta",1);
     fChain->SetBranchStatus("GenPart_phi",1);
     fChain->SetBranchStatus("GenPart_mass",1);
     fChain->SetBranchStatus("GenPart_motherIndex",1);

     //fChain->SetBranchStatus("xsec",1);
     fChain->SetBranchStatus("vtxWeight",1);   // weight to have better agreement between data and MC (added in tree from 06/10/15)

     //added on 23/01/2016
     fChain->SetBranchStatus("SF_BTag",1);
     //variables in sfFriend tree
     fChain->SetBranchStatus("SF_trig1lep",1);
     fChain->SetBranchStatus("SF_trigmetnomu",1);
     fChain->SetBranchStatus("SF_LepTightLoose",1);
     fChain->SetBranchStatus("SF_LepTight",1);
     fChain->SetBranchStatus("SF_NLO",1);
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

   setVarFromConfigFile();
   setSelections();
   setMask();
   //set_SF_NLO_pointers(sf_nlo, ptr_sf_nlo_QCD, ptr_sf_nlo_EWK);

   TVector2 metNoLepTV, ele, lep, met; //lepton and met are used to compute transverse mass

   TLorentzVector l1gen;     // gen level W and l1,l2  (W->(l1 l2)
   TLorentzVector l1reco;

   Double_t mT = 0.0; // transverse mass, computed in loop

   // following indices refer to the leading pair of OS/SF in the list of LepGood. They are initialized with 0 and 1 by default, but can be set with function
   // myGetPairIndexInArray (see functionsForAnalysis.cc for reference). 
   // When myGetPairIndexInArray() is called, the index of "correct" particles will be used. If they are not found (e.g. a pair of OS/SF is mismeasured as 2 mu+), 
   // indices are set as 0 and 1 (and following selection asking lep[0] and lep[1] to be OS or whatever will fail).
   Int_t firstIndex = 0;
   //Int_t secondIndex = 1;

   Int_t firstIndexGen = 0;
   //Int_t secondIndexGen = 1;
   Int_t recoLepFound_flag = 0;
   Int_t genLepFound_flag = 0;
   //Int_t recoGenMatchDR_flag = 0;
   Int_t genTauFound_flag = 0;
   //Int_t W_index = 0; 

   Int_t HLT_passed_flag = 1;          // some computations (for e) require a trigger preselection, while others don't. The former will be done if the flag is set to 1
                                                       // it's set to 1 because if the trigger selection is not applied every event must be considered to be a "good" event having passed all preselections
                                                       // Actually in this code the trigger is not always necessary, but I keep it like this nonetheless.


   Float_t *ptr_nLepLoose = NULL;    // depending on lepton flavour in Z-->ll, it will point to different branches
   Float_t *ptr_nLep10V = NULL;   

   Float_t *ptr_nLepTight = NULL;    // depending on lepton flavour in Z-->ll, it will point to different branches 

   Float_t *ptr_metNoLepPt = NULL;       // only needed for muons, it will point to the branches with the metNoMu_pt, then metNoLepPt = *ptr_metNoLepPt (metNoLepPt defined below)
   //Float_t *ptr_metNoLepEta = NULL; 
   Float_t *ptr_metNoLepPhi = NULL;  

   Int_t *ptr_nRecoLepton = NULL;
   Float_t *ptr_lepton_pt = NULL;
   Float_t *ptr_lepton_eta = NULL;
   Float_t *ptr_lepton_phi = NULL;
   Float_t *ptr_lepton_mass = NULL;

   Float_t nLepLoose = 0.0;               // this variable and the following should be an integer, but in Emanuele's trees they are float, so I keep them as such
   Float_t nLep10V = 0.0;

   Float_t nLepTight = 0.0;               // this variable and the following should be an integer, but in Emanuele's trees they are float, so I keep them as such

   Double_t metNoLepPt = 0.0;        // this variable will be assigned with *ptr_metNoLepPt, where the pointer will point to the branch metNoMu_pt for mu, and with a hand-defined variable for e
   //Double_t metNoLepEta = 0.0;
   Double_t metNoLepPhi = 0.0;   // same story as above

   Int_t nRecoLepton = 0;

   if (fabs(LEP_PDG_ID) == 13) {  // if we have Z -> mumu do stuff...
  
     ptr_nLepLoose = &nMu10V;                      // ask 2 muons
     ptr_nLep10V = &nEle10V;                         // veto on electrons
     ptr_nLepTight = &nMu20T;                                            
     ptr_metNoLepPt = &metNoMu_pt;               // for muons  get this variable from the tree 
     //ptr_metNoLepEta = &metNoMu_eta;               // for muons  get this variable from the tree 
     ptr_metNoLepPhi = &metNoMu_phi;         // for muons  get this variable from the tree
     ptr_nRecoLepton = &nLepGood;
     ptr_lepton_pt = LepGood_pt;
     ptr_lepton_eta = LepGood_eta;
     ptr_lepton_phi = LepGood_phi;
     ptr_lepton_mass = LepGood_mass;

   } else if (fabs(LEP_PDG_ID) == 11) {   // if we have Z -> ee do different stuff...

     ptr_nLepLoose = &nEle10V;                      // ask 2 electrons
     ptr_nLep10V = &nMu10V;                         // veto on muons   
     if (fChain->GetBranch("nEle40T")) ptr_nLepTight = &nEle40T;
     else if (fChain->GetBranch("nEle20T")) ptr_nLepTight = &nEle20T;  // the most recent version for this variable si nEle40T, but older versions use nEle20T

     if (calibEle_flag == 0) {
       ptr_nRecoLepton = &nLepGood;
       ptr_lepton_pt = LepGood_pt;
       ptr_lepton_eta = LepGood_eta;
       ptr_lepton_phi = LepGood_phi;
       ptr_lepton_mass = LepGood_mass;
     } else {
       ptr_nRecoLepton = &nCalibEle;
       ptr_lepton_pt = CalibEle_pt;
       ptr_lepton_eta = CalibEle_eta;
       ptr_lepton_phi = CalibEle_phi;
       ptr_lepton_mass = CalibEle_mass;
     }

   }

   
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

   // Double_t nTotalWeightedEvents = 0.0;     
   // // deciding  what is the event weight
   // Double_t newwgt; 

   Long64_t nentries = fChain->GetEntriesFast();
   cout<<"wlnujetsControlSample::loop()"<<endl;
   cout<<"nentries = "<<nentries<<endl;   

   Long64_t nbytes = 0, nb = 0;

   for (Int_t jentry=0; jentry<nentries; jentry++) {

     Long64_t ientry = LoadTree(jentry);
     if (ientry < 0) break;
     nb = fChain->GetEntry(jentry);   nbytes += nb;
     // if (Cut(ientry) < 0) continue;   

     if (jentry%500000 == 0) cout << jentry << endl;

     UInt_t eventMask = 0; 

     // if (ISDATA_FLAG || unweighted_event_flag) newwgt = 1.0;
     // else {

     //   //newwgt = LUMI * weight * vtxWeight/*/ events_ntot*/;  // starting from 17 November, "events_ntot" substitutes SUMWEIGHT and is already present in the trees. Same for weight, which is now defined as "1000 * xsec * genWeight" (1000*xsec is the cross section in fb, since xsec is in pb.)
     //   // I found out that division by events_ntot was already included in weight definition

     //   if (hasSFfriend_flag != 0) {

     // 	 if (fabs(LEP_PDG_ID) == 13) newwgt = LUMI * weight * vtxWeight * SF_trigmetnomu * SF_LepTight * SF_BTag * SF_NLO;
     // 	 else if (fabs(LEP_PDG_ID) == 11) newwgt = LUMI * weight * vtxWeight * SF_trig1lep * SF_LepTight * SF_BTag * SF_NLO;

     //   } else newwgt = LUMI * weight * vtxWeight * SF_BTag; //SF_BTag is in evVarFriend, not sfFriend

     // }

     newwgt = computeEventWeight();
     nTotalWeightedEvents += newwgt;  // counting events with weights

     if ((nFatJet > 0.5) && (FatJet_pt[0] > 250.) && (fabs(FatJet_eta[0]) < 2.4) && (FatJet_mass[0] > 65.) && (FatJet_mass[0] < 105.) && ((FatJet_tau2[0]/FatJet_tau1[0]) < 0.6) && (metNoMu_pt > 250.)) Vtagged_flag = 1;
     else Vtagged_flag = 0;

     nLepLoose = *ptr_nLepLoose;          
     nLep10V = *ptr_nLep10V;
     nLepTight = *ptr_nLepTight;
     nRecoLepton = *ptr_nRecoLepton;

     // genLepFound_flag is used when analysing WJetsToLNu in MC fo W->munu or W->enu. For other MC samples it's not used.

     if (!ISDATA_FLAG && (GENLEP_TAG != 0)) {

       if (using_wlnujets_MCsample_flag) {

	 //genLepFound_flag = 1;
	 genLepFound_flag = myPartGenWLNuAlgo(nGenPart, GenPart_pdgId, GenPart_motherId, LEP_PDG_ID, firstIndexGen); 
	 if (genLepFound_flag != 0) {
	   eventMask += genLepC.addToMask(1);
	   l1gen.SetPtEtaPhiM(GenPart_pt[firstIndexGen],GenPart_eta[firstIndexGen],GenPart_phi[firstIndexGen],GenPart_mass[firstIndexGen]);
	                                
	 }

       } else if (using_wtaunujets_MCsample_flag) {

	 genTauFound_flag = myPartGenWLNuAlgo(nGenPart, GenPart_pdgId, GenPart_motherId, 15);
	 eventMask += genTauC.addToMask( genTauFound_flag );

       }

     }


     if ( (nRecoLepton >= 1) && fabs(LepGood_pdgId[0]) ==  LEP_PDG_ID) recoLepFound_flag = 1;
     else recoLepFound_flag = 0;

     if (recoLepFound_flag) l1reco.SetPtEtaPhiM(ptr_lepton_pt[0],ptr_lepton_eta[0],ptr_lepton_phi[0],ptr_lepton_mass[0]);

     if (fabs(LEP_PDG_ID) == 13) { 

       if ( HLT_FLAG != 0) {

       	 if ( HLT_MonoJetMetNoMuMHT90 > 0.5 || HLT_MonoJetMetNoMuMHT120 > 0.5 ) HLT_passed_flag = 1; 	 
       	 else HLT_passed_flag = 0;

       }  // end of   if ( HLT_FLAG )

       metNoLepPt = *ptr_metNoLepPt;       
       //metNoLepEta = *ptr_metNoLepEta; 
       metNoLepPhi = *ptr_metNoLepPhi; 
       //metNoLepTV3.SetPtEtaPhi(metNoLepPt,metNoLepEta,metNoLepPhi);   // will use this 3D vector below
       //metNoLepTV.SetMagPhi(metNoLepPt,metNoLepPhi);

     } else if (fabs(LEP_PDG_ID) == 11) { 

       if ( HLT_FLAG != 0 ) {

	 if (HLT_SingleEl == 1) HLT_passed_flag = 1; 	 
	 else HLT_passed_flag = 0;  //continue;

       }  // end of   if ( HLT_FLAG )

       metNoLepTV.SetMagPhi(met_pt,met_phi);
   
       if (recoLepFound_flag) {
	 ele.SetMagPhi(ptr_lepton_pt[0],ptr_lepton_phi[0]);
	 metNoLepTV += ele;
       }

       metNoLepPt = metNoLepTV.Mod();

     }

     // beginning of eventMask building

     
     // genLepC added to mask above if ISDATA_FLAG == false (in order not to repeat here the check) 
     
     eventMask += HLTC.addToMask(HLT_passed_flag); 
     eventMask += jet1C.addToMask(nJetClean30 >= 1 && JetClean_pt[0] > J1PT /*&& fabs(JetClean_eta[0]) < J1ETA*/);
     eventMask += jetMetDphiMinC.addToMask(fabs(dphijm > JMET_DPHI_MIN));
     eventMask += jetNoiseCleaningC.addToMask(JetClean_leadClean[0] > 0.5);
     eventMask += bjetVetoC.addToMask(nBTag15 < 0.5);
     eventMask += lepLooseVetoC.addToMask(nLep10V < 0.5);
     eventMask += tauLooseVetoC.addToMask(nTauClean18V < 0.5);
     eventMask += gammaLooseVetoC.addToMask(nGamma15V < 0.5);     
     eventMask += metC.addToMask(met_pt > 50);    
     eventMask += recoilC.addToMask(metNoLepPt > METNOLEP_START);
     eventMask += metFiltersC.addToMask(cscfilter == 1 && ecalfilter == 1 && hbheFilterNew25ns == 1 && hbheFilterIso == 1 && Flag_eeBadScFilter > 0.5);  
     eventMask += VtagC.addToMask(Vtagged_flag);
     eventMask += noVtagC.addToMask(!Vtagged_flag);
     eventMask += oneLepLooseC.addToMask(nLepLoose > 0.5 && nLepLoose < 1.5);
     if (fabs(LEP_PDG_ID) == 11) eventMask += tightLepC.addToMask(nLepTight < 1.5 && nLepTight > 0.5 && ptr_lepton_pt[0] > LEP1PT && fabs(LepGood_pdgId[0]) == 11);
     else eventMask += tightLepC.addToMask(nLepTight > 0.5 && nLepTight < 1.5);

     // end of eventMask building

     // test matching of reco and gen lep for DY MC 
     if (using_wlnujets_MCsample_flag && (GENLEP_TAG != 0)) {
    
       // now we require a DeltaR cut between gen-level and reco-level leptons (mu or e) from W events to assess that lreco comes from lgen

       if (genLepFound_flag && recoLepFound_flag) {       

	 Double_t DeltaR_lreco_lgen_pair1 = 100.0;
       
	 if(LepGood_pdgId[0] == GenPart_pdgId[firstIndexGen] ) {
	 
	   DeltaR_lreco_lgen_pair1 = l1reco.DeltaR(l1gen);

	 } 
       
	 if (DeltaR_lreco_lgen_pair1 < 0.1) eventMask += recoGenLepMatchC.addToMask(1);
	 else eventMask += recoGenLepMatchC.addToMask(0);

       }

     }

     analysisMask.countEvents(eventMask,newwgt);
     analysisMask_monoJ.countEvents(eventMask,newwgt);
     analysisMask_monoV.countEvents(eventMask,newwgt);

     // transverse mass to fill histograms below
     lep.SetMagPhi(ptr_lepton_pt[0],ptr_lepton_phi[0]);
     met.SetMagPhi(met_pt,met_phi);
     mT = sqrt(2*ptr_lepton_pt[0]*met_pt*(1- TMath::Cos(lep.DeltaPhi(met))));

     if ( ((eventMask & analysisMask_monoJ.globalMask.back()) == analysisMask_monoJ.globalMask.back()) ) {
       
       // this histogram holds the final yields in bins of MET
	 HYieldsMetBin->Fill(metNoLepPt,newwgt);

	 HhtDistribution->Fill(htJet25,newwgt);
	 HtransverseMass->Fill(mT,newwgt);
	 HmetNoLepDistribution->Fill(metNoLepPt,newwgt);
	 HvtxDistribution->Fill(nVert,newwgt);
	 HnjetsDistribution->Fill(nJetClean30,newwgt);
	 Hjet1etaDistribution->Fill(JetClean_eta[0],newwgt);
	 Hjet1ptDistribution->Fill(JetClean_pt[0],newwgt);
	 HjetMetDphiMinDistribution->Fill(dphijm,newwgt);
	 Hlep1ptDistribution->Fill(ptr_lepton_pt[0],newwgt);
	 Hlep1etaDistribution->Fill(ptr_lepton_eta[0],newwgt);
	 if (nJetClean30 >= 2) {
	   Hj1j2dphiDistribution->Fill(dphijj,newwgt);
	   Hjet2etaDistribution->Fill(JetClean_eta[1],newwgt);
	   Hjet2ptDistribution->Fill(JetClean_pt[1],newwgt);
	 }

	 if (using_wlnujets_MCsample_flag) {

	   HYieldsMetBin_qcdRenScaleUp->Fill(metNoLepPt,(newwgt * SF_NLO_QCD_renScaleUp/ SF_NLO_QCD));
	   HYieldsMetBin_qcdRenScaleDown->Fill(metNoLepPt,(newwgt * SF_NLO_QCD_renScaleDown/ SF_NLO_QCD));
	   HYieldsMetBin_qcdFacScaleUp->Fill(metNoLepPt,(newwgt * SF_NLO_QCD_facScaleUp/ SF_NLO_QCD));
	   HYieldsMetBin_qcdFacScaleDown->Fill(metNoLepPt,(newwgt * SF_NLO_QCD_facScaleDown/ SF_NLO_QCD));
	   HYieldsMetBin_qcdPdfUp->Fill(metNoLepPt,(newwgt * SF_NLO_QCD_pdfUp/ SF_NLO_QCD));
	   HYieldsMetBin_qcdPdfDown->Fill(metNoLepPt,(newwgt * SF_NLO_QCD_pdfDown/ SF_NLO_QCD));
	   HYieldsMetBin_ewkUp->Fill(metNoLepPt,(newwgt * SF_NLO_EWK_up/ SF_NLO_EWK));
	   HYieldsMetBin_ewkDown->Fill(metNoLepPt,(newwgt * SF_NLO_EWK_down/ SF_NLO_EWK));

	 }

     } else if (((eventMask & analysisMask_monoV.globalMask.back()) == analysisMask_monoV.globalMask.back())) {

       HYieldsMetBin_monoV->Fill(metNoLepPt,newwgt);

       HhtDistribution_monoV->Fill(htJet25,newwgt);
       HtransverseMass_monoV->Fill(mT,newwgt);
       HmetNoLepDistribution_monoV->Fill(metNoLepPt,newwgt);
       HvtxDistribution_monoV->Fill(nVert,newwgt);
       HnjetsDistribution_monoV->Fill(nJetClean30,newwgt);
       Hjet1etaDistribution_monoV->Fill(FatJet_eta[0],newwgt);
       Hjet1ptDistribution_monoV->Fill(FatJet_pt[0],newwgt);
       Hlep1ptDistribution_monoV->Fill(ptr_lepton_pt[0],newwgt);
       Hlep1etaDistribution_monoV->Fill(ptr_lepton_eta[0],newwgt);
       HprunedMassDistribution_monoV->Fill(FatJet_mass[0],newwgt);
       Htau2OverTau1Distribution_monoV->Fill(FatJet_tau2[0]/FatJet_tau1[0],newwgt);

       if (hasScaledHistograms_flag) {
	 
	 HYieldsMetBin_qcdRenScaleUp_monoV->Fill(metNoMu_pt,(newwgt * SF_NLO_QCD_renScaleUp/ SF_NLO_QCD));
	 HYieldsMetBin_qcdRenScaleDown_monoV->Fill(metNoMu_pt,(newwgt * SF_NLO_QCD_renScaleDown/ SF_NLO_QCD));
	 HYieldsMetBin_qcdFacScaleUp_monoV->Fill(metNoMu_pt,(newwgt * SF_NLO_QCD_facScaleUp/ SF_NLO_QCD));
	 HYieldsMetBin_qcdFacScaleDown_monoV->Fill(metNoMu_pt,(newwgt * SF_NLO_QCD_facScaleDown/ SF_NLO_QCD));
	 HYieldsMetBin_qcdPdfUp_monoV->Fill(metNoMu_pt,(newwgt * SF_NLO_QCD_pdfUp/ SF_NLO_QCD));
	 HYieldsMetBin_qcdPdfDown_monoV->Fill(metNoMu_pt,(newwgt * SF_NLO_QCD_pdfDown/ SF_NLO_QCD));
	 HYieldsMetBin_ewkUp_monoV->Fill(metNoMu_pt,(newwgt * SF_NLO_EWK_up/ SF_NLO_EWK));
	 HYieldsMetBin_ewkDown_monoV->Fill(metNoMu_pt,(newwgt * SF_NLO_EWK_down/ SF_NLO_EWK));

       }

     }

       
   }                        // end of loop on entries

   mySpaces(cout,2);
   selection::printSelectionFlowAndYields(cout, LUMI, nTotalWeightedEvents, &analysisMask);
   mySpaces(cout,2);
   selection::printSelectionFlowAndYields(cout, LUMI, nTotalWeightedEvents, &analysisMask_monoJ);
   mySpaces(cout,2);
   selection::printSelectionFlowAndYields(cout, LUMI, nTotalWeightedEvents, &analysisMask_monoV);

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
   selection::printSelectionFlowAndYields(myfile, LUMI, nTotalWeightedEvents, &analysisMask_monoJ);
   mySpaces(myfile,2);
   selection::printSelectionFlowAndYields(myfile, LUMI, nTotalWeightedEvents, &analysisMask_monoV);
   mySpaces(myfile,2);
   myPrintYieldsMetBinInStream(myfile, HYieldsMetBin, metBinEdgesVector.data(), nMetBins);

   myfile.close();

   // filling with yields and efficiency: I will use efficiency with respect to total and not to previous step, but I could make this choice in the config file

   if ((GENLEP_TAG != 0) && (using_wlnujets_MCsample_flag == 1 || using_wtaunujets_MCsample_flag == 1)) {
     fillRowVector(nTotalWeightedEvents, analysisSelectionManager, analysisMask, yRow, eRow, uncRow,1);
     fillRowVector(nTotalWeightedEvents, analysisSelectionManager_monoJ, analysisMask_monoJ, yRow_monoJ, eRow_monoJ, uncRow_monoJ,1);
     fillRowVector(nTotalWeightedEvents, analysisSelectionManager_monoV, analysisMask_monoV, yRow_monoV, eRow_monoV, uncRow_monoV,1);
   } else {
     fillRowVector(nTotalWeightedEvents, analysisSelectionManager, analysisMask, yRow, eRow, uncRow,0);
     fillRowVector(nTotalWeightedEvents, analysisSelectionManager_monoJ, analysisMask_monoJ, yRow_monoJ, eRow_monoJ, uncRow_monoJ,0);
     fillRowVector(nTotalWeightedEvents, analysisSelectionManager_monoV, analysisMask_monoV, yRow_monoV, eRow_monoV, uncRow_monoV,0);
   }

   // fill last bin with overflow 
   setHistogramLastBinAsOverFlow(hasScaledHistograms_flag);
   if (hasScaledHistograms_flag) createSystematicsHistogram();

   rootFile->Write();

   rootFile->Close();
   delete rootFile;

   //creating a .tex file to build tables with data
   myCreateTexTable(TEX_FNAME, outputFolder, LUMI,nTotalWeightedEvents, anaMasksPtrCollection);

   // end of tex file
   mySpaces(cout,2);

}




