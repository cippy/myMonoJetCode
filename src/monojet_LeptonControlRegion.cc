#define monojet_LeptonControlRegion_cxx
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

#ifdef monojet_LeptonControlRegion_cxx

monojet_LeptonControlRegion::monojet_LeptonControlRegion(TTree *tree) : monojetAna(tree) {
  //cout <<"check in constructor "<<endl;
  //edimarcoTree_v3::Init(tree);
  // suffix = "";
  // uncertainty = "";
  // configFileName = NULL;
  // ISDATA_FLAG = 0;
  // unweighted_event_flag = 0;
  // hasSFfriend_flag = 0;
  //monojetAna::Init(tree);  // could also be just Init(tree)
  // maybe this is useless because the constructor already calls Init

  recoLepFound_flag = 0;
  genLepFound_flag = 0;

  HLT_passed_flag = 1; // set to 1 because when it is not needed, trigger is considered to be automatically passed
  
  // Float_t variables
  nLepLoose = 0.0;               // this variable and the following should be an integer, but in Emanuele's trees they are float, so I keep them as such       
  nLep10V = 0.0;
  nLepTight = 0.0;               // this variable and the following should be an integer, but in Emanuele's trees they are float, so I keep them as such       

  //Double_t variables
  // this variable will be assigned with *ptr_metNoLepPt, where the pointer will point to the branch metNoMu_pt for mu, and with a hand-defined variable for electrons
  metNoLepPt = 0.0;                                                                                                                                  
  //Double_t metNoLepEta = 0.0;                                                                                                                                        
  metNoLepPhi = 0.0;   // same story as above                                                                                                     

  // Int_t variables
  nRecoLepton = 0;


}

#endif

//===============================================

// void monojet_LeptonControlRegion::Init(TTree *tree) {
//   monojetAna::Init(tree);
// } 

//===============================================

void monojet_LeptonControlRegion::setSelections() {

  monojetAna::setSelections();

  if (fabs(LEP_PDG_ID) == 13) {  // if we have Z -> mumu do stuff...

    lepLooseVetoC.set("ele veto","electrons veto");
    if (METNOLEP_START != 0) recoilC.set(Form("recoil > %2.0lf",METNOLEP_START),Form("metNoMu > %2.0lf",METNOLEP_START));

  } else if (fabs(LEP_PDG_ID) == 11) {   // if we have Z -> ee do different stuff...

    lepLooseVetoC.set("muon veto","muons veto");
    if (METNOLEP_START != 0) recoilC.set(Form("recoil > %2.0lf",METNOLEP_START),Form("metNoEle > %2.0lf",METNOLEP_START));
  }

  selection::checkMaskLength();
  // following is directly called in the derived class
  // selection::printActiveSelections(cout); 

}

//===============================================

void monojet_LeptonControlRegion::setHistograms() {

  monojetAna::setHistograms();
  
  Hlep1ptDistribution = new TH1D("Hlep1ptDistribution","",200,0.0,1000.0);
  Hlep1etaDistribution = new TH1D("Hlep1etaDistribution","",100,-5.0,5.0);

  Hlep1ptDistribution_monoV = new TH1D("Hlep1ptDistribution_monoV","",200,0.0,1000.0);
  Hlep1etaDistribution_monoV = new TH1D("Hlep1etaDistribution_monoV","",100,-5.0,5.0);

}

//===============================================

void monojet_LeptonControlRegion::setScaleFactorHistograms() {

  monojetAna::setScaleFactorHistograms();
  
}

//===============================================

void monojet_LeptonControlRegion::setHistogramLastBinAsOverFlow(const Int_t hasScaledHistograms = 0) {

  monojetAna::setHistogramLastBinAsOverFlow(hasScaledHistograms);

  myAddOverflowInLastBin(Hlep1ptDistribution);

  myAddOverflowInLastBin(Hlep1ptDistribution_monoV);

}

//===============================================

void monojet_LeptonControlRegion::setNumberParameterValue(const std::string parameterName, const Double_t value) {

  monojetAna::setNumberParameterValue(parameterName, value);

  if (parameterName == "LEP_PDG_ID") LEP_PDG_ID = (value < 0) ? (-0.5 + value) : (0.5 + value);
  else if (parameterName == "LEP1PT") LEP1PT = value;
  //else if (parameterName == "LEP_ISO_04") LEP_ISO_04 = value;
  else if (parameterName == "LEP1ETA") LEP1ETA = value;
  else if (parameterName == "HLT_LEP1PT") HLT_LEP1PT = value;
  else if (parameterName == "HLT_LEP1ETA") HLT_LEP1ETA = value;

  if (!ISDATA_FLAG) {

    if (parameterName == "GENLEP_TAG") GENLEP_TAG = (value < 0) ? (-0.5 + value) : (0.5 + value);
    else if (parameterName == "GENLEP1PT") GENLEP1PT = value;
    else if (parameterName == "GENLEP1ETA") GENLEP1ETA = value;

  }

}

//===============================================

void monojet_LeptonControlRegion::setControlSampleSpecificParameter() {

  // no common specific parameter among control samples for now
  
}

//===============================================

void monojet_LeptonControlRegion::setVarFromConfigFile() {

  monojetAna::setVarFromConfigFile();
  //monojet_LeptonControlRegion::setControlSampleSpecificParameter();

}

//===============================================

void monojet_LeptonControlRegion::createSystematicsHistogram() {

  monojetAna::createSystematicsHistogram();

}

//===============================================                                                                                                                       

void monojet_LeptonControlRegion::fillEventMask(UInt_t & eventMask) {

  monojetAna::fillEventMask(eventMask);

  //if (HLT_FLAG != 0) eventMask += HLTC.addToMask(HLT_passed_flag);  // trigger added to mask in the specific region
  eventMask += lepLooseVetoC.addToMask(nLep10V < 0.5);     // veto on electrons (if V->mu X) or on electrons (if V->e X)
  eventMask += gammaLooseVetoC.addToMask(nGamma15V < 0.5);
  eventMask += VtagC.addToMask(Vtagged_flag);
  eventMask += noVtagC.addToMask(!Vtagged_flag);
  eventMask += recoilC.addToMask(metNoLepPt > METNOLEP_START);
  eventMask += harderRecoilC.addToMask(metNoLepPt > 250.);


}


//===============================================



