#define monojetAna_cxx
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

#ifdef monojetAna_cxx

//===============================================

monojetAna::monojetAna(TTree *tree) : AnalysisDarkMatter(tree) {
  //cout <<"check in constructor "<<endl;

  // initialize some variables with sensible values. They will be set later depending on config file
  // following variables, if any, are specific to monojet analysis  

  J1PT = 100.0;
  J1ETA = 2.5;

}

//===============================================


void monojetAna::setNumberParameterValue(const string parameterName, const Double_t value) {

  AnalysisDarkMatter::setNumberParameterValue(parameterName, value);

  if (parameterName == "J1PT") J1PT = value;
  else if (parameterName == "J1ETA") J1ETA = value;
  //else if (parameterName == "J2PT") J2PT = value;
  //else if (parameterName == "J2ETA") J2ETA = value;
  //else if (parameterName == "J1J2DPHI") J1J2DPHI = value;

}

//===============================================

void monojetAna::setVarFromConfigFile() {

  AnalysisDarkMatter::setVarFromConfigFile();

}

//===============================================

void monojetAna::setSelections() {

  AnalysisDarkMatter::setSelections();

  jet1C.set("jet1",Form("jet1pt > %3.0lf",J1PT),Form("nJetClean >= 1 && JetClean1_pt > %4.0lf",(Double_t)J1PT));
  jetMetDphiMinC.set("dphiMin(j,MET)",Form("min[dphi(j,MET)] > %1.1lf",JMET_DPHI_MIN),"minimum dphi between jets and MET (using only the first 4 jets)");
  jetNoiseCleaningC.set("jet1 cleaning","noise cleaning","energy fractions (only for jet1): CH > 0.1; NH < 0.8");
  VtagC.set("V tag","V-tagged","nFatJet>0 and leading FatJet(ak8) with: pT>250, |eta|<2.4, pruned mass in [65,105], tau2/tau1 < 0.6, recoil > 250");
  noVtagC.set("no V tag","V-notTagged","if any of V-tagged conditions fails (see VtagC for details)");
  ak8jet1C.set("ak8 jet1","ak8 jet1pt > 250","nFatJetClean > 0.5 && FatJetClean_pt[0] > 250. && fabs(FatJetClean_eta[0]) < 2.4");
  ak8Tau2OverTau1C.set("tau2/tau1","tau2/tau1 < 0.6","FatJetClean_tau2[0]/FatJetClean_tau1[0] < 0.6");
  ak8prunedMassC.set("pruned mass","pruned mass in [65,105]","FatJetClean_prunedMass[0] > 65. && FatJetClean_prunedMass[0] < 105.");
  harderRecoilC.set("hard met cut","recoil > 250");
  if (MET_FILTERS_FLAG != 0) metFiltersC.set("met filters","met filters","cscfilter, ecalfilter, hbheFilterNeww25ns, hbheFilterIso, Flag_eeBadScFilter");
  


  selection::checkMaskLength();

}

//===============================================

void monojetAna::setHistograms() {

  AnalysisDarkMatter::setHistograms();

  // histograms for monoV

  HYieldsMetBin_monoV = new TH1D("HYieldsMetBin_monoV","yields in bins of met; #slash{E}_{T};# of events",nMetBins_monoV,metBinEdgesVector_monoV.data());
  HhtDistribution_monoV = new TH1D("HhtDistribution_monoV","",200,0.0,2000.0);
  HvtxDistribution_monoV = new TH1D("HvtxDistribution_monoV","",40,-0.5,39.5);   
  HnjetsDistribution_monoV = new TH1D("HnjetsDistribution_monoV","njets using nJetClean30",10,-0.5,9.5);   
  Hjet1etaDistribution_monoV = new TH1D("Hjet1etaDistribution_monoV","leading ak8 jets",60,-3.0,3.0);
  HrecoilDistribution_monoV = new TH1D("HrecoilDistribution_monoV","",100,0.0,1000.0);
  Hjet1ptDistribution_monoV = new TH1D("Hjet1ptDistribution_monoV","leading ak8 jets",100,0.0,1000.0); 
  HprunedMassDistribution_monoV = new TH1D("HprunedMassDistribution_monoV","",16,65.0,105.0);
  Htau2OverTau1Distribution_monoV = new TH1D("Htau2OverTau1Distribution_monoV","",20,0.0,1.0);
  
  // saving histograms with bin edges of other histograms used (e.g. content of metBinEdges array ...)
  HmetBinEdges_monoV = new TH1D("HmetBinEdges_monoV","bin edges for met distributions",nMetBins_monoV+1,0.0,nMetBins_monoV+1);
  for (Int_t i = 0; i <= nMetBins_monoV; i++) {
    HmetBinEdges_monoV->SetBinContent(i+1,metBinEdgesVector_monoV[i]);
  }


}

//===============================================

void monojetAna::setScaleFactorHistograms() {

  AnalysisDarkMatter::setScaleFactorHistograms();

  //monoV
  HYieldsMetBin_qcdRenScaleUp_monoV = new TH1D("HYieldsMetBin_qcdRenScaleUp_monoV","yields in bins of met; #slash{E}_{T};# of events",nMetBins_monoV,metBinEdgesVector_monoV.data());
  HYieldsMetBin_qcdRenScaleDown_monoV = new TH1D("HYieldsMetBin_qcdRenScaleDown_monoV","yields in bins of met; #slash{E}_{T};# of events",nMetBins_monoV,metBinEdgesVector_monoV.data());
  HYieldsMetBin_qcdFacScaleUp_monoV = new TH1D("HYieldsMetBin_qcdFacScaleUp_monoV","yields in bins of met; #slash{E}_{T};# of events",nMetBins_monoV,metBinEdgesVector_monoV.data());
  HYieldsMetBin_qcdFacScaleDown_monoV = new TH1D("HYieldsMetBin_qcdFacScaleDown_monoV","yields in bins of met; #slash{E}_{T};# of events",nMetBins_monoV,metBinEdgesVector_monoV.data());
  HYieldsMetBin_qcdPdfUp_monoV = new TH1D("HYieldsMetBin_qcdPdfUp_monoV","yields in bins of met; #slash{E}_{T};# of events",nMetBins_monoV,metBinEdgesVector_monoV.data());
  HYieldsMetBin_qcdPdfDown_monoV = new TH1D("HYieldsMetBin_qcdPdfDown_monoV","yields in bins of met; #slash{E}_{T};# of events",nMetBins_monoV,metBinEdgesVector_monoV.data());
  HYieldsMetBin_ewkUp_monoV = new TH1D("HYieldsMetBin_ewkUp_monoV","yields in bins of met; #slash{E}_{T};# of events",nMetBins_monoV,metBinEdgesVector_monoV.data());
  HYieldsMetBin_ewkDown_monoV = new TH1D("HYieldsMetBin_ewkDown_monoV","yields in bins of met; #slash{E}_{T};# of events",nMetBins_monoV,metBinEdgesVector_monoV.data());
    
  HSyst_qcdRenScale_monoV = new TH1D("HSyst_qcdRenScale_monoV","systematic uncertainty for QCD renormalization scale",nMetBins_monoV,metBinEdgesVector_monoV.data());
  HSyst_qcdFacScale_monoV = new TH1D("HSyst_qcdFacScale_monoV","systematic uncertainty for QCD factorization scale",nMetBins_monoV,metBinEdgesVector_monoV.data());
  HSyst_qcdPdf_monoV = new TH1D("HSyst_qcdPdf_monoV","systematic uncertainty for QCD due to PDF uncertainty",nMetBins_monoV,metBinEdgesVector_monoV.data());
  HSyst_ewk_monoV = new TH1D("HSyst_ewk_monoV","systematic uncertainty for EWK",nMetBins_monoV,metBinEdgesVector_monoV.data());
  HSyst_total_monoV = new TH1D("HSyst_total_monoV","total systematic uncertainty (sum in quadrature of all single systematics)",nMetBins_monoV,metBinEdgesVector_monoV.data());
  

}



//===============================================

void monojetAna::setHistogramLastBinAsOverFlow(const Int_t hasScaledHistograms = 0) {

  AnalysisDarkMatter::setHistogramLastBinAsOverFlow(hasScaledHistograms);

  myAddOverflowInLastBin(HYieldsMetBin_monoV);
  myAddOverflowInLastBin(HhtDistribution_monoV);
  myAddOverflowInLastBin(HrecoilDistribution_monoV);
  myAddOverflowInLastBin(Hjet1ptDistribution_monoV);

  if (hasScaledHistograms) {

    myAddOverflowInLastBin(HYieldsMetBin_qcdRenScaleUp_monoV);
    myAddOverflowInLastBin(HYieldsMetBin_qcdRenScaleDown_monoV);
    myAddOverflowInLastBin(HYieldsMetBin_qcdFacScaleUp_monoV);
    myAddOverflowInLastBin(HYieldsMetBin_qcdFacScaleDown_monoV);
    myAddOverflowInLastBin(HYieldsMetBin_qcdPdfUp_monoV);
    myAddOverflowInLastBin(HYieldsMetBin_qcdPdfDown_monoV);
    myAddOverflowInLastBin(HYieldsMetBin_ewkUp_monoV);
    myAddOverflowInLastBin(HYieldsMetBin_ewkDown_monoV);

  }

}

//===============================================

void monojetAna::createSystematicsHistogram() {

  AnalysisDarkMatter::createSystematicsHistogram();  

    // same for monoV

  // computing systematic uncertainties and saving them as histograms.
  myBuildSystematicsHistogram(HSyst_qcdRenScale_monoV, HYieldsMetBin_monoV, HYieldsMetBin_qcdRenScaleUp_monoV, HYieldsMetBin_qcdRenScaleDown_monoV);
  myBuildSystematicsHistogram(HSyst_qcdFacScale_monoV, HYieldsMetBin_monoV, HYieldsMetBin_qcdFacScaleUp_monoV, HYieldsMetBin_qcdFacScaleDown_monoV);
  myBuildSystematicsHistogram(HSyst_qcdPdf_monoV, HYieldsMetBin_monoV, HYieldsMetBin_qcdPdfUp_monoV, HYieldsMetBin_qcdPdfDown_monoV);
  myBuildSystematicsHistogram(HSyst_ewk_monoV, HYieldsMetBin_monoV, HYieldsMetBin_ewkUp_monoV, HYieldsMetBin_ewkDown_monoV);

  // define an empty histogram to sum uncertainties in a clean way
  TH1D *Htmp = new TH1D("Htmp","",nMetBins_monoV,metBinEdgesVector_monoV.data());
  vector<TH1D*> hptr;
  hptr.push_back(HSyst_qcdRenScale_monoV);
  hptr.push_back(HSyst_qcdFacScale_monoV);
  hptr.push_back(HSyst_qcdPdf_monoV);
  hptr.push_back(HSyst_ewk_monoV);

  HSyst_total_monoV->Multiply(HSyst_total_monoV);  // before adding anything we must get the square of this histogram, because we will add histograms in quadrature    
  for (UInt_t i = 0; i < hptr.size(); i++) {
    Htmp->Multiply(hptr[i],hptr[i]); // square of bin content for each single systematic histogram
    HSyst_total_monoV->Add(Htmp);             // adding the squares
  }

  for (Int_t i = 0; i <= (HSyst_total_monoV->GetNbinsX() + 1); i++) {  // computing square root of each bin's content (from underflow to overflow bin, but they should be empty)
    HSyst_total_monoV->SetBinContent(i, sqrt(HSyst_total_monoV->GetBinContent(i)));
  }

  delete Htmp;



}

//===============================================


void monojetAna::fillEventMask(ULong64_t & eventMask) {

  AnalysisDarkMatter::fillEventMask(eventMask);  

  eventMask += jet1C.addToMask(nJetClean30 >= 1 && JetClean_pt[0] > J1PT /*&& fabs(JetClean_eta[0]) < J1ETA*/);
  eventMask += jetNoiseCleaningC.addToMask(JetClean_leadClean[0] > 0.5);
  eventMask += ak8jet1C.addToMask((nFatJetClean > 0.5) && (FatJetClean_pt[0] > 250.) && (fabs(FatJetClean_eta[0]) < 2.4));
  eventMask += ak8Tau2OverTau1C.addToMask(((FatJetClean_tau2[0]/FatJetClean_tau1[0]) < 0.6));
  eventMask += ak8prunedMassC.addToMask((FatJetClean_prunedMass[0] > 65.) && (FatJetClean_prunedMass[0] < 105.));
  if (MET_FILTERS_FLAG != 0) eventMask += metFiltersC.addToMask(cscfilter == 1 && ecalfilter == 1 && hbheFilterNew25ns == 1 && hbheFilterIso == 1 && Flag_eeBadScFilter > 0.5);
  eventMask += jetMetDphiMinC.addToMask(fabs(dphijm) > JMET_DPHI_MIN);

}

//===============================================


#endif
