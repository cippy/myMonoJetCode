#define vbfHiggsToInv_SignalRegion_cxx
#include "EmanTreeAnalysis.h"
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

#ifdef vbfHiggsToInv_SignalRegion_cxx

vbfHiggsToInv_SignalRegion::vbfHiggsToInv_SignalRegion(TTree *tree) : vbfHiggsToInvAna(tree) {
  //cout <<"check in constructor "<<endl;
  // suffix = "";
  // uncertainty = "";
  // configFileName = NULL;
  // ISDATA_FLAG = 0;
  // unweighted_event_flag = 0;
  // hasSFfriend_flag = 0;
  // Init(tree);

}

#endif

//===============================================

void vbfHiggsToInv_SignalRegion::setSelections() {

  //cout << "CHECK ?!?"<< endl;

  vbfHiggsToInvAna::setSelections();

  //cout << "CHECK VBFSIGREGION"<< endl;
  recoilC.set(Form("recoil > %2.0lf",METNOLEP_START),Form("metNoMu > %4.0lf",METNOLEP_START),"first cut on met");
  muonLooseVetoC.set("muon veto","muons veto");
  electronLooseVetoC.set("ele veto","electrons veto");
  gammaLooseVetoC.set("photon veto","photons veto");
  //cout << "CHECK VBFSIGREGION"<< endl;

  selection::checkMaskLength();
  selection::printActiveSelections(cout); 



}

//===============================================

void vbfHiggsToInv_SignalRegion::setMask() {

  // analysisMask.setName("vbfHiggsToInv signal selection (inclusive)");

  // if (HLT_FLAG != 0) analysisMask.append(HLTC.get2ToId());
  // if (MET_FILTERS_FLAG != 0) analysisMask.append(metFiltersC.get2ToId());
  // analysisMask.append(muonLooseVetoC.get2ToId());
  // analysisMask.append(electronLooseVetoC.get2ToId());
  // if (TAU_VETO_FLAG) analysisMask.append(tauLooseVetoC.get2ToId());
  // analysisMask.append(gammaLooseVetoC.get2ToId());
  // analysisMask.append(bjetVetoC.get2ToId());
  // if (METNOLEP_START != 0) analysisMask.append(recoilC.get2ToId());  
  // analysisMask.append(jet1C.get2ToId());
  // analysisMask.append(jetNoiseCleaningC.get2ToId());
  // analysisMask.append(jetMetDphiMinC.get2ToId());
  // //analysisMask.append(noVtagC.get2ToId());

  // analysisSelectionManager.SetMaskPointer(&analysisMask);

  // if (HLT_FLAG != 0) analysisSelectionManager.append(&HLTC);
  // if (MET_FILTERS_FLAG != 0) analysisSelectionManager.append(&metFiltersC);
  // analysisSelectionManager.append(&muonLooseVetoC);
  // analysisSelectionManager.append(&electronLooseVetoC);
  // if (TAU_VETO_FLAG) analysisSelectionManager.append(&tauLooseVetoC);
  // analysisSelectionManager.append(&gammaLooseVetoC);
  // analysisSelectionManager.append(&bjetVetoC);
  // if (METNOLEP_START != 0) analysisSelectionManager.append(&recoilC); 
  // analysisSelectionManager.append(&jet1C);
  // analysisSelectionManager.append(&jetNoiseCleaningC);
  // analysisSelectionManager.append(&jetMetDphiMinC);
  // //analysisSelectionManager.append(&noVtagC);

  analysisMask.setName("vbfHiggsToInv signal selection (inclusive)");

  if (HLT_FLAG != 0) analysisMask.append(HLTC.get2ToId());
  if (MET_FILTERS_FLAG != 0) analysisMask.append(metFiltersC.get2ToId());
  analysisMask.append(muonLooseVetoC.get2ToId());
  analysisMask.append(electronLooseVetoC.get2ToId());
  analysisMask.append(gammaLooseVetoC.get2ToId());
  if (TAU_VETO_FLAG) analysisMask.append(tauLooseVetoC.get2ToId());
  if (B_VETO_FLAG) analysisMask.append(bjetVetoC.get2ToId());
  analysisMask.append(vbfTaggedJets_jetsPtC.get2ToId());
  analysisMask.append(jetNoiseCleaningC.get2ToId());
  analysisMask.append(vbfTaggedJets_inVMassC.get2ToId());
  analysisMask.append(vbfTaggedJets_deltaEtaC.get2ToId());
  analysisMask.append(jetMetDphiMinC.get2ToId());
  if (METNOLEP_START != 0) analysisMask.append(recoilC.get2ToId());  
  //analysisMask.append(noVtagC.get2ToId());

  analysisSelectionManager.SetMaskPointer(&analysisMask);

  if (HLT_FLAG != 0) analysisSelectionManager.append(&HLTC);
  if (MET_FILTERS_FLAG != 0) analysisSelectionManager.append(&metFiltersC);
  analysisSelectionManager.append(&muonLooseVetoC);
  analysisSelectionManager.append(&electronLooseVetoC);
  analysisSelectionManager.append(&gammaLooseVetoC);
  if (TAU_VETO_FLAG) analysisSelectionManager.append(&tauLooseVetoC);
  if (B_VETO_FLAG) analysisSelectionManager.append(&bjetVetoC);
  analysisSelectionManager.append(&vbfTaggedJets_jetsPtC);
  analysisSelectionManager.append(&jetNoiseCleaningC);
  analysisSelectionManager.append(&vbfTaggedJets_inVMassC);
  analysisSelectionManager.append(&vbfTaggedJets_deltaEtaC);
  analysisSelectionManager.append(&jetMetDphiMinC);
  if (METNOLEP_START != 0) analysisSelectionManager.append(&recoilC); 
  //analysisSelectionManager.append(&noVtagC);

  anaMasksPtrCollection.push_back(&analysisMask);


}

//===============================================

void vbfHiggsToInv_SignalRegion::setHistograms() {

  vbfHiggsToInvAna::setHistograms();

  // following histograms are really needed only when using Zvv or Wlv samples, so they are not instantiated in vbfHiggsToInvAna::setHistograms();

  if (suffix == "ZJetsToNuNu" || suffix == "WJetsToLNu") {
    hasScaledHistograms_flag = 1;
    setScaleFactorHistograms();
  }

  //histograms specific to vbfHiggsToInv selection (but not to control regions) must be set here as in vbfHiggsToInvAna::setHistograms()

}

//===============================================

void vbfHiggsToInv_SignalRegion::setHistogramLastBinAsOverFlow(const Int_t hasScaledHistograms = 0) {

  vbfHiggsToInvAna::setHistogramLastBinAsOverFlow(hasScaledHistograms);

}

void vbfHiggsToInv_SignalRegion::setBranchStatusForAnalysis() {

  vbfHiggsToInvAna::setBranchStatusForAnalysis();

}

//===============================================

void vbfHiggsToInv_SignalRegion::setNumberParameterValue(const std::string parameterName, const Double_t value) {
	 	 
  vbfHiggsToInvAna::setNumberParameterValue(parameterName, value);

  //parameters specific to vbfHiggsToInv selection (but not to control regions) must be set here as in vbfHiggsToInvAna::setNumberParameterValue()

}

//===============================================

void vbfHiggsToInv_SignalRegion::setVarFromConfigFile() {

  vbfHiggsToInvAna::setVarFromConfigFile();

}

//===============================================

Double_t vbfHiggsToInv_SignalRegion::computeEventWeight() {

  if (ISDATA_FLAG || unweighted_event_flag) return 1.0;
  else {
    // sf_nlo_weight = (*ptr_sf_nlo_QCD) * (*ptr_sf_nlo_EWK);
    // return LUMI * weight * vtxWeight * SF_BTag * sf_nlo_weight; //SF_BTag is in evVarFriend, not sfFriend
    Double_t tmp = LUMI * weight * puw * SF_BTag ;
    if (hasSFfriend_flag != 0) tmp = tmp * SF_NLO_QCD * SF_NLO_EWK * SF_trigmetnomu;
    if (suffix == "ZJetsToNuNu" || suffix == "DYJetsToLL") return tmp / 1.23; //SF_BTag is in evVarFriend, not sfFriend
    else if (suffix == "WJetsToLNu") return tmp / 1.21; //SF_BTag is in evVarFriend, not sfFriend
    else return tmp; //SF_BTag is in evVarFriend, not sfFriend

  }

}

//=============================================== 


void vbfHiggsToInv_SignalRegion::fillEventMask(ULong64_t & eventMask) {

  vbfHiggsToInvAna::fillEventMask(eventMask);

  eventMask += muonLooseVetoC.addToMask(nMu10V < 0.5);
  eventMask += electronLooseVetoC.addToMask(nEle10V < 0.5);
  eventMask += gammaLooseVetoC.addToMask(nGamma15V < 0.5);
  if ( HLT_FLAG != 0) eventMask += HLTC.addToMask(HLT_MonoJetMetNoMuMHT90 > 0.5 || HLT_MonoJetMetNoMuMHT120 > 0.5 || HLT_Met170 > 0.5); // will modify this 
  //HLT_* variables are stored as float, so using "== 1" might yield unexpected results
  eventMask += recoilC.addToMask(metNoMu_pt > METNOLEP_START);

}


//===============================================

// void vbfHiggsToInv_SignalRegion::loop(vector< vector<Double_t> *> &yieldsVectorList, vector< vector<Double_t> *> &uncertaintyVectorList, vector< vector<Double_t> *> &efficiencyVectorList)

void vbfHiggsToInv_SignalRegion::loop(vector< Double_t > &yRow, vector< Double_t > &eRow, vector< Double_t > &uncRow, 
				vector< Double_t > &yRow_monoJ, vector< Double_t > &eRow_monoJ, vector< Double_t > &uncRow_monoJ,
				vector< Double_t > &yRow_monoV, vector< Double_t > &eRow_monoV, vector< Double_t > &uncRow_monoV)
{

   if (fChain == 0) return;

   setBranchStatusForAnalysis();
   setVarFromConfigFile();
   setSelections();
   setMask();

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

   // =====================================

   ofstream wfile("weights.txt",ios::out);

   if ( !wfile.is_open() ) {
     cout<<"Error: unable to open file weights.txt !"<<endl;
     exit(EXIT_FAILURE);
   }

   //======================================

   Long64_t nentries = fChain->GetEntriesFast();
   cout<<"vbfHiggsToInv_SignalRegion::loop()"<<endl;
   cout<<"nentries = "<<nentries<<endl;   

   Long64_t nbytes = 0, nb = 0;

   for (Int_t jentry=0; jentry<nentries; jentry++) {

     Long64_t ientry = LoadTree(jentry);
     if (ientry < 0) break;
     nb = fChain->GetEntry(jentry);   nbytes += nb;
     // if (Cut(ientry) < 0) continue;   

     ULong64_t eventMask = 0; 

     if (jentry%500000 == 0) {
       cout << "entry: " << jentry << endl;
     }

     newwgt = computeEventWeight();
     nTotalWeightedEvents += newwgt;  // counting events with weights

     if (suffix == "ZJetsToNuNu") wfile << jentry << " " << newwgt << endl;

     // beginning of eventMask building
     fillEventMask(eventMask);     
     // end of eventMask building

     analysisMask.countEvents(eventMask,newwgt);
     // analysisMask_monoJ.countEvents(eventMask,newwgt);
     // analysisMask_monoV.countEvents(eventMask,newwgt);

     if ( ((eventMask & analysisMask.globalMask.back()) == analysisMask.globalMask.back()) ) {
       
       // this histogram holds the final yields in bins of MET
       HYieldsMetBin->Fill(metNoMu_pt,newwgt);
	 
       HhtDistribution->Fill(htJet25,newwgt);
       HrecoilDistribution->Fill(metNoMu_pt,newwgt);
       HvtxDistribution->Fill(nVert,newwgt);
       HnjetsDistribution->Fill(nJetClean,newwgt);
       Hjet1etaDistribution->Fill(JetClean_eta[0],newwgt);
       Hjet1ptDistribution->Fill(JetClean_pt[0],newwgt);
       HjetMetDphiMinDistribution->Fill(dphijm,newwgt);
       HjetMetDphiMinAllJets->Fill(dphijmAllJets,newwgt);
       Hjet2etaDistribution->Fill(JetClean_eta[1],newwgt);
       Hjet2ptDistribution->Fill(JetClean_pt[1],newwgt);
       HvbfTaggedJets_mT->Fill(vbfJetsMT(),newwgt);
       HvbfTaggedJets_deltaEta->Fill(fabs(JetClean_eta[0]-JetClean_eta[1]),newwgt);
       HvbfTaggedJets_invMass->Fill(vbfJetsInvMass(),newwgt);
       Hjet1CHEF->Fill(Jet_chHEF[0],newwgt);
       if (nJetClean > 1) Hj1j2dphiDistribution->Fill(dphijj,newwgt);

       if (hasScaledHistograms_flag) {

	 HYieldsMetBin_qcdRenScaleUp->Fill(metNoMu_pt,(newwgt * SF_NLO_QCD_renScaleUp/ SF_NLO_QCD));
	 HYieldsMetBin_qcdRenScaleDown->Fill(metNoMu_pt,(newwgt * SF_NLO_QCD_renScaleDown/ SF_NLO_QCD));
	 HYieldsMetBin_qcdFacScaleUp->Fill(metNoMu_pt,(newwgt * SF_NLO_QCD_facScaleUp/ SF_NLO_QCD));
	 HYieldsMetBin_qcdFacScaleDown->Fill(metNoMu_pt,(newwgt * SF_NLO_QCD_facScaleDown/ SF_NLO_QCD));
	 HYieldsMetBin_qcdPdfUp->Fill(metNoMu_pt,(newwgt * SF_NLO_QCD_pdfUp/ SF_NLO_QCD));
	 HYieldsMetBin_qcdPdfDown->Fill(metNoMu_pt,(newwgt * SF_NLO_QCD_pdfDown/ SF_NLO_QCD));
	 HYieldsMetBin_ewkUp->Fill(metNoMu_pt,(newwgt * SF_NLO_EWK_up/ SF_NLO_EWK));
	 HYieldsMetBin_ewkDown->Fill(metNoMu_pt,(newwgt * SF_NLO_EWK_down/ SF_NLO_EWK));

       }


     } 

   }                        // end of loop on entries

   mySpaces(cout,2);
   selection::printSelectionFlowAndYields(cout, LUMI, nTotalWeightedEvents, &analysisMask);
   // mySpaces(cout,2);
   // selection::printSelectionFlowAndYields(cout, LUMI, nTotalWeightedEvents, &analysisMask_monoJ);
   // mySpaces(cout,2);
   // selection::printSelectionFlowAndYields(cout, LUMI, nTotalWeightedEvents, &analysisMask_monoV);

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
   // selection::printSelectionFlowAndYields(myfile, LUMI, nTotalWeightedEvents, &analysisMask_monoJ);
   // mySpaces(myfile,2);
   // selection::printSelectionFlowAndYields(myfile, LUMI, nTotalWeightedEvents, &analysisMask_monoV);
   // mySpaces(myfile,2);
   myPrintYieldsMetBinInStream(myfile, HYieldsMetBin, metBinEdgesVector.data(), nMetBins);

   myfile.close();


   fillRowVector(nTotalWeightedEvents, analysisSelectionManager, analysisMask, yRow, eRow, uncRow,0);
   // fillRowVector(nTotalWeightedEvents, analysisSelectionManager_monoJ, analysisMask_monoJ, yRow_monoJ, eRow_monoJ, uncRow_monoJ,0);
   // fillRowVector(nTotalWeightedEvents, analysisSelectionManager_monoV, analysisMask_monoV, yRow_monoV, eRow_monoV, uncRow_monoV,0);

   // filling last bin with overflow
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

   //cout << " ====  CHECK  ==== " << endl;

}




