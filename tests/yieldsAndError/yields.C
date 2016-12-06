//ROOT header files                                                                                                          
#include <TROOT.h>
#include <TAttFill.h>
#include <TAxis.h>
#include <TCanvas.h>
#include <TChain.h>
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
#include <TObjArray.h>
#include <TObjString.h>
#include <TPad.h>
#include <TPaveStats.h>
#include <TPaveText.h>
#include <TStyle.h>
#include <TString.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
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

// declaring the following static when global variables means that other source codes accessing this one will not see them (they belong to this source)
// using them as if declared with define
static Double_t leadingJetPt  = 80;
static Double_t trailingJetPt = 50;
static Double_t missingEnergy = 150;
static Double_t jetMetDPhi    = 0.5;
static Double_t deltaEtaJJ    = 2.0;
static Double_t Mjj = 450.;


///////////////////////////////////////////////////////////////

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

void yields(const string dirpath = "80X/lumi_24p5fb/SignalRegion/") {
  
  string SRorCR = "";
  if (dirpath.find("SignalRegion") == string::npos) SRorCR = "SR";
  else if (dirpath.find("ControlRegion") == string::npos) SRorCR = "CR";
  else  {
    cout << "WARNING: not signal nor control region. Did you mispell the path?" << endl;
    exit(EXIT_FAILURE);
  }

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
  string path = "/afs/cern.ch/work/m/mciprian/myCode/CMSSW_8_0_19/src/myMonoJetCode/output/vbfHiggsToInv/" + dirpath + "vbfHiggsToInv_" + SRorCR + "_";
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


void entriesAfterSelection(const Double_t lumi = 1.0, const string reg = "sig") {

  TH1::SetDefaultSumw2();

  // reg can be: sig, zmm, zee, wen, wmn
  // lumi is in fb^-1

  cout << endl;
  cout << endl;

  ofstream wfile("weights.txt",ios::out);

  if ( !wfile.is_open() ) {
    cout<<"Error: unable to open file weights.txt !"<<endl;
    exit(EXIT_FAILURE);
  }

  string prefix = "root://eoscms//eos/cms";
  string path = "/store/cmst3/group/susy/emanuele/monox/trees/";
  if (reg != "wen" && reg != "zee") path += "TREES_MET_80X_V4/";
  else path += "TREES_1LEP_80X_V4/";

  vector<TString> allSamples = {"vbfH125=VBF_HToInvisible_M125",
				"ZJetsToNuNu=ZJetsToNuNu_HT100to200 ZJetsToNuNu_HT200to400 ZJetsToNuNu_HT400to600 ZJetsToNuNu_HT600to800 "
				"ZJetsToNuNu_HT800t1200 ZJetsToNuNu_HT1200to2500 ZJetsToNuNu_HT2500toInf",
				"Top=TBar_tWch TTJets TToLeptons_sch T_tWch"
  };

  for (UInt_t i = 0; i < allSamples.size(); i++) {

    TH1F *htest = new TH1F("htest","",1000,0.0,10000.0);
    
    cout << "Check the samples." << endl;

    TObjArray* array = allSamples[i].Tokenize("=");
    TString sampleName = ((TObjString *) array->At(0))->String();
    TString subSamples = ((TObjString *) array->At(1))->String();
    cout << sampleName << "\t" << subSamples << endl;
    TObjArray* array2 = subSamples.Tokenize(" "); 
    vector<TString> subSampleNameVector;

    for (Int_t j = 0; j < array2->GetEntries(); j++) {
      subSampleNameVector.push_back(((TObjString *) array2->At(j))->String());
    }    
    cout << endl;
    cout << "Check ended." << endl;
    cout << endl;

    TChain* chain = new TChain("tree");
    TChain* chFriend = new TChain("mjvars/t");
    TChain* chSfFriend = new TChain("sf/t");

    for(UInt_t i = 0; i < subSampleNameVector.size(); i++) {

      std::string treeRootFile = "";
      std::string friend_treeRootFile = "";
      std::string sf_friend_treeRootFile = "";

      treeRootFile = prefix + path + subSampleNameVector[i] + "_treeProducerDarkMatterMonoJet_tree.root";
      cout << "tree --> " << treeRootFile << endl;
      friend_treeRootFile = prefix + path + "evVarFriend_" + subSampleNameVector[i]+ ".root";
      sf_friend_treeRootFile = prefix + path + "sfFriend_" + subSampleNameVector[i]+ ".root";

      chain->Add(TString(treeRootFile.c_str()));
      chFriend->Add(TString(friend_treeRootFile.c_str()));
      chSfFriend->Add(TString(sf_friend_treeRootFile.c_str()));

    }
    
    //std::cout << "Adding friend to chain ..." << std::endl;
    chain->AddFriend(chFriend);  //adding whole friend chain as friend            
    //std::cout << "Adding friend with scale factors to chain ..." << std::endl;
    chain->AddFriend(chSfFriend);  //adding whole friend chain as friend          
      
    if(!chain || !chFriend || !chSfFriend) {
      std::cout << "Error: chain not created. End of programme" << std::endl;
      exit(EXIT_FAILURE);
    }

    TTreeReader myReader (chain);

    // add new branches uncomment the following
    // TTreeReaderValue<Float_t>  (myReader,""); 
    // TTreeReaderValue<Float_t>  (myReader,""); 

    // # of objects for veto (they are really float in the trees, not int)
    TTreeReaderValue<Float_t>  nMu10V(myReader,"nMu10V"); 
    TTreeReaderValue<Float_t>  nEle10V(myReader,"nEle10V"); 
    TTreeReaderValue<Float_t>  nGamma15V(myReader,"nGamma15V"); 
    TTreeReaderValue<Int_t>  nTauClean18V(myReader,"nTauClean18V"); 
    TTreeReaderValue<Float_t>  nBTag15(myReader,"nBTag15"); 

    // jets and met quantities
    TTreeReaderValue<Float_t>  dphijj(myReader,"dphijj"); 
    TTreeReaderValue<Int_t>    nJetClean(myReader,"nJetClean"); 
    TTreeReaderArray<Float_t>  JetClean_pt(myReader,"JetClean_pt"); 
    TTreeReaderArray<Float_t>  JetClean_eta(myReader,"JetClean_eta"); 
    TTreeReaderArray<Float_t>  Jet_chHEF(myReader,"Jet_chHEF"); 
    TTreeReaderValue<Float_t>  metNoMu_pt(myReader,"metNoMu_pt"); 
    TTreeReaderValue<Float_t>  dphijm(myReader,"dphijm"); 
    TTreeReaderArray<Float_t>  JetClean_leadClean(myReader,"JetClean_leadClean"); 

    // scale factors and weights
    TTreeReaderValue<Float_t>  weight(myReader,"weight"); 
    TTreeReaderValue<Float_t>  puw(myReader,"puw"); 
    TTreeReaderValue<Float_t>  SF_trig1lep(myReader,"SF_trig1lep"); 
    TTreeReaderValue<Float_t>  SF_trigmetnomu(myReader,"SF_trigmetnomu"); 
    TTreeReaderValue<Float_t>  SF_LepTightLoose(myReader,"SF_LepTightLoose"); 
    TTreeReaderValue<Float_t>  SF_LepTight(myReader,"SF_LepTight"); 
    TTreeReaderValue<Float_t>  SF_NLO_QCD(myReader,"SF_NLO_QCD"); 
    TTreeReaderValue<Float_t>  SF_NLO_EWK(myReader,"SF_NLO_EWK"); 
    TTreeReaderValue<Float_t>  SF_BTag(myReader,"SF_BTag"); 

    // met filters
    TTreeReaderValue<Int_t>  Flag_EcalDeadCellTriggerPrimitiveFilter(myReader,"Flag_EcalDeadCellTriggerPrimitiveFilter"); 
    TTreeReaderValue<Int_t>  Flag_HBHENoiseFilter(myReader,"Flag_HBHENoiseFilter"); 
    TTreeReaderValue<Int_t>  Flag_HBHENoiseIsoFilter(myReader,"Flag_HBHENoiseIsoFilter"); 
    TTreeReaderValue<Int_t>  Flag_goodVertices(myReader,"Flag_goodVertices"); 
    TTreeReaderValue<Int_t>  Flag_eeBadScFilter(myReader,"Flag_eeBadScFilter"); 
    TTreeReaderValue<Int_t>  Flag_globalTightHalo2016Filter(myReader,"Flag_globalTightHalo2016Filter"); // for now use the following
    //TTreeReaderValue<Int_t>  Flag_CSCTightHaloFilter(myReader,"Flag_CSCTightHaloFilter"); 

    string selection = "";
    //string selection = "Flag_EcalDeadCellTriggerPrimitiveFilter == 1 && Flag_HBHENoiseFilter == 1 && Flag_HBHENoiseIsoFilter == 1 && Flag_goodVertices == 1 && Flag_eeBadScFilter == 1 && Flag_globalTightHalo2016Filter == 1";
    // implement weights as a multiplication.  
    // string selection = "lumi * weight * puw * SF_NLO_QCD * SF_NLO_EWK";
    // if (reg == "sig") selection += " * SF_trigmetnomu";
    // if (sampleName == "ZJetsToNuNu" || sampleName == "DYJetsToLL") selection += "/ 1.23";
    // if (sampleName == "WJetsToLNu") selection += "/ 1.21";

    //    TH1F *evtAfterSel = new TH1F("evtAfterSel","events passing selection (0 if not passing it)",2,-0.5,1.5);
    // tcanvas only to use TTree::Draw() and get number of entries (no need to save the canvas, but without it the tree would be drawn in a GUI)
    // TCanvas *c = new TCanvas("c","",700,700);
    string tmp(Form("%2.1f",lumi));

    //string theSelection = "";
    //string theSelection = "(Flag_EcalDeadCellTriggerPrimitiveFilter == 1 && Flag_HBHENoiseFilter == 1 && Flag_HBHENoiseIsoFilter == 1 && Flag_goodVertices == 1 && Flag_eeBadScFilter == 1 && Flag_globalTightHalo2016Filter == 1) * (" + tmp + " * weight * puw * SF_NLO_QCD * SF_NLO_EWK)";
    string theSelection = "";
    if (sampleName == "ZJetsToNuNu") theSelection = "(" + tmp + " * weight * puw * SF_NLO_QCD * SF_NLO_EWK * SF_trigmetnomu / 1.23)";
    else theSelection = "(" + tmp + " * weight * puw * SF_NLO_QCD * SF_NLO_EWK * SF_trigmetnomu)";
    cout << "theSelection:\t" << theSelection << endl;

    chain->Draw("metNoMu_pt>>evtAfterSel",theSelection.c_str());
    Long64_t nentries = chain->GetEntries(selection.c_str());
    // delete c;

    cout << "sample:\t" << sampleName << endl;;
    // for (UInt_t k = 0; k < subSampleNameVector.size(); k++) {
    //   cout << " " << subSampleNameVector[k];
    // }
    //cout << endl;

    Double_t wgtEventsSum = 0.0;
    Double_t wgtEventsSum2 = 0.0;

    Int_t unwgtEventsSum = 0;
    Int_t nevt = 0;

    while(myReader.Next()){

      // if (*Flag_EcalDeadCellTriggerPrimitiveFilter != 1 || *Flag_HBHENoiseFilter != 1 || *Flag_HBHENoiseIsoFilter != 1 || *Flag_goodVertices != 1 || *Flag_eeBadScFilter != 1 || *Flag_globalTightHalo2016Filter != 1) continue;

      // if (*Flag_EcalDeadCellTriggerPrimitiveFilter == 1 && *Flag_HBHENoiseFilter == 1 && *Flag_HBHENoiseIsoFilter == 1 && *Flag_goodVertices == 1 && *Flag_eeBadScFilter == 1 && *Flag_globalTightHalo2016Filter == 1) {
      if (1) {
	
	//Int_t wgt = 2;
	Double_t wgt = lumi * *weight * *puw * *SF_NLO_QCD * *SF_NLO_EWK;
	if (reg == "sig" || reg == "zmm" || reg == "wmn") wgt *= *SF_trigmetnomu;  
	if(sampleName == "ZJetsToNuNu") {
	  wgt /= 1.23;
	  wfile << nevt << " " << wgt << endl;
	}
	nevt++;
	// if (sampleName == "ZJetsToNuNu" || sampleName == "DYJetsToLL") wgt /= 1.23; 
	// if (sampleName == "WJetsToLNu") wgt /= 1.21;

	wgtEventsSum += wgt;
	wgtEventsSum2 += wgt * wgt;
	unwgtEventsSum++;
	htest->Fill(*metNoMu_pt,wgt);

      }

    }

    TH1D *evtAfterSel = (TH1D*)gDirectory->Get("evtAfterSel");
    Double_t error = -1.0;
    cout << "evtAfterSel->IntegralAndError(allBins):\t" << evtAfterSel->IntegralAndError(0,evtAfterSel->GetNbinsX()+1,error); 
    cout << "+/-" << error << endl;
    cout << "nentries:\t" << nentries << endl;
    cout << "unwgtEventsSum:\t" << unwgtEventsSum << endl;
    cout << "wgtEventsSum:\t" << wgtEventsSum << " +/- " << sqrt(wgtEventsSum2) << endl;
    cout << "htest->Integral(allBins):" << "\t" << htest->IntegralAndError(0,htest->GetNbinsX()+1,error); 
    cout << " +/- " << error << endl;
    //cout << "htest->GetBinContent(htest->GetNbinsX()+1):" << htest->GetBinContent(htest->GetNbinsX()+1) << endl;
    delete htest;

    cout << endl;
    cout << endl;
    cout << endl;

  }

}
