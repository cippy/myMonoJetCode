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
#include <TH2.h>
#include <TKey.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TMatrixDSym.h>
#include <TMultiGraph.h>
#include <TPad.h>
#include <TPaletteAxis.h>
#include <TPaveStats.h>
#include <TPaveText.h>
#include <TProfile.h>
#include <TTreeIndex.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>
#include <TStyle.h>
#include <TVector3.h>
#include <TVirtualFitter.h>

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

#include <Rtypes.h> // to use kColor

#define NO_MC_WEIGHT 0
#define READ_FROM_LOCAL 1
#define ISTEST 0

using namespace std;

// =====================================================

string getSubDir() {

  string subDir = "genLep/";
  return subDir;

}

//======================================================

string getDirName() {

  string dirName = "/afs/cern.ch/user/m/mciprian/www/vbfHiggsToInv/genStudy/";
  dirName += getSubDir();
  return dirName;

}


//=======================================================

void createDir(const string dirName) {

  // system("cmsenv");
  system(Form("mkdir -p %s",dirName.c_str()));
  // system("currPath=$PWD");
  // system("cd /afs/cern.ch/user/m/mciprian/www/");
  // system("./copyphp.sh");
  // system("cd $currPath");

}

// =====================================================

void fillTH1(TH1F* histo, const Float_t val, const Double_t wgt = 1.0) {

  // embed overflow

  if(val < histo->GetXaxis()->GetBinLowEdge(histo->GetNbinsX()+1))
    histo->Fill(val,wgt);
  else
    histo->Fill(histo->GetXaxis()->GetBinCenter(histo->GetNbinsX()),wgt);

}


//======================================================

void buildChainWithFriend(TChain* chain, TChain* chFriend, TChain* chsfFriend, string sampleName) {
  
  cout << "Creating chain ..." << endl;
  
  vector<string> subSampleNameVector;

  if (sampleName == "WJetsToLNu") {
    if (ISTEST) {
      subSampleNameVector.push_back("WJetsToLNu_HT100to200");
      if (NO_MC_WEIGHT) subSampleNameVector.push_back("WJetsToLNu_HT100to200_ext");
    }
    subSampleNameVector.push_back("WJetsToLNu_HT200to400");
    if (NO_MC_WEIGHT) subSampleNameVector.push_back("WJetsToLNu_HT200to400_ext");
    subSampleNameVector.push_back("WJetsToLNu_HT400to600");
    if (NO_MC_WEIGHT) subSampleNameVector.push_back("WJetsToLNu_HT400to600_ext");
    subSampleNameVector.push_back("WJetsToLNu_HT600to800");
    subSampleNameVector.push_back("WJetsToLNu_HT800to1200");
    if (NO_MC_WEIGHT) subSampleNameVector.push_back("WJetsToLNu_HT800to1200_ext");
    subSampleNameVector.push_back("WJetsToLNu_HT1200to2500");
    if (NO_MC_WEIGHT) subSampleNameVector.push_back("WJetsToLNu_HT1200to2500_ext");
    subSampleNameVector.push_back("WJetsToLNu_HT2500toInf");
    if (NO_MC_WEIGHT) subSampleNameVector.push_back("WJetsToLNu_HT2500toInf_ext");
  }
  
  string treePath = "";
  if (READ_FROM_LOCAL) {
    if (ISTEST) treePath = "/u2/mciprian/TREES_1TIGHTELEP30_80X_4EoP/MC_WJETSHT_LO/";
    else        treePath = "/u2/mciprian/TREES_WJetsToLNu_noSkim/MCTrees/";
  }

  for(UInt_t i = 0; i < subSampleNameVector.size(); i++) {
  
    string treeRootFile = treePath + subSampleNameVector[i] + "/treeProducerDarkMatterMonoJet/tree.root";
    string friend_treeRootFile = treePath + "friends/evVarFriend_" + subSampleNameVector[i]+ ".root";
    string sffriend_treeRootFile = treePath + "friends/sfFriend_" + subSampleNameVector[i]+ ".root";

    chain->Add(TString(treeRootFile.c_str()));
    chFriend->Add(TString(friend_treeRootFile.c_str()));
    chsfFriend->Add(TString(sffriend_treeRootFile.c_str()));

  }

  cout << "Adding friend to chain ..." << endl;
  chain->AddFriend(chFriend);  //adding whole friend chain as friend                                                           
  chain->AddFriend(chsfFriend);  //adding whole friend chain as friend                

  if(!chain || !chFriend || !chsfFriend) {
    cout << "Error: chain not created. End of programme" << endl;
    exit(EXIT_FAILURE);
  }

  cout << "entries in chain      = " << chain->     GetEntries() << endl;
  cout << "entries in chFriend   = " << chFriend->  GetEntries() << endl;
  cout << "entries in chsfFriend = " << chsfFriend->GetEntries() << endl;

  if (chain->GetEntries() != chFriend->GetEntries() || chain->GetEntries() != chsfFriend->GetEntries()) {
    cout << "Error: chain and friends have different number of events. End of programme" << endl;
    exit(EXIT_FAILURE);
  } 

}

// ============================================================

Int_t findGenPartFromWLNu(const Int_t nGenParticles, 
			  TTreeReaderArray<Int_t>& particleId, 
			  TTreeReaderArray<Int_t>& particleMotherId) 
{

  // This function looks whether there are a lepton and a neutrino from a W in                                                                
  //the array. It works when using Emanuele's tree.                                                                                           
  // It also looks for the index of the mother in the list of particles, that is a W                                                          
  // The algorithm looks for a lepton and neutrino (or antineutrino) originating from mother in the list of generated particles.              
  // The function returns 1 if the search is successful, and 0 otherwise                                                                      

  Int_t genPartIndex = -1;  // will become the index of the lepton from  W if search is successful. Otherwise it is a negative value
  Int_t i = 0;      // index of the array                                                                                           
  Int_t Wid = 0;
  Int_t lepCharge = 0; 
  Int_t j = 0;
  Int_t Wfound = 0;

  //first look for a W (assuming there's only 1)    

  // PGD ID W+(W-)     = 24(-24)                            
  // PDG ID e-(e+)     = 11(-11)
  // PDG ID mu-(mu+)   = 13(-13)
  // PDG ID tau-(tau+) = 15(-15)

  Int_t lepPgdId[3] = {11, 13, 15};

  while (!Wfound && (j < nGenParticles)) {

    if (fabs(particleId[j]) == 24) {

      Wid = particleId[j];
      lepCharge =  (Wid > 0) ? -1 : 1;
      Wfound = 1;

    }
    j++;
  }

  if (Wfound) {

    while ((genPartIndex < 0) && (i < nGenParticles)) {

      if (particleMotherId[i] == Wid) {

	if (particleId[i] == (lepCharge * 11)) {
	  genPartIndex = i;
	} else if (particleId[i] == (lepCharge * 13)){
	  genPartIndex = i;
	} else if (particleId[i] == (lepCharge * 15)){
	  genPartIndex = i;
	} 

      }

      i++;

    }

  }

  return genPartIndex;

}

// ===========================================================

void drawDistribution(TH1F* h1, TH1F* h2, TH1F* h3, const string var = "", const string xAxisName = "") {

  gROOT->SetBatch(kTRUE);

  TH1::SetDefaultSumw2(); //all the following histograms will automatically call TH1::Sumw2()  

  TCanvas *canvas = new TCanvas("canvas","",700,700);
  TLegend *leg = new TLegend(0.5,0.7,0.79,0.90);

  canvas->SetGridy();

  h1->SetLineColor(kBlack);
  h2->SetLineColor(kRed);
  h3->SetLineColor(kBlue);

  h1->SetLineWidth(2);
  h2->SetLineWidth(2);
  h3->SetLineWidth(2);

  h1->GetXaxis()->SetTitle(xAxisName.c_str());
  h1->GetYaxis()->SetTitle("a.u");
  h1->GetYaxis()->SetTitleOffset(0.8);
  h1->GetYaxis()->SetTitleSize(0.06);

  h1->SetTitle("");
  h1->SetStats(0);

  h1->Scale(1./h1->Integral());
  h2->Scale(1./h2->Integral());
  h3->Scale(1./h3->Integral());
 
  Double_t maximumYaxisValue = -100000.0;
  Double_t minimumYaxisValue = 100000.0;
  maximumYaxisValue = TMath::Max(h1->GetBinContent(h1->GetMaximumBin()),h2->GetBinContent(h2->GetMaximumBin()));
  maximumYaxisValue = TMath::Max(maximumYaxisValue,h3->GetBinContent(h3->GetMaximumBin()));
  minimumYaxisValue = TMath::Min(h1->GetBinContent(h1->GetMinimumBin()),h2->GetBinContent(h2->GetMinimumBin()));
  minimumYaxisValue = TMath::Min(minimumYaxisValue,h3->GetBinContent(h3->GetMinimumBin()));
  // Double_t diff = maximumYaxisValue - minimumYaxisValue;
  // minimumYaxisValue -= diff * 0.1;
  // maximumYaxisValue += diff * 0.1;

  h1->SetMaximum(maximumYaxisValue * 1.2);
  h1->SetMinimum(0);

  h1->Draw("HE");
  h2->Draw("HE SAME");
  h3->Draw("HE SAME");

  leg->AddEntry(h1,"e","lf");
  leg->AddEntry(h2,"#mu","lf");
  leg->AddEntry(h3,"#tau","lf");
  leg->Draw();
  leg->SetMargin(0.3);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);  // transparent legend
  leg->SetFillColor(0);

  Int_t lastVisibleBin = h1->FindLastBinAbove(h1->GetMinimum());
  if (lastVisibleBin < h2->FindLastBinAbove(h1->GetMinimum())) lastVisibleBin = h2->FindLastBinAbove(h1->GetMinimum());
  if (lastVisibleBin < h3->FindLastBinAbove(h1->GetMinimum())) lastVisibleBin = h3->FindLastBinAbove(h1->GetMinimum());
  h1->GetXaxis()->SetRange(0,lastVisibleBin);

  //  cout << "check" << endl;

  canvas->SaveAs( (getDirName() + var + ".pdf").c_str() );
  canvas->SaveAs( (getDirName() + var + ".png").c_str() );

  // now log scale plot
  h1->SetMaximum( maximumYaxisValue * 1000 );
  h1->SetMinimum( 0.8 * TMath::Max( 0.0001, h1->GetMinimum() ) );
  lastVisibleBin = h1->FindLastBinAbove(h1->GetMinimum());
  if (lastVisibleBin < h2->FindLastBinAbove(h1->GetMinimum())) lastVisibleBin = h2->FindLastBinAbove(h1->GetMinimum());
  if (lastVisibleBin < h3->FindLastBinAbove(h1->GetMinimum())) lastVisibleBin = h3->FindLastBinAbove(h1->GetMinimum());
  h1->GetXaxis()->SetRange(0,lastVisibleBin);

  canvas->SetLogy();
  canvas->SaveAs( (getDirName() + var + "_logy.pdf").c_str() );
  canvas->SaveAs( (getDirName() + var + "_logy.png").c_str() );
  canvas->SetLogy(0);

  delete canvas;
  delete leg;  

}


//============================================================   


void doGenAna(TChain * chain) {


  /*
  Int_t           nGenPart;
  Int_t           GenPart_motherId[24];   //[nGenPart]                                                                                                                 
  Int_t           GenPart_grandmotherId[24];   //[nGenPart]                                                                                                            
  Int_t           GenPart_sourceId[24];   //[nGenPart]                                                                                                                 
  Int_t           GenPart_status[24];   //[nGenPart]                                                                                                                   
  Int_t           GenPart_pdgId[24];   //[nGenPart]                                                                                                                    
  Float_t         GenPart_pt[24];   //[nGenPart]                                                                                                                       
  Float_t         GenPart_eta[24];   //[nGenPart]                                                                                                                      
  Float_t         GenPart_phi[24];   //[nGenPart]                                                                                                                      
  Float_t         GenPart_mass[24];   //[nGenPart]                                                                                                                     
  Int_t           GenPart_motherIndex[24];   //[nGenPart] 
  */

  TH1::SetDefaultSumw2(); //all the following histograms will automatically call TH1::Sumw2() 

  Int_t GenPartIndex = -1;

  TTreeReader myReader(chain);
  
  TTreeReaderValue<Int_t> nGenPart(myReader,"nGenPart");

  // for weights
  TTreeReaderValue<Float_t> weight(myReader,"weight");
  TTreeReaderValue<Float_t> puw(myReader,"puw");
  TTreeReaderValue<Float_t> SF_NLO_QCD(myReader,"SF_NLO_QCD");
  TTreeReaderValue<Float_t> SF_NLO_EWK(myReader,"SF_NLO_EWK");
  TTreeReaderValue<Float_t> SF_BTag(myReader,"SF_BTag");

  TTreeReaderArray<Int_t> GenPart_pdgId(myReader,"GenPart_pdgId"); 
  TTreeReaderArray<Int_t> GenPart_motherId(myReader,"GenPart_motherId"); 
  TTreeReaderArray<Int_t> GenPart_motherIndex(myReader,"GenPart_motherIndex"); 

  TTreeReaderArray<Float_t> GenPart_pt(myReader,"GenPart_pt"); 
  TTreeReaderArray<Float_t> GenPart_eta(myReader,"GenPart_eta"); 
  TTreeReaderArray<Float_t> GenPart_phi(myReader,"GenPart_phi"); 
  TTreeReaderArray<Float_t> GenPart_mass(myReader,"GenPart_mass"); 
  
  Long64_t evtCounter = -1;

  TH1F* hElePt = new TH1F("hElePt","",100,0,1000);
  TH1F* hMuPt = new TH1F("hMuPt","",100,0,1000);
  TH1F* hTauPt = new TH1F("hTauPt","",100,0,1000);

  Double_t wgt = 0.0;

  while (myReader.Next()) {

    evtCounter++;
    if (evtCounter%500000 == 0) cout << "processing event : " << evtCounter << endl;;
   
    GenPartIndex = findGenPartFromWLNu(*nGenPart, GenPart_pdgId, GenPart_motherId);
    // if (GenPartIndex >= 0) {
    //   cout << "lep: PDG_ID = " << GenPart_pdgId[GenPartIndex] << "\tW: PDG_ID = " << GenPart_motherId[GenPartIndex];
    //   cout << "(" << GenPart_pdgId[GenPart_motherIndex[GenPartIndex]] << ")" << endl;  
    // } else {
    //   cout << "event : " << evtCounter << "\tGenPartIndex : " << GenPartIndex << endl;
    // }
    if (GenPartIndex < 0) continue;

    wgt = (*weight) * (*puw) * (*SF_NLO_QCD) * (*SF_NLO_EWK) * (*SF_BTag) / 1.21;
    
    Int_t lepPdgId = fabs(GenPart_pdgId[GenPartIndex]);

    if (lepPdgId == 11) {

      fillTH1(hElePt, GenPart_pt[GenPartIndex], wgt);
      //      hElePt->Fill(GenPart_pt[GenPartIndex],wgt);

    } else if (lepPdgId == 13) {

      fillTH1(hMuPt, GenPart_pt[GenPartIndex], wgt);

    } else if (lepPdgId == 15) {

      fillTH1(hTauPt, GenPart_pt[GenPartIndex], wgt);

    }
 
    
  }  // end of loop

  drawDistribution(hElePt, hMuPt, hTauPt, "pT", "lepton p_{T} [GeV]");

}


//=============================================================


void genLep() {

  TChain* chain = new TChain("tree");
  TChain* chFriend = new TChain("mjvars/t");
  TChain* chsfFriend = new TChain("sf/t");

  createDir(getDirName());

  buildChainWithFriend(chain, chFriend, chsfFriend, "WJetsToLNu");

  doGenAna(chain);

  delete chain;
  delete chFriend;
  delete chsfFriend;


}
