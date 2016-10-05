#include <stdio.h>
#include <stdlib.h>
#include <cstdlib> //as stdlib.h
#include <cstdio>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TList.h>

//#include "edimarcoTree_v7.h"

void skim(const std::string selection = "") {

  TFile* infile = new TFile("tree.root");
  TTree* intree = (TTree*)infile->Get("tree");
  if (!intree) intree = (TTree*)infile->Get("tree");

  intree->SetName("tree");
  TFile* outfile = new TFile("skimmedTree.root", "RECREATE");
  TTree* outtree = intree->CopyTree(selection.c_str());
  outfile->Write();

  infile->Close();
  outfile->Close();

}


void skimWithFriend(const std::string sampleName = "", const std::string path = "/store/cmst3/group/susy/emanuele/monox/trees/TREES_MET_80X_V4/") {

  //TCut c1 = selection.c_str();

  std::string prefix = "root://eoscms//eos/cms";
  std::string infileName = path + sampleName + "_treeProducerDarkMatterMonoJet_tree.root";
  std::string infileFriendName = path  + "evVarFriend_" + sampleName + ".root";
  //std::string infileFriendName = "treeWith_myNlep.root";  // was used to test the skimming of the friend depending on a variable in tree.root
  std::string outfileName = sampleName + "_skim.root";
  std::string outfileFriendName = "evVarFriend_" + sampleName + "_skim.root";

  TFile* infile = new TFile(infileName.c_str());
  if (!infile || infile->IsZombie()) {
    std::cout << "Cannot open file " << infileName << std::endl;
    exit(EXIT_FAILURE);
  }
  TTree* intree = (TTree*)infile->Get("tree");
  if (!intree) intree = (TTree*)infile->Get("tree");
  Long64_t nentries = intree->GetEntries();

  TFile* infileFriend = new TFile(infileFriendName.c_str());
  if (!infileFriend || infileFriend->IsZombie()) {
    std::cout << "Cannot open file " << infileFriendName << std::endl;
    exit(EXIT_FAILURE);
  }
  TTree* infriend = (TTree*)infileFriend->Get("mjvars/t");
  if (!infriend) infriend = (TTree*)infileFriend->Get("mjvars/t");
  // TTree* infriend = (TTree*)infileFriend->Get("tree");          // was used to test the skimming of the friend depending on a variable in tree.root
  // if (!infriend) infriend = (TTree*)infileFriend->Get("tree");  // was used to test the skimming of the friend depending on a variable in tree.root
  Long64_t nentriesFriend = infriend->GetEntries();

  if (nentries != nentriesFriend) {
    
    std::cout << "Warning: input tree and friend have different number of entries" << std::endl;
    std::cout << "End of programme" << std::endl;
    exit(EXIT_FAILURE);

  }

  //intree->SetBranchStatus("*",1);
  Int_t nLepGood = 0;
  //Int_t *LepGood_pdgId = NULL;
  intree->SetBranchAddress("nLepGood",&nLepGood);
  //intree->SetBranchAddress("LepGood_pdgId",LepGood_pdgId);

  intree->SetName("tree");

  TFile* outfile = new TFile(outfileName.c_str(), "RECREATE");
  if (!outfile || outfile->IsZombie()) {
    std::cout << "Cannot open file " << outfileName << std::endl;
    exit(EXIT_FAILURE);
  }
  TTree* outtree = intree->CopyTree("nLepGood >= 2");
  outfile->Write();
  outfile->Close();

  std::cout << "==== check ====" << std::endl;

  TFile* outfileFriend = new TFile(outfileFriendName.c_str(), "RECREATE");
  if (!outfileFriend || outfileFriend->IsZombie()) {
    std::cout << "Cannot open file " << outfileFriendName << std::endl;
    exit(EXIT_FAILURE);
  }
  TDirectory *dir = outfileFriend->mkdir("mjvars");
  dir->cd();
  TTree* outfriend = infriend->CloneTree(0);

  for (Long64_t i=0; i<nentries; i++) {
    intree->GetEntry(i);
    if (nLepGood >= 2) {
      //if (fabs(LepGood_pdgId[0]) == 13 && fabs(LepGood_pdgId[1]) == 13) {
	infriend->GetEntry(i);
	outfriend->Fill();
	//}
    }
  }

  outfileFriend->Write();
  outfileFriend->Close();

  infile->Close();
  infileFriend->Close();

}


void addBranchToTree() {


  // don't understand why, but it also save the tree t to g file
   TFile *f = new TFile("tree.root");
   TTree *t = (TTree*)f->Get("tree");
   TFile *g = new TFile("treeWith_myNlep.root","RECREATE");
   Int_t nLepGood;
   Int_t myNlep;
   TTree *T = t->CloneTree();
   TBranch *bmyNlep = T->Branch("myNlep",&myNlep,"myNlep/I");
   T->SetBranchAddress("nLepGood",&nLepGood);
   Long64_t nentries = T->GetEntries();
   for (Long64_t i=0;i<nentries;i++) {
      t->GetEntry(i);
      myNlep = nLepGood;
      bmyNlep->Fill();
   }
   //T->Print();
   T->Write();
   g->Close();
   delete f;

}

void pruneTree(const std::string path = "/store/cmst3/group/susy/emanuele/monox/trees/TREES_MET_80X_V4/", const std::string prunedPath = "./") {

  std::vector<std::string> subSampleNameVector;
  subSampleNameVector.push_back("VBF_HToInvisible_M125");
  std::string prunedSampleName = "VBF_HToInvisible";

  //TCut c1 = selection.c_str();

  std::string prefix = "root://eoscms//eos/cms";

  TChain* chain = new TChain("tree");
  TChain* chFriend = new TChain("mjvars/t");
  TChain* chSfFriend = new TChain("sf/t");

  for(UInt_t i = 0; i < subSampleNameVector.size(); i++) {

    std::string treeRootFile = "";
    std::string friend_treeRootFile = "";
    std::string sf_friend_treeRootFile = "";

    treeRootFile = prefix + path + subSampleNameVector[i] + "_treeProducerDarkMatterMonoJet_tree.root";
    friend_treeRootFile = prefix + path + "evVarFriend_" + subSampleNameVector[i]+ ".root";
    sf_friend_treeRootFile = prefix + path + "sfFriend_" + subSampleNameVector[i]+ ".root";
   
    chain->Add(TString(treeRootFile.c_str()));
    chFriend->Add(TString(friend_treeRootFile.c_str()));
    chSfFriend->Add(TString(sf_friend_treeRootFile.c_str()));

  }

  std::cout << "Adding friend to chain ..." << std::endl;
  chain->AddFriend(chFriend);  //adding whole friend chain as friend                                                                                            
  std::cout << "Adding friend with scale factors to chain ..." << std::endl;
  chain->AddFriend(chSfFriend);  //adding whole friend chain as friend                                                                                        

  if(!chain) {
    std::cout << "Error: chain not created. End of programme" << std::endl;
    exit(EXIT_FAILURE);
  }

  chain->SetBranchStatus("*",0);

  chain->SetBranchStatus("nMu10V",1);  // # of muons passing loose selection
  chain->SetBranchStatus("nEle10V",1);  // # of electrons passing loose selection for electron veto
  chain->SetBranchStatus("nGamma15V",1);  // # of photons passing loose selection for photon veto
  //chain->SetBranchStatus("nMu20T",1);  // # of muons passing tight selection (pt > 20 + everything else)
  //chain->SetBranchStatus("nTau18V",1);
  chain->SetBranchStatus("nTauClean18V",1);

  chain->SetBranchStatus("dphijj",1);          // dphi between 1st and 2nd jet, 999 if second jet doesn't exist
  chain->SetBranchStatus("nJetClean",1);    // # of jet with pt > 30 & eta < 4.7 and cleaning for against muons misidentified as PFjets   
  chain->SetBranchStatus("JetClean_pt",1);  
  chain->SetBranchStatus("JetClean_eta",1);  
  chain->SetBranchStatus("JetClean_phi",1);  
   
  // chain->SetBranchStatus("nLepGood",1);
  // chain->SetBranchStatus("LepGood_pdgId",1);  // must be 13 for muons ( -13 for mu+), 11 for electrons and 15 for taus
  // chain->SetBranchStatus("LepGood_pt",1);
  // chain->SetBranchStatus("LepGood_eta",1);
  // chain->SetBranchStatus("LepGood_phi",1);   
  // chain->SetBranchStatus("LepGood_mass",1);
  // chain->SetBranchStatus("LepGood_tightId",1);
  // chain->SetBranchStatus("LepGood_relIso04",1);

  chain->SetBranchStatus("metNoMu_pt",1);
  chain->SetBranchStatus("metNoMu_eta",1);
  chain->SetBranchStatus("metNoMu_phi",1);
  chain->SetBranchStatus("htJet25",1);

  chain->SetBranchStatus("nVert",1);  // number of good vertices 

  chain->SetBranchStatus("Flag_EcalDeadCellTriggerPrimitiveFilter",1);
  chain->SetBranchStatus("Flag_HBHENoiseFilter",1);
  chain->SetBranchStatus("Flag_HBHENoiseIsoFilter",1);
  chain->SetBranchStatus("Flag_goodVertices",1);
  chain->SetBranchStatus("Flag_eeBadScFilter",1);
  chain->SetBranchStatus("Flag_globalTightHalo2016Filter",1);
  
  chain->SetBranchStatus("nBTag15",1);  // for b-jet veto
  chain->SetBranchStatus("dphijm",1);          // flag for dphi minimum between met and any of the jets in the event (using only the first four jets
  chain->SetBranchStatus("weight",1);   // modified since 17 November 2015: now it includes the whol weight, e.g. 1000*xsec*genWeight ...
  chain->SetBranchStatus("JetClean_leadClean",1); // has new cleaning on energy fractions (added on 17 November 2015) 

  chain->SetBranchStatus("dphijmAllJets",1);   
  chain->SetBranchStatus("vbfTaggedJet_deltaEta",1);   
  chain->SetBranchStatus("vbfTaggedJet_invMass",1);   
  chain->SetBranchStatus("vbfTaggedJet_leadJetPt",1);   
  chain->SetBranchStatus("vbfTaggedJet_trailJetPt",1);   
  chain->SetBranchStatus("vbfTaggedJet_leadJetEta",1);   
  chain->SetBranchStatus("vbfTaggedJet_trailJetEta",1);   

  chain->SetBranchStatus("nEle40T",1);
  chain->SetBranchStatus("puw",1); // added on 21 Sept 2016, substituting vtxWeight
 
  chain->SetBranchStatus("SF_BTag",1);
  //variables in sfFriend tree
  chain->SetBranchStatus("SF_trig1lep",1);
  chain->SetBranchStatus("SF_trigmetnomu",1);
  chain->SetBranchStatus("SF_LepTightLoose",1);
  chain->SetBranchStatus("SF_LepTightLooseUp",1);
  chain->SetBranchStatus("SF_LepTightLooseDown",1);
  chain->SetBranchStatus("SF_LepTight",1);
  chain->SetBranchStatus("SF_LepTightUp",1);
  chain->SetBranchStatus("SF_LepTightDown",1);
  chain->SetBranchStatus("SF_NLO_QCD",1);
  chain->SetBranchStatus("SF_NLO_QCD_renScaleUp",1);
  chain->SetBranchStatus("SF_NLO_QCD_renScaleDown",1);
  chain->SetBranchStatus("SF_NLO_QCD_facScaleUp",1);
  chain->SetBranchStatus("SF_NLO_QCD_facScaleDown",1);
  chain->SetBranchStatus("SF_NLO_QCD_pdfUp",1);
  chain->SetBranchStatus("SF_NLO_QCD_pdfDown",1);
  chain->SetBranchStatus("SF_NLO_EWK",1);
  chain->SetBranchStatus("SF_NLO_EWK_up",1);
  chain->SetBranchStatus("SF_NLO_EWK_down",1);
  
  std::string outfileName = "";

  // for(UInt_t i = 0; i < subSampleNameVector.size(); i++) {

  //   outfileName = prunedPath + subSampleNameVector[i] + "_pruned.root";

  //   TFile *newfile = new TFile(outfileName,"recreate");
  //   TTree *newtree = chain->CloneTree(0);
  //   newtree->CopyEntries(chain);

  //   newtree->Print();
  //   newfile->Write();
  //   delete newfile;

  // }

  outfileName = prunedPath + prunedSampleName + "_pruned.root";

  TFile *newfile = new TFile(outfileName.c_str(),"recreate");
  TTree *newtree = chain->CloneTree(0);

  //  newtree->CopyEntries(chain);
  
  newtree->Print();
  newfile->Write();
  //  delete newfile;

  // delete chain;
  // delete chFriend;
  // delete chSfFriend;

  std::cout << std::endl; 
  std::cout << "THE END" << std::endl; 
  std::cout << std::endl; 

}


// ==============================


void mergeTree(const std::string path = "/store/cmst3/group/susy/emanuele/monox/trees/TREES_MET_80X_V4/", const std::string mergedPath = "./") {

  /// NOT WORKING 

  std::vector<std::string> subSampleNameVector;
  subSampleNameVector.push_back("VBF_HToInvisible_M125");
  std::string mergedSampleName = "VBF_HToInvisible";

  //TCut c1 = selection.c_str();

  std::string prefix = "root://eoscms//eos/cms";

  TChain* chain = new TChain("tree");
  TChain* chFriend = new TChain("mjvars/t");
  TChain* chSfFriend = new TChain("sf/t");

  for(UInt_t i = 0; i < subSampleNameVector.size(); i++) {

    std::string treeRootFile = "";
    std::string friend_treeRootFile = "";
    std::string sf_friend_treeRootFile = "";

    treeRootFile = prefix + path + subSampleNameVector[i] + "_treeProducerDarkMatterMonoJet_tree.root";
    friend_treeRootFile = prefix + path + "evVarFriend_" + subSampleNameVector[i]+ ".root";
    sf_friend_treeRootFile = prefix + path + "sfFriend_" + subSampleNameVector[i]+ ".root";
   
    chain->Add(TString(treeRootFile.c_str()));
    chFriend->Add(TString(friend_treeRootFile.c_str()));
    chSfFriend->Add(TString(sf_friend_treeRootFile.c_str()));

  }

  // std::cout << "Adding friend to chain ..." << std::endl;
  // chain->AddFriend(chFriend);  //adding whole friend chain as friend                                                                                            
  // std::cout << "Adding friend with scale factors to chain ..." << std::endl;
  // chain->AddFriend(chSfFriend);  //adding whole friend chain as friend                                                                                        

  if(!chain || !chFriend || !chSfFriend) {
    std::cout << "Error: chain not created. End of programme" << std::endl;
    exit(EXIT_FAILURE);
  }

  chain->SetBranchStatus("*",0);
  chFriend->SetBranchStatus("*",0);
  chSfFriend->SetBranchStatus("*",0);

  // variables in chain
  
  // chain->SetBranchStatus("nLepGood",1);
  // chain->SetBranchStatus("LepGood_pdgId",1);  // must be 13 for muons ( -13 for mu+), 11 for electrons and 15 for taus
  // chain->SetBranchStatus("LepGood_pt",1);
  // chain->SetBranchStatus("LepGood_eta",1);
  // chain->SetBranchStatus("LepGood_phi",1);   
  // chain->SetBranchStatus("LepGood_mass",1);
  // chain->SetBranchStatus("LepGood_tightId",1);
  // chain->SetBranchStatus("LepGood_relIso04",1);

  chain->SetBranchStatus("metNoMu_pt",1);
  chain->SetBranchStatus("metNoMu_eta",1);
  chain->SetBranchStatus("metNoMu_phi",1);
  chain->SetBranchStatus("htJet25",1);

  chain->SetBranchStatus("nVert",1);  // number of good vertices 

  chain->SetBranchStatus("Flag_EcalDeadCellTriggerPrimitiveFilter",1);
  chain->SetBranchStatus("Flag_HBHENoiseFilter",1);
  chain->SetBranchStatus("Flag_HBHENoiseIsoFilter",1);
  chain->SetBranchStatus("Flag_goodVertices",1);
  chain->SetBranchStatus("Flag_eeBadScFilter",1);
  chain->SetBranchStatus("Flag_globalTightHalo2016Filter",1);
  
  // variables in evVarFriend
  
  chFriend->SetBranchStatus("nMu10V",1);  // # of muons passing loose selection
  chFriend->SetBranchStatus("nEle10V",1);  // # of electrons passing loose selection for electron veto
  chFriend->SetBranchStatus("nGamma15V",1);  // # of photons passing loose selection for photon veto
  //chFriend->SetBranchStatus("nMu20T",1);  // # of muons passing tight selection (pt > 20 + everything else)
  //chFriend->SetBranchStatus("nTau18V",1);
  chFriend->SetBranchStatus("nTauClean18V",1);

  chFriend->SetBranchStatus("dphijj",1);          // dphi between 1st and 2nd jet, 999 if second jet doesn't exist
  chFriend->SetBranchStatus("nJetClean",1);    // # of jet with pt > 30 & eta < 4.7 and cleaning for against muons misidentified as PFjets   
  chFriend->SetBranchStatus("JetClean_pt",1);  
  chFriend->SetBranchStatus("JetClean_eta",1);  
  chFriend->SetBranchStatus("JetClean_phi",1);  
  chFriend->SetBranchStatus("nBTag15",1);  // for b-jet veto
  chFriend->SetBranchStatus("dphijm",1);          // flag for dphi minimum between met and any of the jets in the event (using only the first four jets
  chFriend->SetBranchStatus("weight",1);   // modified since 17 November 2015: now it includes the whol weight, e.g. 1000*xsec*genWeight ...
  chFriend->SetBranchStatus("JetClean_leadClean",1); // has new cleaning on energy fractions (added on 17 November 2015) 

  chFriend->SetBranchStatus("dphijmAllJets",1);   
  chFriend->SetBranchStatus("vbfTaggedJet_deltaEta",1);   
  chFriend->SetBranchStatus("vbfTaggedJet_invMass",1);   
  chFriend->SetBranchStatus("vbfTaggedJet_leadJetPt",1);   
  chFriend->SetBranchStatus("vbfTaggedJet_trailJetPt",1);   
  chFriend->SetBranchStatus("vbfTaggedJet_leadJetEta",1);   
  chFriend->SetBranchStatus("vbfTaggedJet_trailJetEta",1);   

  chFriend->SetBranchStatus("nEle40T",1);
  chFriend->SetBranchStatus("puw",1); // added on 21 Sept 2016, substituting vtxWeight
 
  chFriend->SetBranchStatus("SF_BTag",1); 

  //variables in sfFriend tree
  chSfFriend->SetBranchStatus("SF_trig1lep",1);
  chSfFriend->SetBranchStatus("SF_trigmetnomu",1);
  chSfFriend->SetBranchStatus("SF_LepTightLoose",1);
  chSfFriend->SetBranchStatus("SF_LepTightLooseUp",1);
  chSfFriend->SetBranchStatus("SF_LepTightLooseDown",1);
  chSfFriend->SetBranchStatus("SF_LepTight",1);
  chSfFriend->SetBranchStatus("SF_LepTightUp",1);
  chSfFriend->SetBranchStatus("SF_LepTightDown",1);
  chSfFriend->SetBranchStatus("SF_NLO_QCD",1);
  chSfFriend->SetBranchStatus("SF_NLO_QCD_renScaleUp",1);
  chSfFriend->SetBranchStatus("SF_NLO_QCD_renScaleDown",1);
  chSfFriend->SetBranchStatus("SF_NLO_QCD_facScaleUp",1);
  chSfFriend->SetBranchStatus("SF_NLO_QCD_facScaleDown",1);
  chSfFriend->SetBranchStatus("SF_NLO_QCD_pdfUp",1);
  chSfFriend->SetBranchStatus("SF_NLO_QCD_pdfDown",1);
  chSfFriend->SetBranchStatus("SF_NLO_EWK",1);
  chSfFriend->SetBranchStatus("SF_NLO_EWK_up",1);
  chSfFriend->SetBranchStatus("SF_NLO_EWK_down",1);
  
  TList *list = new TList;
  list->Add(chain);
  list->Add(chFriend);
  list->Add(chSfFriend);

  std::string outfileName = "";
  outfileName = mergedPath + mergedSampleName + "_merged.root";
  TFile *newfile = new TFile(outfileName.c_str(),"recreate");

  TTree *newtree = TTree::MergeTrees(list);
  newtree->SetName("tree");
  newtree->Write();

  //  delete newfile;

  // delete chain;
  // delete chFriend;
  // delete chSfFriend;

  std::cout << std::endl; 
  std::cout << "THE END" << std::endl; 
  std::cout << std::endl; 

}


// ====================================================================


void pruneAndAddBranchToTree(const std::string path = "/store/cmst3/group/susy/emanuele/monox/trees/TREES_MET_80X_V4/", const std::string outPath = "./", const Int_t Z0_or_W1 = -1) {

  /// NOT WORKING 

  std::vector<std::string> subSampleNameVector;
  subSampleNameVector.push_back("VBF_HToInvisible_M125");
  std::string outSampleName = "VBF_HToInvisible";

  //TCut c1 = selection.c_str();

  std::string prefix = "root://eoscms//eos/cms";

  TChain* chain = new TChain("tree");
  TChain* chFriend = new TChain("mjvars/t");
  TChain* chSfFriend = new TChain("sf/t");

  for(UInt_t i = 0; i < subSampleNameVector.size(); i++) {

    std::string treeRootFile = "";
    std::string friend_treeRootFile = "";
    std::string sf_friend_treeRootFile = "";

    treeRootFile = prefix + path + subSampleNameVector[i] + "_treeProducerDarkMatterMonoJet_tree.root";
    friend_treeRootFile = prefix + path + "evVarFriend_" + subSampleNameVector[i]+ ".root";
    sf_friend_treeRootFile = prefix + path + "sfFriend_" + subSampleNameVector[i]+ ".root";
   
    chain->Add(TString(treeRootFile.c_str()));
    chFriend->Add(TString(friend_treeRootFile.c_str()));
    chSfFriend->Add(TString(sf_friend_treeRootFile.c_str()));

  }

  std::cout << "Adding friend to chain ..." << std::endl;
  chain->AddFriend(chFriend);  //adding whole friend chain as friend                                                                                            
  std::cout << "Adding friend with scale factors to chain ..." << std::endl;
  chain->AddFriend(chSfFriend);  //adding whole friend chain as friend                                                                                        

  if(!chain || !chFriend || !chSfFriend) {
    std::cout << "Error: chain not created. End of programme" << std::endl;
    exit(EXIT_FAILURE);
  }

  chain->SetBranchStatus("*",0);

  // variables in chain
  
  // chain->SetBranchStatus("nLepGood",1);
  // chain->SetBranchStatus("LepGood_pdgId",1);  // must be 13 for muons ( -13 for mu+), 11 for electrons and 15 for taus
  // chain->SetBranchStatus("LepGood_pt",1);
  // chain->SetBranchStatus("LepGood_eta",1);
  // chain->SetBranchStatus("LepGood_phi",1);   
  // chain->SetBranchStatus("LepGood_mass",1);
  // chain->SetBranchStatus("LepGood_tightId",1);
  // chain->SetBranchStatus("LepGood_relIso04",1);

  chain->SetBranchStatus("metNoMu_pt",1);
  chain->SetBranchStatus("metNoMu_eta",1);
  chain->SetBranchStatus("metNoMu_phi",1);
  chain->SetBranchStatus("htJet25",1);

  chain->SetBranchStatus("nVert",1);  // number of good vertices 

  chain->SetBranchStatus("Flag_EcalDeadCellTriggerPrimitiveFilter",1);
  chain->SetBranchStatus("Flag_HBHENoiseFilter",1);
  chain->SetBranchStatus("Flag_HBHENoiseIsoFilter",1);
  chain->SetBranchStatus("Flag_goodVertices",1);
  chain->SetBranchStatus("Flag_eeBadScFilter",1);
  chain->SetBranchStatus("Flag_globalTightHalo2016Filter",1);
  
  // variables in evVarFriend
  
  chain->SetBranchStatus("nMu10V",1);  // # of muons passing loose selection
  chain->SetBranchStatus("nEle10V",1);  // # of electrons passing loose selection for electron veto
  chain->SetBranchStatus("nGamma15V",1);  // # of photons passing loose selection for photon veto
  //chain->SetBranchStatus("nMu20T",1);  // # of muons passing tight selection (pt > 20 + everything else)
  //chain->SetBranchStatus("nTau18V",1);
  chain->SetBranchStatus("nTauClean18V",1);

  chain->SetBranchStatus("dphijj",1);          // dphi between 1st and 2nd jet, 999 if second jet doesn't exist
  chain->SetBranchStatus("nJetClean",1);    // # of jet with pt > 30 & eta < 4.7 and cleaning for against muons misidentified as PFjets   
  chain->SetBranchStatus("JetClean_pt",1);  
  chain->SetBranchStatus("JetClean_eta",1);  
  chain->SetBranchStatus("JetClean_phi",1);  
  chain->SetBranchStatus("nBTag15",1);  // for b-jet veto
  chain->SetBranchStatus("dphijm",1);          // flag for dphi minimum between met and any of the jets in the event (using only the first four jets
  chain->SetBranchStatus("weight",1);   // modified since 17 November 2015: now it includes the whol weight, e.g. 1000*xsec*genWeight ...
  chain->SetBranchStatus("JetClean_leadClean",1); // has new cleaning on energy fractions (added on 17 November 2015) 

  chain->SetBranchStatus("dphijmAllJets",1);   
  chain->SetBranchStatus("vbfTaggedJet_deltaEta",1);   
  chain->SetBranchStatus("vbfTaggedJet_invMass",1);   
  chain->SetBranchStatus("vbfTaggedJet_leadJetPt",1);   
  chain->SetBranchStatus("vbfTaggedJet_trailJetPt",1);   
  chain->SetBranchStatus("vbfTaggedJet_leadJetEta",1);   
  chain->SetBranchStatus("vbfTaggedJet_trailJetEta",1);   

  chain->SetBranchStatus("nEle40T",1);
  chain->SetBranchStatus("puw",1); // added on 21 Sept 2016, substituting vtxWeight
 
  chain->SetBranchStatus("SF_BTag",1); 

  //variables in sfFriend tree
  chain->SetBranchStatus("SF_trig1lep",1);
  chain->SetBranchStatus("SF_trigmetnomu",1);
  chain->SetBranchStatus("SF_LepTightLoose",1);
  chain->SetBranchStatus("SF_LepTightLooseUp",1);
  chain->SetBranchStatus("SF_LepTightLooseDown",1);
  chain->SetBranchStatus("SF_LepTight",1);
  chain->SetBranchStatus("SF_LepTightUp",1);
  chain->SetBranchStatus("SF_LepTightDown",1);
  chain->SetBranchStatus("SF_NLO_QCD",1);
  chain->SetBranchStatus("SF_NLO_QCD_renScaleUp",1);
  chain->SetBranchStatus("SF_NLO_QCD_renScaleDown",1);
  chain->SetBranchStatus("SF_NLO_QCD_facScaleUp",1);
  chain->SetBranchStatus("SF_NLO_QCD_facScaleDown",1);
  chain->SetBranchStatus("SF_NLO_QCD_pdfUp",1);
  chain->SetBranchStatus("SF_NLO_QCD_pdfDown",1);
  chain->SetBranchStatus("SF_NLO_EWK",1);
  chain->SetBranchStatus("SF_NLO_EWK_up",1);
  chain->SetBranchStatus("SF_NLO_EWK_down",1);
  
  Float_t weight, puw, SF_NLO_QCD, SF_NLO_EWK, SF_trigmetnomu; 
  chain->SetBranchAddress("weight",&weight);
  chain->SetBranchAddress("puw",&puw);
  chain->SetBranchAddress("SF_NLO_QCD",&SF_NLO_QCD);
  chain->SetBranchAddress("SF_NLO_EWK",&SF_NLO_EWK);
  chain->SetBranchAddress("SF_trigmetnomu",&SF_trigmetnomu);

  std::string outfileName = "";
  outfileName = outPath + outSampleName + "_out.root";
  TFile *newfile = new TFile(outfileName.c_str(),"recreate");

  TTree *newtree = chain->CloneTree(0);

  Float_t totweight = -999;
  // adding branch
  newtree->Branch("totweight", &totweight, "total weight/F");
  int nentries= chain->GetEntries();
  cout<<nentries<<endl;

  for( int i=0; i < nentries; i++){
    chain->GetEvent(i);
    totweight = weight * puw * SF_NLO_QCD * SF_NLO_EWK * SF_trigmetnomu;
    if (Z0_or_W1 == 0) totweight/= 1.23;
    else if (Z0_or_W1 == 1) totweight/= 1.21;
    newtree->Fill();
  } 

  newfile->Write();

  //  delete newfile;

  // delete chain;
  // delete chFriend;
  // delete chSfFriend;

  std::cout << std::endl; 
  std::cout << "THE END" << std::endl; 
  std::cout << std::endl; 

}
