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

using namespace std;

void doSkim(const string& path = "./",
	    const string& sampleName = "",
	    const string& dirPatternTag = "/",
	    const string& evVarFriend_path = "",
	    const Int_t has_sfFriend = 0,
	    const string& sfFriend_path = "",
	    const string& outpath = "./")
{

  // when using MET trees, directory for evVar friends can be 'friends_SR' or 'friends_VM', so we skim both together (base tree is the same)

  string dirPattern = dirPatternTag + "treeProducerDarkMatterMonoJet" + dirPatternTag;

  string infileName         = path +                                     sampleName + dirPattern + "tree.root";
  string infileFriendName   = path + evVarFriend_path + "evVarFriend_" + sampleName +                  ".root";
  string infileSfFriendName = path + sfFriend_path    + "sfFriend_"    + sampleName +                  ".root";

  string infileFriendName_bis = "";
  string evVarFriend_path_bis = "";

  if (evVarFriend_path == "friends_SR/") evVarFriend_path_bis = "friends_VM/";
  else if (evVarFriend_path == "friends_VM/") evVarFriend_path_bis = "friends_SR/";

  if (evVarFriend_path_bis != "") infileFriendName_bis = path + evVarFriend_path_bis+ "evVarFriend_" + sampleName + ".root";

  // create directories for output
  string createDirCommand = "";
  if (dirPatternTag == "/") createDirCommand = "mkdir -p " + outpath + sampleName + dirPattern; 

  string outfileName         = outpath +                                     sampleName + dirPattern + "tree.root";
  string outfileFriendName   = outpath + evVarFriend_path + "evVarFriend_" + sampleName +                  ".root";
  string outfileSfFriendName = outpath + sfFriend_path    + "sfFriend_"    + sampleName +                  ".root";


  cout << "======================" << endl;
  cout << "======================" << endl;
  cout << sampleName << endl;
  cout << "----------------------" << endl;


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
  Long64_t nentriesFriend = infriend->GetEntries();

  TFile* infileSfFriend = NULL;
  TTree* insffriend = NULL;
  Long64_t nentriesSfFriend = -1;

  if (has_sfFriend) {

    infileSfFriend = new TFile(infileSfFriendName.c_str());
    if (!infileSfFriend || infileSfFriend->IsZombie()) {
      std::cout << "Cannot open file " << infileSfFriendName << std::endl;
      exit(EXIT_FAILURE);
    }
    insffriend = (TTree*)infileSfFriend->Get("sf/t");
    if (!insffriend) insffriend = (TTree*)infileSfFriend->Get("sf/t");
    nentriesSfFriend = insffriend->GetEntries();

  }
  

  if (nentries != nentriesFriend) {
    
    std::cout << "Warning: input tree and evVar friend have different number of entries" << std::endl;
    std::cout << "End of programme" << std::endl;
    exit(EXIT_FAILURE);

  }

  if (has_sfFriend && nentries != nentriesSfFriend) {
    
    std::cout << "Warning: input tree and sf friend have different number of entries" << std::endl;
    std::cout << "End of programme" << std::endl;
    exit(EXIT_FAILURE);

  }

  intree->AddFriend(infriend);  // no need to attach sfFriend, evVarFriend is needed to cut on variables

  //intree->SetBranchStatus("*",1);
  Float_t recoil_pt;
  //Int_t *LepGood_pdgId = NULL;
  intree->SetBranchAddress("recoil_pt",&recoil_pt);
  //intree->SetBranchAddress("LepGood_pdgId",LepGood_pdgId);

  intree->SetName("tree");

  TFile* outfile = new TFile(outfileName.c_str(), "RECREATE");
  if (!outfile || outfile->IsZombie()) {
    std::cout << "Cannot open file " << outfileName << std::endl;
    exit(EXIT_FAILURE);
  }
  TTree* outtree = intree->CopyTree("recoil_pt > 200.0");

  cout << "Entries before skim = " << nentries << endl;
  cout << "Entries after skim = " << outtree->GetEntries() << "\t(" << (Double_t) outtree->GetEntries()/nentries << " % efficiency)" << endl;

  outfile->Write();
  outfile->Close();

  //std::cout << "==== check ====" << std::endl;

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
    if (recoil_pt > 200.0) {
	infriend->GetEntry(i);
	outfriend->Fill();
    }
  }

  outfileFriend->Write();
  outfileFriend->Close();

  TFile* outfileSfFriend = NULL;
  TDirectory *sfdir = NULL;
  TTree* outsffriend = NULL;

  if (has_sfFriend) {

    outfileSfFriend = new TFile(outfileSfFriendName.c_str(), "RECREATE");
    if (!outfileSfFriend || outfileSfFriend->IsZombie()) {
      std::cout << "Cannot open file " << outfileSfFriendName << std::endl;
      exit(EXIT_FAILURE);
    }
    sfdir = outfileSfFriend->mkdir("sf");
    sfdir->cd();
    outsffriend = insffriend->CloneTree(0);
    
    for (Long64_t i=0; i<nentries; i++) {
      intree->GetEntry(i);
      if (recoil_pt > 200.0) {
	insffriend->GetEntry(i);
	outsffriend->Fill();
      }
    }
    
    outfileSfFriend->Write();
    outfileSfFriend->Close();
    
  }

  infile->Close();
  infileFriend->Close();
  if (has_sfFriend) infileSfFriend->Close();



}



//====================================================

// use as:
// root -l -b -q checkTreeAndFriend.C++
// you can store the output inside a file. 
// to use from terminal and change inputs:
// root -l -b -q  'skimTreeAndFriend.C++("/u2/emanuele/","TREES_MET_80X_V4/", "/", "friends_SR/", "friends/","/u2/mciprian/tests/skimmedTrees/)'

void skimTreeAndFriend(const string& path             = "/u2/emanuele/",
		       const string& treeDir          = "TREES_MET_80X_V4/",
		       const string& dirPatternTag    = "/",
		       const string& evVarFriend_path = "friends_SR/", // can be 'friends', 'friends_SR', 'friends_VM' or 'friends_VE'  
		       const string& sfFriend_path    = "friends/",
		       const string& outpath          = "./")
{

  // when using MET trees, directory for evVar friends can be 'friends_SR' or 'friends_VM', so we skim both together (base tree is the same)

  // dirPatternTag says if base tree are to be taken inside /treeProducerDarkMatterMonoJet/ folder or if their name has just 'treeProducerDarkMatterMonoJet'
  // e.g. , can have:
  // <path>/QCD/treeProducerDarkMatterMonoJet/tree.root
  // or
  // <path>/QCD_treeProducerDarkMatterMonoJet_tree.root

  // evVarFriend_path and sfFriend_path tell say if friends are inside a folder (you pass the name) or at the same level of the base trees (default)
  // evVarFriend_path can be "", "friends", "friends_SR", "friends_VM" or "friends_VE"
  // sfFriend_path can be "" or "friends"

  cout << endl;
  cout << "Starting skim ..." << endl;
  cout << endl;

  // create directories
  string createDirCommand = "";
  if (outpath != "./") {
    createDirCommand = "mkdir -p " + outpath;
    system(createDirCommand.c_str());
  }
  createDirCommand = "mkdir -p " + outpath + evVarFriend_path;
  system(createDirCommand.c_str());
  createDirCommand = "mkdir -p " + outpath + sfFriend_path;
  system(createDirCommand.c_str());
  if (evVarFriend_path == "friends_SR/") { 
    createDirCommand = "mkdir -p " + outpath + "friends_VM/";
    system(createDirCommand.c_str());
  } else if (evVarFriend_path == "friends_VM/") {
    createDirCommand = "mkdir -p " + outpath + "friends_SR/";
    system(createDirCommand.c_str());
  } 

  cout << endl;
  cout << "Tree location: " << path + treeDir << endl;
  cout << "evVarFriend location: " << path + treeDir + evVarFriend_path << endl;
  cout << "sfFriend location: " << path + treeDir + sfFriend_path << endl;
  cout << endl;

  vector<string> sampleNameVector;

  /////////////
  // MC samples

  // Zvv QCD
  sampleNameVector.push_back("ZJetsToNuNu_HT100to200");
  sampleNameVector.push_back("ZJetsToNuNu_HT100to200_ext");
  sampleNameVector.push_back("ZJetsToNuNu_HT200to400");
  sampleNameVector.push_back("ZJetsToNuNu_HT200to400_ext");
  sampleNameVector.push_back("ZJetsToNuNu_HT400to600");
  sampleNameVector.push_back("ZJetsToNuNu_HT400to600_ext");
  sampleNameVector.push_back("ZJetsToNuNu_HT600to800");
  sampleNameVector.push_back("ZJetsToNuNu_HT800to1200");
  sampleNameVector.push_back("ZJetsToNuNu_HT800to1200_ext");
  sampleNameVector.push_back("ZJetsToNuNu_HT1200to2500");
  sampleNameVector.push_back("ZJetsToNuNu_HT1200to2500_ext");
  sampleNameVector.push_back("ZJetsToNuNu_HT2500toInf");
  sampleNameVector.push_back("ZJetsToNuNu_HT2500toInf_ext");
  // Wlv QCD
  sampleNameVector.push_back("WJetsToLNu_HT100to200");
  sampleNameVector.push_back("WJetsToLNu_HT100to200_ext");
  sampleNameVector.push_back("WJetsToLNu_HT200to400");
  sampleNameVector.push_back("WJetsToLNu_HT200to400_ext");
  sampleNameVector.push_back("WJetsToLNu_HT400to600");
  sampleNameVector.push_back("WJetsToLNu_HT400to600_ext");
  sampleNameVector.push_back("WJetsToLNu_HT600to800");
  sampleNameVector.push_back("WJetsToLNu_HT800to1200");
  sampleNameVector.push_back("WJetsToLNu_HT800to1200_ext");
  sampleNameVector.push_back("WJetsToLNu_HT1200to2500");
  sampleNameVector.push_back("WJetsToLNu_HT1200to2500_ext");
  sampleNameVector.push_back("WJetsToLNu_HT2500toInf");
  sampleNameVector.push_back("WJetsToLNu_HT2500toInf_ext");
  // Top
  sampleNameVector.push_back("TBar_tWch");
  sampleNameVector.push_back("TTJets");
  sampleNameVector.push_back("TToLeptons_sch");
  sampleNameVector.push_back("T_tWch");
  // Diboson
  sampleNameVector.push_back("ZZ");
  sampleNameVector.push_back("WW");
  sampleNameVector.push_back("WZ");
  // QCD
  sampleNameVector.push_back("QCD_HT200to300_ext");
  sampleNameVector.push_back("QCD_HT300to500");
  sampleNameVector.push_back("QCD_HT300to500_ext");
  sampleNameVector.push_back("QCD_HT500to700");
  sampleNameVector.push_back("QCD_HT500to700_ext");
  sampleNameVector.push_back("QCD_HT700to1000");
  sampleNameVector.push_back("QCD_HT700to1000_ext");
  sampleNameVector.push_back("QCD_HT1000to1500");
  sampleNameVector.push_back("QCD_HT1000to1500_ext");
  sampleNameVector.push_back("QCD_HT1500to2000");
  sampleNameVector.push_back("QCD_HT1500to2000_ext");
  sampleNameVector.push_back("QCD_HT2000toInf");
  sampleNameVector.push_back("QCD_HT2000toInf_ext");
  // Zll QCD
  sampleNameVector.push_back("DYJetsToLL_M50_HT100to200");
  sampleNameVector.push_back("DYJetsToLL_M50_HT100to200_ext");
  sampleNameVector.push_back("DYJetsToLL_M50_HT200to400");
  sampleNameVector.push_back("DYJetsToLL_M50_HT200to400_ext");
  sampleNameVector.push_back("DYJetsToLL_M50_HT400to600");
  sampleNameVector.push_back("DYJetsToLL_M50_HT400to600_ext");
  sampleNameVector.push_back("DYJetsToLL_M50_HT600toInf");
  sampleNameVector.push_back("DYJetsToLL_M50_HT600toInf_ext");
  // Gamma
  sampleNameVector.push_back("GJets_HT40to100");
  sampleNameVector.push_back("GJets_HT100to200");
  sampleNameVector.push_back("GJets_HT400to600");
  sampleNameVector.push_back("GJets_HT600toInf");
  // Zvv EWK
  sampleNameVector.push_back("EWKZToNuNu2Jets");
  // Wlv EWK
  sampleNameVector.push_back("EWKWMinus2Jets");
  sampleNameVector.push_back("EWKWPlus2Jets");
  // Zll EWK
  sampleNameVector.push_back("EWKZToLL2Jets");
  // signals only for SR or VM (not VE)
  if (evVarFriend_path.find("SR") != string::npos || evVarFriend_path.find("VM") != string::npos) {
    // GGH signal
    sampleNameVector.push_back("GluGlu_HToInvisible_M110");
    sampleNameVector.push_back("GluGlu_HToInvisible_M125");
    sampleNameVector.push_back("GluGlu_HToInvisible_M150");
    sampleNameVector.push_back("GluGlu_HToInvisible_M200");
    sampleNameVector.push_back("GluGlu_HToInvisible_M300");
    sampleNameVector.push_back("GluGlu_HToInvisible_M400");
    sampleNameVector.push_back("GluGlu_HToInvisible_M500");
    sampleNameVector.push_back("GluGlu_HToInvisible_M600");
    sampleNameVector.push_back("GluGlu_HToInvisible_M800");
    // VBF H signal
    sampleNameVector.push_back("VBF_HToInvisible_M110");
    sampleNameVector.push_back("VBF_HToInvisible_M125");
    sampleNameVector.push_back("VBF_HToInvisible_M150");
    sampleNameVector.push_back("VBF_HToInvisible_M200");
    sampleNameVector.push_back("VBF_HToInvisible_M300");
    sampleNameVector.push_back("VBF_HToInvisible_M400");
    sampleNameVector.push_back("VBF_HToInvisible_M500");
    sampleNameVector.push_back("VBF_HToInvisible_M600");
    sampleNameVector.push_back("VBF_HToInvisible_M800");
  }

  sampleNameVector.clear();
  sampleNameVector.push_back("ZJetsToNuNu_HT100to200");
  sampleNameVector.push_back("ZJetsToNuNu_HT100to200_ext");

  for (UInt_t i = 0; i < sampleNameVector.size(); i++) {
    doSkim(path+treeDir, sampleNameVector[i], dirPatternTag, evVarFriend_path, 1, sfFriend_path, outpath);
  }

  ////////
  // Data samples

  // vector<string> sampleNameDataVector;

  // if (evVarFriend_path.find("SR") != string::npos || evVarFriend_path.find("VM") != string::npos) {
  //   sampleNameDataVector.push_back("MET_Run2016B_PromptReco_v1_runs_272023_273146");
  //   sampleNameDataVector.push_back("MET_Run2016B_PromptReco_v2_runs_273150_275376");
  //   sampleNameDataVector.push_back("MET_Run2016C_PromptReco_v2_runs_275420_276283");
  //   sampleNameDataVector.push_back("MET_Run2016D_PromptReco_v2_runs_276315_276811");
  //   sampleNameDataVector.push_back("MET_Run2016E_PromptReco_v2_runs_276830_277420");
  //   sampleNameDataVector.push_back("MET_Run2016G_PromptReco_v1_runs_278817_279931");
  // } else if (evVarFriend_path.find("VE") != string::npos) {
  //   sampleNameDataVector.push_back("SingleElectron_Run2016B_PromptReco_v1_runs_272023_273146");
  //   sampleNameDataVector.push_back("SingleElectron_Run2016B_PromptReco_v2_runs_273150_275376");
  //   sampleNameDataVector.push_back("SingleElectron_Run2016C_PromptReco_v2_runs_275420_276283");
  //   sampleNameDataVector.push_back("SingleElectron_Run2016D_PromptReco_v2_runs_276315_276811");
  //   sampleNameDataVector.push_back("SingleElectron_Run2016E_PromptReco_v2_runs_276830_277420");
  //   sampleNameDataVector.push_back("SingleElectron_Run2016F_PromptReco_v1_runs_277820_278808");
  //   sampleNameDataVector.push_back("SingleElectron_Run2016G_PromptReco_v1_runs_278817_279931");
  // }

  // for (UInt_t i = 0; i < sampleNameDataVector.size(); i++) {
  //   doSkim(path+treeDir, sampleNameDataVector[i], dirPatternTag, evVarFriend_path, 0, sfFriend_path, outpath);
  // }



}
