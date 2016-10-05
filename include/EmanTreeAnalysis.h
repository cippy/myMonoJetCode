#ifndef EmanTree_Analysis_h
#define EmanTree_Analysis_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <TVector2.h>
#include <TLorentzVector.h>

#include <iostream>
#include <cstdlib>
#include <vector>

#include "AnalysisDarkMatter.h"

namespace myAnalyzerTEman {


  //====================================================

  class vbfHiggsToInvAna : public AnalysisDarkMatter {
  public:

    vbfHiggsToInvAna(TTree *tree); 
    virtual ~vbfHiggsToInvAna() { std::cout<<"~vbfHiggsToInvAna() called"<<std::endl; }
    
    virtual void setSelections(); 
    virtual void setHistograms();
    virtual void setHistogramLastBinAsOverFlow(const Int_t);
    virtual void setNumberParameterValue(const std::string, const Double_t);
    virtual void setVarFromConfigFile();
    virtual void createSystematicsHistogram();
    virtual void fillEventMask(ULong64_t &); // method to set eventMask event-by-event depending on some selections
    virtual void loop(std::vector< Double_t > &, std::vector< Double_t > &, std::vector< Double_t > &, std::vector< Double_t > &, std::vector< Double_t > &, std::vector< Double_t > &, std::vector< Double_t > &, std::vector< Double_t > &, std::vector< Double_t > &) {};
    virtual Double_t vbfJetsMT();
    virtual Double_t vbfJetsInvMass();
   
    selection vbfTaggedJets_deltaEtaC;  // eta of two jets and their separation
    selection vbfTaggedJets_inVMassC;
    selection vbfTaggedJets_jetsPtC;

    TLorentzVector TLVjet1, TLVjet2; // useful for some computation (invariant mass of two jets computed by vbfJetsInvMass() method)

    TH1D *HjetMetDphiMinAllJets = NULL;
    TH1D *HvbfTaggedJets_mT = NULL;
    TH1D *HvbfTaggedJets_deltaEta = NULL;
    TH1D *HvbfTaggedJets_invMass = NULL;

    Double_t DELTAETAMIN_VBFJETS;
    Double_t INVMASS_VBFJETS;
    Double_t JET1PT_VBFJETS;
    Double_t JET2PT_VBFJETS;



  };

  class vbfHiggsToInv_SignalRegion : public vbfHiggsToInvAna {
  public:

    vbfHiggsToInv_SignalRegion(TTree *tree); 
    virtual ~vbfHiggsToInv_SignalRegion() { std::cout<<"~vbfHiggsToInv_SignalRegion() called"<<std::endl; }
    
    void setSelections(); 
    void setMask();
    void setHistograms();
    void setHistogramLastBinAsOverFlow(const Int_t);
    void setNumberParameterValue(const std::string, const Double_t);
    void setVarFromConfigFile();
    Double_t computeEventWeight();
    //void loop(std::vector< std::vector<Double_t>* > &, std::vector< std::vector<Double_t>* > &, std::vector< std::vector<Double_t>* > &);
    void loop(std::vector< Double_t > &, std::vector< Double_t > &, std::vector< Double_t > &, std::vector< Double_t > &, std::vector< Double_t > &, std::vector< Double_t > &, std::vector< Double_t > &, std::vector< Double_t > &, std::vector< Double_t > &);
    void fillEventMask(ULong64_t &); // method to set eventMask event-by-event depending on some selections

  };

  //====================================================

  class vbfHiggsToInv_LeptonControlRegion : public vbfHiggsToInvAna {
  public:

    vbfHiggsToInv_LeptonControlRegion(TTree *tree);  
    virtual ~vbfHiggsToInv_LeptonControlRegion() { std::cout<<"~vbfHiggsToInv_LeptonControlRegion() called"<<std::endl; }
   
    //virtual void setMask();
    virtual void setHistograms();
    virtual void setScaleFactorHistograms();
    virtual void setHistogramLastBinAsOverFlow(const Int_t);
    virtual void createSystematicsHistogram();
    virtual void setNumberParameterValue(const std::string, const Double_t);  // set some numeric variables if they are in config file
    virtual void setControlSampleSpecificParameter();
    virtual void setVarFromConfigFile();
    virtual void setSelections();
    virtual void fillEventMask(ULong64_t &); // method to set eventMask event-by-event depending on some selections

    selection lepLooseVetoC;  //if Z->mumu, veto on electrons and viceversa  // could use muon or electron veto selection in AnalysisDarkMatter, but then I should define a selection* pointing to the correct one to be used (for now I define a new selection)
   // selection lep1tightIdIso04C;
   // selection twoLepTightC;
   // selection lep1ptC;
   // selection lep2ptC;
   // selection lep1etaC;  
   //selection lep2etaC;
   selection genLepC;  
   selection recoGenLepMatchC;
   selection genTauC;

   Int_t LEP_PDG_ID;           // choose electrons or muons
   Double_t LEP1PT;
   Double_t LEP1ETA;
   Double_t LEP_ISO_04;
   Int_t GENLEP_TAG; // decide whether to use genLep and recoGenMatch selection
   Double_t GENLEP1PT;
   Double_t GENLEP1ETA;
   Double_t HLT_LEP1PT;
   Double_t HLT_LEP1ETA;

   // control samples specific parameters
   char FLAVOUR[10];                   // e.g. "ele", "mu"
   char LL_FLAVOUR[10];             // e.g. "ee", "mumu"
   char CONTROL_SAMPLE[10];   // e.g. "Z-->ee"

   TH1D *Hlep1ptDistribution = NULL;
   TH1D *Hlep1etaDistribution = NULL;

   TH1D *Hlep1ptDistribution_monoV = NULL;
   TH1D *Hlep1etaDistribution_monoV = NULL;

   // some variables used in the code for the lepton control regions (both Z and W)

   Int_t recoLepFound_flag;
   Int_t genLepFound_flag;

   TVector2 metNoLepTV, ele; // to compute metNoEle

   Int_t HLT_passed_flag; // initialize to 1, because sometimes it is not needed, meaning trigger automatically passed   

   Float_t *ptr_nLepLoose = NULL;    // depending on lepton flavour in Z-->ll, it will point to different branches                                                    
   Float_t *ptr_nLep10V = NULL;
   Float_t *ptr_nLepTight = NULL;    // depending on lepton flavour in Z-->ll, it will point to different branches                                 
   Float_t *ptr_metNoLepPt = NULL;  // only needed for muons, it will point to the branches with the metNoMu_pt, then metNoLepPt = *ptr_metNoLepPt (metNoLepPt defined below)                                                                                                                                                             
   //Float_t *ptr_metNoLepEta = NULL;                                                                                                                                 
   Float_t *ptr_metNoLepPhi = NULL;
   Int_t *ptr_nRecoLepton = NULL;
   Float_t *ptr_lepton_pt = NULL;
   Float_t *ptr_lepton_eta = NULL;
   Float_t *ptr_lepton_phi = NULL;
   Float_t *ptr_lepton_mass = NULL;
   Float_t nLepLoose;               // this variable and the following should be an integer, but in Emanuele's trees they are float, so I keep them as such     
   Float_t nLep10V;
   Float_t nLepTight;               // this variable and the following should be an integer, but in Emanuele's trees they are float, so I keep them as such   
   Double_t metNoLepPt;        // this variable will be assigned with *ptr_metNoLepPt, where the pointer will point to the branch metNoMu_pt for mu, and with a hand-defined variable for e                                                                                                                                            
   //Double_t metNoLepEta;                                                                                                                                        
   Double_t metNoLepPhi;   // same story as above                                                                                                               
   Int_t nRecoLepton;

  };


  //====================================================



  class vbfHiggsToInv_zllControlRegion : public vbfHiggsToInv_LeptonControlRegion {
  public:

    vbfHiggsToInv_zllControlRegion(TTree *tree);  
    virtual ~vbfHiggsToInv_zllControlRegion() { std::cout<<"~vbfHiggsToInv_zllControlRegion() called"<<std::endl; }
  
    void setSelections(); 
    void setMask();
    void setHistograms();
    void setScaleFactorHistograms();
    void setHistogramLastBinAsOverFlow(const Int_t);
    void createSystematicsHistogram();
    void setNumberParameterValue(const std::string, const Double_t);  // set some numeric variables if they are in config file
    void setControlSampleSpecificParameter();
    void setVarFromConfigFile();
    Double_t computeEventWeight();
    //void loop(std::vector< Double_t > &, std::vector< Double_t > &, std::vector< Double_t > &);
    void loop(std::vector< Double_t > &, std::vector< Double_t > &, std::vector< Double_t > &, std::vector< Double_t > &, std::vector< Double_t > &, std::vector< Double_t > &, std::vector< Double_t > &, std::vector< Double_t > &, std::vector< Double_t > &);
    void fillEventMask(ULong64_t &); // method to set eventMask event-by-event depending on some selections

    selection oppChargeLeptonsC;
    selection invMassC;
    selection twoLepLooseC;   //if Z->mumu, select loose muons (if Z->ee select electrons)
    selection tightLepC;

    Double_t LEP2PT;
    Double_t LEP2ETA;
    Double_t DILEPMASS_LOW;
    Double_t DILEPMASS_UP;
    Double_t GENLEP2PT;
    Double_t GENLEP2ETA;
    Double_t GEN_ZMASS_LOW;
    Double_t GEN_ZMASS_UP;
    Double_t HLT_LEP2PT;
    Double_t HLT_LEP2ETA;

    // control samples specific parameters
    Double_t invMassBinWidth;  // invariant mass histogram's bin width in GeV
    Int_t NinvMassBins;            // number of bins (computed given the range and invMassBinWidth
    // the following flag is needed to enable search for Z->ll at generator level. For MC samples different from DYJetsToLL I must not require 2 gen leptons from Z
    // unless it is Z->tautau, in which case I start from generated taus and apply selection (tau can produce muon or electron)
    Int_t using_zlljets_MCsample_flag; 
    Int_t using_ztautaujets_MCsample_flag;

    // vbfHiggsToInv histograms    
    TH1D *HinvMass = NULL;
    TH1D *HzptDistribution = NULL;    
    TH1D *Hlep2ptDistribution = NULL;
    TH1D *Hlep2etaDistribution = NULL;
    //following histograms filled using different scale factor for NLO xsec for Z and W to be used for systematic computation in ratio between MET in signal and control region
    TH1D *HYieldsMetBin_LepTightLooseUp = NULL;
    TH1D *HYieldsMetBin_LepTightLooseDown = NULL;
    // syst. uncertainty
    TH1D *HSyst_LepTightLoose = NULL; 

  };


  class vbfHiggsToInv_wlnuControlRegion : public vbfHiggsToInv_LeptonControlRegion {
  public:

    vbfHiggsToInv_wlnuControlRegion(TTree *tree);  
    virtual ~vbfHiggsToInv_wlnuControlRegion() { std::cout<<"~vbfHiggsToInv_wlnuControlRegion() called"<<std::endl; }
  
    void setSelections(); 
    void setMask();
    void setHistograms();
    void setScaleFactorHistograms();
    void setHistogramsForTests();
    void setHistogramLastBinAsOverFlow(const Int_t);
    void createSystematicsHistogram();
    void setNumberParameterValue(const std::string, const Double_t);  // set some numeric variables if they are in config file
    void setControlSampleSpecificParameter();
    void setVarFromConfigFile();
    Double_t computeEventWeight();
    //void loop(std::vector< Double_t > &, std::vector< Double_t > &, std::vector< Double_t > &);
    void loop(std::vector< Double_t > &, std::vector< Double_t > &, std::vector< Double_t > &, std::vector< Double_t > &, std::vector< Double_t > &, std::vector< Double_t > &, std::vector< Double_t > &, std::vector< Double_t > &, std::vector< Double_t > &);
    void fillEventMask(ULong64_t &); // method to set eventMask event-by-event depending on some selections

    // specific selections
    selection oneLepLooseC;   //if W->munu, select loose muon (if W->enu select electron)
    selection tightLepC;
    selection metC; //real met, used only for W->ev 
    selection transverseMassC;

    Double_t MT_UP;

    // control samples specific parameters
    // nothing for now
      
    // variables used in loop() or other methods. 
    Double_t mT; // transverse mass, computed in loop 

    // the following flag is needed to enable search for W->lnu at generator level. For MC samples different from WJetsToLNu I must not require 2 gen leptons from Z
    // unless it is Z->tautau, in which case I start from generated taus and apply selection (tau can produce muon or electron)
    Int_t using_wlnujets_MCsample_flag; 
    Int_t using_wtaunujets_MCsample_flag;
        
    // vbfHiggsToInv histograms
    TH1D *HtransverseMass = NULL;
    TH1D *HptW_mT0to50 = NULL;
    //following histograms filled using different scale factor for NLO xsec for Z and W to be used for systematic computation in ratio between MET in signal and control region
    TH1D *HYieldsMetBin_LepTightUp = NULL;
    TH1D *HYieldsMetBin_LepTightDown = NULL;
    // syst. uncertainty
    TH1D *HSyst_LepTight = NULL; 

    // histograms for tests
    TH1D *HgenTransverseMass;
    TH1D *HgenDphiLepMet;
    TH1D *HtransverseMass_EB;
    TH1D *HtransverseMass_EE;


  };

  //====================================================


  //====================================================
  //====================================================
  // MONOJET
  //====================================================
  //====================================================


  class monojetAna : public AnalysisDarkMatter {
  public:

    monojetAna(TTree *tree); 
    virtual ~monojetAna() { std::cout<<"~monojetAna() called"<<std::endl; }
    
    virtual void setSelections(); 
    virtual void setHistograms();
    virtual void setScaleFactorHistograms();
    virtual void setHistogramLastBinAsOverFlow(const Int_t);
    virtual void setNumberParameterValue(const std::string, const Double_t);
    virtual void setVarFromConfigFile();
    virtual void createSystematicsHistogram();
    virtual void fillEventMask(ULong64_t &); // method to set eventMask event-by-event depending on some selections
    virtual void loop(std::vector< Double_t > &, std::vector< Double_t > &, std::vector< Double_t > &, std::vector< Double_t > &, std::vector< Double_t > &, std::vector< Double_t > &, std::vector< Double_t > &, std::vector< Double_t > &, std::vector< Double_t > &) {};

    selection jet1C;
    selection VtagC;
    selection noVtagC; // opposite of VtagC: monojet and monoV are now exclusive cathegories
    //adding VTag selections separately
    selection ak8jet1C;
    selection ak8Tau2OverTau1C;
    selection ak8prunedMassC;
    selection harderRecoilC; // met > 250 instead of 200

    //Int_t NJETS;
    Double_t J1PT;
    Double_t J1ETA;
    //Double_t J2PT;
    //Double_t J2ETA;
    //Double_t J1J2DPHI;
    
    Int_t Vtagged_flag; // set to 1 or 0 in the event loop depending on the event passing the selection specific for the mono-V cathegory (not taking the common part into account, that one can or cannot be passed)
    
    //monoV histograms
    TH1D *HYieldsMetBin_monoV = NULL;
    TH1D *HhtDistribution_monoV = NULL; 
    TH1D *HvtxDistribution_monoV = NULL;   
    TH1D *HnjetsDistribution_monoV = NULL;   
    TH1D *Hjet1etaDistribution_monoV = NULL;
    TH1D *HrecoilDistribution_monoV = NULL;
    TH1D *Hjet1ptDistribution_monoV = NULL;
    TH1D *HprunedMassDistribution_monoV = NULL;
    TH1D *Htau2OverTau1Distribution_monoV = NULL;
    TH1D *HmetBinEdges_monoV = NULL;
    //following histograms filled using different scale factor for NLO xsec for Z and W to be used for systematic computation in ratio between MET in signal and control region
    TH1D *HYieldsMetBin_qcdRenScaleUp_monoV = NULL;
    TH1D *HYieldsMetBin_qcdRenScaleDown_monoV = NULL;
    TH1D *HYieldsMetBin_qcdFacScaleUp_monoV = NULL;
    TH1D *HYieldsMetBin_qcdFacScaleDown_monoV = NULL;
    TH1D *HYieldsMetBin_qcdPdfUp_monoV = NULL;
    TH1D *HYieldsMetBin_qcdPdfDown_monoV = NULL;
    TH1D *HYieldsMetBin_ewkUp_monoV = NULL;
    TH1D *HYieldsMetBin_ewkDown_monoV = NULL;
    //systematic uncertainties
    TH1D *HSyst_qcdRenScale_monoV = NULL;
    TH1D *HSyst_qcdFacScale_monoV = NULL;
    TH1D *HSyst_qcdPdf_monoV = NULL;
    TH1D *HSyst_ewk_monoV = NULL;
    TH1D *HSyst_total_monoV = NULL;


  };


  //====================================================


  class zlljets_Axe_noSkim_light : public AnalysisDarkMatter /*,public edimarcoTreeFriend*/ {
  public:

    // only computes acceptance and efficiency in order to be faster
    zlljets_Axe_noSkim_light(TTree *tree);
    virtual ~zlljets_Axe_noSkim_light() { std::cout<<"~zlljets_Axe_noSkim_light() called"<<std::endl; }
  
    void loop(const char* configFileName);

  };

  class zlljets_resoResp : public AnalysisDarkMatter /*,public edimarcoTreeFriend*/ {
  public:

    zlljets_resoResp(TTree *tree);
    virtual ~zlljets_resoResp() { std::cout<<"~zlljets_resoResp() called"<<std::endl; }
  
    void loop(const char* configFileName, const Int_t ISDATA_FLAG, const Int_t unweighted_event_flag);

  };


  //====================================================

  class monojet_SignalRegion : public monojetAna {
  public:

    monojet_SignalRegion(TTree *tree); 
    virtual ~monojet_SignalRegion() { std::cout<<"~monojet_SignalRegion() called"<<std::endl; }
  
    void setSelections(); 
    void setMask();
    void setHistograms();
    void setHistogramLastBinAsOverFlow(const Int_t);
    void setNumberParameterValue(const std::string, const Double_t);
    void setVarFromConfigFile();
    Double_t computeEventWeight();
    //void loop(std::vector< std::vector<Double_t>* > &, std::vector< std::vector<Double_t>* > &, std::vector< std::vector<Double_t>* > &);
    void loop(std::vector< Double_t > &, std::vector< Double_t > &, std::vector< Double_t > &, std::vector< Double_t > &, std::vector< Double_t > &, std::vector< Double_t > &, std::vector< Double_t > &, std::vector< Double_t > &, std::vector< Double_t > &);
    void fillEventMask(ULong64_t &); // method to set eventMask event-by-event depending on some selections


  };

  //====================================================

  class monojet_PhotonControlRegion : public monojetAna {
  public:

    monojet_PhotonControlRegion(TTree *tree);  
    virtual ~monojet_PhotonControlRegion() { std::cout<<"~monojet_PhotonControlRegion() called"<<std::endl; }
    
    void setSelections(); 
    void setMask();
    void setHistograms();
    void setHistogramLastBinAsOverFlow(const Int_t);
    void setNumberParameterValue(const std::string, const Double_t);  // set some numeric variables if they are in config file
    void setControlSampleSpecificParameter();
    void setVarFromConfigFile();
    Double_t computeEventWeight();
    //void loop(std::vector< Double_t > &, std::vector< Double_t > &, std::vector< Double_t > &);
    void loop(std::vector< Double_t > &, std::vector< Double_t > &, std::vector< Double_t > &, std::vector< Double_t > &, std::vector< Double_t > &, std::vector< Double_t > &, std::vector< Double_t > &, std::vector< Double_t > &, std::vector< Double_t > &);
    void fillEventMask(ULong64_t &); // method to set eventMask event-by-event depending on some selections


    selection oneGammaLooseC;
    selection tightPhotonC;
    
    Double_t PH1PT;
    Double_t PH1ETA;
    
    Int_t using_gammajets_MCsample_flag;
    
    // control samples specific parameters
    // char FLAVOUR[10];          // e.g. "photon"
    char CONTROL_SAMPLE[10];   // e.g. "gamma"
    
    TH1D *Hphoton1ptDistribution = NULL;
    TH1D *Hphoton1etaDistribution = NULL;
    
    TH1D *Hphoton1ptDistribution_monoV = NULL;
    TH1D *Hphoton1etaDistribution_monoV = NULL;
    
    
  };
  

  //====================================================

  class monojet_LeptonControlRegion : public monojetAna {
  public:

    monojet_LeptonControlRegion(TTree *tree);  
    virtual ~monojet_LeptonControlRegion() { std::cout<<"~monojet_LeptonControlRegion() called"<<std::endl; }
   
    //virtual void setMask();
    virtual void setHistograms();
    virtual void setScaleFactorHistograms();
    virtual void setHistogramLastBinAsOverFlow(const Int_t);
    virtual void createSystematicsHistogram();
    virtual void setNumberParameterValue(const std::string, const Double_t);  // set some numeric variables if they are in config file
    virtual void setControlSampleSpecificParameter();
    virtual void setVarFromConfigFile();
    virtual void setSelections();
    virtual void fillEventMask(ULong64_t &); // method to set eventMask event-by-event depending on some selections

    selection lepLooseVetoC;  //if Z->mumu, veto on electrons and viceversa  // could use muon or electron veto selection in AnalysisDarkMatter, but then I should define a selection* pointing to the correct one to be used (for now I define a new selection)
   // selection lep1tightIdIso04C;
   // selection twoLepTightC;
   // selection lep1ptC;
   // selection lep2ptC;
   // selection lep1etaC;  
   //selection lep2etaC;
   selection genLepC;  
   selection recoGenLepMatchC;
   selection genTauC;

   Int_t LEP_PDG_ID;           // choose electrons or muons
   Double_t LEP1PT;
   Double_t LEP1ETA;
   Double_t LEP_ISO_04;
   Int_t GENLEP_TAG; // decide whether to use genLep and recoGenMatch selection
   Double_t GENLEP1PT;
   Double_t GENLEP1ETA;
   Double_t HLT_LEP1PT;
   Double_t HLT_LEP1ETA;

   // control samples specific parameters
   char FLAVOUR[10];                   // e.g. "ele", "mu"
   char LL_FLAVOUR[10];             // e.g. "ee", "mumu"
   char CONTROL_SAMPLE[10];   // e.g. "Z-->ee"

   TH1D *Hlep1ptDistribution = NULL;
   TH1D *Hlep1etaDistribution = NULL;

   TH1D *Hlep1ptDistribution_monoV = NULL;
   TH1D *Hlep1etaDistribution_monoV = NULL;

   // some variables used in the code for the lepton control regions (both Z and W)

   Int_t recoLepFound_flag;
   Int_t genLepFound_flag;

   TVector2 metNoLepTV, ele; // to compute metNoEle

   Int_t HLT_passed_flag; // initialize to 1, because sometimes it is not needed, meaning trigger automatically passed   

   Float_t *ptr_nLepLoose = NULL;    // depending on lepton flavour in Z-->ll, it will point to different branches                                                    
   Float_t *ptr_nLep10V = NULL;
   Float_t *ptr_nLepTight = NULL;    // depending on lepton flavour in Z-->ll, it will point to different branches                                 
   Float_t *ptr_metNoLepPt = NULL;  // only needed for muons, it will point to the branches with the metNoMu_pt, then metNoLepPt = *ptr_metNoLepPt (metNoLepPt defined below)                                                                                                                                                             
   //Float_t *ptr_metNoLepEta = NULL;                                                                                                                                 
   Float_t *ptr_metNoLepPhi = NULL;
   Int_t *ptr_nRecoLepton = NULL;
   Float_t *ptr_lepton_pt = NULL;
   Float_t *ptr_lepton_eta = NULL;
   Float_t *ptr_lepton_phi = NULL;
   Float_t *ptr_lepton_mass = NULL;
   Float_t nLepLoose;               // this variable and the following should be an integer, but in Emanuele's trees they are float, so I keep them as such     
   Float_t nLep10V;
   Float_t nLepTight;               // this variable and the following should be an integer, but in Emanuele's trees they are float, so I keep them as such   
   Double_t metNoLepPt;        // this variable will be assigned with *ptr_metNoLepPt, where the pointer will point to the branch metNoMu_pt for mu, and with a hand-defined variable for e                                                                                                                                            
   //Double_t metNoLepEta;                                                                                                                                        
   Double_t metNoLepPhi;   // same story as above                                                                                                               
   Int_t nRecoLepton;


  };


  //====================================================



  class zlljetsControlSample : public monojet_LeptonControlRegion {
  public:

    zlljetsControlSample(TTree *tree);  
    virtual ~zlljetsControlSample() { std::cout<<"~zlljetsControlSample() called"<<std::endl; }
  
    void setSelections(); 
    void setMask();
    void setHistograms();
    void setScaleFactorHistograms();
    void setHistogramLastBinAsOverFlow(const Int_t);
    void createSystematicsHistogram();
    void setNumberParameterValue(const std::string, const Double_t);  // set some numeric variables if they are in config file
    void setControlSampleSpecificParameter();
    void setVarFromConfigFile();
    Double_t computeEventWeight();
    //void loop(std::vector< Double_t > &, std::vector< Double_t > &, std::vector< Double_t > &);
    void loop(std::vector< Double_t > &, std::vector< Double_t > &, std::vector< Double_t > &, std::vector< Double_t > &, std::vector< Double_t > &, std::vector< Double_t > &, std::vector< Double_t > &, std::vector< Double_t > &, std::vector< Double_t > &);
    void fillEventMask(ULong64_t &); // method to set eventMask event-by-event depending on some selections

    selection oppChargeLeptonsC;
    selection invMassC;
    selection twoLepLooseC;   //if Z->mumu, select loose muons (if Z->ee select electrons)
    selection tightLepC;

    Double_t LEP2PT;
    Double_t LEP2ETA;
    Double_t DILEPMASS_LOW;
    Double_t DILEPMASS_UP;
    Double_t GENLEP2PT;
    Double_t GENLEP2ETA;
    Double_t GEN_ZMASS_LOW;
    Double_t GEN_ZMASS_UP;
    Double_t HLT_LEP2PT;
    Double_t HLT_LEP2ETA;

    // control samples specific parameters
    Double_t invMassBinWidth;  // invariant mass histogram's bin width in GeV
    Int_t NinvMassBins;            // number of bins (computed given the range and invMassBinWidth
    // the following flag is needed to enable search for Z->ll at generator level. For MC samples different from DYJetsToLL I must not require 2 gen leptons from Z
    // unless it is Z->tautau, in which case I start from generated taus and apply selection (tau can produce muon or electron)
    Int_t using_zlljets_MCsample_flag; 
    Int_t using_ztautaujets_MCsample_flag;

    // monojet histograms    
    TH1D *HinvMass = NULL;
    TH1D *HzptDistribution = NULL;    
    TH1D *Hlep2ptDistribution = NULL;
    TH1D *Hlep2etaDistribution = NULL;
    //following histograms filled using different scale factor for NLO xsec for Z and W to be used for systematic computation in ratio between MET in signal and control region
    TH1D *HYieldsMetBin_LepTightLooseUp = NULL;
    TH1D *HYieldsMetBin_LepTightLooseDown = NULL;
    // syst. uncertainty
    TH1D *HSyst_LepTightLoose = NULL; 

    // monoV histograms
    TH1D *HinvMass_monoV = NULL;
    TH1D *HzptDistribution_monoV = NULL;    
    TH1D *Hlep2ptDistribution_monoV = NULL;
    TH1D *Hlep2etaDistribution_monoV = NULL;
    //following histograms filled using different scale factor for NLO xsec for Z and W to be used for systematic computation in ratio between MET in signal and control region
    TH1D *HYieldsMetBin_LepTightLooseUp_monoV = NULL;
    TH1D *HYieldsMetBin_LepTightLooseDown_monoV = NULL;
    // syst. uncertainty
    TH1D *HSyst_LepTightLoose_monoV = NULL; 

  };

  //====================================================

  class wlnujetsControlSample : public monojet_LeptonControlRegion {
  public:

    wlnujetsControlSample(TTree *tree);  
    virtual ~wlnujetsControlSample() { std::cout<<"~wlnujetsControlSample() called"<<std::endl; }
  
    void setSelections(); 
    void setMask();
    void setHistograms();
    void setScaleFactorHistograms();
    void setHistogramsForTests();
    void setHistogramLastBinAsOverFlow(const Int_t);
    void createSystematicsHistogram();
    void setNumberParameterValue(const std::string, const Double_t);  // set some numeric variables if they are in config file
    void setControlSampleSpecificParameter();
    void setVarFromConfigFile();
    Double_t computeEventWeight();
    //void loop(std::vector< Double_t > &, std::vector< Double_t > &, std::vector< Double_t > &);
    void loop(std::vector< Double_t > &, std::vector< Double_t > &, std::vector< Double_t > &, std::vector< Double_t > &, std::vector< Double_t > &, std::vector< Double_t > &, std::vector< Double_t > &, std::vector< Double_t > &, std::vector< Double_t > &);
    void fillEventMask(ULong64_t &); // method to set eventMask event-by-event depending on some selections

    // specific selections
    selection oneLepLooseC;   //if W->munu, select loose muon (if W->enu select electron)
    selection tightLepC;
    selection metC; //real met, used only for W->ev 
    selection transverseMassC;

    Double_t MT_UP;

    // control samples specific parameters
    // nothing for now
      
    // variables used in loop() or other methods. 
    Double_t mT; // transverse mass, computed in loop 

    // the following flag is needed to enable search for W->lnu at generator level. For MC samples different from WJetsToLNu I must not require 2 gen leptons from Z
    // unless it is Z->tautau, in which case I start from generated taus and apply selection (tau can produce muon or electron)
    Int_t using_wlnujets_MCsample_flag; 
    Int_t using_wtaunujets_MCsample_flag;
        
    // monoJet histograms
    TH1D *HtransverseMass = NULL;
    TH1D *HptW_mT0to50 = NULL;
    //following histograms filled using different scale factor for NLO xsec for Z and W to be used for systematic computation in ratio between MET in signal and control region
    TH1D *HYieldsMetBin_LepTightUp = NULL;
    TH1D *HYieldsMetBin_LepTightDown = NULL;
    // syst. uncertainty
    TH1D *HSyst_LepTight = NULL; 

    // monoV histograms
    TH1D *HtransverseMass_monoV = NULL;
    //following histograms filled using different scale factor for NLO xsec for Z and W to be used for systematic computation in ratio between MET in signal and control region
    TH1D *HYieldsMetBin_LepTightUp_monoV = NULL;
    TH1D *HYieldsMetBin_LepTightDown_monoV = NULL;
    // syst. uncertainty
    TH1D *HSyst_LepTight_monoV = NULL; 

    // histograms for tests
    TH1D *HgenTransverseMass;
    TH1D *HgenDphiLepMet;
    TH1D *HtransverseMass_EB;
    TH1D *HtransverseMass_EE;


  };

  //====================================================

  class zlljets_metResoResp : public AnalysisDarkMatter {
  public:

    zlljets_metResoResp(TTree *tree); 
    virtual ~zlljets_metResoResp() { std::cout<<"~zlljets_metResoResp() called"<<std::endl; }
  
    void loop(std::vector< Double_t > &, std::vector< Double_t > &, std::vector< Double_t > &);
  
  };

  //====================================================

}
#endif




