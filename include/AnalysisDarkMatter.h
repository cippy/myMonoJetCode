#ifndef ANALYSIS_DARK_MATTER_h
#define ANALYSIS_DARK_MATTER_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>

#include <iostream>
#include <cstdlib>
#include <vector>
#include <string>

#include "edimarcoTree_v4.h"
#include "functionsForAnalysis.h"
#include "myClasses.h"

class AnalysisDarkMatter : public edimarcoTree_v4 {
 public:

  AnalysisDarkMatter(TTree *tree);

  virtual ~AnalysisDarkMatter() { std::cout<<"~AnalysisDarkMatter() called"<<std::endl; }
  virtual Int_t    GetEntry(Long64_t entry);
  virtual Long64_t LoadTree(Long64_t entry);
  virtual void Init(TTree *tree);

  //common selections (between signal and control regions) are declared here
  selection HLTC;
  selection metFiltersC;
  selection recoilC; //selection metNoLepC;
  selection jet1C;
  selection jetMetDphiMinC;
  selection jetNoiseCleaningC;
  selection bjetVetoC;
  selection muonLooseVetoC; 
  selection electronLooseVetoC; 
  selection tauLooseVetoC;
  selection gammaLooseVetoC;
  selection VtagC;
  selection noVtagC; // opposite of VtagC: monojet and monoV are now exclusive cathegories

  mask analysisMask;
  selectionManager analysisSelectionManager;
  mask analysisMask_monoJ;
  selectionManager analysisSelectionManager_monoJ;
  mask analysisMask_monoV;
  selectionManager analysisSelectionManager_monoV;

  std::vector<mask*> anaMasksPtrCollection; // list of pointers to masks in the class 

  virtual void setBasicConf(const char* inputSuffix, const std::string inputUncertainty, const char* inputConfigFileName, const Int_t inputIsDataFlag, const Int_t inputUnweightedEeventFlag, const Int_t inputHasSFfriendFlag);
  virtual void setCalibEleFlag();
  virtual void setDirNameSuffix(const std::string);
  virtual void setNumberParameterValue(const std::string, const Double_t);
  virtual void setVarFromConfigFile();
  virtual void setSelections();
  virtual void setHistograms();
  virtual void setScaleFactorHistograms();
  virtual void setHistogramLastBinAsOverFlow(const Int_t);
  virtual void createSystematicsHistogram(); //build histograms with systematic uncertainties
  virtual void fillRowVector(const Double_t, const selectionManager &, const mask &, std::vector<Double_t> &, std::vector<Double_t> &, std::vector<Double_t> &, const Int_t);
  /* virtual void set_SF_NLO_name(const std::string); */
  /* virtual void set_SF_NLO_pointers(const std::string sf_option, Float_t *ptrQCD, Float_t *ptrEWK); */
  //virtual Double_t computeEventWeight() const;   // return weight for the event

  char ROOT_FNAME[100];
  char TXT_FNAME[100];
  char TEX_FNAME[100];

  Double_t LUMI;
  //Int_t NJETS;
  Double_t J1PT;
  Double_t J1ETA;
  //Double_t J2PT;
  //Double_t J2ETA;
  //Double_t J1J2DPHI;
  Int_t TAU_VETO_FLAG;
  Int_t HLT_FLAG;                  // usage depends on specific analysis
  Double_t METNOLEP_START;
  Int_t MET_FILTERS_FLAG;
  Double_t JMET_DPHI_MIN;
  std::string FILENAME_BASE;
  std::string DIRECTORY_TO_SAVE_FILES;
  std::string DIRECTORY_NAME;

  std::string outputFolder;
  std::string dirName_suffix; //used from main to append suffix to directory name without having to change config file

  std::vector<Double_t> metBinEdgesVector;  // filled with values in file named configFileName
  Int_t nMetBins;
  std::vector<Double_t> metBinEdgesVector_monoV;  // filled with values in file named configFileName
  Int_t nMetBins_monoV;

  //obsolete, substituted by selectionManager class
  std::vector<Int_t> selStep;  //set when setting the mask analysisMask 
  //array to store index of step to form selection flow (might want to consider two or more steps together and not separated)
   // in case a step would be present for some sample but not for others (e.g. the RecoGen match done only in Zll MC), the step is referred to as -1 and the corresponding values are set to -1, so that, when printing the table, yields will be filled with " / / " which means " uneffected" (because that step was not done)

  // the following can be set by setBasicConf(), but are initialized in construcor
  std::string suffix;   // it is the sample name (e.g. QCD, ZJetsToNuNu ecc...)
  std::string uncertainty; //uncertainty on yields (from MC, Poisson, X%)
  char* configFileName;
  Int_t ISDATA_FLAG;
  Int_t unweighted_event_flag; //If not 0, a _weq1 suffix is added to directory name
  Int_t hasSFfriend_flag;  //tells if sfFriend are present

  Int_t calibEle_flag; // tells if using calibrated electron properties or traditional ones. If not 0, a _CalibEle suffix is added to directory name
  //std::string sf_nlo; // save option passed to main to decide which scale factor to use for NLO cross section
  Int_t hasScaledHistograms_flag; // sometimes I need to check if the sample I'm running on is a specific one. In order to avoid to check for it every time I need, I compute it once and for all storing results as a flag in tis variable
  /* Double_t sf_nlo_weight; */
  /* Float_t *ptr_sf_nlo_QCD = NULL; */
  /* Float_t *ptr_sf_nlo_EWK = NULL; */

  // Maybe it would be better to left following two variables in analyzers
  Double_t nTotalWeightedEvents;  // counter of total events (with weights if any). Initialized to 0 in constructor
  Double_t newwgt;                           // weight for the event (specific definition depends on sample and on whether data are being analyzed)

  Int_t Vtagged_flag; // set to 1 or 0 in the event loop depending on the event passing the selection specific for the mono-V cathegory (not taking the common part into account, that one can or cannot be passed)

  //histograms for monojet (exclusive, but I don't rename them with *monoJ)

  //root histograms: these are common among all analysis (signal ad control region)
  TH1D *HYieldsMetBin = NULL;
  TH1D *HhtDistribution = NULL; 
  TH1D *HvtxDistribution = NULL;   
  TH1D *HnjetsDistribution = NULL;   
  TH1D *Hj1j2dphiDistribution = NULL;
  TH1D *HjetMetDphiMinDistribution = NULL;
  TH1D *Hjet1etaDistribution = NULL;
  TH1D *Hjet2etaDistribution = NULL;
  TH1D *HmetNoLepDistribution = NULL;
  TH1D *Hjet1ptDistribution = NULL;
  TH1D *Hjet2ptDistribution = NULL;
  TH1D *HmetBinEdges = NULL;
  //following histograms filled using different scale factor for NLO xsec for Z and W to be used for systematic computation in ratio between MET in signal and control region
  TH1D *HYieldsMetBin_qcdRenScaleUp = NULL;
  TH1D *HYieldsMetBin_qcdRenScaleDown = NULL;
  TH1D *HYieldsMetBin_qcdFacScaleUp = NULL;
  TH1D *HYieldsMetBin_qcdFacScaleDown = NULL;
  TH1D *HYieldsMetBin_qcdPdfUp = NULL;
  TH1D *HYieldsMetBin_qcdPdfDown = NULL;
  TH1D *HYieldsMetBin_ewkUp = NULL;
  TH1D *HYieldsMetBin_ewkDown = NULL;
  //systematic uncertainties
  TH1D *HSyst_qcdRenScale = NULL;
  TH1D *HSyst_qcdFacScale = NULL;
  TH1D *HSyst_qcdPdf = NULL;
  TH1D *HSyst_ewk = NULL;
  TH1D *HSyst_total = NULL;

  //monoV histograms
  TH1D *HYieldsMetBin_monoV = NULL;
  TH1D *HhtDistribution_monoV = NULL; 
  TH1D *HvtxDistribution_monoV = NULL;   
  TH1D *HnjetsDistribution_monoV = NULL;   
  TH1D *Hjet1etaDistribution_monoV = NULL;
  TH1D *HmetNoLepDistribution_monoV = NULL;
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

#endif
