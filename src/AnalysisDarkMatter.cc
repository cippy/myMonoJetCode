#define AnalysisDarkMatter_cxx
#include "edimarcoTree_v7.h"
#include "AnalysisDarkMatter.h"
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

#ifdef AnalysisDarkMatter_cxx

//===============================================

AnalysisDarkMatter::AnalysisDarkMatter(TTree *tree) : edimarcoTree_v7(tree) {
  //cout <<"check in constructor "<<endl;
  suffix = "";
  uncertainty = "";
  configFileName = NULL;
  ISDATA_FLAG = 0;
  unweighted_event_flag = 0;
  hasSFfriend_flag = 0;

  calibEle_flag = 0;
  hasScaledHistograms_flag = 0;
  dirName_suffix = ""; //user can add a suffix to directory name from main instead of changing string in config file
  nTotalWeightedEvents = 0.0;
  newwgt = 0.0;

  // initialize some variables with sensible values. They will be set later depending on config file
  LUMI = 1.0;
  B_VETO_FLAG = 1;
  TAU_VETO_FLAG = 1;
  HLT_FLAG = 0;                  // usage depends on specific analysis
  JMET_DPHI_MIN = 0.5;
  MET_FILTERS_FLAG = 1;
  METNOLEP_START = 0.0;
  ENABLE_HISTOGRAM_FOR_TESTS_FLAG = 0;  

  Init(tree);
  

}

//===============================================

Int_t AnalysisDarkMatter::GetEntry(Long64_t entry) {
  return edimarcoTree_v7::GetEntry(entry);
}

//===============================================

Long64_t AnalysisDarkMatter::LoadTree(Long64_t entry) {
  return edimarcoTree_v7::LoadTree(entry);
}

//===============================================

void AnalysisDarkMatter::Init(TTree *tree) {
  edimarcoTree_v7::Init(tree);
} 

//===============================================

void AnalysisDarkMatter::setBasicConf(const char* inputSuffix, const string inputUncertainty, const char* inputConfigFileName, const Int_t inputIsDataFlag, const Int_t inputUnweightedEeventFlag, const Int_t inputHasSFfriendFlag) {

  suffix = inputSuffix;  // it is the sample name (e.g. QCD, ZJetsToNuNu ecc...)
  uncertainty = inputUncertainty; //sample uncertainty (poisson, MC, X%), for now it is just used when creating the table to assign uncertainties to the sum of MC yields, which depends on the uncertainty on the yields for each MC
  configFileName = (char*) inputConfigFileName;
  ISDATA_FLAG = inputIsDataFlag;
  unweighted_event_flag = inputUnweightedEeventFlag;
  hasSFfriend_flag = inputHasSFfriendFlag;

}

//===============================================

void AnalysisDarkMatter::setCalibEleFlag() {

  calibEle_flag = 1;

}

//===============================================

void AnalysisDarkMatter::setDirNameSuffix(const std::string s) {

  dirName_suffix = s;

}

//===============================================

void AnalysisDarkMatter::setNumberParameterValue(const string parameterName, const Double_t value) {

  // variables in file are read as double, if I need int the compiler will truncate. An integer assigned to a double might actually be a little less than the integer in binary notation due to limited number of digits (floats can be even more dangerous). Thus, if a variable is expected to be an int, I sum 0.5 to value before assigning it.
  // e.g. Suppose you want a flag = 1 (or 1.0). Casting flag to Int might yield 0 if the binary expression for 1.0 is less than 1 (ok that's just an example). On the other hand, I can say Int_t val = (Int_t) 1.0 + 0.5 which will yield val = 1. Note however that casting is not necessary, assignment to int will do it automatically.

  //  -->  N.B.  <---
  // adding 0.5 works for positive number, if using negative ones, I must subtract 0.5 (or sum -0.5). For 0 it is irrelevant, e.g. both 0.3 and -0.3 will become 0 if assigned to int, so I don't need to care if 0.0 is actually a positive or integer number
  // flags are expected to be 0 or 1 but let's do the check to be sure

  if (parameterName == "LUMI") LUMI = value;
  //else if (parameterName == "NJETS") NJETS = (value < 0) ? (-0.5 + value) : (0.5 + value);
  else if (parameterName == "B_VETO_FLAG") B_VETO_FLAG = (value < 0) ? (-0.5 + value) : (0.5 + value);
  else if (parameterName == "TAU_VETO_FLAG") TAU_VETO_FLAG = (value < 0) ? (-0.5 + value) : (0.5 + value);
  else if (parameterName == "HLT_FLAG") HLT_FLAG = (value < 0) ? (-0.5 + value) : (0.5 + value);
  else if (parameterName == "METNOLEP_START") METNOLEP_START = value;
  else if (parameterName == "MET_FILTERS_FLAG") MET_FILTERS_FLAG = (value < 0) ? (-0.5 + value) : (0.5 + value);
  else if (parameterName == "JMET_DPHI_MIN") JMET_DPHI_MIN = value;
  else if (parameterName == "ENABLE_HISTOGRAM_FOR_TESTS_FLAG") ENABLE_HISTOGRAM_FOR_TESTS_FLAG = (value < 0) ? (-0.5 + value) : (0.5 + value); // here our convention is broken, since the flag is written with lower case. This is because at first I wanted to set in as an option passed to main as I do for unweighted_event_flag.

}

//===============================================

void AnalysisDarkMatter::setVarFromConfigFile() {
  
  if (configFileName == NULL) {
    cout << "Error: configFileName points to 'NULL'. End of programme" << endl;
    exit(EXIT_FAILURE);
  }

  ifstream inputFile(configFileName);

  if (inputFile.is_open()) {

    Double_t value;
    string name;
    string parameterName;
    string parameterType;

    //mySpaces(cout,2);
    cout <<endl; cout <<endl;
    cout << "Printing content of " << configFileName << " file" << endl;
    //mySpaces(cout,1);
    cout <<endl;

    while (inputFile >> parameterType ) {

      if (parameterType == "NUMBER") {

	inputFile >> parameterName >> value;
	cout << right << setw(20) << parameterName << "  " << left << value << endl;

	setNumberParameterValue(parameterName, value);

      } else if (parameterType == "STRING") {
	 
	inputFile >> parameterName >> name;
	cout << right << setw(20) << parameterName << "  " << left << name << endl;

	if (parameterName == "FILENAME_BASE") {

	  FILENAME_BASE = name; 

	}

	if (parameterName == "DIRECTORY_PATH") {  // name of directory where files are saved

	  DIRECTORY_TO_SAVE_FILES = name;
	  //std::cout << "Files will be saved in '" << name << "' ." <<std::endl;

	} 

	if (parameterName == "DIRECTORY_NAME") {  // name of directory where files are saved

	  DIRECTORY_NAME = name;
	  //std::cout << "Files will be saved in directory named '" << name << "' ." <<std::endl;

	} 

      } else if (parameterType == "ARRAY_NUM") {

	inputFile >> parameterName;

	if (parameterName == "MET_BIN_EDGES") { 

	  cout << right << setw(20) << parameterName << "  ";
	  string stringvalues;
	  getline(inputFile, stringvalues);    // read whole line starting from current position (i.e. without reading ARRAY_NUM)
	  istringstream iss(stringvalues);
	  Double_t num;

	  while(iss >> num) {
	    
	    metBinEdgesVector.push_back(num);
	    cout << metBinEdgesVector.back() << " ";

	  }

	  cout << endl;

	} else if (parameterName == "MET_BIN_EDGES_MONOV") { 

	  cout << right << setw(20) << parameterName << "  ";
	  string stringvalues;
	  getline(inputFile, stringvalues);    // read whole line starting from current position (i.e. without reading ARRAY_NUM)
	  istringstream iss(stringvalues);
	  Double_t num;

	  while(iss >> num) {
	    
	    metBinEdgesVector_monoV.push_back(num);
	    cout << metBinEdgesVector_monoV.back() << " ";

	  }

	  cout << endl;

	}

      }

    }
     
    //mySpaces(cout,2);
    cout <<endl; cout <<endl;

    inputFile.close();
                                                                                                                         
  } else {

    cout << "Error: could not open file " << configFileName << endl;
    exit(EXIT_FAILURE);
    
  }
  
  string working_cmssw_path = "";
  mySetCMSSW_BasePath(working_cmssw_path);
 
  outputFolder = working_cmssw_path + DIRECTORY_TO_SAVE_FILES + DIRECTORY_NAME;
  if (calibEle_flag == 1) outputFolder += "_CalibEle";  //
  if (unweighted_event_flag == 1) outputFolder += "_weq1";
  if (dirName_suffix != "") outputFolder += ("_" + dirName_suffix);  // set to "" if not specified
  outputFolder += "/";

  // in case the config doesn't have the binning. Btw, having the binning in the config allows me to choose it for every analysis without modifying and compiling again the code
  if (metBinEdgesVector.size() == 0) {
    Double_t tmpBinEdges[] = {200.,230.,260.,290.,320.,350.,390.,430.,470.,510.,550.,590.,640.,690.,740.,790.,840.,900.,960.,1020.,1090.,1160.,1250.}; 
    metBinEdgesVector.assign(tmpBinEdges,tmpBinEdges+23);
  }
  if (metBinEdgesVector_monoV.size() == 0) {
    Double_t tmpBinEdges[] = {250.,300.,350.,400.,500.,600.,1000.}; 
    metBinEdgesVector_monoV.assign(tmpBinEdges,tmpBinEdges+7);
  }

  nMetBins = (metBinEdgesVector.size()) - 1;
  nMetBins_monoV = (metBinEdgesVector_monoV.size()) - 1;
  
  cout << "MetBinEdges: [ ";
  for(Int_t i = 0; i <= nMetBins; i++) {
    if (i != nMetBins) cout << metBinEdgesVector[i] << ", ";
    else cout << metBinEdgesVector[i] << "]" << endl;
  }
  cout << "MetBinEdges_monoV: [ ";
  for(Int_t i = 0; i <= nMetBins_monoV; i++) {
    if (i != nMetBins_monoV) cout << metBinEdgesVector_monoV[i] << ", ";
    else cout << metBinEdgesVector_monoV[i] << "]" << endl;
  }
  cout << endl;
  
  if ( !ISDATA_FLAG && unweighted_event_flag) cout << "Warning: no weight applied to events (w = 1)" << endl;  // if MC with unit weight, make user know
  
  if (ISDATA_FLAG) {
    strcpy(ROOT_FNAME,(FILENAME_BASE + "_DATA.root").c_str());
    strcpy(TXT_FNAME,(FILENAME_BASE + "_DATA.txt").c_str());
    strcpy(TEX_FNAME,(FILENAME_BASE + "_DATA.tex").c_str());
  } else {
    strcpy(ROOT_FNAME,(FILENAME_BASE + "_" + suffix + ".root").c_str());
    strcpy(TXT_FNAME,(FILENAME_BASE + "_" + suffix + ".txt").c_str());
    strcpy(TEX_FNAME,(FILENAME_BASE + "_" + suffix + ".tex").c_str());
  }


}

//===============================================

void AnalysisDarkMatter::setSelections() {


  //cout << "CHECK ANALYSISDARKMATTER"<< endl;
  if (B_VETO_FLAG != 0) bjetVetoC.set("bjet veto","b-jets veto");
  //cout << "CHECK ANALYSISDARKMATTER 2"<< endl;
  if (HLT_FLAG != 0) HLTC.set("trigger","trigger");
  //cout << "CHECK ANALYSISDARKMATTER 3"<< endl;
  if (TAU_VETO_FLAG) tauLooseVetoC.set("tau veto","tau veto");
  //cout << "CHECK ANALYSISDARKMATTER 4"<< endl;

  selection::checkMaskLength();

}

//===============================================

void AnalysisDarkMatter::setHistograms() {

  //  cout << "check" << endl;
  // histograms for monojet (exclusive, but I don't rename them as *_monoJ)

  HYieldsMetBin = new TH1D("HYieldsMetBin","yields in bins of met; #slash{E}_{T};# of events",nMetBins,metBinEdgesVector.data());
  HhtDistribution = new TH1D("HhtDistribution","",400,0.0,4000.0);
  HvtxDistribution = new TH1D("HvtxDistribution","",40,-0.5,39.5);   
  HnjetsDistribution = new TH1D("HnjetsDistribution","njets using nJetClean30",10,-0.5,9.5);   
  Hj1j2dphiDistribution = new TH1D("Hj1j2dphiDistribution","",32,0.0,3.2);
  HjetMetDphiMinDistribution = new TH1D("HjetMetDphiMinDistribution","",32,0.0,3.2);
  Hjet1etaDistribution = new TH1D("Hjet1etaDistribution","",100,-5.0,5.0);
  Hjet2etaDistribution = new TH1D("Hjet2etaDistribution","",100,-5.0,5.0);
  HrecoilDistribution = new TH1D("HrecoilDistribution","",100,0.0,1500.0);
  Hjet1ptDistribution = new TH1D("Hjet1ptDistribution","",100,0.0,1500.0); 
  Hjet2ptDistribution = new TH1D("Hjet2ptDistribution","",100,0.0,1500.0);

  // saving histograms with bin edges of other histograms used (e.g. content of metBinEdges array ...)
  HmetBinEdges = new TH1D("HmetBinEdges","bin edges for met distributions",nMetBins+1,0.0,nMetBins+1);
  for (Int_t i = 0; i <= nMetBins; i++) {
    HmetBinEdges->SetBinContent(i+1,metBinEdgesVector[i]);
  }


}

//===============================================

void AnalysisDarkMatter::setScaleFactorHistograms() {

  // all common scaled histograms are defined here. Some more specific to a CR region (e.g. for lepton ID) are defined in the proper .cc file

  HYieldsMetBin_qcdRenScaleUp = new TH1D("HYieldsMetBin_qcdRenScaleUp","yields in bins of met; #slash{E}_{T};# of events",nMetBins,metBinEdgesVector.data());
  HYieldsMetBin_qcdRenScaleDown = new TH1D("HYieldsMetBin_qcdRenScaleDown","yields in bins of met; #slash{E}_{T};# of events",nMetBins,metBinEdgesVector.data());
  HYieldsMetBin_qcdFacScaleUp = new TH1D("HYieldsMetBin_qcdFacScaleUp","yields in bins of met; #slash{E}_{T};# of events",nMetBins,metBinEdgesVector.data());
  HYieldsMetBin_qcdFacScaleDown = new TH1D("HYieldsMetBin_qcdFacScaleDown","yields in bins of met; #slash{E}_{T};# of events",nMetBins,metBinEdgesVector.data());
  HYieldsMetBin_qcdPdfUp = new TH1D("HYieldsMetBin_qcdPdfUp","yields in bins of met; #slash{E}_{T};# of events",nMetBins,metBinEdgesVector.data());
  HYieldsMetBin_qcdPdfDown = new TH1D("HYieldsMetBin_qcdPdfDown","yields in bins of met; #slash{E}_{T};# of events",nMetBins,metBinEdgesVector.data());
  HYieldsMetBin_ewkUp = new TH1D("HYieldsMetBin_ewkUp","yields in bins of met; #slash{E}_{T};# of events",nMetBins,metBinEdgesVector.data());
  HYieldsMetBin_ewkDown = new TH1D("HYieldsMetBin_ewkDown","yields in bins of met; #slash{E}_{T};# of events",nMetBins,metBinEdgesVector.data());

  HSyst_qcdRenScale = new TH1D("HSyst_qcdRenScale","systematic uncertainty for QCD renormalization scale",nMetBins,metBinEdgesVector.data());
  HSyst_qcdFacScale = new TH1D("HSyst_qcdFacScale","systematic uncertainty for QCD factorization scale",nMetBins,metBinEdgesVector.data());
  HSyst_qcdPdf = new TH1D("HSyst_qcdPdf","systematic uncertainty for QCD due to PDF uncertainty",nMetBins,metBinEdgesVector.data());
  HSyst_ewk = new TH1D("HSyst_ewk","systematic uncertainty for EWK",nMetBins,metBinEdgesVector.data());
  HSyst_total = new TH1D("HSyst_total","total systematic uncertainty (sum in quadrature of all single systematics)",nMetBins,metBinEdgesVector.data());


}

//===============================================

void AnalysisDarkMatter::setHistogramLastBinAsOverFlow(const Int_t hasScaledHistograms = 0) {

  myAddOverflowInLastBin(HYieldsMetBin);
  myAddOverflowInLastBin(HhtDistribution);
  myAddOverflowInLastBin(HrecoilDistribution);
  myAddOverflowInLastBin(Hjet1ptDistribution);
  myAddOverflowInLastBin(Hjet2ptDistribution);

  if (hasScaledHistograms) {

    myAddOverflowInLastBin(HYieldsMetBin_qcdRenScaleUp);
    myAddOverflowInLastBin(HYieldsMetBin_qcdRenScaleDown);
    myAddOverflowInLastBin(HYieldsMetBin_qcdFacScaleUp);
    myAddOverflowInLastBin(HYieldsMetBin_qcdFacScaleDown);
    myAddOverflowInLastBin(HYieldsMetBin_qcdPdfUp);
    myAddOverflowInLastBin(HYieldsMetBin_qcdPdfDown);
    myAddOverflowInLastBin(HYieldsMetBin_ewkUp);
    myAddOverflowInLastBin(HYieldsMetBin_ewkDown);

  }

}

//===============================================

void AnalysisDarkMatter::createSystematicsHistogram() {

  // in this method we handle the scaled histogram to all regions. Some histograms specific to a given CR are treated in the same method for those region, which will first call this ,ethod and then add the additional scale factors if any

  //for monojet

  // computing systematic uncertainties and saving them as histograms.
  myBuildSystematicsHistogram(HSyst_qcdRenScale, HYieldsMetBin, HYieldsMetBin_qcdRenScaleUp, HYieldsMetBin_qcdRenScaleDown);
  myBuildSystematicsHistogram(HSyst_qcdFacScale, HYieldsMetBin, HYieldsMetBin_qcdFacScaleUp, HYieldsMetBin_qcdFacScaleDown);
  myBuildSystematicsHistogram(HSyst_qcdPdf, HYieldsMetBin, HYieldsMetBin_qcdPdfUp, HYieldsMetBin_qcdPdfDown);
  myBuildSystematicsHistogram(HSyst_ewk, HYieldsMetBin, HYieldsMetBin_ewkUp, HYieldsMetBin_ewkDown);

  // define an empty histogram to sum uncertainties in a clean way
  TH1D *Htmp = new TH1D("Htmp","",nMetBins,metBinEdgesVector.data());
  vector<TH1D*> hptr;
  hptr.push_back(HSyst_qcdRenScale);
  hptr.push_back(HSyst_qcdFacScale);
  hptr.push_back(HSyst_qcdPdf);
  hptr.push_back(HSyst_ewk);

  // now HSyst should be empty, so the next step is not necessary
  HSyst_total->Multiply(HSyst_total);  // before adding anything we must get the square of this histogram, because we will add histograms in quadrature
     
  for (UInt_t i = 0; i < hptr.size(); i++) {
    Htmp->Multiply(hptr[i],hptr[i]); // square of bin content for each single systematic histogram
    HSyst_total->Add(Htmp);             // adding the squares
  }

  for (Int_t i = 0; i <= (HSyst_total->GetNbinsX() + 1); i++) {  // computing square root of each bin's content (from underflow to overflow bin, but they should be empty)
    HSyst_total->SetBinContent(i, sqrt(HSyst_total->GetBinContent(i)));
  }

  delete Htmp;
  

}

//===============================================

void AnalysisDarkMatter::fillRowVector(const Double_t nTotalWeightedEvents, const selectionManager &selMan, const mask & m, vector<Double_t> & yRow, vector<Double_t> & eRow, vector<Double_t> & uncRow, const Int_t hasPreselection_flag = 0) {

  //  hasPreselection_flag is used in control samples when the first selection(s) generally includes some gen-level cuts (which can be chosen by the user setting the GENLEP_TAG flag: in these cases, the tables will show the number of events after the preselection as the entry point.

  // The previous task exploits the selectionManager object which holds the indices of the selections in the mask which the user wants to show in final tables with yields

  Int_t index_TotalEntryAndPreselection = selMan.getFirstStepIndex() - 1;  // get the index corresponding to the selection before the one that was set as the first selection stored in selectionManager

  // entry point
  if (hasPreselection_flag) {
    yRow.push_back(m.nEvents[index_TotalEntryAndPreselection]); // step before the first step in the selectionManager (the one before metFiltersC for now)
    eRow.push_back(1.0000);
    uncRow.push_back(myGetUncertainty(&m, index_TotalEntryAndPreselection, uncertainty));
  } else {  
    yRow.push_back(nTotalWeightedEvents); // [0] 
    eRow.push_back(1.0000);
    uncRow.push_back(sqrt(nTotalWeightedEvents)); //should use a kind of myGetUncertainty function, but I don't save sum of newwgt^2 so I can't use MC uncertainty
  }

  for(Int_t i = 0; i < selMan.getVectorSize(); i++) {
     
    yRow.push_back(m.nEvents[selMan.getStepIndex(i)]);
    //uncRow.push_back(sqrt(yRow.back()));
    uncRow.push_back(myGetUncertainty(&m, selMan.getStepIndex(i), uncertainty));
    if (i == 0) {
      if (hasPreselection_flag) eRow.push_back(m.nEvents[selMan.getStepIndex(i)]/m.nEvents[index_TotalEntryAndPreselection]);
      else eRow.push_back(m.nEvents[selMan.getStepIndex(i)]/nTotalWeightedEvents);

    } else if( (i != 0) && (m.nEvents[selMan.getStepIndex(i)-1] == 0) ) eRow.push_back(1.0000);
    else eRow.push_back(m.nEvents[selMan.getStepIndex(i)]/m.nEvents[selMan.getStepIndex(i)-1]);
   
  }

}

//===============================================


void AnalysisDarkMatter::fillEventMask(ULong64_t & eventMask) {

  if (B_VETO_FLAG) eventMask += bjetVetoC.addToMask(nBTag15 < 0.5);
  if (TAU_VETO_FLAG) eventMask += tauLooseVetoC.addToMask(nTauClean18V < 0.5);  
  
}

//===============================================

// void AnalysisDarkMatter::fillEventMask(UInt_t & eventMask) {

//   if (B_VETO_FLAG) eventMask += bjetVetoC.addToMask(nBTag15 < 0.5);
//   if (TAU_VETO_FLAG) eventMask += tauLooseVetoC.addToMask(nTauClean18V < 0.5);  
  
// }

// //===============================================


#endif
