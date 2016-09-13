#define vbfHiggsToInvAna_cxx
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

#ifdef vbfHiggsToInvAna_cxx

//===============================================

vbfHiggsToInvAna::vbfHiggsToInvAna(TTree *tree) : AnalysisDarkMatter(tree) {
  //cout <<"check in constructor "<<endl;

  // initialize some variables with sensible values. They will be set later depending on config file
  // following variables, if any, are specific to vbf H->inv analysis  

}

//===============================================


void vbfHiggsToInvAna::setNumberParameterValue(const string parameterName, const Double_t value) {

  AnalysisDarkMatter::setNumberParameterValue(parameterName, value);

}

//===============================================

void vbfHiggsToInvAna::setVarFromConfigFile() {

  AnalysisDarkMatter::setVarFromConfigFile();

}

//===============================================

void vbfHiggsToInvAna::setSelections() {

  AnalysisDarkMatter::setSelections();

  selection::checkMaskLength();

}

//===============================================

void vbfHiggsToInvAna::setHistograms() {

  AnalysisDarkMatter::setHistograms();

}

//===============================================

void vbfHiggsToInvAna::setHistogramLastBinAsOverFlow(const Int_t hasScaledHistograms = 0) {

  AnalysisDarkMatter::setHistogramLastBinAsOverFlow(hasScaledHistograms);

}

//===============================================

void vbfHiggsToInvAna::createSystematicsHistogram() {

  AnalysisDarkMatter::createSystematicsHistogram();  

}

//===============================================


void vbfHiggsToInvAna::fillEventMask(UInt_t & eventMask) {

  AnalysisDarkMatter::fillEventMask(eventMask);  
  
}

//===============================================


#endif
