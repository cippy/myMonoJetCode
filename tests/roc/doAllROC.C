//ROOT header files                                                                                                          
#include <TROOT.h>
#include <TAttFill.h>
#include <TAxis.h>
#include <TCanvas.h>
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
#include <TPad.h>
#include <TPaveStats.h>
#include <TPaveText.h>
#include <TStyle.h>
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

void checkFile(const TFile* f, const string fname) {

  if (!f || !f->IsOpen()) {
    cout<<"*******************************"<<endl;
    cout<<"Error opening file \""<<fname<<"\".\nApplication will be terminated."<<endl;
    cout<<"*******************************"<<endl;
    exit(EXIT_FAILURE);
  }

}

//=====================================


void cloneTH1F(TH1F* h, TH1F* htmp, TFile* f, const string fname, const string TH1name, const string fdir) {

  htmp = (TH1F*) f->Get((fdir+"/"+TH1name).c_str());
  if (!htmp || htmp == NULL) {
    cout << "Error: histogram  ' " << TH1name << " ' not found in file ' " << fname << "'. End of programme." << endl;
    exit(EXIT_FAILURE);
  }
  h = (TH1F*) htmp->Clone();

}

//=======================================

void setHistColor(vector<Int_t> &histColor, const UInt_t nSamples) {

  Int_t colorList[] = {kCyan+1, kBlue+1, kGreen+1, kOrange+7, kRed+1, kGray+1, kCyan+2, kGreen+4, kViolet, kBlack};    
  for (UInt_t i = 0; i < nSamples; i++) {   // now color are assigned in reverse order (the main contribution is the last object in the sample array)             

    histColor.push_back(colorList[i]);

  }

}


//=====================================


void doROC(const string TMVAdir     = "outputTMVA/", 
	   const string outdir      = "/afs/cern.ch/user/m/mciprian/www/vbfHiggsToInv/optimization/roc/", 
	   const string ftype       = "single",   // can be single, double, combo
	   const string method      = "Cuts",     // Cuts, Likelihood, BDT
	   const string varCompared = "detajj"    // at the moment can be detajj or mjj
 	   ) 
{

  gROOT->SetStyle("Plain");  // to have white legend (on my pc it's already white, but in tier2 it appears grey)                                                  
  gStyle->SetHistTopMargin(0.);
  gStyle->SetFillColor(10);

  system(("mkdir -p " + outdir).c_str());
  // system("cd ~/www/");
  // system("./copyphp.sh");
  // system("cd -");

  vector<string> var;
  vector <string> fdir;
  vector<string> histoname;
  vector<string> fname;
  vector<string> leg_text;

  if (ftype == "single") {

    var.push_back("detajj");
    var.push_back("jetmetdphi");
    var.push_back("dphijj");
    var.push_back("mjj");
    var.push_back("mTjj");

    leg_text.push_back("detajj - " + method);
    leg_text.push_back("jetmetdphi - " + method);
    leg_text.push_back("dphijj - " + method);
    leg_text.push_back("mjj - " + method);
    leg_text.push_back("mTjj - " + method);

    fdir.push_back("Method_Cuts/Cuts");
    fdir.push_back("Method_Cuts/Cuts");
    fdir.push_back("Method_Cuts/Cuts");
    fdir.push_back("Method_Cuts/Cuts");
    fdir.push_back("Method_Cuts/Cuts");
    // var.push_back("jetpt1");
    // var.push_back("jetpt2");

    histoname.push_back("MVA_" + method + "_rejBvsS");
    histoname.push_back("MVA_" + method + "_rejBvsS");
    histoname.push_back("MVA_" + method + "_rejBvsS");
    histoname.push_back("MVA_" + method + "_rejBvsS");
    histoname.push_back("MVA_" + method + "_rejBvsS");

  } else if (ftype == "double") {

    if (varCompared == "detajj") {

      var.push_back("detajj");
      var.push_back("detajj_dphijj");
      var.push_back("detajj_jetmetdphi");
      //var.push_back("detajj_njet");
      var.push_back("detajj_mTjj");
      var.push_back("mjj_detajj");

      leg_text.push_back("detajj - Cuts");
      leg_text.push_back("detajj_dphijj - " + method);
      leg_text.push_back("detajj_jetmetdphi - " + method);
      leg_text.push_back("detajj_mTjj - " + method);
      leg_text.push_back("mjj_detajj - " + method);

    } else if (varCompared == "mjj") {

      var.push_back("mjj");
      var.push_back("mjj_detajj");
      var.push_back("mjj_dphijj");
      var.push_back("mjj_jetmetdphi");
      //var.push_back("mjj_njet");
      var.push_back("mjj_mTjj"); 

      leg_text.push_back("mjj - Cuts");
      leg_text.push_back("mjj_detajj - " + method);
      leg_text.push_back("mjj_dphijj - " + method);
      leg_text.push_back("mjj_jetmetdphi - " + method);
      leg_text.push_back("mjj_mTjj - " + method);

    }

    fdir.push_back("Method_Cuts/Cuts");
    fdir.push_back("Method_" + method + "/" + method);
    fdir.push_back("Method_" + method + "/" + method);
    //fdir.push_back("Method_" + method + "/" + method);
    fdir.push_back("Method_" + method + "/" + method);
    fdir.push_back("Method_" + method + "/" + method);

    histoname.push_back("MVA_Cuts_rejBvsS");
    histoname.push_back("MVA_" + method + "_rejBvsS");
    histoname.push_back("MVA_" + method + "_rejBvsS");
    //histoname.push_back("MVA_" + method + "_rejBvsS");
    histoname.push_back("MVA_" + method + "_rejBvsS");
    histoname.push_back("MVA_" + method + "_rejBvsS");


  } else if (ftype == "combo") {

    var.push_back("detajj");
    var.push_back("all");

    leg_text.push_back("detajj - Cuts");
    leg_text.push_back("all - BDTG");

    fdir.push_back("Method_Cuts/Cuts");
    fdir.push_back("Method_BDT/BDTG");

    histoname.push_back("MVA_Cuts_rejBvsS");
    histoname.push_back("MVA_BDTG_rejBvsS");

  } else if (ftype == "best_of_type") {

    var.push_back("detajj");
    var.push_back("detajj_jetmetdphi");
    var.push_back("mjj_detajj");
    var.push_back("detajj_jetmetdphi");
    var.push_back("mjj_mTjj");
    var.push_back("all");

    leg_text.push_back("detajj - Cuts");
    leg_text.push_back("detajj_jetmetdphi - Cuts");
    leg_text.push_back("mjj_detajj - Cuts");
    leg_text.push_back("detajj_jetmetdphi - Likelihood");
    leg_text.push_back("mjj_mTjj - Likelihood");
    leg_text.push_back("all - BDTG");

    fdir.push_back("Method_Cuts/Cuts");
    fdir.push_back("Method_Cuts/Cuts");
    fdir.push_back("Method_Cuts/Cuts");
    fdir.push_back("Method_Likelihood/Likelihood");
    fdir.push_back("Method_Likelihood/Likelihood");
    fdir.push_back("Method_BDT/BDTG");
 
    histoname.push_back("MVA_Cuts_rejBvsS");
    histoname.push_back("MVA_Cuts_rejBvsS");
    histoname.push_back("MVA_Cuts_rejBvsS");
    histoname.push_back("MVA_Likelihood_rejBvsS");
    histoname.push_back("MVA_Likelihood_rejBvsS");
    histoname.push_back("MVA_BDTG_rejBvsS");

    fname.push_back(TMVAdir + "TrainingResult_detajj_cuts.root");
    fname.push_back(TMVAdir + "TrainingResult_detajj_jetmetdphi.root");
    fname.push_back(TMVAdir + "TrainingResult_mjj_detajj.root");
    fname.push_back(TMVAdir + "TrainingResult_detajj_jetmetdphi.root");
    fname.push_back(TMVAdir + "TrainingResult_mjj_mTjj.root");
    fname.push_back(TMVAdir + "TrainingResult_combo.root");

  }

  for (UInt_t i = 0; i < var.size(); i++) {
    if (ftype != "best_of_type") {  // for "best_of_type" the fname is set above by hand (there is no simple way to have it automatically)
      if ( var[i].find("_") != string::npos) fname.push_back(TMVAdir + "TrainingResult_" + var[i] + ".root");
      else if (var[i].find("all") != string::npos) fname.push_back(TMVAdir + "TrainingResult_combo.root");
      else fname.push_back(TMVAdir + "TrainingResult_" + var[i] + "_cuts.root");
    }
  }


  TFile *f = NULL;

  TH1F *htmp = NULL;
  vector<TH1F*> histo;
  vector <TH1F*> histo_sign; // significance as epsilonS/sqrt(epsilonB)

  vector<Int_t> histColor;
  setHistColor(histColor,fname.size());

  Double_t legYmax = 0.1 + 0.08 * var.size();
  TLegend *leg = NULL;
  if (ftype == "single") leg = new TLegend(0.15,0.15,0.6,legYmax);
  else leg = new TLegend(0.15,0.15,0.7,legYmax);

  //  cloneTH1F(histo, htmp, f, fname, histoname, fdir);
  ///////////////

  // loop on files anfd get histograms

  for (UInt_t i = 0; i < fname.size(); i++) {

    f = new TFile(fname[i].c_str(),"READ");
    checkFile(f,fname[i]);

    htmp = (TH1F*)f->Get((fdir[i]+"/"+histoname[i]).c_str());
    if (!htmp || htmp == NULL) {
      cout << "Error: histogram  ' " << histoname[i] << " ' not found in file ' " << fname[i] << "'. End of programme." << endl;
      exit(EXIT_FAILURE);
    }
    histo.push_back((TH1F*) htmp->Clone());

    histo_sign.push_back((TH1F*) htmp->Clone());
    Int_t rebinFactor = 4;
    histo_sign.back()->Rebin(rebinFactor);
    histo_sign.back()->Scale(1./rebinFactor);
    // histo_sign.push_back( ( (TH1F*) ( (TH1F*) htmp->Clone() )->Rebin( rebinFactor, Form("hist_sign%d",i) ) ) );

  }

  ///////////////

  TCanvas* c = new TCanvas("c","",700,700);

  gPad->SetTickx();
  gPad->SetTicky();

  histo[0]->Draw("L");
  histo[0]->SetTitle("");;
  histo[0]->GetYaxis()->SetTitle("1 - #varepsilon _{B}");
  histo[0]->GetYaxis()->SetLabelSize(0.035);
  histo[0]->GetYaxis()->SetTitleSize(0.05);
  histo[0]->GetYaxis()->SetTitleOffset(0.9);
 
  histo[0]->GetXaxis()->SetTitle("#varepsilon_{S}");
  histo[0]->GetXaxis()->SetLabelSize(0.035);
  histo[0]->GetXaxis()->SetTitleSize(0.05);
  histo[0]->GetXaxis()->SetTitleOffset(0.8);

  for (UInt_t i = 0; i < histo.size(); i++) {
    histo[i]->Draw("L SAME");
    histo[i]->SetLineWidth(2);
    histo[i]->SetLineColor(histColor[i]);
    //histo[i]->SetFillColor(histColor[i]);
    leg->AddEntry(histo[i],leg_text[i].c_str(),"lf");
  }

  gStyle->SetStatStyle(0);
  leg->SetFillStyle(0);  // transparent legend    
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();
  leg->SetMargin(0.3);
  leg->SetBorderSize(0);

  string cname = outdir + "roc_" + ftype + "_" + method;
  if (ftype == "double" || ftype == "combo") cname = outdir + "roc_" + ftype + "_" + method + "_" + varCompared;
  if (ftype == "best_of_type") cname = outdir + "roc_bestROCs";
  c->SaveAs((cname + ".pdf").c_str());
  c->SaveAs((cname + ".png").c_str());

  delete c;

  Double_t max_sign = 0.0;
  Double_t min_sign = 10000000.0;
  // now significances
  for (UInt_t h = 0; h < histo_sign.size(); h++ ) {
    //    histo_sign[h] = (TH1F*) histo_sign[h]->Rebin(rebinFactor);
    for (Int_t i = 1; i <= histo_sign[h]->GetNbinsX(); i++) {
      histo_sign[h]->SetBinContent(i, histo_sign[h]->GetBinCenter(i) / sqrt( 1 - histo_sign[h]->GetBinContent(i)) );  // epsilonS/ sqrt(epsilonB)
    }
    if (max_sign  <= histo_sign[h]->GetBinContent( histo_sign[h]->GetMaximumBin() ) ) max_sign = histo_sign[h]->GetBinContent( histo_sign[h]->GetMaximumBin() );
    if (min_sign  >= histo_sign[h]->GetBinContent( histo_sign[h]->GetMinimumBin() ) ) min_sign = histo_sign[h]->GetBinContent( histo_sign[h]->GetMinimumBin() );
  }

  TCanvas* c_sig = new TCanvas("c_sig","",700,700);
  TLegend *leg_sig = new TLegend(*leg);

  gPad->SetTickx();
  gPad->SetTicky();

  histo_sign[0]->Draw("C");
  histo_sign[0]->SetTitle("");
  histo_sign[0]->SetMaximum(max_sign * 1.1);  
  histo_sign[0]->SetMinimum(min_sign * 0.9);  
  histo_sign[0]->GetYaxis()->SetTitle("#varepsilon_{S} / #sqrt{#varepsilon _{B}}");
  histo_sign[0]->GetYaxis()->SetLabelSize(0.035);
  histo_sign[0]->GetYaxis()->SetTitleSize(0.05);
  histo_sign[0]->GetYaxis()->SetTitleOffset(0.9);
 
  histo_sign[0]->GetXaxis()->SetTitle("#varepsilon_{S}");
  histo_sign[0]->GetXaxis()->SetLabelSize(0.035);
  histo_sign[0]->GetXaxis()->SetTitleSize(0.05);
  histo_sign[0]->GetXaxis()->SetTitleOffset(0.8);

  for (UInt_t i = 0; i < histo_sign.size(); i++) {
    histo_sign[i]->Draw("C SAME");
    histo_sign[i]->SetLineWidth(2);
    histo_sign[i]->SetLineColor(histColor[i]);
    //histo_sign[i]->SetFillColor(histColor[i]);
    //leg_sig->AddEntry(histo_sign[i],leg_text[i].c_str(),"lf");
  }

  gStyle->SetStatStyle(0);
  leg_sig->SetFillStyle(0);  // transparent legend    
  leg_sig->SetFillColor(0);
  leg_sig->SetBorderSize(0);
  leg_sig->Draw();
  leg_sig->SetMargin(0.3);
  leg_sig->SetBorderSize(0);

  cname = outdir + "significance_" + ftype + "_" + method;
  if (ftype == "double" || ftype == "combo") cname = outdir + "significance_" + ftype + "_" + method + "_" + varCompared;
  if (ftype == "best_of_type") cname = outdir + "significance_bestROCs";
  c_sig->SaveAs((cname + ".pdf").c_str());
  c_sig->SaveAs((cname + ".png").c_str());

  delete c_sig;


 
}



// ================================

void doAllROC(string TMVAdir     = "outputTMVA", 
	      string outdir  = "baseSel_met175"
	      ) 
{
  
  outdir = "/afs/cern.ch/user/m/mciprian/www/vbfHiggsToInv/optimization/roc/" + outdir + "/";
  TMVAdir = TMVAdir + "/";

  doROC(TMVAdir, outdir, "single", "Cuts", "detajj");
  doROC(TMVAdir, outdir, "double", "Cuts", "detajj");
  doROC(TMVAdir, outdir, "double", "Likelihood", "detajj");
  doROC(TMVAdir, outdir, "double", "Cuts", "mjj");
  doROC(TMVAdir, outdir, "double", "Likelihood", "mjj");
  doROC(TMVAdir, outdir, "combo", "BDT", "detajj");
  doROC(TMVAdir, outdir, "best_of_type", "", "");

}

void doAllROC_multiple() {
  
  doAllROC("outputTMVABasicPreselection","baseSel");
  doAllROC("outputTMVABasicPreselectionQCD","baseSel_QCD_bckg");
  doAllROC("outputTMVABasicPreselectionEWK","baseSel_EWK_bckg");
  doAllROC("outputTMVATightPreselection","tightSel");

}
