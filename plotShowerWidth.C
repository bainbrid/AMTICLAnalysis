#include <algorithm>
#include <iomanip>
#include <stdlib.h>
#include <iostream>
#include <fstream>

#include "TROOT.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TF1.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TLatex.h"
#include "TLine.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TProfile.h"
#include "TLatex.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TGaxis.h"
#include "TString.h"
#include "TExec.h"

#include "TDRStyle.h"

int plotShowerWidth(){

  const bool doWeighted = true;
  const bool useSeed = true;

  SetTdrStyle();
  
  TFile *fin = TFile::Open("D49_FineCalo_AMiters/Reco/pt50_eta21/Debug/RecoTree_PiFromVtx.root");
  fin->cd();
  TTree *RecoTree = (TTree*)gDirectory->Get("RecoTree");

  TCanvas *myc = new TCanvas("myc","myc",1);
  myc->cd();
  gStyle->SetOptStat(0);
  gPad->SetLogy(1);
  
  TH1F *hdR[49];
  for (unsigned iL(0); iL<49;++iL){
    std::ostringstream lname;
    lname << "hdR_" << iL;
    if (doWeighted) hdR[iL] = new TH1F(lname.str().c_str(),";max [E_{LC}/E_{totLayer}#times #Delta #eta (truth) ]",30,0,0.6);
    else hdR[iL] = new TH1F(lname.str().c_str(),";max [#Delta #eta (truth)]",30,0,0.6);
    std::ostringstream lvar,lcut;
    lvar << "max";
    if (doWeighted) lvar << "w";
    lvar << "DeltaR";
    if (useSeed) lvar << "seed";
    lvar << "_" << iL+1 << ">>" <<lname.str();
    lcut << "max";
    if (doWeighted) lcut << "w";
    lcut << "DeltaR";
    if (useSeed) lcut << "seed";
    lcut << "_" << iL+1 << ">0.001";
    RecoTree->Draw(lvar.str().c_str());
  }
  
  TLegend *leg = new TLegend(0.8,0.4,0.9,0.9);
  for (unsigned iL(0); iL<49;++iL){
    if (iL%4!=0) continue;
    hdR[iL]->SetLineColor(iL<28?1:iL<36?2:4);
    hdR[iL]->SetMarkerColor(iL<28?1:iL<36?2:4);
    hdR[iL]->SetMarkerStyle(20+iL/4);
    hdR[iL]->DrawNormalized(iL==0?"P":"Psame");
    std::ostringstream label;
    label << "Layer " << iL+1;
    leg->AddEntry(hdR[iL],label.str().c_str(),"P");
  }
  leg->Draw("same");
  
  TLine *l1 = new TLine(0.15/2.,0,0.15/2,1);
  TLine *l2 = new TLine(0.25/2.,0,0.25/2,1);
  
  l1->SetLineColor(6); l1->SetLineWidth(2); l1->Draw("same");
  l2->SetLineColor(6); l2->SetLineWidth(2); l2->Draw("same");

  TLatex lat;
  lat.SetTextColor(6);
  lat.DrawLatex(0.15/2,0.3,"Search window 0.15/2");
  lat.DrawLatex(0.25/2,0.15,"Search window 0.25/2");
  lat.SetTextColor(1);
  if (useSeed) lat.DrawLatexNDC(0.5,0.94,"Using seed position");
  
  myc->Update();
  if (useSeed){
    if (doWeighted) myc->Print("showerWidthWithSeedEnergyWeighted.png");
    else myc->Print("showerWidthWithSeed.png");
  }
  else {
    if (doWeighted) myc->Print("showerWidthEnergyWeighted.png");
    else myc->Print("showerWidth.png");
  }

  return 0;
}
