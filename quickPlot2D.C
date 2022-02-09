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

#include "globalVars.h"

int quickPlot2D(){

  const bool doAll = true;
  
  
  SetTdrStyle();
  gStyle->SetOptStat(0);

  const unsigned nR = 50;
  const unsigned nV = 9;
  const unsigned nPhi = 6;//doAll ? 1 : 6;
  unsigned phi[6] = {0,10,20,30,40,50};

  TCanvas *myc[nV];

  for (unsigned iV(0); iV<nV; ++iV){
    std::ostringstream lName;
    lName << "myc" << iV;
    myc[iV] = new TCanvas(lName.str().c_str(),lName.str().c_str(),doAll?1000:1500,1000);
    if (!doAll) myc[iV]->Divide(3,2);
  }
  
  
  std::string filePath = "D49_PhiScan/Photons/step3ticl_E500_phi";
  
  TH2F* h2D[nPhi][nV];
  TH2F* h2DAll[nV];

  std::string histName[9] = {
    "all","siThick","radius","E",
    "XYall","XYsiThick","XYradius","XYE",
    "zLayer"
  };
  
  for (unsigned iP(0); iP<nPhi; ++iP){//loop on phi
    std::cout << " - Processing phi " << phi[iP] << std::endl;
    bool isFirst = true;
    for (unsigned iR(0); iR<50; ++iR){
      std::ostringstream lFile;
      lFile << filePath << phi[iP] << "_run" << iR << "_FlatTracksters.root";
      TFile * fTmp = 0;
      fTmp = TFile::Open(lFile.str().c_str());
      if (!fTmp) continue;
      else std::cout << fTmp->GetName() << std::endl;
      fTmp->cd("ticlTree");
      
      for (unsigned iV(0); iV<nV; ++iV){
	//std::cout << " - Processing " << histName[iV] << std::endl;
	TH2F *htmp = (TH2F*)gDirectory->Get(histName[iV].c_str());
	if (!htmp) continue;
	std::ostringstream htitle;
	htitle << histName[iV] << "_phi" << phi[iP];
	if (isFirst) h2D[iP][iV] = (TH2F*)htmp->Clone(htitle.str().c_str());
	else h2D[iP][iV]->Add(htmp);
      }
      if (!isFirst) fTmp->Close();
      isFirst = false;
      
    }
  }
  
  for (unsigned iV(0); iV<nV; ++iV){
    std::cout << " - Processing " << histName[iV] << std::endl;
    std::ostringstream htitle;
    htitle << histName[iV] << "_allPhi";
    h2DAll[iV] = (TH2F*)h2D[0][iV]->Clone(htitle.str().c_str());
    for (unsigned iP(1); iP<nPhi; ++iP){//loop on phi
      h2DAll[iV]->Add(h2D[iP][iV]);
    }

    if (!doAll){
      for (unsigned iP(0); iP<nPhi; ++iP){//loop on phi
	if (iV>3 && iV<8)  h2D[iP][iV]->RebinX(10);
	if (iV>3 && iV<8)  h2D[iP][iV]->RebinY(10);
	
	if (iV<4 && iV>0) h2D[iP][iV]->Divide(h2D[iP][iV],h2D[iP][0],1,1);
	if (iV>4 && iV<8) h2D[iP][iV]->Divide(h2D[iP][iV],h2D[iP][4],1,1);
	
	myc[iV]->cd(iP+1);
	if (iV%4==0 || iV%4==3) gPad->SetLogz(1);
	gPad->SetRightMargin(0.2);
	
	if (iV%4==1) h2D[iP][iV]->GetZaxis()->SetRangeUser(100,400);
	if (iV%4==2) h2D[iP][iV]->GetZaxis()->SetRangeUser(0,1);
	
	h2D[iP][iV]->Draw("colz");
      }//loop on phi
    }
    else{
      myc[iV]->cd();
      if (iV>3 && iV<8)  h2DAll[iV]->RebinX(10);
      if (iV>3 && iV<8)  h2DAll[iV]->RebinY(10);
      
      if (iV<4 && iV>0) h2DAll[iV]->Divide(h2DAll[iV],h2DAll[0],1,1);
      if (iV>4 && iV<8) h2DAll[iV]->Divide(h2DAll[iV],h2DAll[4],1,1);
      if (iV%4==0 || iV%4==3) gPad->SetLogz(1);
      gPad->SetRightMargin(0.2);
      
      if (iV%4==1) h2DAll[iV]->GetZaxis()->SetRangeUser(100,400);
      if (iV%4==2) h2DAll[iV]->GetZaxis()->SetRangeUser(0,1);
      
      h2DAll[iV]->Draw("colz");
    }
    
    
    myc[iV]->Update();
    if (doAll) myc[iV]->Print(("Plot2D_allPhi_"+histName[iV]+".pdf").c_str());
    else myc[iV]->Print(("Plot2D_"+histName[iV]+".pdf").c_str());
  }
  

  return 0;
}





