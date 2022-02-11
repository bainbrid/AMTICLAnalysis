#include <iostream>
#include <sstream>
#include <fstream>
//#include <boost/algorithm/string.hpp>

#include "TCanvas.h"
#include "TPad.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
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

#include "TDRStyle.h"
#include "plotUtils.h"

#include "standardVariables.h"


std::string reg[8] = {
  "CloseByPhotons","CloseByPhotonsFromVtx",
  "CloseByPhotonsWithPU","CloseByPhotonsFromVtxWithPU",
  "ChargedPionsFromVtx","ChargedPionsFromVtxWithPU",
  "ElectronsFromVtx","ElectronsFromVtxWithPU"
};
std::string regShort[8] = {
  "Unconv","FromVtx",
  "WithPU","FromVtxWithPU",
  "PiFromVtx","PiFromVtxWithPU",
  "EleFromVtx","EleFromVtxWithPU"
};


int color[12] = {1,2,3,4,6,7,8,9,kRed+2,kCyan+2,kGray,kOrange};

std::string gIter[9] = {
  "Dummy1",
  "Dummy2",
  "Dummy3",
  "EM1",
  "EM2",
  "EM3",
  "EM",
  "HAD",
  "Sim"
  //  "EM3relax",
  //"EMDef"
};

std::string gIterPi[17] = {
  "Dummy1",
  "Dummy2",
  "Dummy3",
  "EM1",
  "EM2",
  "EM3",
  "HAD1",
  "HAD2",
  "HAD3",
  "TRK1",
  "TRK2",
  "TRK3",
  "TrkEM",
  "EM",
  "Trk",
  "HAD",
  "Sim"
  //  "EM3relax",
  //"EMDef"
};

void makeHistogram1D(const histo & lHist,
		     const int iR,
		     std::vector<std::string> lCuts,
		     std::vector<bool> lInvertCuts,
		     TTree* treeIn,
		     TH1F* & hVar,
		     TCanvas* myc
		     ){

  const std::string lVarShort = lHist.varShort;
  const std::string lVarName = lHist.varName;
  const unsigned nBins = lHist.nB;
  const double binMin = lHist.binMin;
  const double binMax = lHist.binMax;
  const std::vector<double> & binsX = lHist.binsX;
  const std::string lVar = lHist.var;
  const bool setUnderflow = lHist.underflow;
  const bool doFixedBinning = lHist.doFixedBinning;

  std::ostringstream label;
  label << "hVar_"
	<< lVarShort;
  if (doFixedBinning) {
    const double * tmpBins = &binsX[0];
    hVar = new TH1F(label.str().c_str(),lVarName.c_str(),nBins,tmpBins);
    std::cout << " --- Variable binning: " << std::endl;
    std::cout << nBins << " " ;
    for (unsigned iB(0);iB<nBins+1;++iB){
      std::cout << binsX[iB] << ",";
    }
    std::cout << std::endl;
    //exit(1);
  }
  else hVar = new TH1F(label.str().c_str(),lVarName.c_str(),nBins,binMin,binMax);
  hVar->Sumw2();
  
  myc->cd();

  std::ostringstream lCut;
  
  const unsigned nC = lCuts.size();
  for (unsigned iC(0); iC<nC; ++iC){
    if (iC>0) lCut << " && ";
    if (lInvertCuts[iC]) lCut << "!";
    lCut << lCuts[iC];
  }

  if (lVarShort.find("_layer")!=lVarShort.npos){
    size_t startPos = lVarShort.find("_layer")+6;
    std::string layerNum = lVarShort.substr(startPos,lVarShort.size()-startPos);
    //std::cout << " - Check layer: " << layerNum << std::endl;
    if (nC>0) lCut << " && " ;
    lCut << "lc_layer == " << layerNum;
  }

  if (lVarShort.find("nSC")==lVarShort.npos &&
      lVarShort.find("nTS")==lVarShort.npos &&
      lVarShort.find("nCP")==lVarShort.npos){
    if (lCut.str().size()>0) lCut << " && ";
    lCut << "nTS>0";
  }

  if (lVarShort.find("triplets")!=lVarShort.npos && (lVarShort.find("A")!=lVarShort.npos || lVarShort.find("nner")!=lVarShort.npos || lVarShort.find("Beta")!=lVarShort.npos)){
    size_t startPos = lVarShort.find_last_of("_")+1;
    std::string layerNum = lVarShort.substr(startPos,lVarShort.size()-startPos);
    if (lCut.str().size()>0) lCut << " && ";
    lCut << "triplets_layerA_" << layerNum << ">0";
  }
  
  //std::cout << " - Cut string: " << lCut.str() << std::endl;
  treeIn->Draw((lVar+">>"+label.str()).c_str(),lCut.str().c_str(),"PE");
  
  if (setUnderflow){
    hVar->SetBinContent(1,hVar->GetBinContent(0)+hVar->GetBinContent(1));
    hVar->SetBinError(1,sqrt(pow(hVar->GetBinError(0),2)+pow(hVar->GetBinError(1),2)));
    hVar->SetBinContent(0,0);
  }
  //set overflows
  //if (lVarShort.find("dijet_M")==lVarShort.npos){
  hVar->SetBinContent(hVar->GetNbinsX(),hVar->GetBinContent(hVar->GetNbinsX())+hVar->GetBinContent(hVar->GetNbinsX()+1));
  hVar->SetBinError(hVar->GetNbinsX(),sqrt(pow(hVar->GetBinError(hVar->GetNbinsX()),2)+pow(hVar->GetBinError(hVar->GetNbinsX()+1),2)));
  hVar->SetBinContent(hVar->GetNbinsX()+1,0);
  //}
  
  //std::cout << " -- Entries: " << hVar->GetEntries() << ", Integral: " << hVar->Integral(0,hVar->GetNbinsX()+1)
  //	    << std::endl;
  
};

void makePlot(const unsigned iR,
	      const histo & lHist,
	      TH1F* & hVar,
	      TCanvas* & myc,
	      const bool first,
	      TLegend* & leg,
	      const std::string lPlotDir
	      ){
  
  const std::string lVarShort = lHist.varShort;
  const bool logy = lHist.logy;
  
  TLatex lat;
  
  TPad *upper = myc;
  upper->cd();
  gPad->SetGridx(1);
  //gPad->SetGridy(1);
  gPad->SetLogy(logy);
  hVar->SetLineColor(color[0]);
  hVar->SetMarkerColor(color[0]);
  hVar->SetMarkerStyle(20);
  //hVar->SetMinimum();
  //hVar->SetMaximum(200);
  hVar->Draw("PE1");
  //if (iY==0 && first) leg->AddEntry(hVar,aProc[iP].c_str(),"P");
  
  //      sprintf(buf,"N_{tot} = %3.0f",total);
  lat.DrawLatexNDC(0.1,0.96,(reg[iR]).c_str());
  //leg->Draw("same");
  if (lVarShort.find("_layer")!=lVarShort.npos){
    size_t startPos = lVarShort.find("_layer")+6;
    std::string layerNum = lVarShort.substr(startPos,lVarShort.size()-startPos);
    lat.DrawLatexNDC(0.8,0.96,("Layer "+layerNum).c_str());
  }
    
  myc->Update();
  myc->Print((lPlotDir+"/"+lVarShort+"_"+reg[iR]+".pdf").c_str());
  
};

