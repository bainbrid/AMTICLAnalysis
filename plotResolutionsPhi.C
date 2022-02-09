#include <algorithm>
#include <iomanip>
#include <stdlib.h>
#include <iostream>
#include <fstream>

#include "TCanvas.h"
#include "TPad.h"
#include "TFile.h"
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

#include "makeResolution.C"
#include "globalVars.h"


struct MyReso {
  std::string iter;
  std::string pteta;
  double pt;
  double eta;
  double phi;
  double E;
  std::string label;
  double integral;
  double mean;
  double meanerr;
  double rms;
  double rmserr;
  double fitPar0;
  double fitPar0err;
  double fitPar1;
  double fitPar1err;
  double fitPar2;
  double fitPar2err;
  double resHist;
  double resFit;
  void Print(){
    std::cout << iter << " "
	      << pt << " "
	      << eta<< " "
	      << phi<< " "
	      << label
	      << std::endl;
  };

};


void readResolution(std::string aFileName,
		    double ptval, double etaval,
		    double phival,
		    std::vector<MyReso> & vec){

  std::ifstream infile(aFileName);
  //vec.clear();
  
  if (!infile.is_open()){
    std::cout << " -- file " << aFileName << " not found." << std::endl;
    return;
  }

  int counter =0 ;

  while(!infile.eof()){
    MyReso lRes;
    lRes.pt = ptval;
    lRes.eta = etaval;
    lRes.phi = phival;
    infile>>lRes.iter
      //>>lRes.pt>>lRes.eta
	  >>lRes.pteta
	  >>lRes.label>>lRes.integral
	  >>lRes.mean>>lRes.meanerr
	  >>lRes.rms>>lRes.rmserr
	  >>lRes.fitPar0>>lRes.fitPar0err
	  >>lRes.fitPar1>>lRes.fitPar1err
	  >>lRes.fitPar2>>lRes.fitPar2err
	  >>lRes.resHist>>lRes.resFit
      ;
    if (lRes.iter.size()>0) vec.push_back(lRes);
    counter++;
    if (counter>2000) break;
  }
  std::cout << " -- Found " << vec.size() << " (" << counter << ") elements." << std::endl;
  
};

void Meanvseta(const std::vector<MyReso> & vec,
	       std::string label,
	       const double phival,
	       TGraphErrors* & gr,
	       bool useFit,
	       unsigned doMean,
	       bool dovsE){
  
  int nP = 0;
  for (unsigned iE(0); iE<vec.size(); ++iE){//loop on vec elements
    MyReso lRes = vec[iE];
    lRes.E = lRes.pt;//*cosh(lRes.eta/10.);
    if (lRes.label != label || fabs(lRes.phi-phival)>0.01) continue;
    //lRes.Print();
    if (useFit) {
      if (doMean==0){
	gr->SetPoint(nP,lRes.eta,lRes.resFit);
	gr->SetPointError(nP,0.05,lRes.fitPar2err/lRes.fitPar2*lRes.resFit);
      }
      else if (doMean==1){
	gr->SetPoint(nP,lRes.eta,lRes.fitPar1);
	gr->SetPointError(nP,0.05,lRes.fitPar1err);
      }
      else if (doMean==2){
	gr->SetPoint(nP,lRes.eta,lRes.fitPar2);
	gr->SetPointError(nP,0.05,lRes.fitPar2err);
      }
    }
    else {
      if (doMean==0){
	gr->SetPoint(nP,lRes.eta,lRes.resHist);
	gr->SetPointError(nP,0.05,lRes.rmserr/lRes.rms*lRes.resHist);
      }
      else if (doMean==1){
	gr->SetPoint(nP,lRes.eta,lRes.mean);
	gr->SetPointError(nP,0.05,lRes.meanerr);
      }
      else if (doMean==2){
	gr->SetPoint(nP,lRes.eta,lRes.rms);
	gr->SetPointError(nP,0.05,lRes.rmserr);
      }
    }
    nP++;
  }//loop on vec elements
  
  gr->SetTitle(";#eta;#sigma/mean");
  if (doMean==1) gr->SetTitle(";#eta;mean");
  if (doMean==2) gr->SetTitle(";#eta;#sigma");
};


int plotResolutionsPhi(){


  SetTdrStyle();

  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);

  //0 = reso, 1 = mean, 2 = sigma
  const unsigned doMean = 2;

  const unsigned iReg = 0;
  const bool dovsE = true;
  const unsigned nIter = 1;
  const unsigned nI = 3;
  std::string iter[4] = {"Dummy1","Dummy3","EM3","EM"};

  const unsigned nRec = 1;
  unsigned nH[3] = {1,2,3};

  const unsigned nEta = iReg==0 ? 10 : 1;
  double etaval[10] = {17,18,19,20,21,
		       22,23,24,25,26};
  if (iReg!=0) etaval[0] = 21;
  double eta[10] = {1.75,1.85,1.95,2.05,2.15,
		    2.25,2.35,2.45,2.55,2.65};
  if (iReg!=0) eta[0] = 2.1;
  double etaerr[10] = {0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05};
  const unsigned nPT = 1;
  double ptval[12] = {500,5,10,15,20,30,40,50,75,100,150,200};
  const unsigned nV = 10;
  const unsigned nV_allLC = 6;//first 3 - independent of iter
  const unsigned nV_TS = 4;//last 3 -- independent of nRec
  std::string label[nV] = {"EtotOverECP_1RH","EtotOverECP_2RH","EtotOverECP_3RH",
			   "E1LCOverECP_1RH","E1LCOverECP_2RH","E1LCOverECP_3RH",
			   "ETSOverECP","E1LCTSOverECP","ETSmaxEOverECP","SumETSOverECP"
  };

  const unsigned nPhi = 6;
  unsigned phi[nPhi] = {0,10,20,30,40,50};

  int colour[6] = {1,2,3,4,6,7};
  
  TGraphErrors *grFit[nIter][nRec][nPhi][nV];
  TGraphErrors *grHist[nIter][nRec][nPhi][nV];
  for (unsigned iPhi(0); iPhi<nPhi; ++iPhi){//loop on phi  
    for (unsigned iI(0); iI<nIter; ++iI){
      for (unsigned iR(0); iR<nRec; ++iR){
	for (unsigned iV(0); iV<nV; ++iV){//loop on variables
	  grFit[iI][iR][iPhi][iV] = 0;
	  grHist[iI][iR][iPhi][iV] = 0;
	}
      }
    }
  }
    
  for (unsigned iPhi(0); iPhi<nPhi; ++iPhi){//loop on phi  
    for (unsigned iI(0); iI<nIter; ++iI){
      for (unsigned iR(0); iR<nRec; ++iR){
	std::cout << " - " << iter[iI] << std::endl;
	//<< " nRecHits=" << nH[iR] << std::endl;
	
	std::vector<MyReso> vec;
	vec.clear();
	for (unsigned ipt(0); ipt<nPT; ++ipt){  
	  for (unsigned ieta(0); ieta<nEta; ++ieta){
	    std::ostringstream etaStr;
	    etaStr << "_eta" << etaval[ieta] << "_" << etaval[ieta]+1;
	    std::cout << " --- Processing E " << ptval[ipt] << " phi " << phi[iPhi]<< std::endl;
	    
	    std::ostringstream pteta;
	    pteta << "E" << ptval[ipt] << "_phi" << phi[iPhi];
	    
	    std::ostringstream filename;
	    //filename << "D49_All/" << iter[iI] << "/LCPlots/NREC" << nH[iR] << "/resolutions.dat";
	    filename << "D49_PhiScan/Merged/" << iter[iI] << "/" << pteta.str() << etaStr.str() << "/resolutions_" << regShort[iReg] << ".dat";
	    
	    readResolution(filename.str().c_str(),ptval[ipt],eta[ieta],phi[iPhi],vec);
	  }
	}
	
	
	if (vec.size()==0) {
	  std::cout << " -- Vec has no elements" << std::endl;
	  continue;
	}
	
	for (unsigned iV(0); iV<nV; ++iV){//loop on variables
	  std::cout << " - Processing " << label[iV] << std::endl;
	  
	  grFit[iI][iR][iPhi][iV] = new TGraphErrors();
	  std::ostringstream lName;
	  lName << "grFit_" << iter[iI] << "_" << nH[iR] << "_" << label[iV] << "_phi" << phi[iPhi];
	  grFit[iI][iR][iPhi][iV]->SetName(lName.str().c_str());
	  Meanvseta(vec,label[iV],phi[iPhi],grFit[iI][iR][iPhi][iV],true,doMean,dovsE);
	  grHist[iI][iR][iPhi][iV] = new TGraphErrors();
	  lName.str("");
	  lName << "grHist_" << iter[iI] << "_" << nH[iR] << "_" << label[iV] << "_phi" << phi[iPhi];
	  grHist[iI][iR][iPhi][iV]->SetName(lName.str().c_str());
	  Meanvseta(vec,label[iV],phi[iPhi],grHist[iI][iR][iPhi][iV],false,doMean,dovsE);
	  std::cout << " - Npoints: " << " " << grFit[iI][iR][iPhi][iV]->GetN()
		    << " " << grHist[iI][iR][iPhi][iV]->GetN() << " "
		    << std::endl;
	}
      }
    }
  }
  
  TCanvas *mycFit  = new TCanvas("mycFit","mycFit",1);
  mycFit->cd();
  
  TLatex lat;
  char buf[200];
  
  
  for (unsigned iV(0); iV<nV; ++iV){//loop on variables
    //if (iV<(nV-nV_TS)) continue;
    for (unsigned iR(0); iR<nRec; ++iR){
      for (unsigned iI(0); iI<nIter; ++iI){
	if (iV<nV_allLC && iI>0) continue;
	bool isFirst = true;
	double maxY = 0;
	double minY = 10;
	mycFit->Clear();
	mycFit->cd();
	gPad->SetGridy(1);
	gPad->SetLogx(0);
	TLegend *leg = new TLegend(0.15,0.8,0.94,0.92);
	leg->SetFillStyle(0);
	leg->SetNColumns(3);

	for (unsigned iPhi(0); iPhi<nPhi; ++iPhi){
	  std::cout << " -- Processing " << label[iV] << " eta " << phi[iPhi] << std::endl;
	  if (!grFit[iI][iR][iPhi][iV]) continue;
	  if (!grHist[iI][iR][iPhi][iV]) continue;
	  std::cout << " -- Npoint = " << grFit[iI][iR][iPhi][iV]->GetN() << std::endl;
	  if (grFit[iI][iR][iPhi][iV]->GetN()==0) {
	    continue;
	  }
	  double lmax = TMath::MaxElement(grFit[iI][iR][iPhi][iV]->GetN(),grFit[iI][iR][iPhi][iV]->GetY());
	  double lmin = TMath::MinElement(grFit[iI][iR][iPhi][iV]->GetN(),grFit[iI][iR][iPhi][iV]->GetY());
	  if (lmax>maxY) maxY = lmax;
	  if (lmin<minY) minY = lmin;
	  std::cout << " - " << phi[iPhi] << " " << lmin << " " << lmax << std::endl;
	}
	
	for (unsigned iPhi(0); iPhi<nPhi; ++iPhi){
	  
	  if (!grHist[iI][iR][iPhi][iV]) continue;
	  grHist[iI][iR][iPhi][iV]->SetMarkerSize(1.5);
	  grHist[iI][iR][iPhi][iV]->SetMarkerStyle(20+iPhi);
	  grHist[iI][iR][iPhi][iV]->SetMarkerColor(iV<nV_allLC? kGray+2:kGray+iPhi);
	  grHist[iI][iR][iPhi][iV]->SetLineColor(iV<nV_allLC? kGray+2:kGray+iPhi);
	  grHist[iI][iR][iPhi][iV]->SetMaximum(maxY*1.1);
	  grHist[iI][iR][iPhi][iV]->SetMinimum(doMean!=1?0:minY*0.9);
	  //grHist[iI][iR][iPhi][iV]->Draw(isFirst ? "APE":"PEsame");
	  std::ostringstream phiStr;
	  phiStr << phi[iPhi];
	  //leg->AddEntry(grHist[iI][iR][iPhi][iV],(phiStr.str()+" hist").c_str(),"P");
	  //if (grHist[iI][iR][iPhi][iV]->GetN()>0) isFirst = false;
	}

	for (unsigned iPhi(0); iPhi<nPhi; ++iPhi){
	  if (!grFit[iI][iR][iPhi][iV]) continue;
	  grFit[iI][iR][iPhi][iV]->SetMarkerSize(1.5);
	  grFit[iI][iR][iPhi][iV]->SetMarkerStyle(20+iPhi);
	  grFit[iI][iR][iPhi][iV]->SetMarkerColor(colour[iPhi]);
	  grFit[iI][iR][iPhi][iV]->SetLineColor(colour[iPhi]);
	  grFit[iI][iR][iPhi][iV]->SetLineWidth(2);
	  grFit[iI][iR][iPhi][iV]->SetMaximum(maxY*1.2);
	  grFit[iI][iR][iPhi][iV]->SetMinimum(doMean!=1?0:minY*0.9);
	  grFit[iI][iR][iPhi][iV]->Draw(isFirst ? "APLE":"PLEsame");
	  std::ostringstream phiStr;
	  phiStr << phi[iPhi];
	  leg->AddEntry(grFit[iI][iR][iPhi][iV],(phiStr.str()+" fit").c_str(),"P");
	  if (grFit[iI][iR][iPhi][iV]->GetN()>0) isFirst = false;
	}
	
	leg->Draw("same");
	
	sprintf(buf,"%s %s",regShort[iReg].c_str(),label[iV].c_str());
	lat.SetTextColor(1);
	lat.SetTextSize(0.05);
	lat.DrawLatexNDC(0.12,0.95,buf);	
	mycFit->Update();
	
	std::ostringstream lSave;
	if (doMean==0) lSave << "PlotsReso/";
	else if (doMean==1) lSave << "PlotsMean/";
	else lSave << "PlotsSigma/";	
	//lSave << label[iV] << "_nRH" << nH[iR] << "_eta" << phi[iPhi];
	lSave << regShort[iReg] << "_"
	      << label[iV] << "_allphi";
	if (dovsE) lSave << "_vsE";
	mycFit->Print((lSave.str()+".pdf").c_str());
	mycFit->Print((lSave.str()+".png").c_str());
	
      }//loop on iter
    }//loop on nrec
  }//loop on vars
  
  return 0;
}//main
