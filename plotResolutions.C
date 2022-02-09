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
	      << label
	      << std::endl;
  };

};


void readResolution(std::string aFileName,
		    double ptval, double etaval,
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

void MeanvspT(const std::vector<MyReso> & vec,
	      std::string label,
	      const double etaval,
	      TGraphErrors* & gr,
	      bool useFit,
	      unsigned doMean,
	      bool dovsE){
  
  int nP = 0;
  for (unsigned iE(0); iE<vec.size(); ++iE){//loop on vec elements
    MyReso lRes = vec[iE];
    lRes.E = lRes.pt*cosh(lRes.eta/10.);
    if (lRes.label != label || fabs(lRes.eta-etaval)>0.01) continue;
    //lRes.Print();
    if (useFit) {
      if (doMean==0){
	gr->SetPoint(nP,dovsE?lRes.E:lRes.pt,lRes.resFit);
	gr->SetPointError(nP,0,lRes.fitPar2err/lRes.fitPar2*lRes.resFit);
      }
      else if (doMean==1){
	gr->SetPoint(nP,dovsE?lRes.E:lRes.pt,lRes.fitPar1);
	gr->SetPointError(nP,0,lRes.fitPar1err);
      }
      else if (doMean==2){
	gr->SetPoint(nP,dovsE?lRes.E:lRes.pt,lRes.fitPar2);
	gr->SetPointError(nP,0,lRes.fitPar2err);
      }
    }
    else {
      if (doMean==0){
	gr->SetPoint(nP,dovsE?lRes.E:lRes.pt,lRes.resHist);
	gr->SetPointError(nP,0,lRes.rmserr/lRes.rms*lRes.resHist);
      }
      else if (doMean==1){
	gr->SetPoint(nP,dovsE?lRes.E:lRes.pt,lRes.mean);
	gr->SetPointError(nP,0,lRes.meanerr);
      }
      else if (doMean==2){
	gr->SetPoint(nP,dovsE?lRes.E:lRes.pt,lRes.rms);
	gr->SetPointError(nP,0,lRes.rmserr);
      }
    }
    nP++;
  }//loop on vec elements
  
  if (dovsE) {
    gr->SetTitle(";E (GeV);#sigma/mean");
    if (doMean==1) gr->SetTitle(";E (GeV);mean");
    if (doMean==2) gr->SetTitle(";E (GeV);#sigma");
  }
  else {
    gr->SetTitle(";p_{T} (GeV);#sigma/mean");
    if (doMean==1) gr->SetTitle(";p_{T} (GeV);mean");
    if (doMean==2) gr->SetTitle(";p_{T} (GeV);#sigma");
  }
};

void plotResoFitResults(TCanvas* & myc, TGraphErrors* &gr,
			std::string grStr,
			const unsigned iI,
			std::string iter,
			const bool isFirst,
			const int nH, std::string label,
			double & minY,
			double & maxY){

  if (!gr) return;

  myc->cd();
  std::ostringstream lName;
  lName.str("");
  lName << "gr" << grStr << "_" << iter << "_" << nH << "_" << label;
  gr->SetName(lName.str().c_str());

  gr->SetTitle((";#eta;"+grStr).c_str());
  
  if (gr->GetN()>0){
    double lmax = TMath::MaxElement(gr->GetN(),gr->GetY());
    double lmin = TMath::MinElement(gr->GetN(),gr->GetY());
    if (lmax>maxY) maxY = lmax;
    if (lmin<minY) minY = lmin;
  }
  gr->SetMarkerStyle(20+3*iI+nH-1);
  gr->SetMarkerSize(1.5);
  gr->SetMarkerColor(kViolet-(3*iI+nH-1));
  gr->SetLineColor(kViolet-(3*iI+nH-1));
  gr->SetLineWidth(2);
  gr->Draw(isFirst ? "APE":"PEsame");
}


int plotResolutions(){

  bool binInPhi = true;

  SetTdrStyle();

  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);

  //0 = reso, 1 = mean, 2 = sigma
  const unsigned doMean = 0;

  const unsigned iReg = 0;
  const bool dovsE = doMean==0 ? true : false;  
  const unsigned nIter = 1;
  const unsigned nI = 3;
  std::string iter[4] = {"Dummy1","Dummy3","EM3","EM"};

  const unsigned nRec = 1;
  unsigned nH[3] = {1,2,3};

  const unsigned nEta = iReg==0 ? 6 : 1;
  double etaval[6] = {17,19,21,23,25,27};
  if (iReg!=0) etaval[0] = 21;
  double eta[6] = {1.7,1.9,2.1,2.3,2.5,2.7};
  if (iReg!=0) eta[0] = 2.1;
  double etaerr[6] = {0,0,0,0,0,0};
  const unsigned nPT = 12;
  double ptval[12] = {3,5,10,15,20,30,40,50,75,100,150,200};
  const unsigned nV = 10;
  const unsigned nV_allLC = 6;//first 3 - independent of iter
  const unsigned nV_TS = 4;//last 3 -- independent of nRec
  std::string label[nV] = {"EtotOverECP_1RH","EtotOverECP_2RH","EtotOverECP_3RH",
			   "E1LCOverECP_1RH","E1LCOverECP_2RH","E1LCOverECP_3RH",
			   "ETSOverECP","E1LCTSOverECP","ETSmaxEOverECP","SumETSOverECP"
  };

  const unsigned nPhi = binInPhi ? 3 : 1;
  const unsigned step = 20;
  const unsigned start = 0;
  unsigned phiMin[nPhi];
  unsigned phiMax[nPhi];
  if (!binInPhi) {
    phiMin[0] = 0;
    phiMax[0] = 360;
  } else {
    for (unsigned iPhi(0); iPhi<nPhi; ++iPhi){
      phiMin[iPhi] = start + step*iPhi;
      phiMax[iPhi] = phiMin[iPhi] + step;
    }
  }
  
  for (unsigned iPhi(0); iPhi<nPhi; ++iPhi){//loop on phi
    std::ostringstream phiStr;
    phiStr << "_phimod" << phiMin[iPhi] << "_" << phiMax[iPhi];

  TGraphErrors *grFit[nIter][nRec][nEta][nV];
  TGraphErrors *grHist[nIter][nRec][nEta][nV];

  for (unsigned iI(0); iI<nIter; ++iI){
    for (unsigned iR(0); iR<nRec; ++iR){
      for (unsigned iV(0); iV<nV; ++iV){//loop on variables
	for (unsigned ieta(0); ieta<nEta; ++ieta){
	  grFit[iI][iR][ieta][iV] = 0;
	  grHist[iI][iR][ieta][iV] = 0;
	}
      }
    }
  }
  
  for (unsigned iI(0); iI<nIter; ++iI){
    for (unsigned iR(0); iR<nRec; ++iR){
      std::cout << " - " << iter[iI] << std::endl;
      //<< " nRecHits=" << nH[iR] << std::endl;
      
      std::vector<MyReso> vec;
      vec.clear();
      for (unsigned ipt(0); ipt<nPT; ++ipt){  
	for (unsigned ieta(0); ieta<nEta; ++ieta){
	  std::cout << " --- Processing pT " << ptval[ipt] << " eta " << etaval[ieta]<< std::endl;
	  
	  std::ostringstream pteta;
	  pteta << "pt" << ptval[ipt] << "_eta" << etaval[ieta];

	  std::ostringstream filename;
	  //filename << "D49_All/" << iter[iI] << "/LCPlots/NREC" << nH[iR] << "/resolutions.dat";
	  filename << "D49_AllTracksters/Merged/" << iter[iI] << "/" << pteta.str() << phiStr.str() << "/resolutions_" << regShort[iReg] << ".dat";
	  
	  readResolution(filename.str().c_str(),ptval[ipt],etaval[ieta],vec);
	}
      }


      if (vec.size()==0) {
	std::cout << " -- Vec has no elements" << std::endl;
	continue;
      }
      
      for (unsigned iV(0); iV<nV; ++iV){//loop on variables
	std::cout << " - Processing " << label[iV] << std::endl;
	for (unsigned ieta(0); ieta<nEta; ++ieta){
	  std::cout << " -- Processing " << label[iV] << " eta " << etaval[ieta] << std::endl;
	  
	  grFit[iI][iR][ieta][iV] = new TGraphErrors();
	  std::ostringstream lName;
	  lName << "grFit_" << iter[iI] << "_" << nH[iR] << "_" << label[iV] << "_eta" << etaval[ieta];
	  grFit[iI][iR][ieta][iV]->SetName(lName.str().c_str());
	  MeanvspT(vec,label[iV],etaval[ieta],grFit[iI][iR][ieta][iV],true,doMean,dovsE);
	  grHist[iI][iR][ieta][iV] = new TGraphErrors();
	  lName.str("");
	  lName << "grHist_" << iter[iI] << "_" << nH[iR] << "_" << label[iV] << "_eta" << etaval[ieta];
	  grHist[iI][iR][ieta][iV]->SetName(lName.str().c_str());
	  MeanvspT(vec,label[iV],etaval[ieta],grHist[iI][iR][ieta][iV],false,doMean,dovsE);
	  std::cout << " - Npoints: " << " " << grFit[iI][iR][ieta][iV]->GetN()
		    << " " << grHist[iI][iR][ieta][iV]->GetN() << " "
	  << std::endl;
	}
      }
    }
  }
  
  TCanvas *mycReso  = new TCanvas("mycReso","mycReso",1);
  TCanvas *mycN  = new TCanvas("mycN","mycN",1);
  TCanvas *mycC  = new TCanvas("mycC","mycC",1);
  TCanvas *mycS  = new TCanvas("mycS","mycS",1);

  TCanvas *mycFit  = new TCanvas("mycFit","mycFit",1);
  mycFit->cd();
  
  TLatex lat;
  char buf[200];

  double sigmaStoch[nV][nRec][nI][nEta];
  double sigmaStochErr[nV][nRec][nI][nEta];
  double sigmaConst[nV][nRec][nI][nEta];
  double sigmaConstErr[nV][nRec][nI][nEta];
  double sigmaNoise[nV][nRec][nI][nEta];
  double sigmaNoiseErr[nV][nRec][nI][nEta];

  for (unsigned iV(0); iV<nV; ++iV){//loop on variables
    for (unsigned iR(0); iR<nRec; ++iR){
      for (unsigned iI(0); iI<nI; ++iI){
	for (unsigned ieta(0); ieta<nEta; ++ieta){
	  sigmaStoch[iV][iR][iI][ieta] = 0;
	  sigmaStochErr[iV][iR][iI][ieta] = 0;
	  sigmaConst[iV][iR][iI][ieta] = 0;
	  sigmaConstErr[iV][iR][iI][ieta] = 0;
	  sigmaNoise[iV][iR][iI][ieta] = 0;
	  sigmaNoiseErr[iV][iR][iI][ieta] = 0;
	}
      }
    }
  }
  
  for (unsigned iV(0); iV<nV; ++iV){//loop on variables
    //if (iV<(nV-nV_TS)) continue;
    for (unsigned ieta(0); ieta<nEta; ++ieta){
      std::cout << " -- Processing " << label[iV] << " eta " << etaval[ieta] << std::endl;
      for (unsigned iR(0); iR<nRec; ++iR){
	mycFit->Clear();
	mycFit->cd();
	gPad->SetGridy(1);
	gPad->SetLogx(1);
	bool isFirst = true;
	TLegend *leg = new TLegend(0.15,0.8,0.94,0.92);
	leg->SetFillStyle(0);
	leg->SetNColumns(nIter);
	double maxY = 0;
	double minY = 10;
	//std::cout << " - nRecHits=" << nH[iR] << std::endl;
	for (unsigned iI(0); iI<nIter; ++iI){
	  if (iV<nV_allLC && iI>0) continue;
	  //std::cout << " - " << iter[iI] << " " << grFit[iI][iR][ieta][iV]
	  //<< " " << grHist[iI][iR][ieta][iV] << " "
	  //<< std::endl;
	  if (!grFit[iI][iR][ieta][iV]) continue;
	  if (!grHist[iI][iR][ieta][iV]) continue;
	  std::cout << " -- Npoint = " << grHist[iI][iR][ieta][iV]->GetN() << std::endl;
	  if (grHist[iI][iR][ieta][iV]->GetN()==0) {
	    continue;
	  }
	  double lmax = TMath::MaxElement(grHist[iI][iR][ieta][iV]->GetN(),grHist[iI][iR][ieta][iV]->GetY());
	  double lmin = TMath::MinElement(grHist[iI][iR][ieta][iV]->GetN(),grHist[iI][iR][ieta][iV]->GetY());
	  if (lmax>maxY) maxY = lmax;
	  if (lmin<minY) minY = lmin;
	  std::cout << " - " << iter[iI] << " " << lmin << " " << lmax << std::endl;
	}
	for (unsigned iI(0); iI<nIter; ++iI){
	  if (!grHist[iI][iR][ieta][iV]) continue;
	  if (iV<nV_allLC && iI>0) continue;
	  grHist[iI][iR][ieta][iV]->SetMarkerSize(1.5);
	  grHist[iI][iR][ieta][iV]->SetMarkerStyle(20+iI);
	  grHist[iI][iR][ieta][iV]->SetMarkerColor(iV<nV_allLC? kGray+2:kGray+iI);
	  grHist[iI][iR][ieta][iV]->SetLineColor(iV<nV_allLC? kGray+2:kGray+iI);
	  grHist[iI][iR][ieta][iV]->SetMaximum(maxY*1.1);
	  grHist[iI][iR][ieta][iV]->SetMinimum(doMean!=1?0:minY*0.9);
	  grHist[iI][iR][ieta][iV]->Draw(isFirst ? "APE":"PEsame");
	  if (iV>nV_allLC) leg->AddEntry(grHist[iI][iR][ieta][iV],(iter[iI]+" hist").c_str(),"P");
	  else leg->AddEntry(grHist[iI][iR][ieta][iV],"Hist","P");
	  if (grHist[iI][iR][ieta][iV]->GetN()>0) isFirst = false;
	}

	for (unsigned iI(0); iI<nIter; ++iI){
	  if (!grFit[iI][iR][ieta][iV]) continue;
	  if (iV<nV_allLC && iI>0) continue;
	  grFit[iI][iR][ieta][iV]->SetMarkerSize(1.5);
	  grFit[iI][iR][ieta][iV]->SetMarkerStyle(20+iI);
	  grFit[iI][iR][ieta][iV]->SetMarkerColor(kViolet-iI);
	  grFit[iI][iR][ieta][iV]->SetLineColor(kViolet-iI);
	  grFit[iI][iR][ieta][iV]->SetLineWidth(2);
	  grFit[iI][iR][ieta][iV]->SetMaximum(maxY*1.1);
	  grFit[iI][iR][ieta][iV]->SetMinimum(doMean!=1?0:minY*0.9);
	  grFit[iI][iR][ieta][iV]->Draw("PEsame");
	  if (iV>=nV_allLC) leg->AddEntry(grFit[iI][iR][ieta][iV],(iter[iI]+" fit").c_str(),"P");
	  else leg->AddEntry(grFit[iI][iR][ieta][iV],"Fit","P");
	}

	
	leg->Draw("same");
	
	//if (iV<(nV-nV_TS)) sprintf(buf,"%s eta=%3.1f nRecH=%d",label[iV].c_str(),etaval[ieta]/10.,nH[iR]);
	//else
	sprintf(buf,"%s %s eta=%3.1f",regShort[iReg].c_str(),label[iV].c_str(),etaval[ieta]/10.);
	lat.SetTextColor(1);
	lat.SetTextSize(0.05);
	lat.DrawLatexNDC(0.12,0.95,buf);	
	mycFit->Update();

	std::ostringstream lSave;
	if (doMean==0) lSave << "PlotsReso/";
	else if (doMean==1) lSave << "PlotsMean/";
	else lSave << "PlotsSigma/";	
	//lSave << label[iV] << "_nRH" << nH[iR] << "_eta" << etaval[ieta];
	lSave << regShort[iReg] << phiStr.str() << "_"
	      << label[iV] << "_eta" << etaval[ieta];
	if (dovsE) lSave << "_vsE";
	mycFit->Print((lSave.str()+".pdf").c_str());
	mycFit->Print((lSave.str()+".png").c_str());

	if (doMean==0){
	  for (unsigned iI(0); iI<nIter; ++iI){
	    if (iV<nV_allLC && iI>0) continue;
	    if (!grFit[iI][iR][ieta][iV]) continue;
	    TPad *lpad = (TPad*)(mycReso->cd());
	    plotResolution(grFit[iI][iR][ieta][iV],lpad,0,etaval[ieta],
			   dovsE?0.23:0.14,0.01,0.2,
			   sigmaStoch[iV][iR][iI][ieta],
			   sigmaStochErr[iV][iR][iI][ieta],
			   sigmaConst[iV][iR][iI][ieta],
			   sigmaConstErr[iV][iR][iI][ieta],
			   sigmaNoise[iV][iR][iI][ieta],
			   sigmaNoiseErr[iV][iR][iI][ieta],
			   dovsE,false);
	    //if (iV<(nV-nV_TS)) sprintf(buf,"%s eta=%3.1f nRecH=%d",label[iV].c_str(),etaval[ieta]/10.,nH[iR]);
	    //else
	    sprintf(buf,"%s %s eta=%3.1f",regShort[iReg].c_str(),label[iV].c_str(),etaval[ieta]/10.);
	    lat.SetTextSize(0.05);
	    lat.DrawLatexNDC(0.12,0.95,buf);	
	    sprintf(buf,"%s",iter[iI].c_str());
	    if (iV>=nV_allLC)
	      lat.DrawLatexNDC(0.15,0.85,buf);	
	    mycReso->Update();

	    std::ostringstream lSave;
	    lSave << "PlotsReso/"
		  << regShort[iReg] << phiStr.str() << "_"
		  <<"ResoFit_";
	    //lSave << iter[iI] << "_" << label[iV] << "_nRH" << nH[iR] << "_eta" << etaval[ieta];
	    lSave << iter[iI] << "_" << label[iV] << "_eta" << etaval[ieta];
	    if (dovsE) lSave << "_vsE";
	    mycReso->Print((lSave.str()+".pdf").c_str());
	    mycReso->Print((lSave.str()+".png").c_str());
	    
	  }
	}
      }//loop on nRec
    }//loop on ieta
  }//loop on variables

  for (unsigned iV(0); iV<nV; ++iV){//loop on variables

    std::cout << " -- Processing " << label[iV] << " " << phiStr.str() << std::endl;
    if (doMean==0){
      
      //plot versus eta
      if (nEta == 1) continue;
      
      TGraphErrors *grN[nRec][nI];
      TGraphErrors *grC[nRec][nI];
      TGraphErrors *grS[nRec][nI];
      TLegend *leg = new TLegend(0.15,(iV<nV_allLC || iV>=(nV-nV_TS))?0.85:0.75,0.94,0.92);
      leg->SetFillStyle(0);
      std::cout << " -- Plots vs eta: " << label[iV] << std::endl;
      
      bool isFirst = true;
      if (iV==0){
	leg->SetNColumns(3);
	mycN->Clear();
	mycC->Clear();
	mycS->Clear();
	double maxY[3] = {0,0,0};
	double minY[3] = {10,10,10};
	for (unsigned iR(0); iR<1; ++iR){
	  for (unsigned iI(0); iI<3; ++iI){
	    grN[iR][iI] = new TGraphErrors(nEta,eta,sigmaNoise[iV+iI][iR][0],etaerr,sigmaNoiseErr[iV+iI][iR][0]);
	    grC[iR][iI] = new TGraphErrors(nEta,eta,sigmaConst[iV+iI][iR][0],etaerr,sigmaConstErr[iV+iI][iR][0]);
	    grS[iR][iI] = new TGraphErrors(nEta,eta,sigmaStoch[iV+iI][iR][0],etaerr,sigmaStochErr[iV+iI][iR][0]);
	    
	    mycN->cd();
	    std::cout << "mycN " << mycN
		      << " grN " << grN[iR][iI]
		      << " iI " << iI
		      << " iter " << iter[0]
		      << " isFirst " << isFirst
		      << " nH " << nH[iI]
		      << " label " << label[iV+iI]
		      << " minY " << minY[0]
		      << " maxY " << maxY[0]
		      << std::endl;
	    plotResoFitResults(mycN,grN[iR][iI],"n",iI,iter[0],isFirst,nH[iI],label[iV+iI],minY[0],maxY[0]);
	    mycC->cd();
	    plotResoFitResults(mycC,grC[iR][iI],"c",iI,iter[0],isFirst,nH[iI],label[iV+iI],minY[1],maxY[1]);
	    mycS->cd();
	    plotResoFitResults(mycS,grS[iR][iI],"s",iI,iter[0],isFirst,nH[iI],label[iV+iI],minY[2],maxY[2]);
	    std::ostringstream lLegLabel;
	    //lLegLabel << iter[iI] << " ";
	    lLegLabel << "nRH = " << nH[iI];
	    leg->AddEntry(grS[iR][iI],lLegLabel.str().c_str(),"P");

	    std::cout << " Test 1 " << iI << " N "
		      << grN[iR][iI]->GetN() << " S "
		      << grS[iR][iI]->GetN() << " " << sigmaStoch[iV+iI][iR][0][0]
		      << " max " << maxY[2] << " min " << minY[2]
		      << " C "
		      << grC[iR][iI]->GetN()<< " " << sigmaConst[iV+iI][iR][0][0]
		      << " max " << maxY[1] << " min " << minY[1]
		      << std::endl;	    
	    isFirst = false;
	    if (iR==0 && iI==2){
	      mycN->cd();
	      grN[0][0]->SetMaximum(maxY[0]*1.2);
	      grN[0][0]->SetMinimum(minY[0]*0.8);
	      gPad->Modified();
	      mycC->cd();
	      grC[0][0]->SetMaximum(0.018);//maxY[1]*1.2);
	      grC[0][0]->SetMinimum(0.007);//minY[1]*0.8);
	      gPad->Modified();
	      mycS->cd();
	      grS[0][0]->SetMaximum(maxY[2]*1.2);
	      grS[0][0]->SetMinimum(minY[2]*0.8);
	      gPad->Modified();
	    }
	  }//loop on rechits 1,2,3
	}
      }//iV==0
      else if (iV==3){
	leg->SetNColumns(3);
	mycN->Clear();
	mycC->Clear();
	mycS->Clear();
	double maxY[3] = {0,0,0};
	double minY[3] = {10,10,10};
	for (unsigned iR(0); iR<1; ++iR){
	  for (unsigned iI(0); iI<3; ++iI){
	    grN[iR][iI] = new TGraphErrors(nEta,eta,sigmaNoise[iV+iI][iR][0],etaerr,sigmaNoiseErr[iV+iI][iR][0]);
	    grC[iR][iI] = new TGraphErrors(nEta,eta,sigmaConst[iV+iI][iR][0],etaerr,sigmaConstErr[iV+iI][iR][0]);
	    grS[iR][iI] = new TGraphErrors(nEta,eta,sigmaStoch[iV+iI][iR][0],etaerr,sigmaStochErr[iV+iI][iR][0]);
	    
	    mycN->cd();
	    plotResoFitResults(mycN,grN[iR][iI],"n",iI,iter[0],isFirst,nH[iI],label[iV+iI],minY[0],maxY[0]);
	    mycC->cd();
	    plotResoFitResults(mycC,grC[iR][iI],"c",iI,iter[0],isFirst,nH[iI],label[iV+iI],minY[1],maxY[1]);
	    mycS->cd();
	    plotResoFitResults(mycS,grS[iR][iI],"s",iI,iter[0],isFirst,nH[iI],label[iV+iI],minY[2],maxY[2]);
	    std::ostringstream lLegLabel;
	    //lLegLabel << iter[iI] << " ";
	    lLegLabel << "nRH = " << nH[iI];
	    leg->AddEntry(grS[iR][iI],lLegLabel.str().c_str(),"P");

	    std::cout << " Test 2 " << iI << " N "
		      << grN[iR][iI]->GetN() << " S "
		      << grS[iR][iI]->GetN() << " " << sigmaStoch[iV+iI][iR][0][0]
		      << " max " << maxY[2] << " min " << minY[2]
		      << " C "
		      << grC[iR][iI]->GetN()<< " " << sigmaConst[iV+iI][iR][0][0]
		      << " max " << maxY[1] << " min " << minY[1]
		      << std::endl;	    
	    
	    isFirst = false;
	    if (iR==0 && iI==2){
	      mycN->cd();
	      grN[0][0]->SetMaximum(maxY[0]*1.2);
	      grN[0][0]->SetMinimum(minY[0]*0.8);
	      gPad->Modified();
	      mycC->cd();
	      grC[0][0]->SetMaximum(0.018);//maxY[1]*1.2);
	      grC[0][0]->SetMinimum(0.007);//minY[1]*0.8);
	      gPad->Modified();
	      mycS->cd();
	      grS[0][0]->SetMaximum(maxY[2]*1.2);
	      grS[0][0]->SetMinimum(minY[2]*0.8);
	      gPad->Modified();
	    }
	  }
	}
      }//iV 3
      else if (iV>=nV_allLC) {//TS variables
	leg->SetNColumns(nIter);
	mycN->Clear();
	mycC->Clear();
	mycS->Clear();
	double maxY[3] = {0,0,0};
	double minY[3] = {10,10,10};
	
	for (unsigned iR(0); iR<nRec; ++iR){
	  for (unsigned iI(0); iI<nIter; ++iI){
	    
	    grN[iR][iI] = new TGraphErrors(nEta,eta,sigmaNoise[iV][iR][iI],etaerr,sigmaNoiseErr[iV][iR][iI]);
	    grC[iR][iI] = new TGraphErrors(nEta,eta,sigmaConst[iV][iR][iI],etaerr,sigmaConstErr[iV][iR][iI]);
	    grS[iR][iI] = new TGraphErrors(nEta,eta,sigmaStoch[iV][iR][iI],etaerr,sigmaStochErr[iV][iR][iI]);
	    
	    mycN->cd();
	    plotResoFitResults(mycN,grN[iR][iI],"n",iI,iter[iI],isFirst,nH[iR],label[iV],minY[0],maxY[0]);
	    mycC->cd();
	    plotResoFitResults(mycC,grC[iR][iI],"c",iI,iter[iI],isFirst,nH[iR],label[iV],minY[1],maxY[1]);
	    mycS->cd();
	    plotResoFitResults(mycS,grS[iR][iI],"s",iI,iter[iI],isFirst,nH[iR],label[iV],minY[2],maxY[2]);
	    std::ostringstream lLegLabel;
	    lLegLabel << iter[iI] << " ";
	    //if (iV<(nV-nV_TS)) lLegLabel << "nRH = " << nH[iR];
	    leg->AddEntry(grS[iR][iI],lLegLabel.str().c_str(),"P");
	    
	    isFirst = false;
	    if (iR==nRec-1 && iI==nIter-1){
	      mycN->cd();
	      grN[0][0]->SetMaximum(maxY[0]*1.2);
	      grN[0][0]->SetMinimum(minY[0]*0.8);
	      gPad->Modified();
	      mycC->cd();
	      grC[0][0]->SetMaximum(0.018);//maxY[1]*1.2);
	      grC[0][0]->SetMinimum(0.007);//minY[1]*0.8);
	      gPad->Modified();
	      mycS->cd();
	      grS[0][0]->SetMaximum(maxY[2]*1.2);
	      grS[0][0]->SetMinimum(minY[2]*0.8);
	      gPad->Modified();
	    }
	  }//loop on iter
	}//loop on RH
      }//TS variables
      else continue;

      std::cout << " -- Now plotting" << std::endl;
      
      mycN->cd();
      gPad->SetGridy(1);
      leg->Draw("same");
      sprintf(buf,"%s %s",regShort[iReg].c_str(),label[iV].c_str());
      lat.SetTextSize(0.05);
      lat.DrawLatexNDC(0.12,0.95,buf);	
      lat.DrawLatexNDC(0.6,0.7,phiStr.str().c_str());	
      mycN->Update();
      std::ostringstream lSave;
      lSave << "PlotsReso/"
	    << regShort[iReg] << phiStr.str() << "_"
	    <<"FitNoisevsEta_";
      lSave << label[iV];// << "_nRH" << nH[iR];
      mycN->Print((lSave.str()+".pdf").c_str());
      mycN->Print((lSave.str()+".png").c_str());
      
      mycC->cd();
      gPad->SetGridy(1);
      leg->Draw("same");
      lat.DrawLatexNDC(0.12,0.95,buf);	
      lat.DrawLatexNDC(0.6,0.7,phiStr.str().c_str());	
      mycC->Update();
      lSave.str("");
      lSave << "PlotsReso/"
	    << regShort[iReg] << phiStr.str() << "_"
	    << "FitConstvsEta_";
      lSave << label[iV];// << "_nRH" << nH[iR];
      mycC->Print((lSave.str()+".pdf").c_str());
      mycC->Print((lSave.str()+".png").c_str());
      
      mycS->cd();
      gPad->SetGridy(1);
      leg->Draw("same");
      lat.DrawLatexNDC(0.12,0.95,buf);	
      lat.DrawLatexNDC(0.6,0.7,phiStr.str().c_str());	
      mycS->Update();
      lSave.str("");
      lSave << "PlotsReso/"
	    << regShort[iReg] << phiStr.str() << "_"
	    << "FitStochvsEta_";
      lSave << label[iV];// << "_nRH" << nH[iR];
      mycS->Print((lSave.str()+".pdf").c_str());
      mycS->Print((lSave.str()+".png").c_str());
      
    }//doreso
  }//loop on iV

  }//loop on phi bins
  
  return 0;
}//main
