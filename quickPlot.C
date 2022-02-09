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

#include "globalVars.h"

struct Outlier {
  unsigned event;
  unsigned entry;
  unsigned nSC10;
  unsigned nTS10;
  unsigned nSC;
  unsigned nTS;
  unsigned nBlobs;
  double eBlobs;
  double eTotAllLC3;
  double eFrac[10];
  double eFracAtB[10];
  void Print(){
    std::cout << event << " " << entry << " "
	      << nSC10 << " " << nTS10 << " "
	      << nSC << " " << nTS << " "
	      << nBlobs << " "
	      << eBlobs << " " << eTotAllLC3 << " ";
    for (unsigned iC(0); iC<10;++iC){
      std::cout << eFrac[iC] << " "
		<< eFracAtB[iC] << " ";
    }
 
    
    std::cout << std::endl;
  };
};

void getOutliers(std::string aFileName,
		 std::vector<Outlier> & vec){
  
  std::ifstream infile(aFileName);
  vec.clear();
  
  if (!infile.is_open()){
    std::cout << " -- file " << aFileName << " not found." << std::endl;
    return;
  }
  else {
    std::cout << " -- Reading file " << aFileName << "." << std::endl;
  }

  int counter =0 ;

  while(!infile.eof()){
    Outlier lOut;
    lOut.nSC=0;
    infile>>lOut.event>>lOut.entry>>lOut.nSC10>>lOut.nTS10>>lOut.nSC>>lOut.nTS>>lOut.nBlobs>>lOut.eBlobs>>lOut.eTotAllLC3;
    for (unsigned iC(0); iC<10;++iC){
      infile>>lOut.eFrac[iC]>>lOut.eFracAtB[iC];
    }
    if (lOut.nSC>0) vec.push_back(lOut);
    counter++;
    if (counter>10000) break;
  }
  std::cout << " -- Found " << vec.size() << " (" << counter << ") elements." << std::endl;
  
};

void getTrueRZ(TChain *treeIn,
	       TH2F* & hRZtrue,
	       const unsigned entry,
	       const bool isFC,
	       double & truthE,
	       const bool doVsLayer=false,
	       const double deta = 0){
  
  TH1F *hT = new TH1F("hT","",100,0,500);
  treeIn->Draw("cp_energy[0]>>hT","","",1,entry);
  truthE = hT->GetMean();
  hT->Delete();
  
  hT = new TH1F("hT","",100,2.05,2.15);
  treeIn->Draw("cp_eta[0]>>hT","","",1,entry);
  double cpEta = hT->GetMean();
  
  std::ostringstream lname;
  lname << "hRZtrue_" << entry;
  if (isFC) lname << "_FC";
  if (deta > 0) lname << "_detap" << static_cast<unsigned>(abs(deta)*100);
  else if  (deta<0) lname << "_detam" << static_cast<unsigned>(abs(deta)*100);
  if (!doVsLayer) hRZtrue = new TH2F(lname.str().c_str(),";r (cm);z (cm);E_{LC} (GeV)",100,0,200,380,320,510);
  else hRZtrue = new TH2F(lname.str().c_str(),";r (cm);layer;E_{LC} (GeV)",100,0,200,52,0.5,52.5);
  std::cout << " making " << lname.str() << std::endl;
  for (unsigned iL(0); iL<49; ++iL){//loop on layers
    double trueR = sqrt(pow(getZpos(iL)/cos(2*atan(exp(-(cpEta+deta)))),2)-pow(getZpos(iL),2));
    hRZtrue->Fill(trueR,doVsLayer?iL+1:getZpos(iL),1000);
  }
  hT->Delete();
};

void getRZ(TChain *treeIn,
	   TH2F* & hRZ,
	   const unsigned entry,
	   const bool isFC,
	   const bool doAllLC,
	   const bool doSeed,
	   const bool dRmatch,
	   const std::string friendName="",
	   const bool doVsLayer=false,
	   const bool inBlobs=false){
  
  std::ostringstream lname;
  lname << "hRZ_" << entry;
  if (isFC) lname << "_FC";
  if (doAllLC) lname << "_allLC";
  if (friendName.find("PU")!=friendName.npos) lname << "_PU";
  if (inBlobs) lname << "_inBlobs";
  std::string nRHcut;
  if (inBlobs || (friendName.find("PU")!=friendName.npos && dRmatch) ) nRHcut = "3";
  else nRHcut = "1";
  if (!doVsLayer){
    if (doAllLC) hRZ = new TH2F(lname.str().c_str(),(";r (cm);z (cm);E_{LC}^{Si#geq"+nRHcut+"RH} (GeV)").c_str(),100,0,200,380,320,510);
    else hRZ = new TH2F(lname.str().c_str(),";r (cm);z (cm);E_{LC}^{Si#geq3RH} (GeV)",100,0,200,380,320,510);
  } else {
    if (doAllLC) hRZ = new TH2F(lname.str().c_str(),(";r (cm);layer;E_{LC}^{Si#geq"+nRHcut+"RH} (GeV)").c_str(),100,0,200,52,0.5,52.5);
    else hRZ = new TH2F(lname.str().c_str(),";r (cm);layer;E_{LC}^{Si#geq3RH} (GeV)",100,0,200,52,0.5,52.5);
  }
  std::ostringstream ldRmatch;
  
  
  if (doAllLC){
    if (dRmatch) ldRmatch << " && sqrt(pow(TMath::Abs(" << friendName << "all_lc_seedEta-Sim.cp_eta[0]),2)+pow(deltaPhi(" << friendName << "all_lc_seedPhi,Sim.cp_phi[0]),2))<0.4";
    if (inBlobs) ldRmatch << " && TRUTH.lc_in_blob==1";
    if (!doSeed){
      if (doVsLayer) treeIn->Draw((friendName+"all_lc_layer:sqrt(pow("+friendName+"all_lc_y,2)+pow("+friendName+"all_lc_x,2))>>"+lname.str()).c_str(),(friendName+"all_lc_energy*(("+friendName+"all_lc_isSi==0 || ("+friendName+"all_lc_isSi>0 && "+friendName+"all_lc_nrechits>="+nRHcut+"))"+ldRmatch.str()+")").c_str(),"colz",1,entry);
      else treeIn->Draw((friendName+"all_lc_z:sqrt(pow("+friendName+"all_lc_y,2)+pow("+friendName+"all_lc_x,2))>>"+lname.str()).c_str(),(friendName+"all_lc_energy*(("+friendName+"all_lc_isSi==0 || ("+friendName+"all_lc_isSi>0 && "+friendName+"all_lc_nrechits>="+nRHcut+"))"+ldRmatch.str()+")").c_str(),"colz",1,entry);
    }
    else {
      std::ostringstream lVar;
      lVar << "sqrt(pow(getZpos(" << friendName << "all_lc_layer-1)/cos(2*atan(exp(-" << friendName << "all_lc_seedEta))),2)-pow(getZpos(" << friendName << "all_lc_layer-1),2))";
      if (doVsLayer) treeIn->Draw((friendName+"all_lc_layer:"+lVar.str()+">>"+lname.str()).c_str(),(friendName+"all_lc_energy*(("+friendName+"all_lc_isSi==0 || ("+friendName+"all_lc_isSi>0 && "+friendName+"all_lc_nrechits>="+nRHcut+"))"+ldRmatch.str()+")").c_str(),"colz",1,entry);
      else treeIn->Draw((friendName+"all_lc_z:"+lVar.str()+">>"+lname.str()).c_str(),(friendName+"all_lc_energy*(("+friendName+"all_lc_isSi==0 || ("+friendName+"all_lc_isSi>0 && "+friendName+"all_lc_nrechits>="+nRHcut+"))"+ldRmatch.str()+")").c_str(),"colz",1,entry);
    }
  } else {
    if (dRmatch) ldRmatch << " && sqrt(pow(TMath::Abs(" << friendName << "lc_seedEta-" << friendName << "cp_eta[0]),2)+pow(deltaPhi(" << friendName << "lc_seedPhi," << friendName << "cp_phi[0]),2))<0.4";
    if (!doSeed) {
      if (doVsLayer) treeIn->Draw((friendName+"lc_layer:sqrt(pow("+friendName+"lc_y,2)+pow("+friendName+"lc_x,2))>>"+lname.str()).c_str(),(friendName+"lc_energy/"+friendName+"lc_tsMult*(("+friendName+"lc_isSi==0 || ("+friendName+"lc_isSi>0 && "+friendName+"lc_nrechits>2))"+ldRmatch.str()+")").c_str(),"colz",1,entry);
      else treeIn->Draw((friendName+"lc_z:sqrt(pow("+friendName+"lc_y,2)+pow("+friendName+"lc_x,2))>>"+lname.str()).c_str(),(friendName+"lc_energy/"+friendName+"lc_tsMult*(("+friendName+"lc_isSi==0 || ("+friendName+"lc_isSi>0 && "+friendName+"lc_nrechits>2))"+ldRmatch.str()+")").c_str(),"colz",1,entry);
    }
    else {
      std::ostringstream lVar;
      lVar << "sqrt(pow(getZpos(" << friendName << "lc_layer-1)/cos(2*atan(exp(-" << friendName << "lc_seedEta))),2)-pow(getZpos(" << friendName << "lc_layer-1),2))";
      if (doVsLayer) treeIn->Draw((friendName+"lc_layer:"+lVar.str()+">>"+lname.str()).c_str(),(friendName+"lc_energy/"+friendName+"lc_tsMult*(("+friendName+"lc_isSi==0 || ("+friendName+"lc_isSi>0 && "+friendName+"lc_nrechits>2))"+ldRmatch.str()+")").c_str(),"colz",1,entry);
      else treeIn->Draw((friendName+"lc_z:"+lVar.str()+">>"+lname.str()).c_str(),(friendName+"lc_energy/"+friendName+"lc_tsMult*(("+friendName+"lc_isSi==0 || ("+friendName+"lc_isSi>0 && "+friendName+"lc_nrechits>2))"+ldRmatch.str()+")").c_str(),"colz",1,entry);
    }
  }
  
};

void getXZ(TChain *treeIn,
	   TH2F* & hRZ,
	   const unsigned entry,
	   const std::string lCut,
	   const std::string friendName=""
	   ){
  std::ostringstream lname;
  lname << "hXZ_" << entry;
  hRZ = new TH2F(lname.str().c_str(),";x (cm);z (cm);E_{LC}^{#geq1RH} (GeV)",100,-200,200,380,320,510);
  std::ostringstream lVar;
  lVar << "getZpos(" << friendName << "lc_layer-1)*tan(2*atan(exp(-" << friendName << "lc_seedEta)))*cos(" << friendName << "lc_seedPhi)";
  treeIn->Draw((friendName+"lc_z:"+lVar.str()+">>"+lname.str()).c_str(),(friendName+"lc_energy/"+friendName+"lc_tsMult*("+lCut+")").c_str(),"colz",1,entry);
};

void getYZ(TChain *treeIn,
	   TH2F* & hRZ,
	   const unsigned entry,
	   const std::string lCut,
	   const std::string friendName=""
	   ){

  std::ostringstream lname;
  lname << "hYZ_" << entry;
  hRZ = new TH2F(lname.str().c_str(),";y (cm);z (cm);E_{LC}^{#geq1RH} (GeV)",100,-200,200,380,320,510);
  std::ostringstream lVar;
  lVar << "getZpos(" << friendName << "lc_layer-1)*tan(2*atan(exp(-" << friendName << "lc_seedEta)))*sin(" << friendName << "lc_seedPhi)";
  treeIn->Draw((friendName+"lc_z:"+lVar.str()+">>"+lname.str()).c_str(),(friendName+"lc_energy/"+friendName+"lc_tsMult*("+lCut+")").c_str(),"colz",1,entry);
};

void getRZ(TChain *treeIn,
	   TH2F* & hRZ,
	   const unsigned entry,
	   const std::string lCut,
	   const bool isFC,
	   const std::string friendName="",
	   const bool doVsLayer=false
	   ){ 
  std::ostringstream lname;
  lname << "hRZsel_" << entry;
  if (isFC) lname << "_FC";
  if (!doVsLayer) hRZ = new TH2F(lname.str().c_str(),";r (cm);z (cm);E_{LC}^{Si#geq3RH} (GeV)",100,0,200,380,320,510);
  else hRZ = new TH2F(lname.str().c_str(),";r (cm);z (cm);E_{LC}^{Si#geq3RH} (GeV)",100,0,200,52,0.5,52.5);
  std::ostringstream lVar;
  lVar << "sqrt(pow(getZpos(" << friendName << "lc_layer-1)/cos(2*atan(exp(-" << friendName << "lc_seedEta))),2)-pow(getZpos(" << friendName << "lc_layer-1),2))";
  if (!doVsLayer) treeIn->Draw((friendName+"lc_z:"+lVar.str()+">>"+lname.str()).c_str(),(friendName+"lc_energy/"+friendName+"lc_tsMult*("+lCut+(lCut.size()>0?" &&":"")+"("+friendName+"lc_isSi==0 || ("+friendName+"lc_isSi>0 && "+friendName+"lc_nrechits>=3)))").c_str(),"colz",1,entry);
  else treeIn->Draw((friendName+"lc_layer:"+lVar.str()+">>"+lname.str()).c_str(),(friendName+"lc_energy/"+friendName+"lc_tsMult*("+lCut+(lCut.size()>0?" && ":"")+"("+friendName+"lc_isSi==0 || ("+friendName+"lc_isSi>0 && "+friendName+"lc_nrechits>=3)))").c_str(),"colz",1,entry);
};

int plotEvtDisplay(TChain *tree,
		   TChain *treeLC,
		   TCanvas* & myc,
		   TCanvas* & mycT,
		   TCanvas* & mycP,
		   TH2F* & hRZtrue,
		   TH2F* & hRZall,
		   TH2F* & hRZ,
		   TH2F* & hPID,
		   TH2F* & hE,
		   TH1F* & hCat,
		   const std::string lCut,
		   const std::string lCutTS,
		   const std::string lCutSave,
		   TH2F* & hRZs,
		   TH2F* & hXZ,
		   TH2F* & hYZ,
		   std::string plotDir,
		   char* buf,
		   const bool doAllLC,
		   const bool doSeed,
		   const Outlier & outlier,
		   const bool isPU
		   ){

  gStyle->SetMarkerSize(2);//1.2);

  const unsigned entry = outlier.event;
  bool isFC = plotDir.find("FineCalo")!=plotDir.npos;

  bool doCut = lCut.size()>0;

  double truthE = 0;
  myc->cd(1);
  getTrueRZ(tree,hRZtrue,entry,isFC,truthE);  
  getRZ(treeLC,hRZall,entry,isFC,true,doSeed,false);
  
  TLatex lat;
  
  if (!doCut){
    myc->cd(1);
    gPad->SetLogz(1);
    gPad->SetRightMargin(0.18);
    
    TExec *ex1 = new TExec("ex1","gStyle->SetPalette(1);");
    hRZall->GetZaxis()->SetRangeUser(0.01,100);
    hRZall->Draw("colz");
    ex1->Draw();
    hRZall->Draw("colzsame");
    
    hRZtrue->SetLineColor(2);
    hRZtrue->SetFillColor(2);
    hRZtrue->SetFillStyle(1001);
    hRZtrue->Draw("boxsame");
    
    lat.DrawLatexNDC(0.02,0.95,buf);
    lat.DrawLatexNDC(0.15,0.84,"All LCs");	
    
    if (isFC) myc->cd(2);
    else myc->cd(3);

    tree->Draw("eTotLC3>>hELC","","",1,entry);
    TH1F *hELC = (TH1F*)gDirectory->Get("hELC");
    double eLC3 = hELC->GetMean();
    hELC->Delete();

    tree->Draw("eTotAllLC3>>hEAllLC","","",1,entry);
    TH1F *hEAllLC = (TH1F*)gDirectory->Get("hEAllLC");
    double eAllLC3 = hEAllLC->GetMean();
    hEAllLC->Delete();

    
    gPad->SetLogz(1);
    gPad->SetRightMargin(0.18);
    
    getRZ(tree,hRZ,entry,isFC,false,doSeed,true);
    hRZ->GetZaxis()->SetRangeUser(0.01,100);
    
    //TExec *ex2 = new TExec("ex2","MyPalette();");
    //ex2->Draw();
    ex1->Draw();

    
    
    hRZ->Draw("colzsame");
    
    if (doSeed) lat.DrawLatexNDC(0.15,0.78,"Seed pos");	
    else lat.DrawLatexNDC(0.15,0.78,"LC pos");	
    lat.DrawLatexNDC(0.15,0.84,"SimCluster LCs");
    lat.DrawLatexNDC(0.15,0.72,isPU?"WithPU":isFC?"FineCalo":"Default");	
    //hRZtrue->Draw("boxsame");
    
    std::ostringstream counters;
    counters << "nSC10=" << outlier.nSC10 
	     << " nSC=" << outlier.nSC;
    lat.DrawLatexNDC(0.15,0.66,counters.str().c_str());
    counters.str("");
    counters << "nTS10=" << outlier.nTS10 
	     << " nTS=" << outlier.nTS;
    lat.DrawLatexNDC(0.15,0.6,counters.str().c_str());    
    //counters.str("");
    //counters << "eTot LC (#geq 3 RH) = " << eLC3
    //	     << " (" << eAllLC3 << ") GeV";
    //lat.DrawLatexNDC(0.15,0.54,counters.str().c_str());    
    
    if (isFC) myc->cd(4);
    else myc->cd(5);
    gPad->SetRightMargin(0.18);
    
    std::ostringstream lname;
    lname << "hPID_" << entry;
    if (isPU) lname << "_PU";
    if (isFC) lname << "_FC";
    unsigned nX = 12;
    unsigned nY = 4;
    double binsE[11] = {0,1,2,5,10,20,50,100,150,200,55*cosh(2.1)}; 
    double binX[13] = {-250,-220,-200,-120,-100,-16,0,16,30,120,200,220,250};
    double binY[5] = {320,321,330,340,350};
    if (isFC) {
      gPad->SetLogz(1);
      hPID = new TH2F(lname.str().c_str(),";PID;z @ boundary (cm);E (GeV)",nX,binX,nY,binY);
      tree->Draw(("sc_zAtB:sc_pdgid>>"+lname.str()).c_str(),"sc_energyAtB*(sc_energyAtB >= 0)","colz",1,entry);
    }
    else {
      //nX=500;
      nY = 10;
      hPID = new TH2F(lname.str().c_str(),";PID;SC E (GeV); nSC",nX,binX,nY,binsE);
      tree->Draw(("sc_energy:sc_pdgid>>"+lname.str()).c_str(),"","colz",1,entry);
    }
    for (unsigned iY(1);iY<nY+1;++iY){
      hPID->SetBinContent(nX,iY,hPID->GetBinContent(nX,iY)+hPID->GetBinContent(nX+1,iY));
      hPID->SetBinContent(1,iY,hPID->GetBinContent(0,iY)+hPID->GetBinContent(1,iY));
    }
    for (unsigned iX(1);iX<nX+1;++iX){
      //hPID->GetXaxis()->SetBinLabel(iX,pdgBinLabel(iX).c_str());
      hPID->SetBinContent(iX,1,hPID->GetBinContent(iX,0)+hPID->GetBinContent(iX,1));
      hPID->SetBinContent(iX,nY,hPID->GetBinContent(iX,nY)+hPID->GetBinContent(iX,nY+1));
    }
    gStyle->SetPaintTextFormat("3.0f");
    ex1->Draw();
    
    hPID->Draw("colztext89same");
    if (isFC) {
      std::cout << " Integral of energy: " << hPID->Integral()
		<< " truth " << truthE
		<< std::endl;
    }
    lat.DrawLatexNDC(0.15,0.84,isPU?"WithPU":isFC?"FineCalo":"Default");	
    
    myc->cd(6);
    gPad->SetRightMargin(0.18);
    
    gPad->SetLogz(1);
    lname.str("");
    lname << "hE_" << entry;
    if (isPU) lname << "_PU";
    if (isFC) lname << "_FC";
    nX = 9;
    nY = 10;
    double binsID[10] = {0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5};
    hE = new TH2F(lname.str().c_str(),";PDG cat;E_{SC} (GeV);nSC",nX,binsID,nY,binsE);
    if (!isFC) {
      hE->SetMarkerStyle(0);
      hE->SetMarkerColor(1);
    }
    hE->GetYaxis()->SetRangeUser(0,binsE[nY]);
    if (isFC) tree->Draw(("sc_energyAtB_eCut0:sc_pdgcat_eCut0>>"+lname.str()).c_str(),"","colz",1,entry);
    else tree->Draw(("sc_energy_eCut0:sc_pdgcat_eCut0>>"+lname.str()).c_str(),"","text same",1,entry);
    for (unsigned iY(1);iY<nY+1;++iY){
      hE->SetBinContent(nX,iY,hE->GetBinContent(nX,iY)+hE->GetBinContent(nX+1,iY));
      hE->SetBinContent(1,iY,hE->GetBinContent(0,iY)+hE->GetBinContent(1,iY));
    }
    for (unsigned iX(1);iX<nX+1;++iX){
      hE->GetXaxis()->SetBinLabel(iX,pdgBinLabel(iX).c_str());
      hE->SetBinContent(iX,1,hE->GetBinContent(iX,0)+hE->GetBinContent(iX,1));
      hE->SetBinContent(iX,nY,hE->GetBinContent(iX,nY)+hE->GetBinContent(iX,nY+1));
    }
    
    if (isFC) ex1->Draw();
    hE->Draw(isFC && !isPU?"colz same":"text same");
    lat.DrawLatexNDC(0.15,0.95,"colz: FineCalo, text: Default");	
    //lat.DrawLatexNDC(0.15,0.95,"colz: NoPU, text: WithPU");	
    
    myc->Update();
    if (!isFC){
      std::ostringstream lSave;
      lSave << plotDir << "/" << "RZ_entry" << entry;
      //if (doAllLC) lSave << "_allLC";
      if (doSeed) lSave << "_seedInfo";
      
      myc->Print((lSave.str()+".pdf").c_str());
      myc->Print((lSave.str()+".png").c_str());
    }
    //hRZtrue->Delete();
    //hRZ->Delete();

    mycP->cd(!isPU && isFC?1:2);
    gStyle->SetPaintTextFormat("3.3f");

    //gPad->SetLogy(1);
    const unsigned nCat = 9;
    lname.str("");
    lname << "hCat_" << entry;
    if (isFC) lname << "_FC";
    if (isPU) lname << "_PU";
    nX = 9;
    if (isFC) hCat = new TH1F(lname.str().c_str(),";PDG cat;E_{frac} @ boundary",nX,0.5,nX+0.5);
    else hCat = new TH1F(lname.str().c_str(),";PDG cat;E_{frac}",nX,0.5,nX+0.5);
    for (unsigned iC(0); iC<nCat;++iC){
      std::ostringstream lName;
      lName << "eFrac";
      if (isFC) lName << "AtB";
      lName << "_cat" << pdgBinLabel(iC+1) << ">>hA";
      tree->Draw((lName.str()+"(101,0,1.01)").c_str(),"","text",1,entry);
      TH1F *hA = (TH1F*)gDirectory->Get("hA");
      std::cout << pdgBinLabel(iC+1) << " " << hA->GetMean() << std::endl;
      hCat->SetBinContent(iC+1,hA->GetMean());
      hA->Delete();
    }
    std::cout << " Check " << hCat->GetEntries() << std::endl;
    for (unsigned iX(1);iX<nX+1;++iX){
      std::cout << iX << " " << hCat->GetBinContent(iX) << std::endl;
      hCat->GetXaxis()->SetBinLabel(iX,pdgBinLabel(iX).c_str());
    }
    hCat->SetMarkerSize(3);
    hCat->Draw("histtext");
    lat.DrawLatexNDC(0.15,0.84,isPU? "WithPU" : isFC ? "FineCalo":"Default");	  
    lat.DrawLatexNDC(0.02,0.95,buf);
    
      //lat.SetTextColor(1);
    mycP->Update();
    if (!isFC || isPU){
      std::ostringstream lSave;
      lSave << plotDir << "/" << "PDGCat_entry" << entry;
      //if (doAllLC) lSave << "_allLC";
      mycP->Print((lSave.str()+".pdf").c_str());
      mycP->Print((lSave.str()+".png").c_str());
    }

  } else {
    if (!isFC){
      mycT->cd(1);
      gPad->SetLogz(1);
      gPad->SetRightMargin(0.18);
      
      hRZall->GetZaxis()->SetRangeUser(0.01,100);
      hRZall->Draw("colz");
      //ex1->Draw();
      //hRZall->Draw("colzsame");
      hRZtrue->Draw("boxsame");
      
      lat.DrawLatexNDC(0.02,0.95,buf);
      lat.DrawLatexNDC(0.15,0.84,"All LCs");	
      
      mycT->cd(2);
      gPad->SetLogz(1);
      gPad->SetRightMargin(0.18);
      
      
      getRZ(tree,hRZs,entry,lCut,isFC);
      hRZs->GetZaxis()->SetRangeUser(0.01,100);
      hRZs->Draw("colz");
      std::ostringstream counters;
      counters.str("");
      counters << "nSC10=" << outlier.nSC10 
	       << " nSC=" << outlier.nSC;
      lat.DrawLatexNDC(0.15,0.94,"FineCalo");
      lat.DrawLatexNDC(0.15,0.84,counters.str().c_str());
      counters.str("");
      counters << "nTS10=" << outlier.nTS10 
	       << " nTS=" << outlier.nTS;
      lat.DrawLatexNDC(0.15,0.78,counters.str().c_str());
      lat.DrawLatexNDC(0.15,0.72,lCut.c_str());
      
      mycT->cd(4);
      gPad->SetLogz(1);
      gPad->SetRightMargin(0.18);
      
      getXZ(tree,hXZ,entry,lCut);
      hXZ->GetZaxis()->SetRangeUser(0.01,100);
      hXZ->Draw("colz");
      lat.DrawLatexNDC(0.15,0.84,lCut.c_str());
      
      mycT->cd(5);
      gPad->SetLogz(1);
      gPad->SetRightMargin(0.18);
      
      getYZ(tree,hYZ,entry,lCut);
      hYZ->GetZaxis()->SetRangeUser(0.01,100);
      hYZ->Draw("colz");
      lat.DrawLatexNDC(0.15,0.84,lCut.c_str());
    } else {
      mycT->cd(3);
      gPad->SetLogz(1);
      gPad->SetRightMargin(0.18);
      
      getRZ(tree,hRZs,entry,lCut,isFC);
      hRZs->GetZaxis()->SetRangeUser(0.01,100);
      hRZs->Draw("colz");
      std::ostringstream counters;
      counters.str("");
      counters << "nSC10=" << outlier.nSC10 
	       << " nSC=" << outlier.nSC;
      lat.DrawLatexNDC(0.15,0.84,counters.str().c_str());
      counters.str("");
      counters << "nTS10=" << outlier.nTS10 
	       << " nTS=" << outlier.nTS;
      lat.DrawLatexNDC(0.15,0.78,counters.str().c_str());
      lat.DrawLatexNDC(0.15,0.72,lCut.c_str());
      lat.DrawLatexNDC(0.15,0.94,isFC?"FineCalo":"Default");
    }
    mycT->cd(6);
    //gPad->SetLogy(1);
    //gPad->SetLogx(1);
    if (isFC && !isPU) tree->Draw("ts_energy>>hA(100,0,100)",lCutTS.c_str(),"",1,entry);
    else tree->Draw("ts_energy>>hB(100,0,100)",lCutTS.c_str(),"same",1,entry);
    if (isFC && !isPU){
      TH1F *hA = (TH1F*)gDirectory->Get("hA");
      hA->SetBinContent(100,hA->GetBinContent(100)+hA->GetBinContent(101));
      hA->SetLineColor(2);
      hA->GetYaxis()->SetRangeUser(0.1,hA->GetMaximum()*2);
      hA->SetTitle(";TS E (GeV);nTS");
    } else {
      TH1F *hB = (TH1F*)gDirectory->Get("hB");
      hB->SetBinContent(100,hB->GetBinContent(100)+hB->GetBinContent(101));
    }
    lat.SetTextColor(isFC?2:kBlue);
    lat.DrawLatexNDC(isFC?0.15:0.4,0.95,isFC ? "FineCalo":"Default");	  
    //lat.DrawLatexNDC(isFC?0.15:0.4,0.95,isPU ? "WithPU":"NoPU");	  
    
    lat.SetTextColor(1);
    
    if (!isFC){
      mycT->Update();
      std::ostringstream lSave;
      lSave << plotDir << "/" << "XYRZ_entry" << entry;
      //if (doAllLC) lSave << "_allLC";
      if (doSeed) lSave << "_seedInfo_" << lCutSave;
      
      mycT->Print((lSave.str()+".pdf").c_str());
      mycT->Print((lSave.str()+".png").c_str());
    }
  }
  
  return 0; 
};

TH1F* makePlot(TChain *tree,
	       TCanvas* & myc,
	       std::string plotDir,
	       char* buf,
	       unsigned index,
	       std::string var,
	       std::string varShort,
	       std::string cut,
	       std::string cutShort,
	       std::string title,
	       const unsigned nBins,
	       double binMin,
	       double binMax,
	       const bool doLogy,
	       std::string drawOption){
  
  myc->cd();
  //gStyle->SetOptStat("eMRuo");
  gStyle->SetOptStat(0);
  gPad->SetLogy(doLogy);
  
  std::ostringstream label;
  label << "myhist" << index;
  TH1F *hTmp = (TH1F*)gDirectory->Get(label.str().c_str());
  if (hTmp) hTmp->Delete();

  gROOT->cd();
  
  TH1F *myhist = new TH1F(label.str().c_str(),title.c_str(),nBins,binMin,binMax);
  
  tree->Draw((var+">>"+label.str()).c_str(),cut.c_str());
  if (doLogy) myhist->SetMinimum(0.1);
  //add overflows
  myhist->SetBinContent(nBins,myhist->GetBinContent(nBins)+myhist->GetBinContent(nBins+1));
  myhist->SetBinContent(1,myhist->GetBinContent(0)+myhist->GetBinContent(1));
  
  myhist->Draw(drawOption.c_str());
  
  TLatex lat;
  lat.DrawLatexNDC(0.01,0.95,buf);
  if (cut.size()>0){
    lat.SetTextSize(0.03);
    lat.DrawLatexNDC(0.5,0.95,cut.c_str());
  }
  std::ostringstream histstat;
  histstat << "Nsel=" << myhist->Integral();
  lat.SetTextSize(0.05);
  lat.DrawLatexNDC(0.01,0.02,histstat.str().c_str());
  
  myc->Update();
  myc->Print((plotDir+"/"+varShort+cutShort+".pdf").c_str());
  myc->Print((plotDir+"/"+varShort+cutShort+".png").c_str());
  
  
  return myhist;
};
TH1F* makePlot(TChain *tree,
	       TCanvas* & myc,
	       std::string plotDir,
	       char* buf,
	       std::string var,
	       std::string varShort,
	       std::string cut,
	       std::string cutShort,
	       std::string title,
	       const unsigned nBins,
	       double binMin,
	       double binMax,
	       const bool doLogy,
	       std::string drawOption){

  return makePlot(tree,
	       myc,
	       plotDir,
	       buf,
	       0,
	       var,
	       varShort,
	       cut,
	       cutShort,
	       title,
	       nBins,
	       binMin,
	       binMax,
	       doLogy,
	       drawOption);
};


int makePlot2D(TChain *tree,
	       TCanvas* & myc,
	       std::string plotDir,
	       char* buf,
	       std::string var,
	       std::string varShort,
	       std::string cut,
	       std::string cutShort,
	       std::string title,
	       const unsigned nBinsX,
	       double binXMin,
	       double binXMax,
	       const unsigned nBinsY,
	       double binYMin,
	       double binYMax,
	       const bool doLogz,
	       std::string drawOption){
  
  myc->cd();
  //gStyle->SetOptStat("eMRuo");
  gStyle->SetOptStat(0);
  gPad->SetLogz(doLogz);
  
  
  TH2F *hTmp = (TH2F*)gDirectory->Get("myhist2D");
  if (hTmp) hTmp->Delete();

  gROOT->cd();
  
  TH2F *myhist = new TH2F("myhist2D",title.c_str(),nBinsX,binXMin,binXMax,nBinsY,binYMin,binYMax);
  
  tree->Draw((var+">>myhist2D").c_str(),cut.c_str());

  //add overflows
  for (unsigned iY(1);iY<nBinsY+1;++iY){
    myhist->SetBinContent(nBinsX,iY,myhist->GetBinContent(nBinsX,iY)+myhist->GetBinContent(nBinsX+1,iY));
    myhist->SetBinContent(1,iY,myhist->GetBinContent(0,iY)+myhist->GetBinContent(1,iY));
  }
  for (unsigned iX(1);iX<nBinsX+1;++iX){
    //myhist->GetXaxis()->SetBinLabel(iX,pdgBinLabel(iX).c_str());
    myhist->SetBinContent(iX,1,myhist->GetBinContent(iX,0)+myhist->GetBinContent(iX,1));
    myhist->SetBinContent(iX,nBinsY,myhist->GetBinContent(iX,nBinsY)+myhist->GetBinContent(iX,nBinsY+1));
  }
  
  myhist->Draw(drawOption.c_str());
  
  TLatex lat;
  lat.DrawLatexNDC(0.01,0.95,buf);
  if (cut.size()>0){
    lat.SetTextSize(0.03);
    lat.DrawLatexNDC(0.5,0.95,cut.c_str());
  }
  std::ostringstream histstat;
  histstat << "Nsel=" << myhist->Integral();
  lat.DrawLatexNDC(0.01,0.01,histstat.str().c_str());
  
  myc->Update();
  myc->Print((plotDir+"/"+varShort+cutShort+".pdf").c_str());
  myc->Print((plotDir+"/"+varShort+cutShort+".png").c_str());
  
  
  return 0;
};

int plotSelection(std::string filePath,
		  std::string cutStr,
		  std::string cutShort,
		  const bool doBoth){
  
  SetTdrStyle();

  if (cutShort.size()>0) cutShort = "_"+cutShort;

  //std::string filePath = "D49_AllTracksters";
  const unsigned iReg = 4;
  const unsigned ptval = 50;
  const unsigned etaval = 21;
 
  std::string iter = "Sim";
  std::ostringstream pteta;
  pteta << "pt" << ptval << "_eta" << etaval;

  std::ostringstream lPlotDir;
  if (!doBoth) lPlotDir << filePath << "/Truth/" << iter << "/" << pteta.str() << "/Plots/";
  else lPlotDir << "DefAndFC/Truth/" << iter << "/" << pteta.str() << "/Plots/";
  
  if (system(("mkdir -p "+lPlotDir.str()).c_str())){
    std::cout << " -- Cannot create output dir..." << lPlotDir.str() << std::endl;
    system(("echo \"-- return value \"$?"));
    return 1;
  }

  const unsigned nRuns = 10;

  TChain *tree = new TChain(("ticlTree/TSTree_"+iter).c_str());  
  TChain *treeLC = new TChain("ticlTree/treeLC");  
  //if (filePath.find("FineCalo")!=filePath.npos) tree->AddFile((filePath+"/"+reg[iReg]+"/FineCalo/step3ticl_"+pteta.str()+"_FlatTracksters.root").c_str());
  if (doBoth) {
    for (unsigned iR(0); iR<nRuns; ++iR){
      std::ostringstream infile;
      infile << "D49_FineCalo/" << reg[iReg]
	     <<  "/FineCalo/step3ticl_" << pteta.str()
	     << "_run" << iR
	     << "_FlatTracksters.root";
      tree->AddFile(infile.str().c_str());
      treeLC->AddFile(infile.str().c_str());
    }
    //tree->AddFile(("D49_FineCalo/FineCalo/step3ticl_"+pteta.str()+"_FlatTracksters.root").c_str());
    //treeLC->AddFile(("D49_FineCalo/FineCalo/step3ticl_"+pteta.str()+"_FlatTracksters.root").c_str());
    tree->AddFriend(treeLC);
    TChain *tree2 = new TChain(("ticlTree/TSTree_"+iter).c_str());  
    tree2->AddFile(("D49_DefSC/"+reg[iReg]+"/step3ticl_"+pteta.str()+"_FlatTracksters.root").c_str());
    tree->AddFriend(tree2,"Def");
    TChain *treeLC2 = new TChain("ticlTree/treeLC");  
    treeLC2->AddFile(("D49_DefSC/"+reg[iReg]+"/step3ticl_"+pteta.str()+"_FlatTracksters.root").c_str());
    tree->AddFriend(treeLC2,"Def");
  } else {
    if (filePath.find("FineCalo")!=filePath.npos) tree->AddFile((filePath+"/FineCalo/step3ticl_"+pteta.str()+"_FlatTracksters.root").c_str());
    else tree->AddFile((filePath+"/"+reg[iReg]+"/step3ticl_"+pteta.str()+"_FlatTracksters.root").c_str());
    if (filePath.find("FineCalo")!=filePath.npos) treeLC->AddFile((filePath+"/FineCalo/step3ticl_"+pteta.str()+"_FlatTracksters.root").c_str());
    else treeLC->AddFile((filePath+"/"+reg[iReg]+"/step3ticl_"+pteta.str()+"_FlatTracksters.root").c_str());
    tree->AddFriend(treeLC);
  }
  
  
  TChain *friendTree = new TChain("TruthTree");
  if (doBoth) {
    TChain *friendTree2 = new TChain("TruthTree");
    friendTree->AddFile(("D49_FineCalo/Truth/"+iter+"/"+pteta.str()+"/Debug/TruthTree_"+regShort[iReg]+".root").c_str());
    friendTree2->AddFile(("D49_DefSC/Truth/"+iter+"/"+pteta.str()+"/Debug/TruthTree_"+regShort[iReg]+".root").c_str());
    tree->AddFriend(friendTree);
    tree->AddFriend(friendTree2,"Def");
  } else {
    friendTree->AddFile((filePath+"/Truth/"+iter+"/"+pteta.str()+"/Debug/TruthTree_"+regShort[iReg]+".root").c_str());
    tree->AddFriend(friendTree);
  }
  
  if (!tree) {
    std::cout << " problem getting tree with friends." << std::endl;
    return 1;
  }
  
  TCanvas *myc = new TCanvas("myc","myc",1);
  
  char buf[500];
  //  sprintf(buf,"%s %s %s",regShort[iReg].c_str(),iter.c_str(),pteta.str().c_str());
  sprintf(buf,"%s %s",regShort[iReg].c_str(),pteta.str().c_str());


  if (doBoth){

    /*makePlot(tree,myc,lPlotDir.str(),buf,
	     "Def.nSC_eCut0","nSC_Def",
	     cutStr,cutShort,
	     ";nSC;Events",
	     50,0,50,1,
	     "hist");
    
    makePlot(tree,myc,lPlotDir.str(),buf,
	     "Def.nSC_eCut10","nSC10_Def",
	     cutStr,cutShort,
	     ";nSC (E_{SC} > 0.10 #times E_{CP});Events",
	     10,0,10,1,
	     "histtext");
    
    makePlot(tree,myc,lPlotDir.str(),buf,
	     "nSC_eCut0","nSC",
	     cutStr,cutShort,
	     ";nSC;Events",
	     50,0,50,1,
	     "hist");
    
    makePlot(tree,myc,lPlotDir.str(),buf,
	     "nSC_eCut10","nSC10",
	     cutStr,cutShort,
	     ";nSC (E_{SC} > 0.10 #times E_{CP});Events",
	     10,0,10,1,
	     "histtext");
    
    
    makePlot(tree,myc,lPlotDir.str(),buf,
	     "eTotSCAtB/cp_energy[0]","eTotSCAtBoverECP",
	     cutStr,cutShort,
	     ";E^{tot}_{SC} @ boundary/E_{CP};Events",
	     100,0.5,1.5,1,
	     "hist");
    
    makePlot(tree,myc,lPlotDir.str(),buf,
	     "eTotSC/cp_energy[0]","eTotSCoverECP_FC",
	     cutStr,cutShort,
	     ";E^{tot}_{SC} @ vtx/E_{CP};Events",
	     100,0.5,1.5,1,
	     "hist");
    
    makePlot(tree,myc,lPlotDir.str(),buf,
	     "Def.eTotSC/cp_energy[0]","eTotSCoverECP_Def",
	     cutStr,cutShort,
	     ";E^{tot}_{SC} @ vtx/E_{CP};Events",
	     100,0.5,1.5,1,
	     "hist");
    
    makePlot(tree,myc,lPlotDir.str(),buf,
	     "Def.eTotSC_noPi/cp_energy[0]","eTotSCNoPioverECP_Def",
	     cutStr,cutShort,
	     ";E^{tot}_{SC,no #pi} @ vtx/E_{CP};Events",
	     100,0.5,1.5,1,
	     "hist");
    
    
    makePlot(tree,myc,lPlotDir.str(),buf,
	     "sc_zAtB","sc_zAtB",
	     cutStr,cutShort,
	     ";z (cm) @ boundary (SC);Events",
	     200,300,500,1,
	     "hist");
    
    makePlot(tree,myc,lPlotDir.str(),buf,
	     "sc_zAtB_eCut0","sc_zAtB_eCut0",
	     cutStr,cutShort,
	     ";z (cm) @ boundary (SC);Events",
	     200,300,500,1,
	     "hist");
    
    makePlot(tree,myc,lPlotDir.str(),buf,
	     "eTotSC/Def.eTotSC_noPi","eTotSC_FCoverDefNoPi",
	     cutStr,cutShort,
	     ";E^{tot}_{SC} FC /E^{tot}_{SC,noPi} Def ;Events",
	     100,0,2,1,
	     "hist");
    
    makePlot(tree,myc,lPlotDir.str(),buf,
	     "(eTotSC-eTotSCAtB)/eTotSC","eTotSC_atB",
	     cutStr,cutShort,
	     ";E^{tot}_{SC} (AtVtx-AtB) /E^{tot}_{SC} AtVtx (FC) ;Events",
	     100,-0.1,0.3,1,
	     "hist");
    
    makePlot(tree,myc,lPlotDir.str(),buf,
	     "eTotAllLC/Def.eTotAllLC","eTotAllLC_FCoverDef",
	     cutStr,cutShort,
	     ";E^{tot}_{AllLC} FC /E^{tot}_{AllLC} Def ;Events",
	     100,0.5,1.5,1,
	     "hist");
    
    makePlot(tree,myc,lPlotDir.str(),buf,
	     "eTotAllLC/cp_energy[0]","eTotAllLC_FCoverCP",
	     cutStr,cutShort,
	     ";E^{tot}_{AllLC} FC /E_{CP} ;Events",
	     100,0.5,1.5,1,
	     "hist");
    
    makePlot(tree,myc,lPlotDir.str(),buf,
	     "eTotAllLC/eTotSCAtB","eTotAllLCoverSCAtB",
	     cutStr,cutShort,
	     ";E^{tot}_{AllLC} /E^{tot}_{SC} @ boundary ;Events",
	     100,0.5,1.5,1,
	     "hist");
    
    makePlot(tree,myc,lPlotDir.str(),buf,
	     "eTotLC/Def.eTotLC","eTotLC_FCoverDef",
	     cutStr,cutShort,
	     ";E^{tot}_{LC} (FC/Def);Events",
	     100,0.5,1.5,1,
	     "hist");
    
    makePlot(tree,myc,lPlotDir.str(),buf,
	     "eTotLC3/Def.eTotLC3","eTotLC3_FCoverDef",
	     cutStr,cutShort,
	     ";E^{tot}_{LC3} (FC/Def);Events",
	     100,0.5,1.5,1,
	     "hist");
    
    makePlot(tree,myc,lPlotDir.str(),buf,
	     "(eTotLC-Def.eTotLC)/eTotAllLC","eTotLC_FCminusDef",
	     cutStr,cutShort,
	     ";E^{tot}_{LC} (FC-Def)/E^{tot}_{AllLC};Events",
	     120,-0.2,0.4,1,
	     "hist");
    
    makePlot(tree,myc,lPlotDir.str(),buf,
	     "(eTotLC3-Def.eTotLC3)/eTotAllLC3","eTotLC3_FCminusDef",
	     cutStr,cutShort,
	     ";E^{tot}_{LC3} (FC-Def)/E^{tot}_{AllLC3};Events",
	     120,-0.2,0.4,1,
	     "hist");
    
    makePlot(tree,myc,lPlotDir.str(),buf,
	     "eTotLC/eTotAllLC","eTotLC_FCoverAll",
	     cutStr,cutShort,
	     ";E^{tot}_{LC} (FC)/E^{tot}_{AllLC};Events",
	     100,0.5,1.5,1,
	     "hist");
    */
    makePlot(tree,myc,lPlotDir.str(),buf,
	     "eTotLC3/eTotAllLC3","eTotLC3_FCoverAll",
	     cutStr,cutShort,
	     ";E^{tot}_{LC3} (FC)/E^{tot}_{AllLC3};Events",
	     100,0.5,1.5,1,
	     "hist");
    
    makePlot(tree,myc,lPlotDir.str(),buf,
	     "(eTotAllLC-eTotAllLC3)/eTotAllLC","eTotLC_RH3overAll",
	     cutStr,cutShort,
	     ";(E^{tot}_{AllLC}-E^{tot}_{AllLC3})/E^{tot}_{AllLC};Events",
	     120,0,1.01,1,
	     "hist");
    
    makePlot(tree,myc,lPlotDir.str(),buf,
	     "(eTotAllLC-eTotAllLC2)/eTotAllLC","eTotLC_RH2overAll",
	     cutStr,cutShort,
	     ";(E^{tot}_{AllLC}-E^{tot}_{AllLC2})/E^{tot}_{AllLC};Events",
	     120,0,1.01,1,
	     "hist");
    
    /*makePlot(tree,myc,lPlotDir.str(),buf,
	     "Def.eTotLC/Def.eTotAllLC","eTotLC_DefoverAll",
	     cutStr,cutShort,
	     ";E^{tot}_{LC} (Def)/E^{tot}_{AllLC};Events",
	     100,0.5,1.5,1,
	     "hist");
    
    makePlot(tree,myc,lPlotDir.str(),buf,
	     "Def.eTotLC3/Def.eTotAllLC3","eTotLC3_DefoverAll",
	     cutStr,cutShort,
	     ";E^{tot}_{LC3} (Def)/E^{tot}_{AllLC3};Events",
	     100,0.5,1.5,1,
	     "hist");*/
    /*
    const unsigned nCat = 10;
    std::ostringstream lname;
    for (unsigned iC(0); iC<nCat;++iC){
      lname.str("");
      lname << "eFracAtB_cat" <<  pdgBinLabel(iC+1);
      makePlot(tree,myc,lPlotDir.str(),buf,
	       lname.str(),lname.str()+"_FC",
	       cutStr,cutShort,
	       ";"+lname.str()+";Events",
	       101,0.,1.01,1,
	       "hist");
      lname.str("");
      lname << "eFrac_cat" <<  pdgBinLabel(iC+1);
      makePlot(tree,myc,lPlotDir.str(),buf,
	       "Def."+lname.str(),lname.str()+"_Def",
	       cutStr,cutShort,
	       ";"+lname.str()+";Events",
	       101,0.,1.01,1,
	       "hist");
    }
    
    makePlot(tree,myc,lPlotDir.str(),buf,
	     "(eFracAtB_catgamma+eFracAtB_cate)-(Def.eFrac_catgamma+Def.eFrac_cate)",
	     "fracEM_FCminusDef",
	     cutStr,cutShort,
	     ";eFrac(gamma+e) FC-Def;Events",
	     100,-1,1,1,
	     "hist");
 
    makePlot(tree,myc,lPlotDir.str(),buf,
	     "(eFracAtB_catmu)-(Def.eFrac_catmu)",
	     "fracMu_FCminusDef",
	     cutStr,cutShort,
	     ";eFrac(mu) FC-Def;Events",
	     100,-1,1,1,
	     "hist");

    makePlot(tree,myc,lPlotDir.str(),buf,
	     "(eFracAtB_catpi0+eFracAtB_catn+eFracAtB_catneHad)-(Def.eFrac_catpi0+Def.eFrac_catn+Def.eFrac_catneHad)",
	     "fracNEHAD_FCminusDef",
	     cutStr,cutShort,
	     ";eFrac(neHad) FC-Def;Events",
	     100,-1,1,1,
	     "hist");

    makePlot(tree,myc,lPlotDir.str(),buf,
	     "(eFracAtB_catpipm+eFracAtB_catp+eFracAtB_catchHad)-(Def.eFrac_catpipm+Def.eFrac_catp+Def.eFrac_catchHad)",
	     "fracCHHAD_FCminusDef",
	     cutStr,cutShort,
	     ";eFrac(chHad) FC-Def;Events",
	     100,-1,1,1,
	     "hist");
    */
    makePlot(tree,myc,lPlotDir.str(),buf,
	     "nBlobs","nBlobs",
	     cutStr,cutShort,
	     ";nBlobs ;Events",
	     10,0,10,
	     1,"histtext");

    makePlot(tree,myc,lPlotDir.str(),buf,
	     "nAllLC3/nAllLC","nAllLC3over1",
	     cutStr,cutShort,
	     ";nAllLC3/nAllLC ;Events",
	     100,0,1.01,
	     1,"hist");

    makePlot(tree,myc,lPlotDir.str(),buf,
	     "nAllLC2/nAllLC","nAllLC2over1",
	     cutStr,cutShort,
	     ";nAllLC2/nAllLC ;Events",
	     100,0,1.01,
	     1,"hist");

    makePlot(tree,myc,lPlotDir.str(),buf,
	     "EBlobs/eTotAllLC","eBlobsOverEtot",
	     cutStr,cutShort,
	     ";E_{Blobs} / E^{tot}_{AllLC} ;Events",
	     100,0,1.1,
	     1,"hist");

    makePlot(tree,myc,lPlotDir.str(),buf,
	     "EBlobs/eTotAllLC3","eBlobsOverEtot3",
	     cutStr,cutShort,
	     ";E_{Blobs} / E^{tot}_{AllLC3} ;Events",
	     100,0,1.1,
	     1,"hist");

    makePlot(tree,myc,lPlotDir.str(),buf,
	     "EBlobs/eTotSCAtB","eBlobsOverESCAtB",
	     cutStr,cutShort,
	     ";E_{Blobs} / E_{SC}^{tot} @ boundary ;Events",
	     100,0,1.2,
	     1,"hist");

    makePlot(tree,myc,lPlotDir.str(),buf,
	     "minMiss","minMiss",
	     cutStr,cutShort,
	     ";minMiss ;Events",
	     20,0,20,
	     1,"hist");

    makePlot(tree,myc,lPlotDir.str(),buf,
	     "maxMiss","maxMiss",
	     cutStr,cutShort,
	     ";maxMiss ;Events",
	     40,0,40,
	     1,"hist");

    makePlot(tree,myc,lPlotDir.str(),buf,
	     "maxEfracBlobs","maxEfracBlobs",
	     cutStr,cutShort,
	     "; maxEfracBlobs;Events",
	     100,0,1.01,
	     1,"hist");

    makePlot(tree,myc,lPlotDir.str(),buf,
	     "minEfracBlobs","minEfracBlobs",
	     cutStr,cutShort,
	     ";minEfracBlobs ;Events",
	     100,0,1.01,
	     1,"hist");

    makePlot(tree,myc,lPlotDir.str(),buf,
	     "maxLengthBlobs","maxLengthBlobs",
	     cutStr,cutShort,
	     ";maxLengthBlobs ;Events",
	     50,0,50,
	     1,"hist");

    makePlot(tree,myc,lPlotDir.str(),buf,
	     "minLengthBlobs","minLengthBlobs",
	     cutStr,cutShort,
	     ";minLengthBlobs ;Events",
	     50,0,50,
	     1,"hist");

    makePlot(tree,myc,lPlotDir.str(),buf,
	     "energy_blobs/eTotAllLC","energy_blobs",
	     cutStr,cutShort,
	     ";energy_blobs/E^{tot}_{AllLC} ;Events",
	     100,0,1.01,
	     1,"hist");

    makePlot(tree,myc,lPlotDir.str(),buf,
	     "length_blobs","length_blobs",
	     cutStr,cutShort,
	     "; length_blobs;Events",
	     50,0,50,
	     1,"hist");

    makePlot(tree,myc,lPlotDir.str(),buf,
	     "eFrac_blobs","eFrac_blobs",
	     cutStr,cutShort,
	     "; eFrac_blobs;Events",
	     100,0,1.01,
	     1,"hist");

    makePlot(tree,myc,lPlotDir.str(),buf,
	     "firstL_blobs","firstL_blobs",
	     cutStr,cutShort,
	     ";firstL_blobs ;Events",
	     50,0,50,
	     1,"hist");

    makePlot(tree,myc,lPlotDir.str(),buf,
	     "nLC_blobs","nLC_blobs",
	     cutStr,cutShort,
	     ";nLC_blobs ;Events",
	     100,0,300,
	     1,"hist");

    makePlot(tree,myc,lPlotDir.str(),buf,
	     "nLC_blobs/length_blobs","nLCOverLength_blobs",
	     cutStr,cutShort,
	     ";nLC / length blobs ;Events",
	     100,0,20,
	     1,"hist");

    makePlot(tree,myc,lPlotDir.str(),buf,
	     "nMiss_blobs","nMiss_blobs",
	     cutStr,cutShort,
	     ";nMiss_blobs ;Events",
	     30,0,30,
	     1,"hist");
    
    makePlot(tree,myc,lPlotDir.str(),buf,
	     "ts_energy/eTotAllLC3","eTSOverEAll3",
	     cutStr,cutShort,
	     ";E_{TS}/E_{AllLC3}^{tot};Events",
	     100,0,2,1,
	     "hist");
    
  } else {
  
    makePlot(tree,myc,lPlotDir.str(),buf,
	     "nSC_eCut0","nSC",
	     cutStr,cutShort,
	     ";nSC;Events",
	     20,0,20,1,
	     "histtext");
    
  makePlot(tree,myc,lPlotDir.str(),buf,
	   "nSC_eCut1","nSC_eCut1",
	   cutStr,cutShort,
	   ";nSC (E_{SC} > 0.01 #times E_{CP});Events",
	   20,0,20,1,
	   "histtext");

  makePlot(tree,myc,lPlotDir.str(),buf,
	   "nSC_eCut5","nSC_eCut5",
	   cutStr,cutShort,
	   ";nSC (E_{SC} > 0.05 #times E_{CP});Events",
	   12,0,12,1,
	   "histtext");

  makePlot(tree,myc,lPlotDir.str(),buf,
	   "nSC_eCut10","nSC_eCut10",
	   cutStr,cutShort,
	   ";nSC (E_{SC} > 0.10 #times E_{CP});Events",
	   10,0,10,1,
	   "histtext");
  
  makePlot(tree,myc,lPlotDir.str(),buf,
	   "sc_pdgid_eCut10","sc_pdgid_eCut10",
	   cutStr,cutShort,
	   ";PDGID (E_{SC} > 0.10 #times E_{CP});Events",
	   500,-250,250,1,
	   "histtext");

  makePlot(tree,myc,lPlotDir.str(),buf,
	   "sc_pdgid_eCut0","sc_pdgid_eCut0",
	   cutStr,cutShort,
	   ";PDGID (SC);Events",
	   500,-250,250,1,
	   "histtext");

  makePlot(tree,myc,lPlotDir.str(),buf,
	   "sc_dRtoCP_eCut0","sc_dRtoCP_eCut0",
	   cutStr,cutShort,
	   ";#Delta R (CP,SC);Events",
	   100,0,2,1,
	   "histtext");

  makePlot(tree,myc,lPlotDir.str(),buf,
	   "sc_zAtB_eCut0","sc_zAtB_eCut0",
	   cutStr,cutShort,
	   ";z (cm) @ boundary (SC);Events",
	   500,0,500,1,
	   "histtext");

  makePlot(tree,myc,lPlotDir.str(),buf,
	   "nTS_eCut10","nTS_eCut10",
	   cutStr,cutShort,
	   ";nTS (E_{TS} > 0.10 #times E_{CP});Events",
	   10,0,10,1,
	   "histtext");
  
  makePlot(tree,myc,lPlotDir.str(),buf,
	   "ts_nBlobs","ts_nBlobs",
	     cutStr,cutShort,
	   ";nBlobs;Events",
	   10,0,10,1,
	   "histtext");
  
  makePlot(tree,myc,lPlotDir.str(),buf,
	   "eTotSCAtB/cp_energy[0]","eTotSCAtBoverECP",
	   cutStr,cutShort,
	   ";E^{tot}_{SC} @ boundary/E_{CP};Events",
	   100,0.5,1.5,1,
	   "hist");

  makePlot(tree,myc,lPlotDir.str(),buf,
	   "eTotSC/cp_energy[0]","eTotSCoverECP",
	   cutStr,cutShort,
	   ";E^{tot}_{SC} @ vtx/E_{CP};Events",
	   100,0.5,1.5,1,
	   "hist");

  /*
  makePlot(tree,myc,lPlotDir.str(),buf,
	   "eTotSC_eCut10/cp_energy[0]","eTotSCeCut10overECP_nSC1",
	   cutStr,cutShort,
	   "nSC_eCut10==1",
	   ";E^{tot}_{SC10}/E_{CP};Events",
	   100,0,2,1,
	   "hist");
  
  makePlot(tree,myc,lPlotDir.str(),buf,
	   "ts_energy/sc_energy","eTSOverESC",
	   cutStr,cutShort,
	   "",
	   ";E_{TS}/E_{SC};Events",
	   100,0,2,1,
	   "hist");

  makePlot(tree,myc,lPlotDir.str(),buf,
	   "eTotSC/cp_energy[0]","eTotSCoverECP",
	   "",
	   ";E^{tot}_{SC}/E_{CP};Events",
	   100,0.5,1.5,1,
	   "hist");

  makePlot(tree,myc,lPlotDir.str(),buf,
	   "eTotSCAtB/cp_energy[0]","eTotSCAtBoverECP",
	   "",
	   ";E^{tot}_{SC} @ boundary/E_{CP};Events",
	   100,0.5,1.5,1,
	   "hist");
  
  makePlot(tree,myc,lPlotDir.str(),buf,
	   "eTotSC_noNeutron/cp_energy[0]","eTotSCNoNeutronoverECP",
	   "",
	   ";E^{tot}_{SC} (no neutron)/E_{CP};Events",
	   100,0.5,1.5,1,
	   "hist");
  
  makePlot(tree,myc,lPlotDir.str(),buf,
	   "eTotSC_noProton/cp_energy[0]","eTotSCNoProtonoverECP",
	   "",
	   ";E^{tot}_{SC} (no proton)/E_{CP};Events",
	   100,0.5,1.5,1,
	   "hist");
  
  makePlot(tree,myc,lPlotDir.str(),buf,
	   "eTotSC_noNucleon/cp_energy[0]","eTotSCNoNucleonoverECP",
	   "",
	   ";E^{tot}_{SC} (no nucleon)/E_{CP};Events",
	   100,0.5,1.5,1,
	   "hist");
  
  
  makePlot(tree,myc,lPlotDir.str(),buf,
	   "eTotSC_noPi/cp_energy[0]","eTotSCNoPioverECP",
	   "",
	   ";E^{tot}_{SC} (#pi removed)/E_{CP};Events",
	   100,0.5,1.5,1,
	   "hist");
  */  
  }
  
  return 0;
};

int getDisplay(std::string plotDir,
	       const Outlier & outlier,
	       TChain* tree,
	       TChain* treeLC,
	       TCanvas* & myc,
	       TCanvas* & mycT,
	       TCanvas* & mycP,
	       char* buf,
	       TH2F * & hRZtrue,
	       TH2F * & hRZall,
	       TH2F * & hRZ,
	       TH2F * & hPID,
	       TH2F * & hE,
	       TH1F * & hCat,
	       const std::string lCut,
	       const std::string lCutTS,
	       const std::string lCutSave,
	       TH2F * & hRZs,
	       TH2F * & hXZ,
	       TH2F * & hYZ,
	       const bool isPU
	       ){
  
  SetTdrStyle();
  gStyle->SetOptStat(0);
  const unsigned entry = outlier.event;
  
  
  return plotEvtDisplay(tree,
			treeLC,
			myc,mycT,mycP,
			hRZtrue,hRZall,
			hRZ,hPID,hE,hCat,
			lCut,lCutTS,lCutSave,
			hRZs,hXZ,hYZ,
			plotDir,
			buf,
			false,true,
			outlier,
			isPU);
}

void getDisplays(){

  std::map<int,std::pair<int,std::string> > pdgMap = getPDGMap("pdgList.dat");

  const unsigned nF = 2;
  //std::string filePath[2] = {"D49_FineCalo","D49_DefSC"};
  std::string filePath[2] = {"D49_FineCalo","D49_FineCalo"};

  const unsigned iReg[nF] = {4,4};
  const unsigned ptval = 50;
  const unsigned etaval = 21;

  std::string iter = "Sim";
  std::ostringstream pteta;
  pteta << "pt" << ptval << "_eta" << etaval;

   
  const unsigned nC = 6;
  const std::string lCut[6] = {
    "",
    "ts_photon_proba[lc_TSidx]+ts_ele_proba[lc_TSidx]>0.5",
    "ts_mu_proba[lc_TSidx]>0.5",
    "ts_pi0_proba[lc_TSidx]>0.5",
    "ts_chHad_proba[lc_TSidx]>0.5",
    "ts_neHad_proba[lc_TSidx]>0.5"
  };
  const std::string lCutTS[6] = {
    "",
    "ts_photon_proba+ts_ele_proba>0.5",
    "ts_mu_proba>0.5",
    "ts_pi0_proba>0.5",
    "ts_chHad_proba>0.5",
    "ts_neHad_proba>0.5"
  };
  const std::string lCutSave[6] = {
    "",
    "EM",
    "mu",
    "pi0",
    "chHad",
    "neHad"
  };
  
  TCanvas *mycTest = new TCanvas("mycTest","mycTest",1);
  TCanvas *mycT = new TCanvas("mycT","mycT",1);
  TCanvas *myc = new TCanvas("mycD","mycD",1);
  TCanvas *mycP = new TCanvas("mycP","mycP",1);
  myc->Divide(3,2);
  mycP->Divide(2,1);
  mycT->Divide(3,2);

  TH2F* hRZtrue[nF];
  TH2F* hRZall[nF];
  TH2F* hRZ[nF];
  TH2F* hPID[nF];
  TH2F* hE[nF];
  TH1F* hCat[nF];
  TH2F* hRZs[nF];
  TH2F* hXZ = 0;
  TH2F* hYZ = 0;

  
  std::vector<Outlier> evtNumbers[nF];
  std::ostringstream lPlotDir[nF];
  TChain *tree[nF];
  TChain *treeLC[nF];
  
  for (unsigned iF(0);iF<nF;++iF){
    hRZtrue[iF] = 0;
    hRZall[iF] = 0;
    hRZ[iF] = 0;
    hPID[iF] = 0;
    hE[iF] = 0;
    hCat[iF] = 0;
    hRZs[iF] = 0;
    
    bool isFC = filePath[iF].find("FineCalo")!=filePath[iF].npos;
    bool isPU = iReg[iF]==5;
      
    getOutliers(filePath[iF]+"/Truth/Sim/pt50_eta21/Debug/eventsToPrint_"+regShort[iReg[iF]]+".dat",
		evtNumbers[iF]);

    lPlotDir[iF] << filePath[iF] << "/Truth/" << iter << "/" << pteta.str() << "/Displays/";
    
    if (system(("mkdir -p "+lPlotDir[iF].str()).c_str())){
      std::cout << " -- Cannot create output dir..." << lPlotDir[iF].str() << std::endl;
      system(("echo \"-- return value \"$?"));
      return;
    }
    
    tree[iF] = new TChain(("ticlTree/TSTree_"+iter).c_str());  
    treeLC[iF] = new TChain("ticlTree/treeLC");  

    if (isFC){
      tree[iF]->AddFile((filePath[iF]+"/"+reg[iReg[iF]]+"/FineCalo/step3ticl_"+pteta.str()+"_run0_FlatTracksters.root").c_str());
      tree[iF]->AddFile((filePath[iF]+"/"+reg[iReg[iF]]+"/FineCalo/step3ticl_"+pteta.str()+"_run1_FlatTracksters.root").c_str());
      tree[iF]->AddFile((filePath[iF]+"/"+reg[iReg[iF]]+"/FineCalo/step3ticl_"+pteta.str()+"_run2_FlatTracksters.root").c_str());
      tree[iF]->AddFile((filePath[iF]+"/"+reg[iReg[iF]]+"/FineCalo/step3ticl_"+pteta.str()+"_run3_FlatTracksters.root").c_str());
      tree[iF]->AddFile((filePath[iF]+"/"+reg[iReg[iF]]+"/FineCalo/step3ticl_"+pteta.str()+"_run4_FlatTracksters.root").c_str());
      tree[iF]->AddFile((filePath[iF]+"/"+reg[iReg[iF]]+"/FineCalo/step3ticl_"+pteta.str()+"_run5_FlatTracksters.root").c_str());
      tree[iF]->AddFile((filePath[iF]+"/"+reg[iReg[iF]]+"/FineCalo/step3ticl_"+pteta.str()+"_run6_FlatTracksters.root").c_str());
      tree[iF]->AddFile((filePath[iF]+"/"+reg[iReg[iF]]+"/FineCalo/step3ticl_"+pteta.str()+"_run7_FlatTracksters.root").c_str());
      tree[iF]->AddFile((filePath[iF]+"/"+reg[iReg[iF]]+"/FineCalo/step3ticl_"+pteta.str()+"_run8_FlatTracksters.root").c_str());
      tree[iF]->AddFile((filePath[iF]+"/"+reg[iReg[iF]]+"/FineCalo/step3ticl_"+pteta.str()+"_run9_FlatTracksters.root").c_str());
    }
    //else if (isFC) tree[iF]->AddFile((filePath[iF]+"/FineCalo/step3ticl_"+pteta.str()+"_FlatTracksters.root").c_str());
    else tree[iF]->AddFile((filePath[iF]+"/"+reg[iReg[iF]]+"/step3ticl_"+pteta.str()+"_FlatTracksters.root").c_str());
    if (!tree[iF]) return;
    if (isFC){
      treeLC[iF]->AddFile((filePath[iF]+"/"+reg[iReg[iF]]+"/FineCalo/step3ticl_"+pteta.str()+"_run0_FlatTracksters.root").c_str());
      treeLC[iF]->AddFile((filePath[iF]+"/"+reg[iReg[iF]]+"/FineCalo/step3ticl_"+pteta.str()+"_run1_FlatTracksters.root").c_str());
      treeLC[iF]->AddFile((filePath[iF]+"/"+reg[iReg[iF]]+"/FineCalo/step3ticl_"+pteta.str()+"_run2_FlatTracksters.root").c_str());
      treeLC[iF]->AddFile((filePath[iF]+"/"+reg[iReg[iF]]+"/FineCalo/step3ticl_"+pteta.str()+"_run3_FlatTracksters.root").c_str());
      treeLC[iF]->AddFile((filePath[iF]+"/"+reg[iReg[iF]]+"/FineCalo/step3ticl_"+pteta.str()+"_run4_FlatTracksters.root").c_str());
      treeLC[iF]->AddFile((filePath[iF]+"/"+reg[iReg[iF]]+"/FineCalo/step3ticl_"+pteta.str()+"_run5_FlatTracksters.root").c_str());
      treeLC[iF]->AddFile((filePath[iF]+"/"+reg[iReg[iF]]+"/FineCalo/step3ticl_"+pteta.str()+"_run6_FlatTracksters.root").c_str());
      treeLC[iF]->AddFile((filePath[iF]+"/"+reg[iReg[iF]]+"/FineCalo/step3ticl_"+pteta.str()+"_run7_FlatTracksters.root").c_str());
      treeLC[iF]->AddFile((filePath[iF]+"/"+reg[iReg[iF]]+"/FineCalo/step3ticl_"+pteta.str()+"_run8_FlatTracksters.root").c_str());
      treeLC[iF]->AddFile((filePath[iF]+"/"+reg[iReg[iF]]+"/FineCalo/step3ticl_"+pteta.str()+"_run9_FlatTracksters.root").c_str());
    }
    //else if (isFC) treeLC[iF]->AddFile((filePath[iF]+"/FineCalo/step3ticl_"+pteta.str()+"_FlatTracksters.root").c_str());
    else treeLC[iF]->AddFile((filePath[iF]+"/"+reg[iReg[iF]]+"/step3ticl_"+pteta.str()+"_FlatTracksters.root").c_str());
    if (!treeLC[iF]) return;
    
    TChain *friendTree = new TChain("TruthTree");
    friendTree->AddFile((filePath[iF]+"/Truth/"+iter+"/"+pteta.str()+"/Debug/TruthTree_"+regShort[iReg[iF]]+".root").c_str());
    
    tree[iF]->AddFriend(friendTree);
    
  }
  
  if (evtNumbers[0].size() != evtNumbers[1].size()) {
    std::cout << " -- problem with event lists: "
	      <<evtNumbers[0].size()
	      << " " << evtNumbers[1].size()
	      << std::endl; 
    return;
  }
  
  const unsigned nE = evtNumbers[0].size();

  unsigned nSel = 0;
  
  for (unsigned iE(0);iE<nE;++iE){
    //if (nSel > 5) continue;
    //    if (!(evtNumbers[0][iE].nSC10==1 && evtNumbers[1][iE].nSC10==1 &&
    //	  evtNumbers[0][iE].nSC < 200 && evtNumbers[1][iE].nSC==1)
    //	){
    //    if ( (evtNumbers[0][iE].eFrac[0]+evtNumbers[0][iE].eFrac[1]-evtNumbers[1][iE].eFrac[0]-evtNumbers[1][iE].eFrac[1] > 0.8) ||
    //	 (evtNumbers[0][iE].eFrac[4]+evtNumbers[0][iE].eFrac[6]+evtNumbers[0][iE].eFrac[8]-evtNumbers[1][iE].eFrac[4]-evtNumbers[1][iE].eFrac[6]-evtNumbers[1][iE].eFrac[8] > 0.2) ||
    //	 (evtNumbers[0][iE].eFrac[3]+evtNumbers[0][iE].eFrac[5]+evtNumbers[0][iE].eFrac[7]-evtNumbers[1][iE].eFrac[3]-evtNumbers[1][iE].eFrac[5]-evtNumbers[1][iE].eFrac[7] > 0.5)
    //	 )
    //if (evtNumbers[0][iE].nBlobs == 3)
    //if (iE==2 || iE==990 || iE == 10 || iE==549 || iE==899 || iE == 306)
    //    if (iE==303 || iE==493 || iE==101 || iE==532  || iE==588  || iE==652  ||
    //	iE==868  || iE==164 || iE==50 ||
    //	iE==2 || iE==990 || iE == 10 || iE==549 || iE==899 || iE == 306)
	//evtNumbers[0][iE].eBlobs/evtNumbers[0][iE].eTotAllLC3 < 0.98)
    if (iE==101 || iE==102 || iE==124 || iE == 191 || iE == 202 || iE == 306
	|| iE == 127 || iE == 271
	|| iE == 109 || iE == 133 || iE == 130 || iE == 215
	|| iE == 222 || iE == 242 || iE == 273 || iE == 262)
    {
	nSel++;
	std::cout << " -- entry " << iE << std::endl;

	mycTest->cd();
	std::map<int,std::pair<double,double> > lAllPid;
	for (unsigned iF(0);iF<nF;++iF){
	  tree[iF]->Draw("abs(sc_pdgid)>>hpid(4000,0,4000)","sc_energy*(sc_energy>0)","",1,evtNumbers[iF][iE].event);
	  TH1F *hpid = (TH1F*)gDirectory->Get("hpid");
	  for (unsigned iB(1); iB<hpid->GetNbinsX()+1;++iB){
	    if (hpid->GetBinContent(iB)>0){
	      int absPid =  hpid->GetXaxis()->GetBinLowEdge(iB);
	      double E1 = iF==0 ? hpid->GetBinContent(iB) : 0;
	      double E2 = iF==1 ? hpid->GetBinContent(iB) : 0;
	      std::pair<std::map<int,std::pair<double,double> >::iterator,bool> lInsert = lAllPid.insert(std::pair<int,std::pair<double,double> >(absPid,std::pair<double,double>(E1,E2)));
	      if(!lInsert.second){
		if (iF==0) lInsert.first->second.first += E1;
		else if (iF==1) lInsert.first->second.second += E2;
	      }
	    }
	  }
	  hpid->Delete();
	}
	std::map<int,std::pair<double,double> >::iterator lIter = lAllPid.begin();
	for (; lIter != lAllPid.end();++lIter){
	  std::cout << std::setprecision(4) //<< lIter->first << " & "
		    << "$" << pdgMap.find(lIter->first)->second.second << "$ & "
		    << lIter->second.first << " & "
		    << lIter->second.second << " \\\\"
		    << std::endl;
	}
	
	
	for (unsigned iC(0); iC<nC;++iC){
	  
	  for (unsigned iF(0);iF<nF;++iF){
	    bool isPU = iReg[iF]==5;
	    //evtNumbers[iF][iE].Print();
	    char buf[500];
	    //  sprintf(buf,"%s %s %s",regShort[iReg[iF]].c_str(),iter.c_str(),pteta.str().c_str());
	    sprintf(buf,"%s %s evt#%d",regShort[iReg[iF]].c_str(),pteta.str().c_str(),evtNumbers[iF][iE].event);
	    
	    getDisplay(lPlotDir[iF].str(),evtNumbers[iF][iE],
		       tree[iF],treeLC[iF],
		       myc,mycT,mycP,buf,
		       hRZtrue[iF],hRZall[iF],hRZ[iF],hPID[iF],hE[iF],hCat[iF],
		       lCut[iC],lCutTS[iC],lCutSave[iC],
		       hRZs[iF],hXZ,hYZ,
		       isPU);
	  }
	}
      }
  }
  
  std::cout << " -- Found " << nSel << " events with decaying pions." << std::endl;
  
};

int quickPlot(){

  getDisplays();

  return 0;
  
  const unsigned nF = 2;
  std::string filePath[2] = {"D49_FineCalo","D49_DefSC"};


  const unsigned nS = 1;
  std::string cutStr[6] = {
    "",
    "nBlobs==0",
    "nBlobs==1",
    "nBlobs==2",
    "nBlobs==3",
    "nBlobs==4"
  };
  std::string cutShort[6] = {
    "NoSel",
    "nB0",
    "nB1",
    "nB2",
    "nB3",
    "nB4"
  };

  const bool doBoth = true;
  
  for (unsigned iF(0);iF < (doBoth?1:nF) ;++iF){
    for (unsigned iS(0);iS<nS;++iS){
      plotSelection(filePath[iF],cutStr[iS],cutShort[iS],doBoth);
    }
  }
  
  return 0;
}
