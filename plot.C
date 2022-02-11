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

//#include "quickPlot.C"
#include "globalVars.h"

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
  std::cout << "SAVE15" << std::endl;
  myc->Print((plotDir+"/"+varShort+cutShort+".pdf").c_str());
  //@@myc->Print((plotDir+"/"+varShort+cutShort+".png").c_str());


  return myhist;
};

class myPlot1D{
public:
  unsigned nPlots;
  std::vector<std::string> varNum;
  std::vector<std::string> varDen;
  std::vector<std::string> leglabel;
  std::string varSave;
  std::string histTitle;
  unsigned nBins;
  double binMin;
  double binMax;
  unsigned logy;
  std::string drawOpt;
  bool doReso;
  void print() {
    std::cout << "myPlot1D:   " << std::endl
	      << " nPlots:    " << nPlots << std::endl
	      << " histTitle: " << histTitle << std::endl
	      << " nBins:     " << nBins << std::endl
	      << " binMin:    " << binMin << std::endl
	      << " binMax:    " << binMax << std::endl
	      << " logy:      " << logy << std::endl
	      << " drawOpt:   " << drawOpt << std::endl;
  };
};


std::vector<myPlot1D> makeAllPlots(){

  std::string fileName = "plotList2.dat";

  std::vector<myPlot1D> lVec;
  std::ifstream infile(fileName);
  if (!infile.is_open()){
    std::cout << " -- file " << fileName << " not found." << std::endl;
    return lVec;
  }

  // nPlots { varNum varDen leglabel } varSave hisTitle nBins binMin binMax logy drawOpt doReso
  
  unsigned counter = 0;
  while(!infile.eof()){
    myPlot1D aPlot;
    aPlot.nPlots = 0;
    infile>>aPlot.nPlots;
    std::string tmpStr;
    for (unsigned i(0); i<aPlot.nPlots;++i){
      infile>>tmpStr;
      aPlot.varNum.push_back(tmpStr);
      infile>>tmpStr;
      aPlot.varDen.push_back(tmpStr);
      if (aPlot.nPlots>1){
	infile>>tmpStr;
	aPlot.leglabel.push_back(tmpStr);
      }
    }
    infile>>aPlot.varSave>>aPlot.histTitle>>aPlot.nBins>>aPlot.binMin>>aPlot.binMax>>aPlot.logy>>aPlot.drawOpt>>aPlot.doReso;
    if (aPlot.nPlots>0 && aPlot.varNum[0].find("#")==aPlot.varNum[0].npos) lVec.push_back(aPlot);
    counter++;
    if (counter > 10000) break;
  }
  
  std::cout << " -- Found " << lVec.size() << " elements in file " << fileName << std::endl;
  return lVec;
};


void overlayPlots(const std::vector<TH1F *> & htmp,
		  const std::vector<std::string> & leglabel,
		  TCanvas* & myc,
		  std::string plotDir,
		  char* buf,
		  std::string varShort,
		  std::string cut,
		  std::string cutShort,
		  const bool doLogy,
		  std::string drawOption){

  myc->cd();
  //gStyle->SetOptStat("eMRuo");
  gStyle->SetOptStat(0);
  gPad->SetLogy(doLogy);
  const unsigned nP = htmp.size();

  TLegend *leg = new TLegend(0.86,0.82,1.,0.92);
  for (unsigned iP(0); iP<nP; ++iP){
    std::cout << " Overlaying " <<  htmp[iP]->GetName() << " " << htmp[iP]->Integral() << std::endl;
    htmp[iP]->SetLineColor(myColor(iP));
    htmp[iP]->SetMarkerColor(myColor(iP));
    htmp[iP]->SetMarkerStyle(myMarker(iP));
    htmp[iP]->Draw(iP==0?drawOption.c_str() : (drawOption+"same").c_str());
    std::ostringstream histstat;
    histstat << ", Nsel=" << htmp[iP]->Integral();
    leg->AddEntry(htmp[iP],(leglabel[iP]+histstat.str()).c_str(),"P");
  }

  leg->Draw("same");
  TLatex lat;
  lat.DrawLatexNDC(0.01,0.95,buf);
  if (cut.size()>0){
    lat.SetTextSize(0.03);
    lat.DrawLatexNDC(0.5,0.95,cut.c_str());
  }
  
  myc->Update();
  std::cout << "SAVE11" << std::endl;
  myc->Print((plotDir+"/"+varShort+cutShort+".pdf").c_str());
  //@@myc->Print((plotDir+"/"+varShort+cutShort+".png").c_str());
  
  
};


//@@int plotEvtDisplayPU(TChain *treeLC,
//@@		     std::string iter,
//@@		     const unsigned iC,
//@@		     TCanvas* & myc,
//@@		     TCanvas* & mycAll,
//@@		     TH2F* & hRZtrue,
//@@		     TH2F* & hRZall,
//@@		     TH2F* & hRZ,
//@@		     TH2F* & hRZalldR,
//@@		     TH2F* & hRZdR,
//@@		     std::string plotDir,
//@@		     char* buf,
//@@		     const bool doAllLC,
//@@		     const bool doSeed,
//@@		     const Outlier & outlier,
//@@		     const bool isPU,
//@@		     const bool doVsLayer,
//@@		     const std::string puLabel
//@@		     ){
//@@
//@@  const bool doAMiters = (puLabel=="AMiters");
//@@  gStyle->SetMarkerSize(2);//1.2);
//@@  
//@@  const unsigned entry = outlier.event;
//@@  const bool isFC = true;
//@@
//@@  TLatex lat;
//@@  TExec *ex1 = new TExec("ex1","gStyle->SetPalette(1);");
//@@  
//@@  bool isSim = iter=="Sim";
//@@  if (isSim){
//@@    if (isPU) mycAll->cd(1);
//@@    else myc->cd(1);
//@@    double truthE = 0;
//@@    getTrueRZ(treeLC,hRZtrue,entry,isFC,truthE,doVsLayer);  
//@@    getRZ(treeLC,hRZall,entry,isFC,true,doSeed,false,isPU?"PU.":"",doVsLayer);
//@@    if (isPU) getRZ(treeLC,hRZalldR,entry,isFC,true,doSeed,true,isPU?"PU.":"",doVsLayer);
//@@    
//@@  
//@@    if (isPU) mycAll->cd(1);
//@@    else myc->cd(1);
//@@    gPad->SetLogz(1);
//@@    gPad->SetRightMargin(0.18);
//@@    
//@@    hRZall->GetZaxis()->SetRangeUser(0.01,100);
//@@    hRZall->Draw("colz");
//@@    ex1->Draw();
//@@    hRZall->Draw("colzsame");
//@@    
//@@    hRZtrue->SetLineColor(2);
//@@    hRZtrue->SetFillColor(2);
//@@    hRZtrue->SetFillStyle(1001);
//@@    hRZtrue->Draw("boxsame");
//@@    
//@@    lat.DrawLatexNDC(0.02,0.95,buf);
//@@    lat.DrawLatexNDC(0.15,0.84,"All LCs");	
//@@    if (isPU){
//@@      myc->cd(2);
//@@      gPad->SetLogz(1);
//@@      gPad->SetRightMargin(0.18);
//@@      
//@@      hRZalldR->GetZaxis()->SetRangeUser(0.01,100);
//@@      hRZalldR->Draw("colz");
//@@      ex1->Draw();
//@@      hRZalldR->Draw("colzsame");
//@@      
//@@      hRZtrue->Draw("boxsame");
//@@      
//@@      lat.DrawLatexNDC(0.02,0.95,buf);
//@@      lat.DrawLatexNDC(0.15,0.84,"All LCs #DeltaR(CP)<0.4");	
//@@    }
//@@
//@@    if (!isPU){
//@@      myc->cd(4);
//@@      gPad->SetLogz(1);
//@@      gPad->SetRightMargin(0.18);
//@@      bool inBlobs = true;
//@@
//@@      TH2F *hRZT[6];
//@@      double dreta[4] = {0.15/2,-0.15/2,0.25/2,-0.25/2};
//@@      for (unsigned ieta(0); ieta<4;++ieta){
//@@	hRZT[ieta]=0;
//@@	getTrueRZ(treeLC,hRZT[ieta],entry,isFC,truthE,doVsLayer,dreta[ieta]);
//@@      }
//@@      
//@@      getRZ(treeLC,hRZ,entry,isFC,true,doSeed,false,isPU?"PU.":"",doVsLayer,inBlobs);
//@@      hRZ->GetZaxis()->SetRangeUser(0.01,100);
//@@      //TExec *ex2 = new TExec("ex2","MyPalette();");
//@@      //ex2->Draw();
//@@
//@@      ex1->Draw();
//@@      hRZ->Draw("colzsame");
//@@
//@@      for (unsigned ieta(0); ieta<4;++ieta){
//@@	hRZT[ieta]->SetLineColor(kRed+ieta/2);
//@@	hRZT[ieta]->SetFillColor(kRed+ieta/2);
//@@	hRZT[ieta]->SetFillStyle(1001);
//@@	hRZT[ieta]->Draw("boxsame");
//@@      }
//@@      hRZ->Draw("colzsame");
//@@
//@@      
//@@      if (doSeed) lat.DrawLatexNDC(0.15,0.78,"Seed pos");	
//@@      else lat.DrawLatexNDC(0.15,0.78,"LC pos");	
//@@      lat.DrawLatexNDC(0.15,0.84,(iter+" LCs in Blobs").c_str());
//@@      lat.DrawLatexNDC(0.15,0.72,isPU?puLabel.c_str():isFC?"FineCalo":"Default");	 
//@@    }
//@@  }
//@@  if (!isPU) myc->cd(8+iC);
//@@  else mycAll->cd(iC==0?2:6+iC);
//@@  
//@@  gPad->SetLogz(1);
//@@  gPad->SetRightMargin(0.18);
//@@  
//@@  getRZ(treeLC,hRZ,entry,isFC,false,doSeed,false,isPU?"PU"+iter+".":iter+".",doVsLayer);
//@@  hRZ->GetZaxis()->SetRangeUser(0.01,100);
//@@  
//@@  //TExec *ex2 = new TExec("ex2","MyPalette();");
//@@  //ex2->Draw();
//@@  ex1->Draw();
//@@  
//@@  
//@@  
//@@  hRZ->Draw("colzsame");
//@@  
//@@  if (doSeed) lat.DrawLatexNDC(0.15,0.78,"Seed pos");	
//@@  else lat.DrawLatexNDC(0.15,0.78,"LC pos");	
//@@  lat.DrawLatexNDC(0.15,0.84,(iter+" LCs").c_str());
//@@  lat.DrawLatexNDC(0.15,0.72,isPU?puLabel.c_str():isFC?"FineCalo":"Default");	
//@@  //hRZtrue->Draw("boxsame");
//@@    
//@@  /*std::ostringstream counters;
//@@  counters << "nSC10=" << outlier.nSC10 
//@@	   << " nSC=" << outlier.nSC;
//@@  lat.DrawLatexNDC(0.15,0.66,counters.str().c_str());
//@@  counters.str("");
//@@  counters << "nTS10=" << outlier.nTS10 
//@@	   << " nTS=" << outlier.nTS;
//@@  lat.DrawLatexNDC(0.15,0.6,counters.str().c_str());    
//@@  //counters.str("");
//@@  //counters << "eTot LC (#geq 3 RH) = " << eLC3
//@@  //	     << " (" << eAllLC3 << ") GeV";
//@@  //lat.DrawLatexNDC(0.15,0.54,counters.str().c_str());    
//@@  */
//@@
//@@  if (!isPU){
//@@    myc->cd(3);
//@@    treeLC->Draw("RECO.eTotLC3>>hELC","","",1,entry);
//@@    TH1F *hELC = (TH1F*)gDirectory->Get("hELC");
//@@    double eLC3 = hELC->GetMean();
//@@    hELC->Delete();
//@@    
//@@    if (doAMiters) treeLC->Draw("PURECO.sumts_energy_EMAM+PURECO.sumts_energy_HADAM>>hELC","","",1,entry);
//@@    else treeLC->Draw("RECO.sumts_energy_EMAM+RECO.sumts_energy_HADAM>>hELC","","",1,entry);
//@@    
//@@    hELC = (TH1F*)gDirectory->Get("hELC");
//@@    double eLC3AM = hELC->GetMean();
//@@    hELC->Delete();
//@@
//@@    treeLC->Draw("RECO.eTotAllLC3>>hEAllLC","","",1,entry);
//@@    TH1F *hEAllLC = (TH1F*)gDirectory->Get("hEAllLC");
//@@    double eAllLC3 = hEAllLC->GetMean();
//@@    hEAllLC->Delete();
//@@    
//@@    treeLC->Draw("TRUTH.EBlobs>>hEblobs","","",1,entry);
//@@    TH1F *hEblobs = (TH1F*)gDirectory->Get("hEblobs");
//@@    double eBlobs = hEblobs->GetMean();
//@@    hEblobs->Delete();
//@@    
//@@    treeLC->Draw("TRUTH.eTotSCAtB>>hESC","","",1,entry);
//@@    TH1F *hESC = (TH1F*)gDirectory->Get("hESC");
//@@    double eSC = hESC->GetMean();
//@@    hESC->Delete();
//@@    
//@@    TH1F *hT = new TH1F("hT","",100,0,500);
//@@    treeLC->Draw("cp_energy[0]>>hT","","",1,entry);
//@@    double truthE = hT->GetMean();
//@@    hT->Delete();
//@@
//@@    treeLC->Draw("RECO.nTS_Trk>>hTrkTS","","",1,entry);
//@@    TH1F *hTrkTS = (TH1F*)gDirectory->Get("hTrkTS");
//@@    double nTrkTS = hTrkTS->GetMean();
//@@    hTrkTS->Delete();
//@@
//@@    lat.SetTextSize(0.07);
//@@    sprintf(buf,"E_{CP}=%3.3f GeV",truthE);
//@@    lat.DrawLatexNDC(0.15,0.9,buf);	
//@@    sprintf(buf,"E_{SC@B}=%3.3f GeV",eSC);
//@@    lat.DrawLatexNDC(0.15,0.83,buf);	
//@@    sprintf(buf,"E_{allLC3}=%3.3f GeV",eAllLC3);
//@@    lat.DrawLatexNDC(0.15,0.76,buf);	
//@@    sprintf(buf,"E_{Blobs}=%3.3f GeV",eBlobs);
//@@    lat.DrawLatexNDC(0.15,0.69,buf);	
//@@    sprintf(buf,"E_{RecoLC3}=%3.3f GeV",eLC3);
//@@    lat.DrawLatexNDC(0.15,0.62,buf);	
//@@    sprintf(buf,"%%Blobs=%3.2f, %%Reco=%3.2f",eBlobs/eAllLC3*100,eLC3/eAllLC3*100);
//@@    lat.DrawLatexNDC(0.1,0.55,buf);	
//@@    sprintf(buf,"Has Trk seed: %3.1f",nTrkTS);
//@@    lat.DrawLatexNDC(0.1,0.48,buf);	
//@@
//@@    sprintf(buf,"E_{AMLC3}=%3.3f GeV",eLC3AM);
//@@    lat.DrawLatexNDC(0.15,0.38,buf);	
//@@    sprintf(buf,"%%RecoAM=%3.2f",eLC3AM/eAllLC3*100);
//@@    lat.DrawLatexNDC(0.1,0.3,buf);
//@@
//@@    lat.SetTextSize(0.05);
//@@
//@@    myc->cd(5);
//@@    gPad->SetRightMargin(0.18);
//@@    gPad->SetLogz(1);
//@@    std::ostringstream lname;
//@@    lname.str("");
//@@    lname << "hE_" << entry;
//@@    if (isFC) lname << "_FC";
//@@    int nX = 9;
//@@    int nY = 10;
//@@    double binsID[10] = {0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5};
//@@    double binsE[11] = {0,1,2,5,10,20,50,100,150,200,55*cosh(2.1)}; 
//@@    TH2F *hE = new TH2F(lname.str().c_str(),";PDG cat;E_{SC} (GeV);nSC",nX,binsID,nY,binsE);
//@@    hE->GetYaxis()->SetRangeUser(0,binsE[nY]);
//@@    treeLC->Draw(("sc_energyAtB_eCut0:sc_pdgcat_eCut0>>"+lname.str()).c_str(),"","colz",1,entry);
//@@    for (unsigned iY(1);iY<nY+1;++iY){
//@@      hE->SetBinContent(nX,iY,hE->GetBinContent(nX,iY)+hE->GetBinContent(nX+1,iY));
//@@      hE->SetBinContent(1,iY,hE->GetBinContent(0,iY)+hE->GetBinContent(1,iY));
//@@    }
//@@    for (unsigned iX(1);iX<nX+1;++iX){
//@@      hE->GetXaxis()->SetBinLabel(iX,pdgBinLabel(iX).c_str());
//@@      hE->SetBinContent(iX,1,hE->GetBinContent(iX,0)+hE->GetBinContent(iX,1));
//@@      hE->SetBinContent(iX,nY,hE->GetBinContent(iX,nY)+hE->GetBinContent(iX,nY+1));
//@@    }
//@@    
//@@    ex1->Draw();
//@@    hE->Draw("colz same");
//@@    
//@@  }
//@@  
//@@
//@@
//@@  if (isPU){
//@@    myc->cd(15+iC);
//@@    gPad->SetLogz(1);
//@@    gPad->SetRightMargin(0.18);
//@@    
//@@    getRZ(treeLC,hRZdR,entry,isFC,false,doSeed,true,isPU?"PU"+iter+".":iter+".",doVsLayer);
//@@    hRZdR->GetZaxis()->SetRangeUser(0.01,100);
//@@    
//@@    ex1->Draw();
//@@  
//@@    hRZdR->Draw("colzsame");
//@@    
//@@    if (doSeed) lat.DrawLatexNDC(0.15,0.78,"Seed pos");	
//@@    else lat.DrawLatexNDC(0.15,0.78,"LC pos");	
//@@    lat.DrawLatexNDC(0.15,0.84,(iter+" LCs #DeltaR(CP)<0.4").c_str());
//@@    lat.DrawLatexNDC(0.15,0.72,isPU?puLabel.c_str():isFC?"FineCalo":"Default");	
//@@  }
//@@  
//@@  myc->Update();
//@@  mycAll->Update();
//@@  if (isPU && iC==6){
//@@    std::ostringstream lSave;
//@@    lSave << plotDir << "/" << "R";
//@@    if (doVsLayer) lSave << "L";
//@@    else lSave << "Z";
//@@    lSave << "_entry" << entry;
//@@    //if (doAllLC) lSave << "_allLC";
//@@    if (doSeed) lSave << "_seedInfo";
//@@    lSave << "_allIters";// << iter;
//@@    std::cout << "SAVE12" << std::endl;
//@@    myc->Print((lSave.str()+".pdf").c_str());
//@@    //@@myc->Print((lSave.str()+".png").c_str());
//@@    lSave << "_allLCs";// << iter;
//@@    std::cout << "SAVE13" << std::endl;
//@@    mycAll->Print((lSave.str()+".pdf").c_str());
//@@    //@@mycAll->Print((lSave.str()+".png").c_str());
//@@  }
//@@
//@@  return 0; 
//@@};

int getTree(std::string filePath,
	    std::string pteta,
	    TChain *& treeLC){
  
  //@@const unsigned nRuns = 10;
  const unsigned nRuns = 1;
  //@@const unsigned nF = 2;
  const unsigned nF = 1;
  //@@const bool doAMiters = true;
  const bool doAMiters = false;
  //@@const bool isDebug = false;
  const bool isDebug = true;
  //std::string regs[nF] = {"/ChargedPionsFromVtx","/ChargedPionsFromVtxWithPU"};
  //@@std::string regs[2] = {"/ElectronsFromVtx","/ElectronsFromVtx"};
  std::string regs[1] = {"/CloseByPhotons"};
  //@@unsigned iRegion = 6;
  unsigned iRegion = 0;
  if (doAMiters && iRegion==4) regs[1] = "_AMiters/ChargedPionsFromVtx";
  
  //@@const unsigned nIter = 7;
  const unsigned nIter = 5;
  std::string iters[2][nIter] = {
    //{"Sim","TrkEM","EM","Trk","HAD","EMAM","HADAM"},
    //{"Sim","TrkEM","EM","Trk","HAD","EMAM","HADAM"}
    //@@{"Sim","TrkEM","EM","Trk","HAD","Dummy1","EM3"},
    //@@{"Sim","TrkEM","EM","Trk","HAD","Dummy1","EM3"}
    {"TrkEM","EM","Trk","HAD","Sim"},
    {"TrkEM","EM","Trk","HAD","Sim"}
  };
  if (!doAMiters){
    for (unsigned iC(0); iC<nIter;++iC){
      iters[1][iC] = iters[0][iC];
    }
  }
  //std::string filePath = "D49_AllTracksters";
  

  TChain *tree[2][nIter+1];
  for (unsigned iF(0); iF<nF;++iF){
    if (iF!=0) tree[iF][nIter] = new TChain("ticlTree/treeLC");
    for (unsigned iR(0); iR<nRuns; ++iR){
      std::ostringstream infile;
      infile << filePath << regs[iF] << "/";
      if (filePath.find("FineCalo")!=filePath.npos) infile <<  "FineCalo/";
      //@@infile << "step3ticl_" << pteta
      infile << "step3_" << pteta
	     << "_run" << iR
	     << "_FlatTracksters.root";
      std::cout << "FILE for ticlTree/treeLC: " << infile.str() << std::endl;
      if (iF==0) treeLC->AddFile(infile.str().c_str());
      else tree[iF][nIter]->AddFile(infile.str().c_str());
    }
    if (iF==1) treeLC->AddFriend(tree[iF][nIter],"PU");
  }
  for (unsigned iC(0); iC<nIter;++iC){
    for (unsigned iF(0); iF<nF;++iF){
      tree[iF][iC] = new TChain(("ticlTree/TSTree_"+iters[iF][iC]).c_str());
      for (unsigned iR(0); iR<nRuns; ++iR){
	std::ostringstream infile;
	infile << filePath << regs[iF] << "/";
	if (filePath.find("FineCalo")!=filePath.npos) infile <<  "FineCalo/";
	//@@infile << "step3ticl_" << pteta
	infile << "step3_" << pteta
	       << "_run" << iR
	       << "_FlatTracksters.root";
	std::cout << "FILE for ticlTree/TSTree: " << infile.str() << std::endl;
	tree[iF][iC]->AddFile(infile.str().c_str());
      }
      if (iF==0) treeLC->AddFriend(tree[iF][iC],iters[iF][iC].c_str());
      else treeLC->AddFriend(tree[iF][iC],("PU"+iters[iF][iC]).c_str());
    }
  }

//  TChain *friendTree = new TChain("TruthTree");
////@@  TChain *friendTree2 = new TChain("TruthTree");
//  std::string tmp = filePath+"/Truth/"+iters[0][0]+"/"+pteta+"/"+(isDebug?"Debug/":"")+"/TruthTree_"+regShort[iRegion]+".root";
//  std::cout << "FILE for TruthTree: " << tmp << std::endl;
//  friendTree->AddFile(tmp.c_str());
////@@  if (!doAMiters) friendTree2->AddFile((filePath+"/Truth/"+iters[0][0]+"/"+pteta+"/"+(isDebug?"Debug/":"")+"/TruthTree_"+regShort[iRegion+1]+".root").c_str());
//  
//  treeLC->AddFriend(friendTree,"TRUTH");
////@@  if (!doAMiters) treeLC->AddFriend(friendTree2,"PUTRUTH");

  TChain *friendTree3 = new TChain("RecoTree");
//@@  TChain *friendTree4 = new TChain("RecoTree");
  std::string tmp3 = filePath+"/Reco/"+pteta+"/"+(isDebug?"Debug/":"")+"/RecoTree_"+regShort[iRegion]+".root";
  std::cout << "FILE for RecoTree: " << tmp3 << std::endl;
  friendTree3->AddFile(tmp3.c_str());
//@@  if (!doAMiters) friendTree4->AddFile((filePath+"/Reco/"+pteta+"/"+(isDebug?"Debug/":"")+"/RecoTree_"+regShort[iRegion+1]+".root").c_str());
//@@  else if (iRegion==4) friendTree4->AddFile((filePath+"_AMiters/Reco/"+pteta+"/"+(isDebug?"Debug/":"")+"/RecoTree_"+regShort[iRegion]+".root").c_str());
//@@  else if (iRegion==6) friendTree4->AddFile((filePath+"/Reco/"+pteta+"/"+(isDebug?"Debug/":"")+"/RecoTree_"+regShort[iRegion]+".root").c_str());
  
  treeLC->AddFriend(friendTree3,"RECO");
//@@  treeLC->AddFriend(friendTree4,"PURECO");

//@@  TChain *friendTree5 = new TChain("RecoTree");
//@@  TChain *friendTree6 = new TChain("RecoTree");
//@@  friendTree5->AddFile((filePath+"/Reco/"+pteta+"/"+(isDebug?"Debug/":"")+"/RecoTreeAM_"+regShort[iRegion]+".root").c_str());
//@@  if (!doAMiters) friendTree6->AddFile((filePath+"/Reco/"+pteta+"/"+(isDebug?"Debug/":"")+"/RecoTreeAM_"+regShort[iRegion+1]+".root").c_str());
//@@  else if (iRegion==4) friendTree6->AddFile((filePath+"_AMiters/Reco/"+pteta+"/"+(isDebug?"Debug/":"")+"/RecoTreeAM_"+regShort[iRegion]+".root").c_str());
//@@  else if (iRegion==6) friendTree6->AddFile((filePath+"/Reco/"+pteta+"/"+(isDebug?"Debug/":"")+"/RecoTreeAM_"+regShort[iRegion]+".root").c_str());
//@@  
//@@  treeLC->AddFriend(friendTree5,"RECOAM");
//@@  treeLC->AddFriend(friendTree6,"PURECOAM");
  
  
//@@  std::cout << "TEST2: " << treeLC << " " << bool(bool(treeLC==NULL)==true) << std::endl;
//@@  if ((bool(treeLC==NULL)==true)) {
//@@    std::cout << "HERE2" << std::endl;
//@@    std::cout << " problem getting tree with friends." << std::endl;
//@@    return 1;
//@@  }

  return 0;

};

void getResolutions(TChain * treeLC,
		    std::string pteta,
		    std::string plotDir,
		    std::string cutShort,
		    TH1F* & histE,
		    std::string labelE){

  TCanvas *mycE = new TCanvas("mycE","mycE",1);
  TLatex lat;
  char buf[200];
	       
  mycE->cd();
  gStyle->SetOptStat("eMRuo");
  gPad->SetLogy(1);
  histE->SetLineColor(1);
  histE->SetMarkerColor(1);
  histE->SetMarkerStyle(20);
  histE->Draw("PE");
    
  gStyle->SetOptStat("eMR");
  gStyle->SetOptFit(1111);
  
  histE->Fit("gaus","same");
  TF1 *fitFunc = (TF1*)histE->GetFunction("gaus");
  if (!fitFunc) return;

  lat.SetTextSize(0.05);
  lat.SetTextColor(1);
  sprintf(buf,"Hist #sigma/E=%3.3f",histE->GetRMS()/histE->GetMean());
  lat.DrawLatexNDC(0.15,0.85,buf);	
  
  sprintf(buf,"Fit #sigma/E=%3.3f",fitFunc->GetParameter(2)/fitFunc->GetParameter(1));
  lat.DrawLatexNDC(0.15,0.78,buf);	

  sprintf(buf,"%s %s",regShort[4].c_str(),pteta.c_str());
  lat.DrawLatexNDC(0.01,0.95,buf);	
  
  mycE->Update();
  std::cout << "SAVE14" << std::endl;
  mycE->Print((plotDir+"/Reso_"+labelE+".pdf").c_str());
  //@@mycE->Print((plotDir+"/Reso_"+labelE+".png").c_str());
  
};

int plotSelectionPU(TChain * treeLC,
		    std::string pteta,
		    std::string filePath,
		    std::string cutStr,
		    std::string cutShort){
  
  SetTdrStyle();

  //@@const unsigned nIter = 7;
  const unsigned nIter = 5;
  //@@std::string iters[nIter] = {"Sim","TrkEM","EM","Trk","HAD","EMAM","HADAM"};
  std::string iters[nIter] = {"TrkEM","EM","Trk","Had","Sim"};

  std::ostringstream lPlotDir;
  lPlotDir << filePath<<"/CompaPU/" << pteta << "/Plots/";//_nL10/";
  
  if (system(("mkdir -p "+lPlotDir.str()).c_str())){
    std::cout << " -- Cannot create output dir..." << lPlotDir.str() << std::endl;
    system(("echo \"-- return value \"$?"));
    return 1;
  }

  if (cutShort.size()>0) cutShort = "_"+cutShort;
  TCanvas *myc = new TCanvas("myc","myc",1);
  
  char buf[500];
  
  //  sprintf(buf,"%s %s %s",regShort[iReg].c_str(),iter.c_str(),pteta.str().c_str());
  //@@sprintf(buf,"%s %s",regShort[4].c_str(),pteta.c_str());
  sprintf(buf,"%s %s",regShort[0].c_str(),pteta.c_str());

  std::vector<myPlot1D> lVec = makeAllPlots();
  for (unsigned iV(0); iV<lVec.size(); ++iV){
    myPlot1D aP = lVec[iV];
    aP.print();
    std::vector<TH1F *> htmp;
    htmp.resize(aP.nPlots,0);
    for (unsigned iP(0); iP<aP.nPlots;++iP){
      
      htmp[iP] = makePlot(treeLC,myc,lPlotDir.str(),buf,
			  iP,aP.varNum[iP]+"/"+aP.varDen[iP],aP.varSave,
			  cutStr+(cutStr.size()>0?" && ":"")+aP.varDen[iP]+">0",cutShort,
			  aP.histTitle,
			  aP.nBins,aP.binMin,aP.binMax,aP.logy,
			  aP.drawOpt);
      //@@if (aP.doReso) getResolutions(treeLC,pteta,lPlotDir.str(),cutShort,htmp[iP],aP.varSave);
    }
    if (aP.nPlots>1){
      overlayPlots(htmp,aP.leglabel,myc,lPlotDir.str(),buf,
		   aP.varSave,cutStr+(cutStr.size()>0?" && ":"")+aP.varDen[0]+">0",
		   cutShort,
		   aP.logy,"PE");
    }
    
  }
  
  return 0;
};

//@@int getDisplayPU(std::string plotDir,
//@@		 const std::string iter,
//@@		 const unsigned iC,
//@@		 const Outlier & outlier,
//@@		 TChain* treeLC,
//@@		 TCanvas* & myc,
//@@		 TCanvas* & mycAll,
//@@		 char* buf,
//@@		 TH2F * & hRZtrue,
//@@		 TH2F * & hRZall,
//@@		 TH2F * & hRZ,
//@@		 TH2F * & hRZalldR,
//@@		 TH2F * & hRZdR,
//@@		 const bool isPU,
//@@		 const bool doVsLayer,
//@@		 const std::string puLabel
//@@		 ){
//@@  
//@@  SetTdrStyle();
//@@  gStyle->SetOptStat(0);
//@@  const unsigned entry = outlier.event;
//@@
//@@  plotEvtDisplayPU(treeLC,
//@@		   iter,iC,
//@@		   myc,mycAll,
//@@		   hRZtrue,hRZall,
//@@		   hRZ,
//@@		   hRZalldR,hRZdR,
//@@		   plotDir,
//@@		   buf,
//@@		   false,true,
//@@		   outlier,
//@@		   isPU,
//@@		   doVsLayer,
//@@		   puLabel);
//@@  return 0;
//@@};


//@@void getDisplaysPU(TChain *treeLC,
//@@		   std::string filePath,
//@@		   std::string pteta,
//@@		   const bool doVsLayer,
//@@		   const std::string puLabel
//@@		   ){
//@@  
//@@  std::map<int,std::pair<int,std::string> > pdgMap = getPDGMap("pdgList.dat");
//@@  
//@@  
//@@  TCanvas *myc = new TCanvas("mycD","mycD",1);
//@@  TCanvas *mycAll = new TCanvas("mycA","mycA",1);
//@@  const unsigned nIter = 7;
//@@  std::string iters[2][nIter] = {
//@@    {"Sim","TrkEM","EM","Trk","HAD","EMAM","HADAM"},
//@@    {"Sim","TrkEM","EM","Trk","HAD","EMAM","HADAM"}
//@@  };
//@@
//@@  const bool doAMiters = (puLabel=="AMiters");
//@@  if (!doAMiters){
//@@    for (unsigned iC(0); iC<nIter;++iC){
//@@      iters[1][iC] = iters[0][iC];
//@@    }
//@@  }
//@@  
//@@  //@@const unsigned nF = 2;
//@@  const unsigned nF = 1;
//@@  
//@@  TH2F* hRZtrue[nF];
//@@  TH2F* hRZall[nF];
//@@  TH2F* hRZalldR[nF];
//@@  TH2F* hRZ[nF];
//@@  TH2F* hRZdR[nF];
//@@  
//@@  std::vector<Outlier> evtNumbers[nF];
//@@  std::ostringstream lPlotDir[nF];
//@@  //@@unsigned iReg[2] = {4,5};
//@@  unsigned iReg[1] = {0};
//@@  
//@@  for (unsigned iF(0);iF<nF;++iF){
//@@    hRZtrue[iF] = 0;
//@@    hRZall[iF] = 0;
//@@    hRZ[iF] = 0;
//@@    hRZalldR[iF] = 0;
//@@    hRZdR[iF] = 0;
//@@    
//@@    bool isFC = filePath.find("FineCalo")!=filePath.npos;
//@@    bool isPU = iReg[iF]==5;
//@@    
//@@    getOutliers(filePath+"/Truth/Sim/pt50_eta21/Debug/eventsToPrint_"+regShort[iReg[iF]]+".dat",
//@@		evtNumbers[iF]);
//@@    
//@@    lPlotDir[iF] << filePath << "/Reco/" << pteta;
//@@    lPlotDir[iF] << "/Displays/";
//@@    if (isPU) lPlotDir[iF] << puLabel << "/";
//@@    
//@@    if (system(("mkdir -p "+lPlotDir[iF].str()).c_str())){
//@@      std::cout << " -- Cannot create output dir..." << lPlotDir[iF].str() << std::endl;
//@@      system(("echo \"-- return value \"$?"));
//@@      return;
//@@    }
//@@    
//@@  }
//@@  
//@@  if (evtNumbers[0].size() != evtNumbers[1].size()) {
//@@    std::cout << " -- problem with event lists: "
//@@	      <<evtNumbers[0].size()
//@@	      << " " << evtNumbers[1].size()
//@@	      << std::endl; 
//@@    return;
//@@  }
//@@  
//@@  const unsigned nE = evtNumbers[0].size();
//@@
//@@  unsigned nSel = 0;
//@@  
//@@  for (unsigned iE(0);iE<nE;++iE){
//@@    //if (nSel > 5) continue;
//@@    //    if (!(evtNumbers[0][iE].nSC10==1 && evtNumbers[1][iE].nSC10==1 &&
//@@    //	  evtNumbers[0][iE].nSC < 200 && evtNumbers[1][iE].nSC==1)
//@@    //	){
//@@    //    if ( (evtNumbers[0][iE].eFrac[0]+evtNumbers[0][iE].eFrac[1]-evtNumbers[1][iE].eFrac[0]-evtNumbers[1][iE].eFrac[1] > 0.8) ||
//@@    //	 (evtNumbers[0][iE].eFrac[4]+evtNumbers[0][iE].eFrac[6]+evtNumbers[0][iE].eFrac[8]-evtNumbers[1][iE].eFrac[4]-evtNumbers[1][iE].eFrac[6]-evtNumbers[1][iE].eFrac[8] > 0.2) ||
//@@    //	 (evtNumbers[0][iE].eFrac[3]+evtNumbers[0][iE].eFrac[5]+evtNumbers[0][iE].eFrac[7]-evtNumbers[1][iE].eFrac[3]-evtNumbers[1][iE].eFrac[5]-evtNumbers[1][iE].eFrac[7] > 0.5)
//@@    //	 )
//@@    //if (evtNumbers[0][iE].nBlobs == 3)
//@@    //if (iE==2){// || iE==990 || iE == 10 || iE==549 || iE==899 || iE == 306)
//@@    if (iE==101 || iE==202 || iE==306 || iE==102  || iE==130  || iE==215  ||
//@@	iE==222  || iE==262
//@@      //|| iE==2 || iE==990 || iE == 10 || iE==549 || iE==899 || iE == 306 ||
//@@      //iE==960 || iE==630
//@@	//|| evtNumbers[0][iE].eBlobs/evtNumbers[0][iE].eTotAllLC3 < 0.95
//@@	)
//@@      {
//@@	nSel++;
//@@	std::cout << " -- entry " << iE << std::endl;
//@@
//@@	myc->Clear();
//@@	myc->Divide(7,3);
//@@	mycAll->Clear();
//@@	mycAll->Divide(6,2);
//@@	for (unsigned iC(0); iC<nIter;++iC){
//@@	  for (unsigned iF(0);iF<nF;++iF){
//@@	    bool isPU = iReg[iF]==5;
//@@	    //evtNumbers[iF][iE].Print();
//@@	    char buf[500];
//@@	    //  sprintf(buf,"%s %s %s",regShort[iReg[iF]].c_str(),iter.c_str(),pteta.c_str());
//@@	    sprintf(buf,"%s %s evt#%d",regShort[iReg[iF]].c_str(),pteta.c_str(),evtNumbers[iF][iE].event);
//@@	    
//@@	    getDisplayPU(lPlotDir[iF].str(),iters[iF][iC],iC,
//@@			 evtNumbers[iF][iE],
//@@			 treeLC,
//@@			 myc,mycAll,buf,
//@@			 hRZtrue[iF],hRZall[iF],hRZ[iF],
//@@			 hRZalldR[iF],hRZdR[iF],
//@@			 isPU,doVsLayer,
//@@			 puLabel);
//@@	    
//@@	  }
//@@	}
//@@      }
//@@  }
//@@  
//@@  std::cout << " -- Found " << nSel << " events with decaying pions." << std::endl;
//@@  
//@@};

int plot(){

  gROOT->SetBatch(1);
  
  //std::string filePath = "D49_FineCalo_AMiters";
  //@@std::string filePath = "D49_DefSC";
  std::string filePath = "../TiCLTreeProducer/D49_FineCalo";

  //@@const unsigned nPT = 12;
  const unsigned nPT = 1;
  //@@double ptval[12] = {3,5,10,15,20,30,40,50,75,100,150,200};
  double ptval[1] = {50};
  const unsigned etaval = 21;

  for (unsigned ipt(0); ipt<nPT; ++ipt){  
    
    std::ostringstream pteta;
    pteta << "pt" << ptval[ipt] << "_eta" << etaval;
    
    TChain *treeLC = new TChain("ticlTree/treeLC");  
    //@@std::cout << "TEST: " << treeLC << std::endl;
    
    getTree(filePath,pteta.str(),treeLC);
    //@@return 0;
    
    //bool doVsLayer = true;
    //std::string puLabel = "WithPU";
    //std::string puLabel = "AMiters";
    //getDisplaysPU(treeLC,filePath,pteta.str(),doVsLayer,puLabel);
    //return 0;
    
    const unsigned nS = 1;
    std::string cutStr[5] = {
      "",
      "TRUTH.eFracAtB_catpipm>=0.99",
      "TRUTH.eFracAtB_catpipm<0.99",
      "(TRUTH.eFracAtB_cate+TRUTH.eFracAtB_catgamma)>0.5",
      "(TRUTH.eFracAtB_catn+TRUTH.eFracAtB_catneHad)>0.1"
      //"TRUTH.nSC_eCut10==0",
      //"TRUTH.nSC_eCut10==1",
      //"TRUTH.nSC_eCut10==2",
      //"TRUTH.nSC_eCut10==3",
      //"TRUTH.nSC_eCut10>=4"
    };
    std::string cutShort[5] = {
      "NoSel",
      "Enriched_pipm",
      "Enriched_notallpipm",
      "Enriched_egamma",
      "Enriched_neHad"
      //    "nSC0",
    };
    
    
    for (unsigned iS(0);iS<nS;++iS){
      plotSelectionPU(treeLC,pteta.str(),filePath,cutStr[iS],cutShort[iS]);
    }

  }//loop over pT vals
  
  return 0;
}
