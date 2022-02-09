#include <algorithm>
#include <iomanip>
#include <stdlib.h>
#include <iostream>
#include <fstream>

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

void merge(std::string filePath,
	   const unsigned iR,
	   const std::string aIter,
	   const std::string pteta,
	   const std::string phiStr,
	   const bool doDebug = false,
	   const unsigned phiMin = 0,
	   const unsigned phiMax = 360,
	   const float etaMin = 1.5,
	   const float etaMax = 3){
  
  SetTdrStyle();
  
  std::string partName = iR<2? "#gamma" : iR==4 || iR==5 ? "#pi^{#pm}" : "e^{#pm}";
  
  //double zposCEH = 364;//cm

  
  std::ostringstream lPlotDir;
  lPlotDir << filePath << "/Merged/" << aIter << "/" << pteta << phiStr << "/";
  if (doDebug) lPlotDir << "Debug";
  
  if (system(("mkdir -p "+lPlotDir.str()).c_str())){
    std::cout << " -- Cannot create output dir..." << lPlotDir.str() << std::endl;
    system(("echo \"-- return value \"$?"));
    return;
  }
  
  std::ofstream fileOut;
  fileOut.open((lPlotDir.str()+"/events_"+regShort[iR]+".dat").c_str());
  std::ofstream fileOutRes;
  
  fileOutRes.open((lPlotDir.str()+"/resolutions_"+regShort[iR]+".dat").c_str());
  
  TFile *fout = TFile::Open((lPlotDir.str()+"/Histos_all_"+regShort[iR]+".root").c_str(),"RECREATE");
  if (!fout) return;
  

  const unsigned nL = 49;//layers
  const unsigned nLEE = 30;//layers
  const unsigned nV1D = 23;
  const unsigned nV2D = 8;
  const unsigned nV = nV1D+nV2D;

  TCanvas *myc[nV];
  for (unsigned iL(0); iL<nV; ++iL){//loop on variables
    std::ostringstream label;
    label << "myc_" << iL;
    myc[iL] = new TCanvas(label.str().c_str(),label.str().c_str(),1);
  }


  int color[12] = {1,2,3,4,6,7,8,9,kRed+2,kCyan+2,kGray,kOrange};


  TLatex lat;
  char buf[200];
  fout->cd();
  TH1F *hist[nV1D];
  TH2F *hist2D[nV2D];
  std::string label[nV] = {"EtotOverECP_1RH","E1LCOverECP_1RH",
			   "EtotOverECP_2RH","E1LCOverECP_2RH",
			   "EtotOverECP_3RH","E1LCOverECP_3RH",
			   "ETSOverEtot","ETSOverECP","E1LCTSOverECP",
			   "nTS",
			   "ETSmaxEOverECP","SumETSOverECP",
			   "minDR","minDEta","minDPhi",
			   "ETSmaxEOverEtot",
			   "TSfirstLayer","TSlastLayer","TSlength",
			   "TSmaxEfirstLayer","TSmaxElastLayer","TSmaxElength",
			   "cp_missingEnergyFraction",
  			   "DRvsEfrac","DEtavsEfrac","DPhivsEfrac",
			   "TSidxmindRvsTSidxmaxE",
			   "mindRLCCP_1RH","mindRLCCP_2RH","mindRLCCP_3RH",
			   "mindRLCCPTS"
  };
  
  unsigned nBins[nV1D];
  for (unsigned iV(0); iV<nV1D; ++iV){//loop on variables
    nBins[iV] = 60;
  };
  nBins[9] = 10;//nTS
  
  double binMin[nV1D] = {
    0.5,0.5,0.5,0.5,0.5,0.5,
    0.,0.5,0.5,
    0,0.5,0.5,0,0,0,
    0,0,0,0,0,0,0,0
  };
  double binMax[nV1D] = {
    1.5,1.5,1.5,1.5,1.5,1.5,
    1.05,1.5,1.5,
    10,1.5,1.5,0.3,0.3,0.3,
    1.05,60,60,60,60,60,60,1.1
  };

  unsigned nBinsX[nV2D];
  for (unsigned iV(0); iV<nV2D; ++iV){//loop on variables
    nBinsX[iV] = 50;
    if (iV>=4) nBinsX[iV] = nL;
  };
  if (iR<4) nBinsX[7] = nLEE;

  double binXmin[nV2D] = {0,0,0,-1,1,1,1,1};
  double binXmax[nV2D] = {1.1,1.1,1.1,10,nL+1,nL+1,nL+1,nL+1};
  if (iR<4) binXmax[7] = nLEE+1;
  unsigned nBinsY[nV2D];
  for (unsigned iV(0); iV<nV2D; ++iV){//loop on variables
    nBinsY[iV] = 50;
  };
  double binYmin[nV2D] = {0,0,0,-1,0,0,0,0};
  double binYmax[nV2D] = {0.5,0.5,0.5,10,0.1,0.1,0.1,0.1};

  nBinsX[3] = binXmax[3]-binXmin[3];
  nBinsY[3] = binYmax[3]-binYmin[3];

  std::string axisLabel[nV] = {
    ";E_{tot}^{#geq 1RH}/E_{CP};"+partName,
    ";E_{tot}^{#geq 1RH, 1LC}/E_{CP};"+partName,
    ";E_{tot}^{#geq 2RH}/E_{CP};"+partName,
    ";E_{tot}^{#geq 2RH, 1LC}/E_{CP};"+partName,
    ";E_{tot}^{#geq 3RH}/E_{CP};"+partName,
    ";E_{tot}^{#geq 3RH, 1LC}/E_{CP};"+partName,
    ";E_{TS}^{min#Delta R}/E_{tot};"+partName,
    ";E_{TS}^{min#Delta R}/E_{CP};"+partName,
    ";E_{TS}^{min#Delta R,1LC}/E_{CP};"+partName,
    ";nTS;"+partName,
    ";E_{TS}^{maxE}/E_{CP};"+partName,
    ";#SigmaE_{TS}/E_{CP};"+partName,
    ";min#DeltaR(CP,TS-LC);"+partName,
    ";min#Delta#eta(CP,TS-LC);"+partName,
    ";min#Delta#phi(CP,TS-LC);"+partName,
    ";E_{TS}^{maxE}/E_{tot};"+partName,
    ";TS first layer;Tracksters",
    ";TS last layer;Tracksters",
    ";TS length;Tracksters",
    ";TS^{maxE} first layer;"+partName,
    ";TS^{maxE} last layer;"+partName,
    ";TS^{maxE} length;"+partName,
    ";|E_{CP}-E_{LCrechitsMatched}|/E_{CP};"+partName,    
    ";E_{TS}/E_{CP};#DeltaR(CP,TS-LC);Tracksters",
    ";E_{TS}/E_{CP};#Delta#eta(CP,TS-LC);Tracksters",
    ";E_{TS}/E_{CP};#Delta#phi(CP,TS-LC);Tracksters",
    ";TSidx (maxE); TSidx (min#Delta R);"+partName,
    ";HGC layer;min #DeltaR (LC #geq 1RH,CP);"+partName,
    ";HGC layer;min #DeltaR (LC #geq 2RH,CP);"+partName,
    ";HGC layer;min #DeltaR (LC #geq 3RH,CP);"+partName,
    ";HGC layer;TS min #DeltaR (LC,CP);"+partName
  };
  
  const bool dology[nV] = {
    1,1,1,1,1,1,
    1,1,1,1,
    1,1,1,1,1,1,
    1,1,1,1,1,1,1,
    0,0,0,0,0,0,0,0
  };
  double maxY[nV];
  double minY[nV];
  
  for (unsigned iV(0); iV<nV; ++iV){//loop on variables
    maxY[iV] = 0;
    minY[iV] = 1000;
  }

  TChain *treeLC = new TChain("ticlTree/treeLC");
  TChain *tree = new TChain(("ticlTree/TSTree_"+aIter).c_str());
  
  //E[iR] = ptval[ipt]*cosh(etaval[ieta]/10.);

  //tree->AddFile((filePath+"/"+reg[iR]+"/step3ticl_"+pteta+"_FlatTracksters.root").c_str());
  
  //treeLC->AddFile((filePath+"/"+reg[iR]+"/step3ticl_"+pteta+"_FlatTracksters.root").c_str());
  tree->AddFile((filePath+"/Photons/step3ticl_"+pteta+"_FlatTracksters.root").c_str());
  
  treeLC->AddFile((filePath+"/Photons/step3ticl_"+pteta+"_FlatTracksters.root").c_str());
  
  
  if (!tree) return;
  if (!treeLC) return;

  
  unsigned nPhotons = tree->GetEntries();
  std::cout << " -- Number of particles: " << nPhotons << std::endl;
  
  if (nPhotons != treeLC->GetEntries()) {
    std::cout << " - Mismatch LC tree " << treeLC->GetEntries() << " TS tree " << nPhotons << std::endl;
    return;
  }
      
  for (unsigned iV(0); iV<nV; ++iV){//loop on variables
    
    fout->cd();
    std::ostringstream histLabel;
    histLabel << "hist_" << label[iV];
    if (iV<nV1D) hist[iV] = new TH1F(histLabel.str().c_str(),axisLabel[iV].c_str(),nBins[iV],binMin[iV],binMax[iV]);
    else hist2D[iV-nV1D] = new TH2F(histLabel.str().c_str(),axisLabel[iV].c_str(),nBinsX[iV-nV1D],binXmin[iV-nV1D],binXmax[iV-nV1D],nBinsY[iV-nV1D],binYmin[iV-nV1D],binYmax[iV-nV1D]);
  }
  
  TBranch        *b_nAllLC;
  TBranch        *b_all_lc_energy;
  TBranch        *b_all_lc_z;
  TBranch        *b_all_lc_seedEta;
  TBranch        *b_all_lc_seedPhi;
  TBranch        *b_all_lc_layer;
  TBranch        *b_all_lc_nrechits;
  TBranch        *b_all_lc_mult;
  TBranch        *b_nCP;
  TBranch        *b_cp_eta;
  TBranch        *b_cp_phi;
  TBranch        *b_cp_energy;
  TBranch        *b_cp_missingEnergyFraction;
  TBranch        *b_nTS;
  TBranch        *b_ts_firstLayer;
  TBranch        *b_ts_lastLayer;
  TBranch        *b_ts_energy;
  TBranch        *b_ts_eta_fromLC;
  TBranch        *b_ts_phi_fromLC;
  TBranch        *b_lc_TSidx;
  TBranch        *b_nLC;
  TBranch        *b_lc_energy;
  TBranch        *b_lc_seedEta;
  TBranch        *b_lc_seedPhi;
  TBranch        *b_lc_layer;
  TBranch        *b_lc_nrechits;
  TBranch        *b_lc_mult;

  Int_t           nCP = 0;
  vector<double>        *cp_energy = 0;
  vector<double>        *cp_missingEnergyFraction = 0;
  vector<double>        *cp_eta = 0;
  vector<double>        *cp_phi = 0;
  Int_t           nAllLC = 0;
  Int_t           nTS = 0;
  vector<double>        *ts_firstLayer = 0;
  vector<double>        *ts_lastLayer = 0;
  vector<double>        *ts_energy = 0;
  vector<double>        *ts_eta = 0;
  vector<double>        *ts_phi = 0;
  vector<double>  *all_lc_energy = 0;
  vector<double>  *all_lc_z = 0;
  vector<double>  *all_lc_seedEta = 0;
  vector<double>  *all_lc_seedPhi = 0;
  vector<int>     *all_lc_layer = 0;
  vector<int>     *all_lc_nrechits = 0;
  vector<int>     *all_lc_mult = 0;
  Int_t           nLC = 0;
  vector<double>  *lc_energy = 0;
  vector<double>  *lc_seedEta = 0;
  vector<double>  *lc_seedPhi = 0;
  vector<int>     *lc_layer = 0;
  vector<int>     *lc_nrechits = 0;
  vector<int>     *lc_mult = 0;
  vector<int>     *lc_TSidx = 0;

  tree->SetBranchAddress("nTS", &nTS, &b_nTS);
  tree->SetBranchAddress("ts_firstLayer", &ts_firstLayer, &b_ts_firstLayer);
  tree->SetBranchAddress("ts_lastLayer", &ts_lastLayer, &b_ts_lastLayer);
  tree->SetBranchAddress("ts_energy", &ts_energy, &b_ts_energy);
  tree->SetBranchAddress("ts_eta_fromLC", &ts_eta, &b_ts_eta_fromLC);
  tree->SetBranchAddress("ts_phi_fromLC", &ts_phi, &b_ts_phi_fromLC);
  tree->SetBranchAddress("nCP", &nCP, &b_nCP);
  tree->SetBranchAddress("cp_energy", &cp_energy, &b_cp_energy);
  tree->SetBranchAddress("cp_missingEnergyFraction", &cp_missingEnergyFraction, &b_cp_missingEnergyFraction);
  tree->SetBranchAddress("cp_eta", &cp_eta, &b_cp_eta);
  tree->SetBranchAddress("cp_phi", &cp_phi, &b_cp_phi);
  tree->SetBranchAddress("nLC", &nLC, &b_nLC);
  tree->SetBranchAddress("lc_TSidx", &lc_TSidx, &b_lc_TSidx);
  tree->SetBranchAddress("lc_energy", &lc_energy, &b_lc_energy);
  tree->SetBranchAddress("lc_seedEta", &lc_seedEta, &b_lc_seedEta);
  tree->SetBranchAddress("lc_seedPhi", &lc_seedPhi, &b_lc_seedPhi);
  tree->SetBranchAddress("lc_layer", &lc_layer, &b_lc_layer);
  tree->SetBranchAddress("lc_nrechits", &lc_nrechits, &b_lc_nrechits);
  tree->SetBranchAddress("lc_mult", &lc_mult, &b_lc_mult);
  treeLC->SetBranchAddress("nAllLC", &nAllLC, &b_nAllLC);
  treeLC->SetBranchAddress("all_lc_energy", &all_lc_energy, &b_all_lc_energy);
  treeLC->SetBranchAddress("all_lc_z", &all_lc_z, &b_all_lc_z);
  treeLC->SetBranchAddress("all_lc_seedEta", &all_lc_seedEta, &b_all_lc_seedEta);
  treeLC->SetBranchAddress("all_lc_seedPhi", &all_lc_seedPhi, &b_all_lc_seedPhi);
  treeLC->SetBranchAddress("all_lc_layer", &all_lc_layer, &b_all_lc_layer);
  treeLC->SetBranchAddress("all_lc_nrechits", &all_lc_nrechits, &b_all_lc_nrechits);
  treeLC->SetBranchAddress("all_lc_mult", &all_lc_mult, &b_all_lc_mult);
  
  const int nEntries = tree->GetEntries();
  
  double zPos[nL];
  for (unsigned iL(0); iL<nL; ++iL){//loop on layers
    zPos[iL] = 0;
  }
  
  for (unsigned iE(0); iE < nEntries; ++iE){//loop on entries
    Long64_t ientry = treeLC->LoadTree(iE); 
    if (iE%100==0)
      std::cout << " -- Processing entry " << iE  << " " << ientry << std::endl;

    tree->GetEntry(iE);
    treeLC->GetEntry(iE);

    double newPhi = (*cp_phi)[0];
    if (newPhi<0) newPhi +=2*TMath::Pi();
    newPhi = newPhi*360./(2*TMath::Pi());
    
    if ((static_cast<int>(newPhi+0.5)%60) < phiMin ||
	(static_cast<int>(newPhi+0.5)%60) >= phiMax) continue;

    if ((*cp_eta)[0]<etaMin || (*cp_eta)[0]>=etaMax) continue;
    
    //vs nRH
    double Etot[3] = {0,0,0};
    double Etot1LC[3] = {0,0,0};

    //std::cout << " PASS phi [" << phiMin << "," << phiMax << "[ Evt " << iE << " found " << nAllLC << " LC. CP eta,phi = " << (*cp_eta)[0] << "," << newPhi << " " << static_cast<int>(newPhi+0.5) << " " << static_cast<int>(newPhi+0.5)%60 << std::endl;

    for (unsigned iL(0); iL<nL; ++iL){//loop on layers
      //std::cout << " -- Processing layer " << iL << std::endl;

      
      for (int iLC(0); iLC<nAllLC; ++iLC){
	if (doDebug && (*all_lc_layer)[iLC] == iL+1){
	  if (zPos[iL] <1) zPos[iL] = (*all_lc_z)[iLC];
	}
      }

      for (unsigned iRH(0); iRH<3; ++iRH){
	double minDRCP = 1000;
	double minE = 0;
	bool foundLC = false;
	
	for (int iLC(0); iLC<nAllLC; ++iLC){
	  
	  if ((*all_lc_layer)[iLC] == iL+1 && (*all_lc_nrechits)[iLC]>=iRH+1){
	    Etot[iRH] += (*all_lc_energy)[iLC];
	    double dR = sqrt(pow(TMath::Abs((*all_lc_seedEta)[iLC]-(*cp_eta)[0]),2)+pow(deltaPhi((*all_lc_seedPhi)[iLC],(*cp_phi)[0]),2));
	    if (dR<minDRCP) {
	      minDRCP = dR;
	      minE = (*all_lc_energy)[iLC];
	      foundLC = true;
	    }
	  
	  }
	}
	if (foundLC){
	  Etot1LC[iRH]+=minE;
	  //if (minDRCP>0.01) std::cout << " -- Entry " << iE << " layer " << iL+1 << " minDR = " << minDRCP << " CP eta-phi = " << cp_eta << " " << cp_phi << std::endl;
	  hist2D[4+iRH]->Fill(iL+1,minDRCP);
	}
      }//loop on RH
    }//loop on layers
    //std::cout << " --- End loop on layers " << std::endl;	

    for (unsigned iRH(0); iRH<3; ++iRH){
      hist[2*iRH]->Fill(Etot[iRH]/(*cp_energy)[0]);
      hist[2*iRH+1]->Fill(Etot1LC[iRH]/(*cp_energy)[0]);
    }

    //loop on TS
    double EsumTS = 0;
    double ETS_mindR = 0;
    double ETS_Emax = 0;
    double mindr = 10;
    double maxE = 0;
    double mindeta = 10;
    double mindphi = 10;
    int idxMindR = -1;
    int idxMaxE = -1;
    
    if (nCP!=1){
      std::cout << " -- Did not find 1 CP but " << nCP << " passing event." << std::endl;
      continue;
    }
    
    for (int its(0); its<nTS; ++its){
      double E = (*ts_energy)[its];
      EsumTS += E;
      double dr = sqrt(pow(TMath::Abs((*cp_eta)[0]-(*ts_eta)[its]),2)+pow(deltaPhi((*cp_phi)[0],(*ts_phi)[its]),2));
      double dphi = deltaPhi((*ts_phi)[its],(*cp_phi)[0]);
      double deta = TMath::Abs((*ts_eta)[its]-(*cp_eta)[0]);
      hist2D[0]->Fill((*ts_energy)[its]/(*cp_energy)[0],dr);
      hist2D[1]->Fill((*ts_energy)[its]/(*cp_energy)[0],deta);
      hist2D[2]->Fill((*ts_energy)[its]/(*cp_energy)[0],dphi);
      hist[16]->Fill((*ts_firstLayer)[its]);
      hist[17]->Fill((*ts_lastLayer)[its]);
      hist[18]->Fill((*ts_lastLayer)[its]-(*ts_firstLayer)[its]+1);
      if (dr<mindr){
	mindr=dr;
	mindphi=dphi;
	mindeta=deta;
	idxMindR = its;
      }
      if (E>maxE){
	maxE = E;
	idxMaxE = its;
      }
    }//loop on TS
    //std::cout << " --- End loop on TS " << std::endl;	

    double Etot1LCTS = 0;
    //for closest dR TS
    for (unsigned iL(0); iL<nL; ++iL){//loop on layers
      //skip CE-H for photons
      if ((iR<4 || iR>5) && iL>=nLEE) continue;
      
      double minDRTS = 1000;
      double minETS = 0;
      bool foundLCTS = false;
      
      for (int iLC(0); iLC<nLC; ++iLC){
	if ((*lc_layer)[iLC] == iL+1 && (*lc_TSidx)[iLC]==idxMindR){
	  double dR = sqrt(pow(TMath::Abs((*lc_seedEta)[iLC]-(*cp_eta)[0]),2)+pow(deltaPhi((*lc_seedPhi)[iLC],(*cp_phi)[0]),2));
	  if (dR<minDRTS) {
	    minDRTS = dR;
	    minETS = (*lc_energy)[iLC];
	    foundLCTS = true;
	  }
	}
      }
      if (foundLCTS){
	Etot1LCTS+=minETS;
	hist2D[7]->Fill(iL+1,minDRTS);
      }
    }
    //std::cout << " --- Fill histos " << std::endl;	

    hist2D[3]->Fill(idxMaxE,idxMindR);
    //Etot dependent on iter...
    unsigned iRH = 2;
    if (aIter.find("1")!=aIter.npos || aIter.find("Sim")!=aIter.npos) iRH=0;
    else if (aIter.find("2")!=aIter.npos) iRH=1;
    if (idxMindR>=0){
      hist[6]->Fill((*ts_energy)[idxMindR]/Etot[iRH]);
      hist[7]->Fill((*ts_energy)[idxMindR]/(*cp_energy)[0]);
    }
    hist[8]->Fill(Etot1LCTS/(*cp_energy)[0]);
    hist[9]->Fill(nTS);
    if (idxMaxE>=0){
      hist[10]->Fill((*ts_energy)[idxMaxE]/(*cp_energy)[0]);
      hist[15]->Fill((*ts_energy)[idxMaxE]/Etot[iRH]);
      hist[19]->Fill((*ts_firstLayer)[idxMaxE]);
      hist[20]->Fill((*ts_lastLayer)[idxMaxE]);
      hist[21]->Fill((*ts_lastLayer)[idxMaxE]-(*ts_firstLayer)[idxMaxE]+1);
    }
    hist[11]->Fill(EsumTS/(*cp_energy)[0]);
    hist[12]->Fill(mindr);
    hist[13]->Fill(mindeta);
    hist[14]->Fill(mindphi);
    hist[22]->Fill((*cp_missingEnergyFraction)[0]);
    
    if (nTS==0) fileOut << (*cp_energy)[0] << " " << (*cp_energy)[0]/cosh((*cp_eta)[0]) << " " << (*cp_eta)[0] << " " << iE << " " << ientry << std::endl;
    
    //std::cout << " --- End of event " << std::endl;	
  }//loop on entries

  fileOut.close();
  if (doDebug) {
    for (unsigned iL(0); iL<nL; ++iL){//loop on layers
      std::cout << "else if (iL==" << iL << ") return " << zPos[iL] << ";//cm" << std::endl;
    }	
    //return;
  }
  
  for (unsigned iE(0); iE<nV1D; ++iE){

    //set overflows
    hist[iE]->SetBinContent(hist[iE]->GetNbinsX(),hist[iE]->GetBinContent(hist[iE]->GetNbinsX())+hist[iE]->GetBinContent(hist[iE]->GetNbinsX()+1));
    hist[iE]->SetBinError(hist[iE]->GetNbinsX(),sqrt(pow(hist[iE]->GetBinError(hist[iE]->GetNbinsX()),2)+pow(hist[iE]->GetBinError(hist[iE]->GetNbinsX()+1),2)));
    hist[iE]->SetBinContent(hist[iE]->GetNbinsX()+1,0);
    hist[iE]->SetBinError(hist[iE]->GetNbinsX()+1,0);


    myc[iE]->cd();
    gStyle->SetOptStat("eMRuo");
    gPad->SetLogy(dology[iE]);
    hist[iE]->SetLineColor(1);
    hist[iE]->SetMarkerColor(1);
    hist[iE]->SetMarkerStyle(20);
    //hist[iE]->SetMaximum(dology[iE]? maxY[iE]*10 : maxY[iE]*1.2);
    //hist[iE]->SetMinimum(dology[iE]? 0.1 : minY[iE]*0.8);
    //plot only for one eta
    hist[iE]->Draw("PE");
    
    gStyle->SetOptStat("eMR");

    if (label[iE].find("Over")!=label[iE].npos){
      gStyle->SetOptFit(1111);
      
      hist[iE]->Fit("gaus","same");
      TF1 *fitFunc = (TF1*)hist[iE]->GetFunction("gaus");
      if (!fitFunc) continue;
      fileOutRes << aIter << " "
		 << pteta << " "
		 << label[iE] << " "
		 << hist[iE]->Integral(0,hist[iE]->GetNbinsX()+1) << " "
		 << hist[iE]->GetMean() << " " << hist[iE]->GetMeanError() << " "
		 << hist[iE]->GetRMS() << " " << hist[iE]->GetRMSError() << " "
		 << fitFunc->GetParameter(0) << " " << fitFunc->GetParError(0) << " "
		 << fitFunc->GetParameter(1) << " " << fitFunc->GetParError(1) << " "
		 << fitFunc->GetParameter(2) << " " << fitFunc->GetParError(2) << " "
		 << hist[iE]->GetRMS()/hist[iE]->GetMean() << " "
		 << fitFunc->GetParameter(2)/fitFunc->GetParameter(1)
		 << std::endl;
    }
    std::string nRHStr;
    if (iE<2) nRHStr = "n_{RH} #geq 1";
    else if (iE==2 || iE==3) nRHStr = "n_{RH} #geq 2";
    else if (iE==4 || iE==5) nRHStr = "n_{RH} #geq 3";
    else nRHStr = "";
    sprintf(buf,"%s %s %s %s",regShort[iR].c_str(),aIter.c_str(),pteta.c_str(),nRHStr.c_str());
    lat.SetTextSize(0.05);
    lat.SetTextColor(1);
    lat.DrawLatexNDC(0.12,0.95,buf);	
    
    myc[iE]->Update();
    myc[iE]->Print((lPlotDir.str()+"/"+label[iE]+"_"+regShort[iR]+".pdf").c_str());
    myc[iE]->Print((lPlotDir.str()+"/"+label[iE]+"_"+regShort[iR]+".png").c_str());
  }
  
  for (unsigned iE(0); iE<nV2D; ++iE){
    myc[nV1D+iE]->cd();
    gPad->SetLogz(1);
    gPad->SetRightMargin(0.15);
    hist2D[iE]->SetStats(0);
    hist2D[iE]->Draw("colz");

    std::string nRHStr;
    if (iE==4) nRHStr = "n_{RH} #geq 1";
    else if (iE==5) nRHStr = "n_{RH} #geq 2";
    else if (iE==6) nRHStr = "n_{RH} #geq 3";
    else nRHStr = "";
    sprintf(buf,"%s %s %s %s",regShort[iR].c_str(),aIter.c_str(),pteta.c_str(),nRHStr.c_str());
    lat.SetTextSize(0.05);
    lat.SetTextColor(1);
    lat.DrawLatexNDC(0.12,0.95,buf);	
  
    myc[nV1D+iE]->Update();
    myc[nV1D+iE]->Print((lPlotDir.str()+"/"+label[iE+nV1D]+"_"+regShort[iR]+".pdf").c_str());
    myc[nV1D+iE]->Print((lPlotDir.str()+"/"+label[iE+nV1D]+"_"+regShort[iR]+".png").c_str());

  }
  fout->Write();
  //cleanup
  fileOutRes.close();
  for (unsigned iE(0); iE<nV2D; ++iE){
    hist2D[iE]->Delete();
  }
  for (unsigned iE(0); iE<nV1D; ++iE){
    hist[iE]->Delete();
  }


  tree->Reset();
  treeLC->Reset();
  tree->Delete();
  treeLC->Delete();
  
}//merge

int mergeTracksters(){
  
  bool doDebug = false;
  bool binInPhi = false;
  bool binInEta = false;
  bool useEPhi = false;
  
  std::string filePath = useEPhi ? "D49_PhiScan" : "D49_DefSC";//D49_AllTracksters";
  //std::string filePath = "D49_DevMaxMiss";
  
  const unsigned nR = 4;

  const unsigned regIdx[4] = {0,1,4,6};
  //std::string iter[11] = {"Dummy1","Dummy2","Dummy3","EM1","EM2","EM3","TrkEM","EM","Trk","HAD","Sim"};
  
  const unsigned nPhi = binInPhi ? 3 : 1;
  const unsigned nEtaBin = binInEta ? 10 : 1;
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
  float etaMin[nEtaBin];
  float etaMax[nEtaBin];
  if (!binInEta) {
    etaMin[0] = 1.5;
    etaMax[0] = 3;
  } else {
    for (unsigned iEta(0); iEta<nEtaBin; ++iEta){
      etaMin[iEta] = 1.7 + 0.1*iEta;
      etaMax[iEta] = etaMin[iEta] + 0.1;
    }
  }

  for (unsigned iR(3); iR<nR; ++iR){

    const unsigned iReg = regIdx[iR];
    std::cout << " -- Processing region" << reg[iReg] << std::endl;
    
    double ptval[12] = {3,5,10,15,20,30,40,50,75,100,150,200};
    const unsigned nPT = doDebug ? 1 : iReg==4 ? 7 : useEPhi ? 1 : 12;
    if (iReg==4){
      ptval[1] = 10;
      ptval[2] = 20;
      ptval[3] = 50;
      ptval[4] = 100;
      ptval[5] = 150;
      ptval[6] = 200;
    }

    if (useEPhi) ptval[0] = 500;
    
    double etaval[6] = {17,19,21,23,25,27};
    const unsigned nEta = iReg==0? 6 : 1;
    if (iReg!=0) etaval[0] = 21;

    if (useEPhi) {
      for (unsigned ieta(0); ieta<6;++ieta){
	etaval[ieta] = 10*ieta;
      }
    }
    //std::string iter[17] = {"Dummy1","Dummy2","Dummy3","TRK1","TRK2","TRK3","EM1","EM2","EM3","HAD1","HAD2","HAD3","TrkEM","EM","Trk","HAD","Sim"};
    //std::string iter[9] = {"TrkEM","EM","Trk","HAD","Sim","HAD3","HAD3a","HAD3b","HAD3c"};
    std::string iter[11] = {"TrkEM","EM","Trk","HAD","Sim","Dummy1","Dummy2","Dummy3","EM1","EM2","EM3"};
    //std::string iter[8] = {"Dummy1","Dummy2","Dummy3","EM1","EM2","EM3","EM","Sim"};
    
    const unsigned nI = doDebug ? 1 : iReg==0 ? 4 : iReg==1 ? 8: iReg==6? 11 : 9;
    if (iReg==1) {
      iter[0] = "Dummy1";
      iter[1] = "Dummy3";
      iter[2] = "EM3";
      //iter[3] = "EM3a";
      //iter[4] = "EM3b";
      //iter[5] = "EM3c";
      iter[3] = "EM";
      iter[4] = "TrkEM";
      iter[5] = "Trk";
      iter[6] = "HAD";
      iter[7] = "Sim";
    }
    else if (iReg==0) {
      iter[0] = "Dummy1";
      iter[1] = "Dummy3";
      iter[2] = "EM3";
      iter[3] = "EM";
      iter[4] = "Sim";
    }
    
    if (doDebug) {
      ptval[0] = 20;
      etaval[0] = 21;
      iter[0] = "EM3";
    }
    
    for (unsigned iI(0); iI<nI; ++iI){
      std::cout << " -- Processing iter " << iter[iI] << std::endl;
      for (unsigned ipt(0); ipt<nPT; ++ipt){  
	for (unsigned ieta(0); ieta<nEta; ++ieta){
	  if (doDebug && ieta>0) continue;
	  if (!useEPhi) std::cout << " --- Processing pT " << ptval[ipt] << " eta " << etaval[ieta]<< std::endl;
	  else std::cout << " --- Processing E " << ptval[ipt] << " phi " << etaval[ieta]<< std::endl;
	  for (unsigned iPhi(0); iPhi<nPhi; ++iPhi){
	    for (unsigned iEta(0); iEta<nEtaBin; ++iEta){
	      
	      std::ostringstream pteta;
	      if (!useEPhi) pteta << "pt" << ptval[ipt] << "_eta" << etaval[ieta];
	      else pteta << "E" << ptval[ipt] << "_phi" << etaval[ieta];
	      std::ostringstream phiStr;
	      if (binInPhi) phiStr << "_phimod" << phiMin[iPhi] << "_" << phiMax[iPhi];
	      else if (binInEta) phiStr << "_eta" << etaMin[iEta]*10 << "_" << etaMax[iEta]*10;
	      merge(filePath,iReg,iter[iI],pteta.str(),phiStr.str(),doDebug,phiMin[iPhi],phiMax[iPhi],etaMin[iEta],etaMax[iEta]);
	    }//loop on eta bin
	  }//loop on phi bin
	}//loop on eta/phi
      }//loop on pT/E
    }//loop on iterations
    
  }//loop on reg
  return 0;
};
