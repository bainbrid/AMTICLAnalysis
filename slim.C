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

struct myTSsort{
  double E;
  unsigned iter;
  unsigned idx;  
};

bool compareTS(const myTSsort & lhs,const myTSsort & rhs){
  return lhs.E > rhs.E;
};


void fillRecoTree(std::string filePath,
		  const unsigned iR,
		  const std::string pteta,
		  const bool doDebug = false,
		  const bool doAMsum = false){

  SetTdrStyle();

  std::string partName = iR<2? "#gamma" : iR==4 || iR==5 ? "#pi^{#pm}" : "e^{#pm}";

  const bool doAM = filePath.find("AMiters")!=filePath.npos;
  
  std::ostringstream lPlotDir;
  lPlotDir << filePath << "/Reco/" << pteta << "/";
  if (doDebug) lPlotDir << "Debug";
  
  if (system(("mkdir -p "+lPlotDir.str()).c_str())){
    std::cout << " -- Cannot create output dir..." << lPlotDir.str() << std::endl;
    system(("echo \"-- return value \"$?"));
    return;
  }

  std::ofstream fileOut;
  fileOut.open((lPlotDir.str()+"/eventsToPrint_"+regShort[iR]+".dat").c_str());
  
  TFile *fout = 0;
  if (doAMsum) fout = TFile::Open((lPlotDir.str()+"/RecoTreeAM_"+regShort[iR]+".root").c_str(),"RECREATE");
  else  fout = TFile::Open((lPlotDir.str()+"/RecoTree_"+regShort[iR]+".root").c_str(),"RECREATE");
  if (!fout) return;
  

  const unsigned nL = 49;//layers
  const unsigned nLEE = 28;//layers

  fout->cd();
  TTree *newTree = new TTree("RecoTree","Tree with reco info");
  double eTotTS;
  double nTotTS;
  double eTS1;
  double eTS2;
  double eTS3;
  double deltaRTS12;
  double deltaRTS13;
  double deltaRTS23;
  double deltaEtaTS12;
  double deltaEtaTS13;
  double deltaEtaTS23;
  double deltaPhiTS12;
  double deltaPhiTS13;
  double deltaPhiTS23;
  double eTotLC3;
  double eTotSimTSLC3;
  double eTotTSdR;
  double nTotTSdR;
  double eTotLC3dR;
  double eTotAllLC3;
  double eTotSimLC3;
  double nAllLC3;
  double eTotAllLC2;
  double nAllLC2;
  double eTotAllLC1;
  double nAllLC1;
  unsigned allTS_firstL;
  unsigned allTS_lastL;
  double eTotLC3_ECAL;
  double eTotLC3_FHCAL;
  double eTotLC3_BHCAL;
  double eTotAllLC3_ECAL;
  double eTotAllLC3_FHCAL;
  double eTotAllLC3_BHCAL;
  double eTotSimLC3_ECAL;
  double eTotSimLC3_FHCAL;
  double eTotSimLC3_BHCAL;
  double eTotSimTSLC3_ECAL;
  double eTotSimTSLC3_FHCAL;
  double eTotSimTSLC3_BHCAL;

  //@@ const unsigned nIter = 7;
  const unsigned nIter = 1;
  //std::string iters[nIter] = {"Sim","TrkEM","EM","Trk","HAD","EMAM","HADAM"};
  //@@std::string iters[nIter] = {"TrkEM","EM"};//@@,"Trk","HAD","Sim"};
  std::string iters[nIter] = {"EM3"};//@@,"Trk","HAD","Sim"};
  //std::string AMiters[nIter] = {"Sim","TrkEMAM","EMAM","TrkAM","HADAM"};
  std::string AMiters[nIter];
  for (unsigned iI(0); iI<nIter;++iI){
    AMiters[iI] = iters[iI];
  }

  //@@bool skipInSum[nIter] = {true,false,false,false,false,true,true};
  //bool skipInSum[nIter] = {true,false,false,true,true,true,true};
  bool skipInSum[nIter] = {false};//,false,false};//@@,false,false,true};
  if (doAMsum) {
    skipInSum[1] = true;
    skipInSum[2] = true;
    skipInSum[3] = true;
    skipInSum[4] = true;
    skipInSum[5] = true;
    skipInSum[6] = false;
  }
  unsigned nSel[nIter];
  double sel_energy[nIter];
  double sumts_energy[nIter];
  unsigned sumts_firstL[nIter];
  unsigned sumts_lastL[nIter];
  unsigned nSel_dRcut[nIter];
  double sel_energy_dRcut[nIter];
  double maxwDeltaR[nL];
  double maxwDeltaRseed[nL];
  double maxDeltaR[nL];
  double maxDeltaRseed[nL];
  
  double dR_TrkEM_EM;
  double dR_TrkEM_Trk;
  double dR_EM_Trk;

  double deltaRCP12;
  double deltaEtaCP12;
  double deltaPhiCP12;

  double deltaXCP12;
  double deltaYCP12;
  double deltaXYCP12;
  
  newTree->Branch("nTotTS",&nTotTS);
  newTree->Branch("eTotTS",&eTotTS);
  newTree->Branch("eTS1",&eTS1);
  newTree->Branch("eTS2",&eTS2);
  newTree->Branch("eTS3",&eTS3);
  newTree->Branch("deltaRTS12",&deltaRTS12);
  newTree->Branch("deltaRTS13",&deltaRTS13);
  newTree->Branch("deltaRTS23",&deltaRTS23);
  newTree->Branch("deltaEtaTS12",&deltaEtaTS12);
  newTree->Branch("deltaEtaTS13",&deltaEtaTS13);
  newTree->Branch("deltaEtaTS23",&deltaEtaTS23);
  newTree->Branch("deltaPhiTS12",&deltaPhiTS12);
  newTree->Branch("deltaPhiTS13",&deltaPhiTS13);
  newTree->Branch("deltaPhiTS23",&deltaPhiTS23);
  newTree->Branch("eTotLC3",&eTotLC3);
  newTree->Branch("eTotLC3_ECAL",&eTotLC3_ECAL);
  newTree->Branch("eTotLC3_FHCAL",&eTotLC3_FHCAL);
  newTree->Branch("eTotLC3_BHCAL",&eTotLC3_BHCAL);
  newTree->Branch("eTotSimTSLC3",&eTotSimTSLC3);
  newTree->Branch("eTotSimTSLC3_ECAL",&eTotSimTSLC3_ECAL);
  newTree->Branch("eTotSimTSLC3_FHCAL",&eTotSimTSLC3_FHCAL);
  newTree->Branch("eTotSimTSLC3_BHCAL",&eTotSimTSLC3_BHCAL);
  newTree->Branch("eTotSimLC3",&eTotSimLC3);
  newTree->Branch("eTotSimLC3_ECAL",&eTotSimLC3_ECAL);
  newTree->Branch("eTotSimLC3_FHCAL",&eTotSimLC3_FHCAL);
  newTree->Branch("eTotSimLC3_BHCAL",&eTotSimLC3_BHCAL);
  newTree->Branch("eTotAllLC3_ECAL",&eTotAllLC3_ECAL);
  newTree->Branch("eTotAllLC3_FHCAL",&eTotAllLC3_FHCAL);
  newTree->Branch("eTotAllLC3_BHCAL",&eTotAllLC3_BHCAL);
  newTree->Branch("nTotTSdR",&nTotTSdR);
  newTree->Branch("eTotTSdR",&eTotTSdR);
  newTree->Branch("eTotLC3dR",&eTotLC3dR);
  newTree->Branch("eTotAllLC1",&eTotAllLC1);
  newTree->Branch("eTotAllLC2",&eTotAllLC2);
  newTree->Branch("eTotAllLC3",&eTotAllLC3);
  newTree->Branch("nAllLC1",&nAllLC1);
  newTree->Branch("nAllLC2",&nAllLC2);
  newTree->Branch("nAllLC3",&nAllLC3);
  newTree->Branch("allTS_lastL",&allTS_lastL);
  newTree->Branch("allTS_firstL",&allTS_firstL);

  newTree->Branch("dR_TrkEM_EM",&dR_TrkEM_EM);
  newTree->Branch("dR_TrkEM_Trk",&dR_TrkEM_Trk);
  newTree->Branch("dR_EM_Trk",&dR_EM_Trk);

  newTree->Branch("deltaRCP12",&deltaRCP12);
  newTree->Branch("deltaEtaCP12",&deltaEtaCP12);
  newTree->Branch("deltaPhiCP12",&deltaPhiCP12);

  newTree->Branch("deltaXCP12",&deltaXCP12);
  newTree->Branch("deltaYCP12",&deltaYCP12);
  newTree->Branch("deltaXYCP12",&deltaXYCP12);

//@@  for (unsigned iL(0); iL<nL;++iL){
//@@    std::ostringstream lName;
//@@    lName << "maxwDeltaR_" << iL+1;
//@@    newTree->Branch(lName.str().c_str(),&maxwDeltaR[iL]);
//@@    lName.str("");
//@@    lName << "maxwDeltaRseed_" << iL+1;
//@@    newTree->Branch(lName.str().c_str(),&maxwDeltaRseed[iL]);
//@@    lName.str("");
//@@    lName << "maxDeltaR_" << iL+1;
//@@    newTree->Branch(lName.str().c_str(),&maxDeltaR[iL]);
//@@    lName.str("");
//@@    lName << "maxDeltaRseed_" << iL+1;
//@@    newTree->Branch(lName.str().c_str(),&maxDeltaRseed[iL]);
//@@  }
  
  TChain *tree[nIter];
  TChain *treeLC = new TChain("ticlTree/treeLC");
  
  for (unsigned iC(0); iC<nIter;++iC){
    tree[iC] = 0;
    std::ostringstream lName;
    lName << "nTS_" << iters[iC];
    newTree->Branch(lName.str().c_str(),&nSel[iC]);
    lName.str("");
    lName << "sumlc_energy_" << iters[iC];
    newTree->Branch(lName.str().c_str(),&sel_energy[iC]); 
    lName.str("");
    lName << "sumts_energy_" << iters[iC];
    newTree->Branch(lName.str().c_str(),&sumts_energy[iC]); 
    lName.str("");
    lName << "sumts_lastL_" << iters[iC];
    newTree->Branch(lName.str().c_str(),&sumts_lastL[iC]); 
    lName.str("");
    lName << "sumts_firstL_" << iters[iC];
    newTree->Branch(lName.str().c_str(),&sumts_firstL[iC]); 
    lName.str("");
    lName << "nTS_dRcut_" << iters[iC];
    newTree->Branch(lName.str().c_str(),&nSel_dRcut[iC]);
    lName.str("");
    lName << "sumlc_energy_dRcut_" << iters[iC];
    newTree->Branch(lName.str().c_str(),&sel_energy_dRcut[iC]); 
    if (!doAM) tree[iC] = new TChain(("ticlTree/TSTree_"+iters[iC]).c_str());
    else tree[iC] = new TChain(("ticlTree/TSTree_"+AMiters[iC]).c_str());

    //@@const unsigned nRuns = 10;
    const unsigned nRuns = 4;
    for (unsigned iRun(0); iRun<nRuns; ++iRun){
      std::ostringstream label;
      label << filePath << "/" << reg[iR] << "/";
      if (filePath.find("FineCalo")!=filePath.npos) label << "FineCalo/current/";
      //@@label << "step3ticl_" << pteta << "_run" << iRun << "_FlatTracksters.root";
      label << "step3ticl_" << pteta << "_run" << iRun << "_FlatTracksters.root";

      std::cout << "FILE for iter '" << iters[iC] << "': " << label.str() << std::endl;
      tree[iC]->AddFile(label.str().c_str());
      if (iC==0) treeLC->AddFile(label.str().c_str());
    }
    if (!tree[iC]) return;
  }//loop on iters
  
  
  //E[iR] = ptval[ipt]*cosh(etaval[ieta]/10.);
  
  if (!treeLC) return;
  
  unsigned nPhotons = treeLC->GetEntries();
  std::cout << " -- Number of particles: " << nPhotons << std::endl;

  TBranch        *b_nAllLC;
  TBranch        *b_all_lc_energy;
  TBranch        *b_all_lc_z;
  TBranch        *b_all_lc_seedEta;
  TBranch        *b_all_lc_seedPhi;
  TBranch        *b_all_lc_layer;
  TBranch        *b_all_lc_nrechits;
  TBranch        *b_all_lc_isSi;
  TBranch        *b_all_lc_mult;

  Int_t           nAllLC = 0;
  vector<double>  *all_lc_energy = 0;
  vector<double>  *all_lc_z = 0;
  vector<double>  *all_lc_seedEta = 0;
  vector<double>  *all_lc_seedPhi = 0;
  vector<int>     *all_lc_layer = 0;
  vector<int>     *all_lc_nrechits = 0;
  vector<int>     *all_lc_isSi = 0;
  vector<int>     *all_lc_mult = 0;

  treeLC->SetBranchAddress("nAllLC", &nAllLC, &b_nAllLC);
  treeLC->SetBranchAddress("all_lc_energy", &all_lc_energy, &b_all_lc_energy);
  treeLC->SetBranchAddress("all_lc_z", &all_lc_z, &b_all_lc_z);
  treeLC->SetBranchAddress("all_lc_seedEta", &all_lc_seedEta, &b_all_lc_seedEta);
  treeLC->SetBranchAddress("all_lc_seedPhi", &all_lc_seedPhi, &b_all_lc_seedPhi);
  treeLC->SetBranchAddress("all_lc_layer", &all_lc_layer, &b_all_lc_layer);
  treeLC->SetBranchAddress("all_lc_nrechits", &all_lc_nrechits, &b_all_lc_nrechits);
  treeLC->SetBranchAddress("all_lc_isSi", &all_lc_isSi, &b_all_lc_isSi);
  treeLC->SetBranchAddress("all_lc_mult", &all_lc_mult, &b_all_lc_mult);

  TBranch        *b_nCP = 0;
  TBranch        *b_cp_eta = 0;
  TBranch        *b_cp_phi = 0;
  TBranch        *b_cp_energy = 0;

  Int_t           nCP = 0;
  vector<double>        *cp_energy = 0;
  vector<double>        *cp_eta = 0;
  vector<double>        *cp_phi = 0;

  tree[0]->SetBranchAddress("nCP", &nCP, &b_nCP);
  tree[0]->SetBranchAddress("cp_energy", &cp_energy, &b_cp_energy);
  tree[0]->SetBranchAddress("cp_eta", &cp_eta, &b_cp_eta);
  tree[0]->SetBranchAddress("cp_phi", &cp_phi, &b_cp_phi);

  TBranch        *b_nTS[nIter];
  TBranch        *b_ts_energy[nIter];
  TBranch        *b_ts_eta[nIter];
  TBranch        *b_ts_phi[nIter];
  TBranch        *b_ts_firstLayer[nIter];
  TBranch        *b_ts_lastLayer[nIter];
  TBranch        *b_nLC[nIter];
  TBranch        *b_lc_idx[nIter];
  TBranch        *b_lc_TSidx[nIter];
  TBranch        *b_lc_energy[nIter];
  TBranch        *b_lc_seedEta[nIter];
  TBranch        *b_lc_eta[nIter];
  TBranch        *b_lc_seedPhi[nIter];
  TBranch        *b_lc_layer[nIter];
  TBranch        *b_lc_tsMult[nIter];
  TBranch        *b_lc_nrechits[nIter];
  TBranch        *b_lc_isSi[nIter];

  Int_t           nTS[nIter];
  vector<double>        *ts_energy[nIter];
  vector<double>        *ts_eta[nIter];
  vector<double>        *ts_phi[nIter];
  vector<double>        *ts_firstLayer[nIter];
  vector<double>        *ts_lastLayer[nIter];
  Int_t           nLC[nIter];
  vector<double>  *lc_energy[nIter];
  vector<double>  *lc_eta[nIter];
  vector<double>  *lc_seedEta[nIter];
  vector<double>  *lc_seedPhi[nIter];
  vector<int>     *lc_layer[nIter];
  vector<int>     *lc_idx[nIter];
  vector<double>     *lc_tsMult[nIter];
  vector<int>     *lc_TSidx[nIter];
  vector<int>     *lc_nrechits[nIter];
  vector<int>     *lc_isSi[nIter];
    
    
  unsigned nEvtsProb[nIter];
  for (unsigned iC(0); iC<nIter;++iC){
    nEvtsProb[iC] = 0;
    
    b_nTS[iC] = 0;
    b_ts_energy[iC] = 0;
    b_ts_eta[iC] = 0;
    b_ts_phi[iC] = 0;
    b_ts_firstLayer[iC] = 0;
    b_ts_lastLayer[iC] = 0;
    b_nLC[iC] = 0;
    b_lc_idx[iC] = 0;
    b_lc_TSidx[iC] = 0;
    b_lc_energy[iC] = 0;
    b_lc_eta[iC] = 0;
    b_lc_seedEta[iC] = 0;
    b_lc_seedPhi[iC] = 0;
    b_lc_layer[iC] = 0;
    b_lc_tsMult[iC] = 0;
    b_lc_nrechits[iC] = 0;
    b_lc_isSi[iC] = 0;
    
    nTS[iC] = 0;
    ts_energy[iC] = 0;
    ts_eta[iC] = 0;
    ts_phi[iC] = 0;
    ts_firstLayer[iC] = 0;
    ts_lastLayer[iC] = 0;
    nLC[iC] = 0;
    lc_energy[iC] = 0;
    lc_eta[iC] = 0;
    lc_seedEta[iC] = 0;
    lc_seedPhi[iC] = 0;
    lc_layer[iC] = 0;
    lc_idx[iC] = 0;
    lc_tsMult[iC] = 0;
    lc_TSidx[iC] = 0;
    lc_nrechits[iC] = 0;
    lc_isSi[iC] = 0;
    
    tree[iC]->SetBranchAddress("nTS", &nTS[iC], &b_nTS[iC]);
    tree[iC]->SetBranchAddress("ts_energy", &ts_energy[iC], &b_ts_energy[iC]);
    tree[iC]->SetBranchAddress("ts_eta_fromLC", &ts_eta[iC], &b_ts_eta[iC]);
    tree[iC]->SetBranchAddress("ts_phi_fromLC", &ts_phi[iC], &b_ts_phi[iC]);
    tree[iC]->SetBranchAddress("ts_firstLayer", &ts_firstLayer[iC], &b_ts_firstLayer[iC]);
    tree[iC]->SetBranchAddress("ts_lastLayer", &ts_lastLayer[iC], &b_ts_lastLayer[iC]);
    tree[iC]->SetBranchAddress("nLC", &nLC[iC], &b_nLC[iC]);
    tree[iC]->SetBranchAddress("lc_idx", &lc_idx[iC], &b_lc_idx[iC]);
    tree[iC]->SetBranchAddress("lc_TSidx", &lc_TSidx[iC], &b_lc_TSidx[iC]);
    tree[iC]->SetBranchAddress("lc_energy", &lc_energy[iC], &b_lc_energy[iC]);
    tree[iC]->SetBranchAddress("lc_eta", &lc_eta[iC], &b_lc_eta[iC]);
    tree[iC]->SetBranchAddress("lc_seedEta", &lc_seedEta[iC], &b_lc_seedEta[iC]);
    tree[iC]->SetBranchAddress("lc_seedPhi", &lc_seedPhi[iC], &b_lc_seedPhi[iC]);
    tree[iC]->SetBranchAddress("lc_layer", &lc_layer[iC], &b_lc_layer[iC]);
    tree[iC]->SetBranchAddress("lc_tsMult", &lc_tsMult[iC], &b_lc_tsMult[iC]);
    tree[iC]->SetBranchAddress("lc_nrechits", &lc_nrechits[iC], &b_lc_nrechits[iC]);
    tree[iC]->SetBranchAddress("lc_isSi", &lc_isSi[iC], &b_lc_isSi[iC]);
  }
  
  const int nEntries = treeLC->GetEntries();

  
  for (unsigned iE(0); iE < nEntries; ++iE){//loop on entries
    Long64_t ientry = treeLC->LoadTree(iE); 
    if (iE%100==0)
      std::cout << " -- Processing entry " << iE  << " " << ientry << std::endl;
    treeLC->GetEntry(iE);
    
    eTotAllLC1 = 0;
    eTotAllLC2 = 0;
    eTotAllLC3 = 0;
    nAllLC1 = 0;
    nAllLC2 = 0;
    nAllLC3 = 0;
    eTotAllLC3_ECAL = 0;
    eTotAllLC3_FHCAL = 0;
    eTotAllLC3_BHCAL = 0;
    
    for (int iLC(0); iLC<nAllLC; ++iLC){
      eTotAllLC1 += (*all_lc_energy)[iLC];
      nAllLC1++;
      if (((*all_lc_isSi)[iLC]>0 && (*all_lc_nrechits)[iLC]>1) || ((*all_lc_isSi)[iLC]==0)){
	eTotAllLC2 += (*all_lc_energy)[iLC];
	nAllLC2++;
      }
      if ((*all_lc_isSi)[iLC]>0 && (*all_lc_nrechits)[iLC]<3) continue;
      eTotAllLC3 += (*all_lc_energy)[iLC];
      unsigned lLay = (*all_lc_layer)[iLC];
      if (lLay<=nLEE) eTotAllLC3_ECAL += (*all_lc_energy)[iLC];
      else if ((*all_lc_isSi)[iLC]>0)  eTotAllLC3_FHCAL += (*all_lc_energy)[iLC];
      else   eTotAllLC3_BHCAL += (*all_lc_energy)[iLC];
      nAllLC3++;
    }
    
    nTotTS = 0;
    eTotTS = 0;
    eTS1 = 0;
    eTS2 = 0;
    eTS3 = 0;
    deltaRTS12 = -1;
    deltaRTS13 = -1;
    deltaRTS23 = -1;
    deltaEtaTS12 = -1;
    deltaEtaTS13 = -1;
    deltaEtaTS23 = -1;
    deltaPhiTS12 = -1;
    deltaPhiTS13 = -1;
    deltaPhiTS23 = -1;
    eTotLC3 = 0;
    eTotLC3_ECAL = 0;
    eTotLC3_FHCAL = 0;
    eTotLC3_BHCAL = 0;

    eTotSimTSLC3 = 0;
    eTotSimTSLC3_ECAL = 0;
    eTotSimTSLC3_FHCAL = 0;
    eTotSimTSLC3_BHCAL = 0;
    eTotSimLC3 = 0;
    eTotSimLC3_ECAL = 0;
    eTotSimLC3_FHCAL = 0;
    eTotSimLC3_BHCAL = 0;

    nTotTSdR = 0;
    eTotTSdR = 0;
    eTotLC3dR = 0;

    dR_TrkEM_EM = 100;
    dR_TrkEM_Trk = 100;
    dR_EM_Trk = 100;

    allTS_lastL = 0;
    allTS_firstL = 60;
    for (unsigned iL(0); iL<nL;++iL){
      maxwDeltaR[iL] = 0;
      maxwDeltaRseed[iL] = 0;
      maxDeltaR[iL] = 0;
      maxDeltaRseed[iL] = 0;
    }

    deltaPhiCP12 = -1;
    deltaEtaCP12 = -1;
    deltaRCP12 = -1;
    deltaXCP12 = -1;
    deltaYCP12 = -1;
    deltaXYCP12 = -1;
    if ( nCP > 1 ) {
      deltaPhiCP12 = abs(deltaPhi((*cp_phi)[0],(*cp_phi)[1]));
      deltaEtaCP12 = abs(deltaEta((*cp_eta)[0],(*cp_eta)[1]));
      deltaRCP12 = deltaR((*cp_eta)[0],(*cp_eta)[1],(*cp_phi)[0],(*cp_phi)[1]);
//      std::cout << "DELTA: " << std::endl
//		<< " dP: " << deltaPhi((*cp_phi)[0],(*cp_phi)[1]) << std::endl
//		<< " dE: " << deltaEta((*cp_eta)[0],(*cp_eta)[1]) << std::endl
//		<< " dR: " << deltaR((*cp_eta)[0],(*cp_eta)[1],(*cp_phi)[0],(*cp_phi)[1]) << std::endl;

      // photons are "PV pointing", uncharged, massless 
      // no need to extrapolate momentum vector to ECAL surface
      // simply convert eta / phi 
      float tmp_theta = (180./3.141)*2.*std::atan(std::exp(-2.44));
      float cp1_theta = 2.*std::atan(std::exp(-1.*(*cp_eta)[0]));
      float cp2_theta = 2.*std::atan(std::exp(-1.*(*cp_eta)[1]));
      float cp1_phi = (*cp_phi)[0];
      float cp2_phi = (*cp_phi)[1];
      float cp1_z = 320.;
      float cp2_z = 320.;
      float cp1_r = std::tan(cp1_theta)*cp1_z;
      float cp2_r = std::tan(cp2_theta)*cp2_z;
      float cp1_x = std::cos(cp1_phi)*cp1_r;
      float cp2_x = std::cos(cp2_phi)*cp2_r;
      float cp1_y = std::sin(cp1_phi)*cp1_r;
      float cp2_y = std::sin(cp2_phi)*cp2_r;
      float cp12_dx = cp1_x-cp2_x;
      float cp12_dy = cp1_y-cp2_y;
      float cp12_dr = std::sqrt(cp12_dx*cp12_dx+cp12_dy*cp12_dy);
      deltaXCP12 = cp12_dx;
      deltaYCP12 = cp12_dy;
      deltaXYCP12 = cp12_dr;

//      std::cout << "DELTA: " << std::endl
//		<< " cp1_theta: " << cp1_theta << std::endl
//		<< " cp2_theta: " << cp2_theta << std::endl
//		<< " cp1_phi:   " << cp1_phi << std::endl
//		<< " cp2_phi:   " << cp2_phi << std::endl
//		<< " cp1_z:     " << cp1_z << std::endl
//		<< " cp2_z:     " << cp2_z << std::endl
//		<< " cp1_x:     " << cp1_x << std::endl
//		<< " cp2_x:     " << cp2_x << std::endl
//		<< " cp1_y:     " << cp1_y << std::endl
//		<< " cp2_y:     " << cp2_y << std::endl
//		<< " cp12_dx:   " << cp12_dx << std::endl
//		<< " cp12_dy:   " << cp12_dy << std::endl
//		<< " cp12_dr:   " << cp12_dr << std::endl;
    }
    
    std::map<int,int> lSimLC;
    std::vector<myTSsort> sumTSeVec;
    for (unsigned iC(0); iC<nIter;++iC){

      //std::cout << iters[iC] << std::endl;
      tree[iC]->GetEntry(iE);
      
      //sim iter only, get max dR to truth
      if (iC==0){
	for (unsigned iL(0); iL<nL;++iL){
	  double trueR = getZpos(iL)*tan(2*atan(exp(-(*cp_eta)[0])));
	  double eTotLayer = 0;
	  //std::cout << iL << " " << trueR << std::endl;
	  for (int iLC(0); iLC<nLC[0]; ++iLC){
	    unsigned lLay = (*lc_layer[0])[iLC];
	    if (lLay != iL+1) continue;
	    if ((*lc_isSi[0])[iLC]>0 && (*lc_nrechits[0])[iLC]<3) continue;
	    eTotLayer += (*lc_energy[0])[iLC];
	  }
	  if (eTotLayer==0) continue;
	  for (int iLC(0); iLC<nLC[0]; ++iLC){
	    unsigned lLay = (*lc_layer[0])[iLC];
	    if (lLay != iL+1) continue;
	    if ((*lc_isSi[0])[iLC]>0 && (*lc_nrechits[0])[iLC]<3) continue;
	 
	    double lcRseed = getZpos(iL)*tan(2*atan(exp(-(*lc_seedEta[0])[iLC])));
	    double diffTruth = abs(lcRseed-trueR)*(*lc_energy[0])[iLC]/eTotLayer;
	    if (diffTruth> maxDeltaRseed[iL])  maxDeltaRseed[iL] = diffTruth;
	    double lcR = getZpos(iL)*tan(2*atan(exp(-(*lc_eta[0])[iLC])));
	    diffTruth = abs(lcRseed-trueR)*(*lc_energy[0])[iLC]/eTotLayer;
	    if (diffTruth> maxDeltaR[iL])  maxDeltaR[iL] = diffTruth;
	  }

	  maxwDeltaR[iL] = getDEta((*cp_eta)[0],getZpos(iL),maxDeltaR[iL]);
	  maxwDeltaRseed[iL] = getDEta((*cp_eta)[0],getZpos(iL),maxDeltaRseed[iL]);
	}
      }
    
      
      sel_energy[iC] = 0;
      sel_energy_dRcut[iC] = 0;
      
      double eProb = 0;
      unsigned nLCprob = 0;
      for (int iLC(0); iLC<nLC[iC]; ++iLC){
	//for Sim iteration: just do 3 RH min in Si....
	unsigned lLay = (*lc_layer[iC])[iLC];
	if (iC==0 && (*lc_isSi[iC])[iLC]>0 && (*lc_nrechits[iC])[iLC]<3) continue;

	double tsMult = (*lc_tsMult[iC])[iLC];
	if (iC==0) {
	  lSimLC[(*lc_idx[iC])[iLC]] = 1;
	  eTotSimLC3 += (*lc_energy[iC])[iLC]/tsMult;
	  if (lLay<=nLEE) eTotSimLC3_ECAL += (*lc_energy[iC])[iLC]/tsMult;
	  else if ((*lc_isSi[iC])[iLC]>0)  eTotSimLC3_FHCAL += (*lc_energy[iC])[iLC]/tsMult;
	  else   eTotSimLC3_BHCAL += (*lc_energy[iC])[iLC]/tsMult;	  
	}
	else if (!skipInSum[iC]){
	  if (lSimLC.find((*lc_idx[iC])[iLC])!=lSimLC.end()){
	    eTotSimTSLC3 += (*lc_energy[iC])[iLC]/tsMult;
	    if (lLay<=nLEE) eTotSimTSLC3_ECAL += (*lc_energy[iC])[iLC]/tsMult;
	    else if ((*lc_isSi[iC])[iLC]>0)  eTotSimTSLC3_FHCAL += (*lc_energy[iC])[iLC]/tsMult;
	    else   eTotSimTSLC3_BHCAL += (*lc_energy[iC])[iLC]/tsMult;
	  }
	}
	//else if (iC!=0 && (*lc_nrechits[iC])[iLC]<3) {
	//std::cout << " Arf, less than 3 RH !" << std::endl;
	//}
	if (iC>0 && tsMult>nTS[iC]) {
	  eProb += (*lc_energy[iC])[iLC];
	  nLCprob++;
	}

	if (tsMult>0){
	  if (!skipInSum[iC]) {
	    eTotLC3 += (*lc_energy[iC])[iLC]/tsMult;
	    if (lLay<=nLEE) eTotLC3_ECAL += (*lc_energy[iC])[iLC]/tsMult;
	    else if ((*lc_isSi[iC])[iLC]>0)  eTotLC3_FHCAL += (*lc_energy[iC])[iLC]/tsMult;
	    else   eTotLC3_BHCAL += (*lc_energy[iC])[iLC]/tsMult;
	  }
	  sel_energy[iC] += (*lc_energy[iC])[iLC]/tsMult;

	  double dRtoCP = sqrt(pow(TMath::Abs((*lc_seedEta[iC])[iLC]-(*cp_eta)[0]),2)+pow(deltaPhi((*lc_seedPhi[iC])[iLC],(*cp_phi)[0]),2));
	  if (dRtoCP < 0.4) {
	    if (!skipInSum[iC]) eTotLC3dR += (*lc_energy[iC])[iLC]/tsMult;
	    sel_energy_dRcut[iC] += (*lc_energy[iC])[iLC]/tsMult;
	  }
	  
	}
	else {
	  std::cout << " -- found 0 multiplicity LC: E = " << (*lc_energy[iC])[iLC]
		    << std::endl;
	}
	//	if ((*lc_nrechits)[iLC]<3) continue;
	//else eTotLC3 += (*lc_energy)[iLC];
      }//loop on LCs
      if (nLCprob>0){
	if (iters[iC] != "Trk") std::cout << " -- Inconsistency: entry " << iE << " iter " << iters[iC] << " nTS " << nTS[iC] << " nLCprob = " << nLCprob << " eProb = " << eProb << std::endl;
	nEvtsProb[iC]++;
      }

      
      nSel[iC] = 0;
      nSel_dRcut[iC] = 0;
      sumts_energy[iC] = 0;
      sumts_firstL[iC] = 0;
      sumts_lastL[iC] = 0;
      unsigned minL = 60;
      unsigned maxL = 0;
      
      for (int its(0); its<nTS[iC]; ++its){
	nSel[iC]++;
	sumts_energy[iC] += (*ts_energy[iC])[its];
	if ((*ts_lastLayer[iC])[its]>maxL) maxL = (*ts_lastLayer[iC])[its];
	if ((*ts_firstLayer[iC])[its]<minL) minL = (*ts_firstLayer[iC])[its];
	if ((*ts_lastLayer[iC])[its]>allTS_lastL) allTS_lastL = (*ts_lastLayer[iC])[its];
	if ((*ts_firstLayer[iC])[its]<allTS_firstL) allTS_firstL = (*ts_firstLayer[iC])[its];

	if (!skipInSum[iC]) {
	  nTotTS++;
	  eTotTS += (*ts_energy[iC])[its];
	  myTSsort lele;
	  lele.E = (*ts_energy[iC])[its];
	  lele.iter = iC;
	  lele.idx = its;
	  sumTSeVec.push_back(lele);
	}
	double dRtoCP = sqrt(pow(TMath::Abs((*ts_eta[iC])[its]-(*cp_eta)[0]),2)+pow(deltaPhi((*ts_phi[iC])[its],(*cp_phi)[0]),2));
	if (dRtoCP < 0.4) {
	  nSel_dRcut[iC]++;
	  if (!skipInSum[iC]) {
	    nTotTSdR++;
	    eTotTSdR += (*ts_energy[iC])[its];
	  }
	}
      }
      sumts_firstL[iC] = minL;
      sumts_lastL[iC] = maxL;
      //std::cout << iC << " " << nSel[iC] << " " << sel_energy[iC] << std::endl;
    }//loop on iterations 
    
    if (sumTSeVec.size()>0){
      std::sort(sumTSeVec.begin(),sumTSeVec.end(),compareTS);
      eTS1 = sumTSeVec[0].E;
      if (sumTSeVec.size()>1){
	eTS2 = sumTSeVec[1].E;
	deltaRTS12 = deltaR((*ts_eta[sumTSeVec[0].iter])[sumTSeVec[0].idx],
			    (*ts_eta[sumTSeVec[1].iter])[sumTSeVec[1].idx],
			    (*ts_phi[sumTSeVec[0].iter])[sumTSeVec[0].idx],
			    (*ts_phi[sumTSeVec[1].iter])[sumTSeVec[1].idx]);
	deltaEtaTS12 = deltaEta((*ts_eta[sumTSeVec[0].iter])[sumTSeVec[0].idx],
				(*ts_eta[sumTSeVec[1].iter])[sumTSeVec[1].idx]);
	deltaPhiTS12 = deltaPhi((*ts_phi[sumTSeVec[0].iter])[sumTSeVec[0].idx],
				(*ts_phi[sumTSeVec[1].iter])[sumTSeVec[1].idx]);
	if (sumTSeVec.size()>2){
	  eTS3 = sumTSeVec[2].E;
	  deltaRTS13 = deltaR((*ts_eta[sumTSeVec[0].iter])[sumTSeVec[0].idx],
			      (*ts_eta[sumTSeVec[2].iter])[sumTSeVec[2].idx],
			      (*ts_phi[sumTSeVec[0].iter])[sumTSeVec[0].idx],
			      (*ts_phi[sumTSeVec[2].iter])[sumTSeVec[2].idx]);
	  deltaEtaTS13 = deltaEta((*ts_eta[sumTSeVec[0].iter])[sumTSeVec[0].idx],
				  (*ts_eta[sumTSeVec[2].iter])[sumTSeVec[2].idx]);
	  deltaPhiTS13 = deltaPhi((*ts_phi[sumTSeVec[0].iter])[sumTSeVec[0].idx],
				  (*ts_phi[sumTSeVec[2].iter])[sumTSeVec[2].idx]);
	  deltaRTS23 = deltaR((*ts_eta[sumTSeVec[2].iter])[sumTSeVec[2].idx],
			      (*ts_eta[sumTSeVec[1].iter])[sumTSeVec[1].idx],
			      (*ts_phi[sumTSeVec[2].iter])[sumTSeVec[2].idx],
			      (*ts_phi[sumTSeVec[1].iter])[sumTSeVec[1].idx]);
	  deltaEtaTS23 = deltaEta((*ts_eta[sumTSeVec[2].iter])[sumTSeVec[2].idx],
				  (*ts_eta[sumTSeVec[1].iter])[sumTSeVec[1].idx]);
	  deltaPhiTS23 = deltaPhi((*ts_phi[sumTSeVec[2].iter])[sumTSeVec[2].idx],
				  (*ts_phi[sumTSeVec[1].iter])[sumTSeVec[1].idx]);
	  
	}
	
      }
    }

    
    newTree->Fill();
    
    fileOut << iE << " " << ientry << " "
	    << eTotTS << " "
	    << eTotLC3 << " dR4 "
	    << eTotTSdR << " "
	    << eTotLC3dR << " ";
    for (unsigned iC(0); iC<nIter;++iC){
      fileOut << iters[iC] << " "
	      << nSel[iC] << " "
	      << sel_energy[iC] << " dR4 "
	      << nSel_dRcut[iC] << " "
	      << sel_energy_dRcut[iC];
      if (iC<nIter-1) fileOut << " ";
    }
    
    fileOut << std::endl;

    if (fabs(eTotTS-eTotLC3)>0.01*eTotTS){
      std::cout << "-- Problematic event: "
		<< iE << " " << ientry << " "
		<< eTotTS << " "
		<< eTotLC3;
      for (unsigned iC(0); iC<nIter;++iC){
	std::cout << iters[iC] << " nTS="
		<< nSel[iC] << " ETS="
		<< sumts_energy[iC] << " ELS="
		<< sel_energy[iC];
      }
      std::cout << std::endl;
    }	  
    
    //std::cout << " --- Fill histos " << std::endl;	
    
    //std::cout << " --- End of event " << std::endl;	
    
  }//loop on entries

  
  fout->Write();
  std::cout << " -- Tree " << fout->GetName() << " saved." << std::endl;

  std::cout << " -- Found " ;
  for (unsigned iC(0); iC<nIter;++iC){
    std::cout << iters[iC] << " "
	      << nEvtsProb[iC] << " ";
  }
  std::cout << " problematic events." << std::endl;
  
  treeLC->Delete();
  for (unsigned iC(0); iC<nIter;++iC){
    tree[iC]->Delete();
  }
  fileOut.close();
  
}//fill

//@@int exploreTSReco(){
int slim(){
  
  //@@bool doDebug = false;
  bool doDebug = true;

  //std::string filePath = "D49_AllTracksters";
  //std::string filePath = "D49_FineCalo_AMiters";
  //std::string filePath = "D49_DevMaxMiss";
  //@@std::string filePath = "D49_DefSC";
  std::string filePath = "../TiCLTreeProducer/D49_FineCalo";
  
  //@@const unsigned regIdx[5] = {0,1,4,5,6};
  const unsigned regIdx[1] = {0};
  //@@for (unsigned iR=4; iR<5; ++iR){
  for (unsigned iR=0; iR<1; ++iR){
      
    const unsigned iReg = regIdx[iR];
    std::cout << " -- Processing region " << reg[iReg] << std::endl;
    
    //@@ double ptval[12] = {3,5,10,15,20,30,40,50,75,100,150,200};
    double ptval[1] = {200};//@@{50};
    const unsigned nPT = doDebug ? 1 : iReg==4 ? 7 : 12;
    if (iReg==4 || iReg==5){
      ptval[1] = 10;
      ptval[2] = 20;
      ptval[3] = 50;
      ptval[4] = 100;
      ptval[5] = 150;
      ptval[6] = 200;
    }
    //@@double etaval[6] = {17,19,21,23,25,27};
    double etaval[1] = {21};
    //@@const unsigned nEta = iReg==0? 6 : 1;
    const unsigned nEta = iReg==0? 1 : 1;
    if (iReg!=0) etaval[0] = 21;
    
    if (doDebug) {
      ptval[0] = 200;//@@50;
      etaval[0] = 21;
    }
    
    for (unsigned ipt(0); ipt<nPT; ++ipt){  
      for (unsigned ieta(0); ieta<nEta; ++ieta){
	std::cout << " --- Processing pT " << ptval[ipt] << " eta " << etaval[ieta]<< std::endl;
	
	std::ostringstream pteta;
	pteta << "pt" << ptval[ipt] << "_eta" << etaval[ieta];
	fillRecoTree(filePath,iReg,pteta.str(),doDebug,false);
	//@@fillRecoTree(filePath,iReg,pteta.str(),doDebug,true); // doAMsum
      }//loop on eta
    }//loop on pT
    
  }//loop on reg
  return 0;
};


