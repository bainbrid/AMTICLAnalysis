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

void fillTruthTree(std::string filePath,
		   const unsigned iR,
		   const std::string aIter,
		   const std::string pteta,
		   const bool doDebug = false){

  std::map<int,std::pair<int,std::string> > pdgMap = getPDGMap("pdgList.dat");

  SetTdrStyle();

  std::string partName = iR<2? "#gamma" : iR==4 || iR==5 ? "#pi^{#pm}" : "e^{#pm}";
  
  std::ostringstream lPlotDir;
  lPlotDir << filePath << "/Truth/" << aIter << "/" << pteta << "/";
  if (doDebug) lPlotDir << "Debug";
  
  if (system(("mkdir -p "+lPlotDir.str()).c_str())){
    std::cout << " -- Cannot create output dir..." << lPlotDir.str() << std::endl;
    system(("echo \"-- return value \"$?"));
    return;
  }


  std::ofstream fileOut;
  fileOut.open((lPlotDir.str()+"/eventsToPrint_"+regShort[iR]+".dat").c_str());

  
  TFile *fout = TFile::Open((lPlotDir.str()+"/TruthTree_"+regShort[iR]+".root").c_str(),"RECREATE");
  if (!fout) return;
  

  const unsigned nL = 52;//layers
  const unsigned nLEE = 28;//layers

  fout->cd();
  TTree *newTree = new TTree("TruthTree","Tree with truth info for PR");
  const unsigned nC = 5;
  const double eCut[nC] = {0,0.01,0.03,0.05,0.1};
  unsigned nSC_eCut[nC];
  const unsigned nCat = 10;
  double eFrac[nCat];
  double eFracAtB[nCat];
  double eTotSC;
  double eTotLC;
  double eTotLC3;
  double eTotAllLC;
  double eTotAllLC3;
  double nAllLC3;
  double eTotAllLC2;
  double nAllLC2;
  double eTotSCAtB;
  double eTotSC_noPi;
  double eTotSC_eCut[nC];
  double eTotSCAtB_eCut[nC];
  vector<double> sc_dRtoCP_eCut[nC];
  vector<double> sc_energy_eCut[nC];
  vector<double> sc_energyAtB_eCut[nC];
  vector<double> sc_zAtB_eCut[nC];
  vector<double> sc_eta_eCut[nC];
  vector<double> sc_phi_eCut[nC];
  vector<double> sc_etaAtB_eCut[nC];
  vector<double> sc_phiAtB_eCut[nC];
  vector<int> sc_pdgid_eCut[nC];
  vector<int> sc_pdgcat_eCut[nC];
  vector<int> sc_tsIdx_eCut[nC];
  vector<int> sc_idx_eCut[nC];
  vector<double> sc_tsminDR_eCut[nC];

  unsigned nTS_eCut10;
  vector<double> ts_energy_eCut10;

  unsigned nBlobs;
  double EBlobs;
  double EBlobs_ECAL;
  double EBlobs_FHCAL;
  double EBlobs_BHCAL;
  double EBlobs_HoverE;
  unsigned minMiss;
  unsigned maxMiss;
  double maxEfracBlobs;
  double minEfracBlobs;
  unsigned maxLengthBlobs;
  unsigned minLengthBlobs;

  vector<double> energy_blobs;
  vector<double> eFrac_blobs;
  vector<unsigned> length_blobs;
  vector<unsigned> firstL_blobs;
  vector<unsigned> nLC_blobs;
  vector<unsigned> nMiss_blobs;
  vector<unsigned> lc_in_blob;

  
  newTree->Branch("eTotSC",&eTotSC);
  newTree->Branch("eTotLC",&eTotLC);
  newTree->Branch("eTotLC3",&eTotLC3);
  newTree->Branch("eTotAllLC",&eTotAllLC);
  newTree->Branch("eTotAllLC3",&eTotAllLC3);
  newTree->Branch("nAllLC3",&nAllLC3);
  newTree->Branch("eTotAllLC2",&eTotAllLC2);
  newTree->Branch("nAllLC2",&nAllLC2);
  newTree->Branch("eTotSCAtB",&eTotSCAtB);
  newTree->Branch("eTotSC_noPi",&eTotSC_noPi);
  newTree->Branch("nTS_eCut10",&nTS_eCut10);
  newTree->Branch("ts_energy_eCut10",&ts_energy_eCut10);
  newTree->Branch("nBlobs",&nBlobs);
  newTree->Branch("EBlobs",&EBlobs);
  newTree->Branch("EBlobs_ECAL",&EBlobs_ECAL);
  newTree->Branch("EBlobs_FHCAL",&EBlobs_FHCAL);
  newTree->Branch("EBlobs_BHCAL",&EBlobs_BHCAL);
  newTree->Branch("EBlobs_HoverE",&EBlobs_HoverE);
  newTree->Branch("minMiss",&minMiss);
  newTree->Branch("maxMiss",&maxMiss);
  newTree->Branch("maxEfracBlobs",&maxEfracBlobs);
  newTree->Branch("minEfracBlobs",&minEfracBlobs);
  newTree->Branch("maxLengthBlobs",&maxLengthBlobs);
  newTree->Branch("minLengthBlobs",&minLengthBlobs);
  newTree->Branch("energy_blobs",&energy_blobs);
  newTree->Branch("length_blobs",&length_blobs);
  newTree->Branch("eFrac_blobs",&eFrac_blobs);
  newTree->Branch("firstL_blobs",&firstL_blobs);
  newTree->Branch("nLC_blobs",&nLC_blobs);
  newTree->Branch("lc_in_blob",&lc_in_blob);
  newTree->Branch("nMiss_blobs",&nMiss_blobs);
  
  for (unsigned iC(0); iC<nCat;++iC){
    std::ostringstream lName;
    lName << "eFrac_cat" << pdgBinLabel(iC+1);
    newTree->Branch(lName.str().c_str(),&eFrac[iC]);
    lName.str("");
    lName << "eFracAtB_cat" << pdgBinLabel(iC+1);
    newTree->Branch(lName.str().c_str(),&eFracAtB[iC]);
  }


  for (unsigned iC(0); iC<nC;++iC){
    std::ostringstream lName;
    lName << "nSC_eCut" << eCut[iC]*100;
    newTree->Branch(lName.str().c_str(),&nSC_eCut[iC]);
    lName.str("");
    lName << "sc_dRtoCP_eCut" << eCut[iC]*100;
    newTree->Branch(lName.str().c_str(),&sc_dRtoCP_eCut[iC]);
    lName.str("");
    lName << "sc_energy_eCut" << eCut[iC]*100;
    newTree->Branch(lName.str().c_str(),&sc_energy_eCut[iC]);
    lName.str("");
    lName << "sc_energyAtB_eCut" << eCut[iC]*100;
    newTree->Branch(lName.str().c_str(),&sc_energyAtB_eCut[iC]);
    lName.str("");
    lName << "sc_zAtB_eCut" << eCut[iC]*100;
    newTree->Branch(lName.str().c_str(),&sc_zAtB_eCut[iC]);
    lName.str("");
    lName << "sc_eta_eCut" << eCut[iC]*100;
    newTree->Branch(lName.str().c_str(),&sc_eta_eCut[iC]);
    lName.str("");
    lName << "sc_phi_eCut" << eCut[iC]*100;
    newTree->Branch(lName.str().c_str(),&sc_phi_eCut[iC]);
    lName.str("");
    lName << "sc_etaAtB_eCut" << eCut[iC]*100;
    newTree->Branch(lName.str().c_str(),&sc_etaAtB_eCut[iC]);
    lName.str("");
    lName << "sc_phiAtB_eCut" << eCut[iC]*100;
    newTree->Branch(lName.str().c_str(),&sc_phiAtB_eCut[iC]);
    lName.str("");
    lName << "sc_pdgid_eCut" << eCut[iC]*100;
    newTree->Branch(lName.str().c_str(),&sc_pdgid_eCut[iC]);
    lName.str("");
    lName << "sc_pdgcat_eCut" << eCut[iC]*100;
    newTree->Branch(lName.str().c_str(),&sc_pdgcat_eCut[iC]);
    lName.str("");
    lName << "sc_idx_eCut" << eCut[iC]*100;
    newTree->Branch(lName.str().c_str(),&sc_idx_eCut[iC]);
    lName.str("");
    lName << "sc_tsIdx_eCut" << eCut[iC]*100;
    newTree->Branch(lName.str().c_str(),&sc_tsIdx_eCut[iC]);
    lName.str("");
    lName << "sc_tsminDR_eCut" << eCut[iC]*100;
    newTree->Branch(lName.str().c_str(),&sc_tsminDR_eCut[iC]);
    lName.str("");
    lName << "eTotSC_eCut" << eCut[iC]*100;
    newTree->Branch(lName.str().c_str(),&eTotSC_eCut[iC]);
    lName.str("");
    lName << "eTotSCAtB_eCut" << eCut[iC]*100;
    newTree->Branch(lName.str().c_str(),&eTotSCAtB_eCut[iC]);
  }



  
  TChain *tree = new TChain(("ticlTree/TSTree_"+aIter).c_str());
  TChain *treeLC = new TChain("ticlTree/treeLC");
  
  //E[iR] = ptval[ipt]*cosh(etaval[ieta]/10.);

  const unsigned nRuns = 10;
  for (unsigned iRun(0); iRun<nRuns; ++iRun){
    std::ostringstream label;
    label << filePath << "/" << reg[iR] << "/";
    
    if (filePath.find("FineCalo")!=filePath.npos) label << "FineCalo/";
    label << "step3ticl_" << pteta << "_run" << iRun << "_FlatTracksters.root";
    tree->AddFile(label.str().c_str());
    treeLC->AddFile(label.str().c_str());
  }
  
  if (!tree) return;
  if (!treeLC) return;
  
  unsigned nPhotons = tree->GetEntries();
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
 
  TBranch        *b_nCP;
  TBranch        *b_cp_eta;
  TBranch        *b_cp_phi;
  TBranch        *b_cp_energy;
  TBranch        *b_nSC;
  TBranch        *b_sc_eta;
  TBranch        *b_sc_phi;
  TBranch        *b_sc_energy;
  TBranch        *b_sc_etaAtB;
  TBranch        *b_sc_phiAtB;
  TBranch        *b_sc_energyAtB;
  TBranch        *b_sc_zAtB;
  TBranch        *b_sc_pdgid;
  TBranch        *b_nTS;
  TBranch        *b_ts_firstLayer;
  TBranch        *b_ts_lastLayer;
  TBranch        *b_ts_energy;
  TBranch        *b_ts_eta_fromLC;
  TBranch        *b_ts_phi_fromLC;
  TBranch        *b_nLC;
  TBranch        *b_lc_idx;
  TBranch        *b_lc_TSidx;
  TBranch        *b_lc_energy;
  TBranch        *b_lc_seedEta;
  TBranch        *b_lc_seedPhi;
  TBranch        *b_lc_layer;
  TBranch        *b_lc_nrechits;
  TBranch        *b_lc_isSi;
  TBranch        *b_lc_mult;
  TBranch        *b_lc_tsMult;
  TBranch        *b_lc_nSC;
  
  TBranch        *b_ts_photon_proba;   //!
  TBranch        *b_ts_ele_proba;   //!
  TBranch        *b_ts_mu_proba;   //!
  TBranch        *b_ts_pi0_proba;   //!
  TBranch        *b_ts_chHad_proba;   //!
  TBranch        *b_ts_neHad_proba;   //!
  TBranch        *b_ts_ambg_proba;   //!
  TBranch        *b_ts_unkwn_proba;   //!
  std::vector<double>        *ts_photon_proba = 0;
  std::vector<double>        *ts_ele_proba = 0;
  std::vector<double>        *ts_mu_proba = 0;
  std::vector<double>        *ts_pi0_proba = 0;
  std::vector<double>        *ts_chHad_proba = 0;
  std::vector<double>        *ts_neHad_proba = 0;
  std::vector<double>        *ts_ambg_proba = 0;
  std::vector<double>        *ts_unkwn_proba = 0;
   
  Int_t           nCP = 0;
  vector<double>        *cp_energy = 0;
  vector<double>        *cp_eta = 0;
  vector<double>        *cp_phi = 0;
  Int_t           nSC = 0;
  vector<double>        *sc_energy = 0;
  vector<double>        *sc_eta = 0;
  vector<double>        *sc_phi = 0;
  vector<double>        *sc_energyAtB = 0;
  vector<double>        *sc_etaAtB = 0;
  vector<double>        *sc_phiAtB = 0;
  vector<double>        *sc_zAtB = 0;
  vector<int>           *sc_pdgid = 0;
  Int_t           nTS = 0;
  vector<int>        *ts_firstLayer = 0;
  vector<int>        *ts_lastLayer = 0;
  vector<double>        *ts_energy = 0;
  vector<double>        *ts_eta = 0;
  vector<double>        *ts_phi = 0;
  Int_t           nLC = 0;
  vector<double>  *lc_energy = 0;
  vector<double>  *lc_seedEta = 0;
  vector<double>  *lc_seedPhi = 0;
  vector<int>     *lc_layer = 0;
  vector<int>     *lc_nrechits = 0;
  vector<int>     *lc_isSi = 0;
  vector<int>     *lc_mult = 0;
  vector<int>     *lc_idx = 0;
  vector<double>     *lc_tsMult = 0;
  vector<int>     *lc_TSidx = 0;
  vector<int>     *lc_nSC = 0;

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

  tree->SetBranchAddress("nTS", &nTS, &b_nTS);
  tree->SetBranchAddress("ts_firstLayer", &ts_firstLayer, &b_ts_firstLayer);
  tree->SetBranchAddress("ts_lastLayer", &ts_lastLayer, &b_ts_lastLayer);
  tree->SetBranchAddress("ts_energy", &ts_energy, &b_ts_energy);
  tree->SetBranchAddress("ts_eta_fromLC", &ts_eta, &b_ts_eta_fromLC);
  tree->SetBranchAddress("ts_phi_fromLC", &ts_phi, &b_ts_phi_fromLC);
  tree->SetBranchAddress("ts_photon_proba", &ts_photon_proba, &b_ts_photon_proba);
  tree->SetBranchAddress("ts_ele_proba", &ts_ele_proba, &b_ts_ele_proba);
  tree->SetBranchAddress("ts_mu_proba", &ts_mu_proba, &b_ts_mu_proba);
  tree->SetBranchAddress("ts_pi0_proba", &ts_pi0_proba, &b_ts_pi0_proba);
  tree->SetBranchAddress("ts_chHad_proba", &ts_chHad_proba, &b_ts_chHad_proba);
  tree->SetBranchAddress("ts_neHad_proba", &ts_neHad_proba, &b_ts_neHad_proba);
  tree->SetBranchAddress("ts_ambg_proba", &ts_ambg_proba, &b_ts_ambg_proba);
  tree->SetBranchAddress("ts_unkwn_proba", &ts_unkwn_proba, &b_ts_unkwn_proba);

  
  tree->SetBranchAddress("nCP", &nCP, &b_nCP);
  tree->SetBranchAddress("cp_energy", &cp_energy, &b_cp_energy);
  tree->SetBranchAddress("cp_eta", &cp_eta, &b_cp_eta);
  tree->SetBranchAddress("cp_phi", &cp_phi, &b_cp_phi);
  tree->SetBranchAddress("nSC", &nSC, &b_nSC);
  tree->SetBranchAddress("sc_energy", &sc_energy, &b_sc_energy);
  tree->SetBranchAddress("sc_eta", &sc_eta, &b_sc_eta);
  tree->SetBranchAddress("sc_phi", &sc_phi, &b_sc_phi);
  tree->SetBranchAddress("sc_energyAtB", &sc_energyAtB, &b_sc_energyAtB);
  tree->SetBranchAddress("sc_etaAtB", &sc_etaAtB, &b_sc_etaAtB);
  tree->SetBranchAddress("sc_phiAtB", &sc_phiAtB, &b_sc_phiAtB);
  tree->SetBranchAddress("sc_zAtB", &sc_zAtB, &b_sc_zAtB);
  tree->SetBranchAddress("sc_pdgid", &sc_pdgid, &b_sc_pdgid);
  tree->SetBranchAddress("nLC", &nLC, &b_nLC);
  tree->SetBranchAddress("lc_idx", &lc_idx, &b_lc_idx);
  tree->SetBranchAddress("lc_TSidx", &lc_TSidx, &b_lc_TSidx);
  tree->SetBranchAddress("lc_energy", &lc_energy, &b_lc_energy);
  tree->SetBranchAddress("lc_seedEta", &lc_seedEta, &b_lc_seedEta);
  tree->SetBranchAddress("lc_seedPhi", &lc_seedPhi, &b_lc_seedPhi);
  tree->SetBranchAddress("lc_layer", &lc_layer, &b_lc_layer);
  tree->SetBranchAddress("lc_nrechits", &lc_nrechits, &b_lc_nrechits);
  tree->SetBranchAddress("lc_isSi", &lc_isSi, &b_lc_isSi);
  tree->SetBranchAddress("lc_mult", &lc_mult, &b_lc_mult);
  tree->SetBranchAddress("lc_tsMult", &lc_tsMult, &b_lc_tsMult);
  tree->SetBranchAddress("lc_nSC", &lc_nSC, &b_lc_nSC);
  
  const int nEntries = tree->GetEntries();
  
  for (unsigned iE(0); iE < nEntries; ++iE){//loop on entries
    Long64_t ientry = tree->LoadTree(iE); 
    if (iE%100==0)
      std::cout << " -- Processing entry " << iE  << " " << ientry << std::endl;

    nSC = 0;
    
    tree->GetEntry(iE);
    treeLC->GetEntry(iE);

    //std::cout << " Evt " << iE << " found " << nAllLC << " LC and " << nTS << " TS. TS first layer = " << (*ts_firstLayer)[0] << " CP eta,phi = " << (*cp_eta)[0] << "," << (*cp_phi)[0] << std::endl;

    if (nCP!=1){
      std::cout << " -- Did not find 1 CP but " << nCP << " passing event." << std::endl;
      //continue;
    }


    /*    if (nSC!=nTS){
      std::cout << " -- Different: nSC = " << nSC << " nTS = " << nTS << std::endl;
      std::sort((*sc_energy).begin(),(*sc_energy).end(), std::greater<double>());
      std::sort((*ts_energy).begin(),(*ts_energy).end(), std::greater<double>());
      for (int isc(0); isc<nSC; ++isc){
	std::cout << isc << " SC " << (*sc_energy)[isc]
		  << " eta " << (*sc_eta)[isc]
		  << " phi " << (*sc_phi)[isc]
		  << " " << (*sc_pdgid)[isc];
	if (isc<nTS) std::cout << " TS " << (*ts_energy)[isc]
			       << " eta " << (*ts_eta)[isc]
			       << " phi " << (*ts_phi)[isc]
			       << " g " << (*ts_photon_proba)[isc]
			       << " e " << (*ts_ele_proba)[isc]
			       << " mu " << (*ts_mu_proba)[isc]
			       << " pi0 " << (*ts_pi0_proba)[isc]
			       << " chHad " << (*ts_chHad_proba)[isc]
			       << " neHad " << (*ts_neHad_proba)[isc]
			       << " ambg " << (*ts_ambg_proba)[isc]
			       << " unk " << (*ts_unkwn_proba)[isc]
			       << std::endl;
	else std::cout << std::endl;
      }
      }*/
    
    eTotLC = 0;
    eTotLC3 = 0;
    //std::map<unsigned,double>checkMult;
    //std::pair<std::map<unsigned,double>::iterator,bool> checkMultInsert;
    
    
    for (int iLC(0); iLC<nLC; ++iLC){
      if ((*lc_tsMult)[iLC]>0) eTotLC += (*lc_energy)[iLC]/(*lc_tsMult)[iLC];
      //checkMultInsert = checkMult.insert(std::pair<unsigned,double>((*lc_idx)[iLC],1./(*lc_tsMult)[iLC]));
      //if (!checkMultInsert.second) checkMultInsert.first->second += 1./(*lc_tsMult)[iLC];

      if ((*lc_isSi)[iLC]>0 && (*lc_nrechits)[iLC]<3) continue;
      if ((*lc_tsMult)[iLC]>0) eTotLC3 += (*lc_energy)[iLC]/(*lc_tsMult)[iLC];
      else eTotLC3 += (*lc_energy)[iLC];
    }

    //std::cout << " Check map" << std::endl;
    //std::map<unsigned,double>::iterator checkMultIter = checkMult.begin();
    //for (;checkMultIter!=checkMult.end();++checkMultIter){
    //if (checkMultIter->second>1.01) std::cout << " event " << ientry << " lc " << checkMultIter->first << " sumFrac " << checkMultIter->second << std::endl;
    //}
    
    std::map<int,int>lMapLC;
    lMapLC.clear();
    eTotAllLC = 0;
    eTotAllLC2 = 0;
    eTotAllLC3 = 0;
    nAllLC2 = 0;
    nAllLC3 = 0;

    for (int iLC(0); iLC<nAllLC; ++iLC){
      eTotAllLC += (*all_lc_energy)[iLC];
      if ( ((*all_lc_isSi)[iLC]>0 && (*all_lc_nrechits)[iLC]>1) || (*all_lc_isSi)[iLC]==0){
	eTotAllLC2 += (*all_lc_energy)[iLC];
	nAllLC2++;
      }
      if ((*all_lc_isSi)[iLC]>0 && (*all_lc_nrechits)[iLC]<3) continue;
      eTotAllLC3 += (*all_lc_energy)[iLC];
      nAllLC3++;
      std::pair<std::map<int,int>::iterator,bool> isInserted =  lMapLC.insert(std::pair<int,int>((*all_lc_layer)[iLC],1));
      if (!isInserted.second) isInserted.first->second += 1; 
    }
    
    std::vector<unsigned> new_lc_mult;
    for (unsigned iL(0); iL<49;++iL){
      std::map<int,int>::iterator lEle = lMapLC.find(iL+1);
      if (lEle !=lMapLC.end()) new_lc_mult.push_back(lEle->second);
      else new_lc_mult.push_back(0);
    }

    //std::cout << " -- Filled " <<  new_lc_mult.size() << " layers." << std::endl;
    
    //loop on SC
    eTotSC = 0;
    eTotSCAtB = 0;
    eTotSC_noPi = 0;
    for (unsigned iC(0); iC<nCat;++iC){
      eFrac[iC] = 0;
      eFracAtB[iC] = 0;
    }
    for (unsigned iC(0); iC<nC;++iC){
      nSC_eCut[iC] = 0;
      eTotSC_eCut[iC] = 0;
      eTotSCAtB_eCut[iC] = 0;
      sc_dRtoCP_eCut[iC].clear();
      sc_energy_eCut[iC].clear();
      sc_eta_eCut[iC].clear();
      sc_phi_eCut[iC].clear();
      sc_energyAtB_eCut[iC].clear();
      sc_etaAtB_eCut[iC].clear();
      sc_phiAtB_eCut[iC].clear();
      sc_zAtB_eCut[iC].clear();
      sc_pdgid_eCut[iC].clear();
      sc_pdgcat_eCut[iC].clear();
      sc_idx_eCut[iC].clear();
      sc_tsIdx_eCut[iC].clear();
      sc_tsminDR_eCut[iC].clear();
    }


    //std::cout << " -- Initialised SC. nSC = " << nSC << std::endl;

    for (int isc(0); isc<nSC; ++isc){
      int absPid = abs((*sc_pdgid)[isc]);

      //cut off particles re-entering the calo...
      //if ((*sc_zAtB)[isc]<0) continue;
      //cut knock-off nucleons not at front-face...
      //if ( (absPid == 2212 || absPid == 2112 ) && (*sc_zAtB)[isc]>321) continue;
      //std::cout << " sc " << isc << " " << absPid << std::flush;
      std::map<int,std::pair<int,std::string> >::iterator lElement;
      lElement = pdgMap.find(absPid);
     
      unsigned pdgCat = lElement!=pdgMap.end() ? getPDGcat((*lElement)): 10;
      //std::cout << " " << pdgCat << std::endl;
      
      if (pdgCat > 9) std::cout << " !!! unknow : " << absPid << std::endl;
      eFrac[pdgCat-1] += (*sc_energy)[isc];
      eFracAtB[pdgCat-1] += (*sc_energyAtB)[isc];
      
      
      double dRtoCP = sqrt(pow(TMath::Abs((*sc_eta)[isc]-(*cp_eta)[0]),2)+pow(deltaPhi((*sc_phi)[isc],(*cp_phi)[0]),2));
      eTotSC += (*sc_energy)[isc];
      eTotSCAtB += (*sc_energyAtB)[isc];
      
      
      for (unsigned iC(0); iC<nC;++iC){
	if ((*sc_energy)[isc]/(*cp_energy)[0]>eCut[iC]){
	  eTotSCAtB_eCut[iC] += (*sc_energyAtB)[isc];
	  sc_energyAtB_eCut[iC].push_back((*sc_energyAtB)[isc]);
	  sc_etaAtB_eCut[iC].push_back((*sc_etaAtB)[isc]);
	  sc_phiAtB_eCut[iC].push_back((*sc_phiAtB)[isc]);
	  sc_zAtB_eCut[iC].push_back((*sc_zAtB)[isc]);
	  sc_pdgcat_eCut[iC].push_back(pdgCat);

	  sc_dRtoCP_eCut[iC].push_back(dRtoCP);
	  sc_idx_eCut[iC].push_back(isc);
	  nSC_eCut[iC]++;
	  eTotSC_eCut[iC] += (*sc_energy)[isc];
	  sc_energy_eCut[iC].push_back((*sc_energy)[isc]);
	  sc_eta_eCut[iC].push_back((*sc_eta)[isc]);
	  sc_phi_eCut[iC].push_back((*sc_phi)[isc]);
	  sc_pdgid_eCut[iC].push_back((*sc_pdgid)[isc]);

	  
	  //loop on TS
	  double mindr = 10;
	  int idxMindR = -1;
	  for (int its(0); its<nTS; ++its){
	    double dr = sqrt(pow(TMath::Abs((*sc_eta)[isc]-(*ts_eta)[its]),2)+pow(deltaPhi((*sc_phi)[isc],(*ts_phi)[its]),2));
	    if (dr<mindr){
	      mindr=dr;
	      idxMindR = its;
	    }
	  }//loop on ts

	  /*int tmpPID = 0;
	  if ((*ts_photon_proba)[idxMindR]>0.5) tmpPID=22; 
	  else if ((*ts_ele_proba)[idxMindR]>0.5) tmpPID=11; 
	  else if ((*ts_mu_proba)[idxMindR]>0.5) tmpPID=13; 
	  else if ((*ts_pi0_proba)[idxMindR]>0.5) tmpPID=110; 
	  else if ((*ts_chHad_proba)[idxMindR]>0.5) tmpPID=211; 
	  else if ((*ts_neHad_proba)[idxMindR]>0.5) tmpPID=130; 
	  */
	  
	  sc_tsIdx_eCut[iC].push_back(idxMindR);
	  sc_tsminDR_eCut[iC].push_back(mindr);
	}
      }//loop on ecut
    }//loop on SC

    for (unsigned iC(0); iC<nCat;++iC){
      eFrac[iC] = eFrac[iC]/eTotSC;
      if (eTotSCAtB>0) eFracAtB[iC] = eFracAtB[iC]/eTotSCAtB;
      else eFracAtB[iC] = 0;
    }

    //remove duplicated energy from pion ...
    eTotSC_noPi = eTotSC;
    if (eTotSC>1.5*(*cp_energy)[0]){
      for (int isc(0); isc<nSC; ++isc){
	if ((*sc_pdgid)[isc]==211 && (*sc_energy)[isc] > 0.5*(*cp_energy)[0]){
	  eTotSC_noPi -= (*sc_energy)[isc];
	  break;
	}
      }
    }

    //std::cout << " -- Done with SC. Moving to TS" << std::endl;
    
    //separate in blobs / missing layers
    nTS_eCut10 = 0;
    ts_energy_eCut10.clear();

    nBlobs = 0;
    EBlobs = 0;
    EBlobs_ECAL = 0;
    EBlobs_FHCAL = 0;
    EBlobs_BHCAL = 0;
    EBlobs_HoverE = 0;
    minMiss = 0;
    maxMiss = 0;
    maxEfracBlobs = 0;
    minEfracBlobs = 0;
    maxLengthBlobs = 0;
    minLengthBlobs = 0;
    energy_blobs.clear();
    eFrac_blobs.clear();
    length_blobs.clear();
    firstL_blobs.clear();
    nLC_blobs.clear();
    lc_in_blob.clear();
    nMiss_blobs.clear();
    

    //std::cout << " Found " << nTS << " TS " << std::endl;
    for (int its(0); its<nTS; ++its){
      //discard too low energy/too split content....
      if ((*ts_energy)[its]/(*cp_energy)[0]<0.1) continue;
      nTS_eCut10++;
      ts_energy_eCut10.push_back((*ts_energy)[its]);
    }

    
    std::vector<unsigned int> uniqueLayerIds;
    uniqueLayerIds.reserve(60);
    //loop on lcs
    for (int iLC(0); iLC<nAllLC; ++iLC){
      lc_in_blob.push_back(0);
      if ((*all_lc_isSi)[iLC]>0 && (*all_lc_nrechits)[iLC]<3) continue;
      unsigned layerId = (*all_lc_layer)[iLC];
      if (new_lc_mult[layerId-1]>0) uniqueLayerIds.push_back(layerId);
    }

    lMapLC.clear();
    
    std::sort(uniqueLayerIds.begin(), uniqueLayerIds.end());
    uniqueLayerIds.erase(std::unique(uniqueLayerIds.begin(), uniqueLayerIds.end()), uniqueLayerIds.end());
    
    std::vector<unsigned int> idxInVec;
    
    if (uniqueLayerIds.size()==0) {
      std::cout << " -- Problem, no LC ?? All: " << nAllLC << " nRH>2: " << nAllLC3 << std::endl;
      //return;
      newTree->Fill();
      
      fileOut << iE << " " << ientry << " "
	      << nSC_eCut[nC-1] << " "
	      << nTS_eCut10<< " "
	      << nSC_eCut[0] << " "
	      << nTS << " "
	      << nBlobs << " "
	      << EBlobs << " "
	      << eTotAllLC3 << " ";
      for (unsigned iC(0); iC<nCat;++iC){
	fileOut << eFrac[iC] << " "
		<< eFracAtB[iC];
	if (iC<nCat-1) fileOut << " ";
      }
      
      fileOut << std::endl;
      continue;
    }
    unsigned int j = uniqueLayerIds[0];
    unsigned int indexInVec = 0;
    
    //decompose candidate into continuous blobs and missing chunks
    std::vector<unsigned int> nContinuous;
    std::vector<unsigned int> nMissing;
    
    unsigned int tmpMiss = 0;
    unsigned int tmpCont = 0;
    for (const auto &layer : uniqueLayerIds) {
      if (layer > j) {
	nContinuous.push_back(tmpCont);
	idxInVec.push_back(indexInVec-tmpCont);
	tmpCont = 0;
	while (layer != j){
	  tmpMiss++;
	  j++;
	}
      }
      if (tmpCont==0 && tmpMiss>0){
	nMissing.push_back(tmpMiss);
	tmpMiss = 0;
      }
      tmpCont++;
      indexInVec++;
      j++;
    }
    //push last one....
    nContinuous.push_back(tmpCont);
    idxInVec.push_back(indexInVec-tmpCont);
    nBlobs = nContinuous.size();
    if (nMissing.size()>0){
      minMiss = *(std::min_element(nMissing.begin(), nMissing.end()));
      maxMiss = *(std::max_element(nMissing.begin(), nMissing.end()));
    }
    else {
      minMiss = 0;
      maxMiss = 0;
    }
    //loop on blobs
    if (nBlobs==1){
      EBlobs = eTotAllLC3;
      double EH = 0;
      for (int iLC(0); iLC<nAllLC; ++iLC){
	if ((*all_lc_isSi)[iLC]>0 && (*all_lc_nrechits)[iLC]<3) continue;
	if ((*all_lc_layer)[iLC]>nLEE && (*all_lc_layer)[iLC]<=nLEE+5) EH += (*all_lc_energy)[iLC];
	if ((*all_lc_layer)[iLC]<=nLEE) EBlobs_ECAL += (*all_lc_energy)[iLC];
	else if ((*all_lc_isSi)[iLC]>0) EBlobs_FHCAL += (*all_lc_energy)[iLC];
	else EBlobs_BHCAL += (*all_lc_energy)[iLC];
	lc_in_blob[iLC] = 1;
      }
      if (EBlobs_ECAL>0) EBlobs_HoverE = EH/EBlobs_ECAL;
      else EBlobs_HoverE = 100;
      maxEfracBlobs = 1;
      minEfracBlobs = 1;
      maxLengthBlobs = nContinuous[0];
      minLengthBlobs = nContinuous[0];
      energy_blobs.push_back(eTotAllLC3);
      nLC_blobs.push_back(nAllLC3);
      eFrac_blobs.push_back(1);
      length_blobs.push_back(nContinuous[0]);
      firstL_blobs.push_back(uniqueLayerIds[0]);
      nMiss_blobs.push_back(0);
    } else if (nBlobs>0) {
      EBlobs = 0;
      EBlobs_ECAL = 0;
      EBlobs_FHCAL = 0;
      EBlobs_BHCAL = 0;
      EBlobs_HoverE = 0;
      double EH = 0;
      double maxEfrac = 0;
      double minEfrac = 1;
      double maxLength = 0;
      double minLength = 60;
      unsigned nSel = 0;
      for (unsigned iC(0); iC<nBlobs;++iC){
	unsigned length = nContinuous[iC];
	if (length < 5) continue;
	nSel++;
	double Eblob = 0;
	unsigned firstL=uniqueLayerIds[idxInVec[iC]];
	unsigned nLCblob = 0;
	for (int iLC(0); iLC<nAllLC; ++iLC){
	  if ((*all_lc_isSi)[iLC]>0 && (*all_lc_nrechits)[iLC]<3) continue;
	  if ((*all_lc_layer)[iLC]<firstL ||
	      (*all_lc_layer)[iLC]>=firstL+length) continue;
	  Eblob += (*all_lc_energy)[iLC];
	  lc_in_blob[iLC] = 1;
	  nLCblob++;
	  if ((*all_lc_layer)[iLC]>nLEE && (*all_lc_layer)[iLC]<=nLEE+5) EH += (*all_lc_energy)[iLC];
	  if ((*all_lc_layer)[iLC]<=nLEE) EBlobs_ECAL += (*all_lc_energy)[iLC];
	  else if ((*all_lc_isSi)[iLC]>0) EBlobs_FHCAL += (*all_lc_energy)[iLC];
	  else EBlobs_BHCAL += (*all_lc_energy)[iLC];
	}
	double eFrac = Eblob/eTotAllLC3;
	if (eFrac > maxEfrac) maxEfrac = eFrac;
	if (eFrac < minEfrac) minEfrac = eFrac;
	EBlobs += Eblob;
	if (length>maxLength) maxLength = length;
	if (length<minLength) minLength = length;
	energy_blobs.push_back(Eblob);
	nLC_blobs.push_back(nLCblob);
	eFrac_blobs.push_back(eFrac);
	length_blobs.push_back(length);
	firstL_blobs.push_back(firstL);
	if (iC<(nBlobs-1)) nMiss_blobs.push_back(nMissing[iC]);
	else nMiss_blobs.push_back(0);
      }
      if (EBlobs_ECAL>0) EBlobs_HoverE = EH/EBlobs_ECAL;
      else EBlobs_HoverE = 100;
      //redefine n after selection...
      nBlobs = nSel;
      minEfracBlobs = minEfrac;
      maxEfracBlobs = maxEfrac;
      minLengthBlobs = minLength;
      maxLengthBlobs = maxLength;
      //if (nBlobs==2 && minEfracBlobs>0.4) std::cout << " Check: " << ientry << " " << eFrac_blobs[0] << " " << eFrac_blobs[1] << std::endl;
    }


    
    newTree->Fill();

    fileOut << iE << " " << ientry << " "
	    << nSC_eCut[nC-1] << " "
	    << nTS_eCut10<< " "
	    << nSC_eCut[0] << " "
	    << nTS << " "
	    << nBlobs << " "
	    << EBlobs << " "
	    << eTotAllLC3 << " ";
    for (unsigned iC(0); iC<nCat;++iC){
      fileOut << eFrac[iC] << " "
	      << eFracAtB[iC];
      if (iC<nCat-1) fileOut << " ";
    }
    
    fileOut << std::endl;
    
    //std::cout << " --- Fill histos " << std::endl;	

    //std::cout << " --- End of event " << std::endl;	
  }//loop on entries
  
  fout->Write();
  
  tree->Delete();
  fileOut.close();

}//fill

int exploreTruthInfo(){
  
  bool doDebug = false;

  //std::string filePath = "D49_AllTracksters";
  //std::string filePath = "D49_FineCalo_AMiters";
  //std::string filePath = "D49_DevMaxMiss";
  std::string filePath = "D49_DefSC";
  
  const unsigned nR = 5;

  const unsigned regIdx[5] = {0,1,4,5,6};
  //std::string iter[11] = {"Dummy1","Dummy2","Dummy3","EM1","EM2","EM3","TrkEM","EM","Trk","HAD","Sim"};
  
  
  for (unsigned iR(nR-1); iR<nR; ++iR){

    const unsigned iReg = regIdx[iR];
    std::cout << " -- Processing region " << reg[iReg] << std::endl;
    
    double ptval[12] = {3,5,10,15,20,30,40,50,75,100,150,200};
    const unsigned nPT = doDebug ? 1 : iReg==4 ? 7 : 12;
    if (iReg==4 || iReg==5){
      ptval[1] = 10;
      ptval[2] = 20;
      ptval[3] = 50;
      ptval[4] = 100;
      ptval[5] = 150;
      ptval[6] = 200;
    }
    double etaval[6] = {17,19,21,23,25,27};
    const unsigned nEta = iReg==0? 6 : 1;
    if (iReg!=0) etaval[0] = 21;
    
    if (doDebug) {
      ptval[0] = 50;
      etaval[0] = 21;
    }
    
    //std::string iter[17] = {"Dummy1","Dummy2","Dummy3","TRK1","TRK2","TRK3","EM1","EM2","EM3","HAD1","HAD2","HAD3","TrkEM","EM","Trk","HAD","Sim"};
    std::string iter[5] = {"Sim","TrkEM","EM","Trk","HAD"};
    const unsigned nI = 1;
    
    for (unsigned iI(0); iI<nI; ++iI){
      std::cout << " -- Processing iter " << iter[iI] << std::endl;
      for (unsigned ipt(0); ipt<nPT; ++ipt){  
	for (unsigned ieta(0); ieta<nEta; ++ieta){
	  std::cout << " --- Processing pT " << ptval[ipt] << " eta " << etaval[ieta]<< std::endl;
	  
	  std::ostringstream pteta;
	  pteta << "pt" << ptval[ipt] << "_eta" << etaval[ieta];
	  fillTruthTree(filePath,iReg,iter[iI],pteta.str(),doDebug);
	}//loop on eta
      }//loop on pT
    }

  }//loop on reg
  return 0;
};


