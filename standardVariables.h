#include <iostream>
#include <sstream>
#include <fstream>
#include <map>

struct histo{
  std::string var;
  std::string varShort;
  std::string varName;
  unsigned nB;
  double binMin;
  double binMax;
  std::vector<double> binsX;
  bool logy;
  bool underflow;
  bool doFixedBinning;
};

static std::map<std::string,histo> histMap;

void print(const histo & lHist){
  std::cout << " --- Adding histo: " << lHist.varShort << " " << lHist.var << std::endl;
};

void addVariable(std::string var,
		 std::string varShort,
		 std::string varName,
		 unsigned nB,		 
		 double binMin,
		 double binMax,
		 const std::vector<double> & binsX,
		 bool logy,
		 bool underflow,
		 bool doFixedBinning=false){
  histo tmp;
  tmp.var = var;
  tmp.varShort = varShort;
  tmp.varName = varName;
  tmp.nB = nB;
  tmp.binMin = binMin;
  tmp.binMax = binMax;
  tmp.binsX = binsX;
  tmp.logy = logy;
  tmp.underflow = underflow;
  tmp.doFixedBinning = doFixedBinning;
  histMap[varShort] = tmp;
  //print(tmp);
};

void fillMap(const bool doTriplets = true, const bool lowStat = false){
  

  //std::vector<double> detaBins{0,1,2,3,4,5,5.5,6,7.2,7.7,9.4};
  //addVariable(iT==0?"dijet_dEta":"lMjj_dijet_deta","dijet_dEta",";#Delta#eta_{jj};Events",10,1,8,detaBins,1,0,true);


  std::vector<double> dummyX;
  
  addVariable("nTS","nTS", ";number of reco TS;Events",53,-1.5,51.5,dummyX,1,1);
  addVariable("ts_nLC","ts_nLC", ";number of reco LC per TS;Photons",103,-1.5,101.5,dummyX,1,1);
  addVariable("nCP","nCP", ";number of CP;Events",5,-1.5,3.5,dummyX,1,1);
  addVariable("cp_nSC","cp_nSC", ";number of SC per CP;Photons",13,-1.5,11.5,dummyX,1,1);
  addVariable("abs(cp_pdgid)","cp_abspdgid", ";|PDG id| of CP;Photons",400,0,400,dummyX,1,1);
  addVariable("nSC","nSC", ";number of SC;Events",18,-1.5,16.5,dummyX,1,1);
  addVariable("sc_CPidx","sc_CPidx", ";SC idx of CP;Photons",10,-1.5,8.5,dummyX,1,1);
  addVariable("abs(sc_pdgid)","sc_abspdgid", ";|PDG id| of SC;Photons",400,0,400,dummyX,1,1);

  //addVariable("trackster","TSidx", ";TS index;Photons",23,-1.5,21.5,dummyX,1,1);
  addVariable("ts_energy/cp_energy","TSEoverCPE", ";E_{TS}/E_{CP};Photons",100,0,2,dummyX,1,1);
  addVariable("ts_regEnergy/cp_energy","TSRegEoverCPE", ";E^{reg}_{TS}/E_{CP};Photons",100,0,2,dummyX,1,1);
  addVariable("ts_emEnergy/ts_energy","TSEMoverTSE", ";E^{EM}_{TS}/E_{TS};Photons",100,0.6,1.2,dummyX,1,1);
  addVariable("TMath::Abs(ts_eta_PCA-cp_eta)","TSCPEtaPCADiff",";|#eta_{TS}^{PCA}-#eta_{CP}|;Photons",40,0,0.4,dummyX,1,1);
  addVariable("TMath::Abs(ts_eta_fromLC-cp_eta)","TSCPEtaLCDiff",";|#eta_{TS}^{LC}-#eta_{CP}|;Photons",40,0,0.4,dummyX,1,1);
  addVariable("deltaPhi(ts_phi_PCA,cp_phi)","TSCPPhiPCADiff",";|#phi_{TS}^{PCA}-#phi_{CP}|;Photons",40,0,0.4,dummyX,1,1);
  addVariable("deltaPhi(ts_phi_fromLC,cp_phi)","TSCPPhiLCDiff",";|#phi_{TS}^{LC}-#phi_{CP}|;Photons",40,0,0.4,dummyX,1,1);
  addVariable("sqrt(pow(TMath::Abs(cp_eta-ts_eta_PCA),2)+pow(deltaPhi(cp_phi,ts_phi_PCA),2))","TSCPDeltaRPCA", ";#DeltaR(CP,TS-PCA);Photons",50,0,1,dummyX,1,1);
  addVariable("sqrt(pow(TMath::Abs(cp_eta-ts_eta_fromLC),2)+pow(deltaPhi(cp_phi,ts_phi_fromLC),2))","TSCPDeltaRLC", ";#DeltaR(CP,TS-LC);Photons",50,0,0.1,dummyX,1,1);
  addVariable("ts_photon_proba","TSphotonProba",";#gamma proba;Photons",100,0,1,dummyX,1,1);
  addVariable("ts_ele_proba","TSeleProba",";e proba;Photons",100,0,1,dummyX,1,1);
  addVariable("ts_mu_proba","TSmuProba",";#mu proba;Photons",100,0,1,dummyX,1,1);
  addVariable("ts_pi0_proba","TSpi0Proba",";#pi^{0} proba;Photons",100,0,1,dummyX,1,1);
  addVariable("ts_neHad_proba","TSneHadProba",";neutral Had proba;Photons",100,0,1,dummyX,1,1);
  addVariable("ts_chHad_proba","TSchHadProba",";charged Had proba;Photons",100,0,1,dummyX,1,1);
  addVariable("ts_ambg_proba","TSambgProba",";ambiguous proba;Photons",100,0,1,dummyX,1,1);
  addVariable("ts_unkwn_proba","TSunkwnProba",";unknown proba;Photons",100,0,1,dummyX,1,1);
  addVariable("ts_firstLayer","TSfirstLayer",";first layer;Photons",30,-0.5,29.5,dummyX,1,1);
  addVariable("ts_lastLayer","TSlastLayer",";last layer;Photons",30,-0.5,29.5,dummyX,1,1);
  addVariable("ts_lastLayer-ts_firstLayer+1","TSlength",";TS length (Nlayers);Photons",30,-0.5,29.5,dummyX,1,1);
  addVariable("ts_outInHopsPerformed","TSoutInHopsPerformed",";outInHopsPerformed;Photons",12,-1.5,10.5,dummyX,1,1);
  addVariable("cp_missingEnergyFraction","MissEtFrac", ";|E_{CP}-E_{LCrechitsMatched}|/E_{CP};Photons",110,0,1.1,dummyX,1,1);

  addVariable("ts_sigma1","TSsigma1", ";#sigma^{1}_{TS} (cm);Photons",100,0,50,dummyX,1,1);
  addVariable("ts_sigma2","TSsigma2", ";#sigma^{2}_{TS} (cm);Photons",100,0,20,dummyX,1,1);
  addVariable("ts_sigma3","TSsigma3", ";#sigma^{3}_{TS} (cm);Photons",100,0,20,dummyX,1,1);
  addVariable("ts_BCx","TSBaryCenterX", ";x^{BC}_{TS} (cm);Photons",100,-300,300,dummyX,1,1);
  addVariable("ts_BCy","TSBaryCenterY", ";y^{BC}_{TS} (cm);Photons",100,-300,300,dummyX,1,1);
  addVariable("ts_BCz","TSBaryCenterZ", ";z^{BC}_{TS} (cm);Photons",80,320,400,dummyX,1,1);

  /*

  const unsigned nL = 29;
  for (unsigned iL(0); iL<nL; ++iL){//loop on layers
    std::ostringstream label;
    label << "_layer" << iL+1;
    
    addVariable("lc_energy/ts_energy",("LCEoverTSE"+label.str()).c_str(), ";E_{LC}/E_{TS};Photons",100,0,0.5,dummyX,1,1);
    addVariable("TMath::Abs(lc_eta-ts_eta_PCA)",("LCTSEtaDiff"+label.str()).c_str(), ";|#eta_{LC}-#eta_{TS}^{PCA}|;Photons",100,0,0.7,dummyX,1,1);
    addVariable("deltaPhi(lc_phi,ts_phi_PCA)",("LCTSPhiDiff"+label.str()).c_str(), ";|#phi_{LC}-#phi_{TS}^{PCA}|;Photons",100,0,0.7,dummyX,1,1);
    addVariable("sqrt(pow(TMath::Abs(lc_eta-ts_eta_PCA),2)+pow(deltaPhi(lc_phi,ts_phi_PCA),2))",("LCTSDeltaR"+label.str()).c_str(), ";#DeltaR(LC,TS-PCA);Photons",100,0,0.7,dummyX,1,1);

    addVariable("TMath::Abs(lc_eta-cp_eta)",("LCCPEtaDiff"+label.str()).c_str(), ";|#eta_{LC}-#eta_{CP}|;Photons",100,0,0.7,dummyX,1,1);
    addVariable("deltaPhi(lc_phi,cp_phi)",("LCCPPhiDiff"+label.str()).c_str(), ";|#phi_{LC}-#phi_{CP}|;Photons",100,0,0.7,dummyX,1,1);
    addVariable("sqrt(pow(TMath::Abs(lc_eta-cp_eta),2)+pow(deltaPhi(lc_phi,cp_phi),2))",("LCCPDeltaR"+label.str()).c_str(), ";#DeltaR(LC,CP);Photons",100,0,0.7,dummyX,1,1);

    
    //addVariable("lc_layer",("LCLayer"+label.str()).c_str(), ";LC layer;Photons",nL,1,nL+1,dummyX,1,1);
    //addVariable("lc_algo",("LCAlgo"+label.str()).c_str(), ";LC algo;Photons",10,0,10,dummyX,1,1);
    addVariable("lc_nrechits",("LCNrechits"+label.str()).c_str(), ";N_{RecHits}^{LC};Photons",81,-0.5,80.5,dummyX,1,1);
    addVariable("lc_tsMult",("LCTSmult"+label.str()).c_str(), ";N_{TS}^{LC};Photons",21,-0.5,20.5,dummyX,1,1);

    if (iL<28 && doTriplets){
      label.str("");
      label << "_" <<iL+1;
      
      addVariable(("nTriplets"+label.str()).c_str(),("nTriplets"+label.str()).c_str(),";nTriplets;Photons",10,-0.5,9.5,dummyX,1,1);
      std::ostringstream lTmp;
      lTmp << iL+1 << "-triplets_layerA" << label.str();
      addVariable(lTmp.str().c_str(),("triplets_deltalayerA"+label.str()).c_str(),";triplets L_{B}-L_{A};Photons",5,-0.5,4.5,dummyX,1,1);
      lTmp.str("");
      lTmp << "triplets_layerC" << label.str() << "-" << iL+1;
      addVariable(lTmp.str().c_str(),("triplets_deltalayerC"+label.str()).c_str(),";triplets L_{C}-L_{B};Photons",5,-0.5,4.5,dummyX,1,1);

      addVariable(("(triplets_energyA"+label.str()+"-triplets_energyB"+label.str()+")/triplets_energyB"+label.str()).c_str(),("triplets_deltaenergyA"+label.str()).c_str(),";triplets (E_{A}-E_{B})/E_{B};Photons",200,-10,10,dummyX,1,1);
      addVariable(("(triplets_energyC"+label.str()+"-triplets_energyB"+label.str()+")/triplets_energyB"+label.str()).c_str(),("triplets_deltaenergyC"+label.str()).c_str(),";triplets (E_{C}-E_{B})/E_{B};Photons",200,-10,10,dummyX,1,1);
      addVariable(("triplets_cosBeta"+label.str()).c_str(),("triplets_cosBeta"+label.str()).c_str(),";triplets cos#beta;Photons",101,0,1.01,dummyX,1,1);
      addVariable(("triplets_cosAlphaInner"+label.str()).c_str(),("triplets_cosAlphaInner"+label.str()).c_str(),";triplets cos#alpha(inner);Photons",101,0,1.01,dummyX,1,1);
      addVariable(("triplets_cosAlphaOuter"+label.str()).c_str(),("triplets_cosAlphaOuter"+label.str()).c_str(),";triplets cos#alpha(outer);Photons",101,0,1.01,dummyX,1,1);
      
      addVariable(("triplets_inner_in_links"+label.str()).c_str(),("triplets_inner_in_links"+label.str()).c_str(),";triplets inner-in-links;Photons",5,-0.5,4.5,dummyX,1,1);
      addVariable(("triplets_inner_out_links"+label.str()).c_str(),("triplets_inner_out_links"+label.str()).c_str(),";triplets inner-out-links;Photons",5,-0.5,4.5,dummyX,1,1);
      addVariable(("triplets_outer_in_links"+label.str()).c_str(),("triplets_outer_in_links"+label.str()).c_str(),";triplets outer-in-links;Photons",5,-0.5,4.5,dummyX,1,1);
      addVariable(("triplets_outer_out_links"+label.str()).c_str(),("triplets_outer_out_links"+label.str()).c_str(),";triplets outer-out-links;Photons",5,-0.5,4.5,dummyX,1,1);
    }
    
    if (iL<28){
      std::ostringstream mult;
      mult << "lc_mult[" << iL << "]";
      addVariable(mult.str().c_str(),("LCmult"+label.str()).c_str(), ";N_{LC};Photons",21,-0.5,20.5,dummyX,1,1);


      
    }
  }//loop on layers
  */
};

