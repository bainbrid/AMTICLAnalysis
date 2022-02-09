#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <map>

#include "TCanvas.h"
#include "TPad.h"
#include "TStyle.h"
#include "Math/Vector3D.h"

int myColor(const unsigned idx){
  if (idx<4) return idx+1;
  else if (idx<8) return idx+2;
  else if (idx<13) return kRed+idx-10;
  else if (idx<18) return kCyan+idx-15;
}

int myMarker(const unsigned idx){
  return idx+20;
};

std::map<int,std::pair<int,std::string> > getPDGMap(const std::string aFileName){

  std::map<int,std::pair<int,std::string> > lMap;
  std::ifstream infile(aFileName);
  //vec.clear();
  
  if (!infile.is_open()){
    std::cout << " -- file " << aFileName << " not found." << std::endl;
    return lMap;
  }

  unsigned counter = 0;
  while(!infile.eof()){
    int id = 0;
    int ch = 0;
    std::string name;
    infile>>id>>ch>>name;
    if (id>0) lMap[id] = std::pair<int,std::string>(ch,name);
    counter++;
    if (counter > 10000) break;
  }
  
  std::cout << " -- Found " << lMap.size() << " elements in file " << aFileName << std::endl;
  return lMap;
  
};

void MyPalette()
{
  static int myPalette[10] = {1,2,3,4,6,7,kOrange,kGray,8,9};
  gStyle->SetPalette(10,myPalette);
};

unsigned getPDGcat(const std::pair<int,std::pair<int,std::string> > & aEle){
  if (aEle.first==22) return 1;
  else if (aEle.first==11) return 2;
  else if (aEle.first==13) return 3;
  else if (aEle.first==111) return 4;
  else if (aEle.first==211) return 5;
  else if (aEle.first==2112) return 6;
  else if (aEle.first==2212) return 7;
  else if (aEle.second.first==0) return 8;
  else if (aEle.second.first!=0) return 9;

  else return 10;
};

unsigned getPDGcat(int pdgid){
  std::map<int,std::pair<int,std::string> > pdgMap = getPDGMap("pdgList.dat");
  unsigned pdgCat = getPDGcat(*(pdgMap.find(abs(pdgid))));
  return pdgCat;
};

std::string pdgBinLabel(const unsigned iB){
  if (iB==1) return "gamma";
  else if (iB==2) return "e";
  else if (iB==3) return "mu";
  else if (iB==4) return "pi0";
  else if (iB==5) return "pipm";
  else if (iB==6) return "n";
  else if (iB==7) return "p";
  else if (iB==8) return "neHad";
  else if (iB==9) return "chHad";
  else return "unknown";
}


double deltaPhi(const double phi1, const double phi2){
  double delta_phi = fabs(phi2 - phi1);
  if (delta_phi > M_PI){ 
    delta_phi = 2*M_PI-delta_phi;
  } 
  return delta_phi;
};

double deltaEta(const double eta1, const double eta2){
  return TMath::Abs(eta1-eta2);
};

double deltaR(const double eta1, const double eta2, const double phi1, const double phi2){
  return sqrt(pow(deltaEta(eta1,eta2),2)+pow(deltaPhi(phi1,phi2),2));
};

std::pair<double,double> getXY(const double & z,
			       const double & eta,
			       const double & phi){

  double theta = 2*atan(exp(-eta));
  double x=z*tan(theta)*cos(phi);
  double y=z*tan(theta)*sin(phi);
  return std::pair<double,double>(x,y);
  
};

double cosTheta(const ROOT::Math::XYZVector & AC,
		const ROOT::Math::XYZVector & AB){
  
  return sqrt( AC.Dot(AB)*AC.Dot(AB) / (AC.Mag2()*AB.Mag2()) );
  
};


TPad* plot_ratio(TCanvas *canv, bool up){
  canv->SetFillColor      (0);
  canv->SetBorderMode     (0);
  canv->SetBorderSize     (10);
  // Set margins to reasonable defaults
  canv->SetLeftMargin     (0.18);
  canv->SetLeftMargin     (0.17);
  canv->SetRightMargin    (0.05);
  canv->SetTopMargin      (0.05);
  canv->SetBottomMargin   (0.18);
  // Setup a frame which makes sense
  canv->SetFrameFillStyle (0);
  canv->SetFrameLineStyle (0);
  canv->SetFrameBorderMode(0);
  canv->SetFrameBorderSize(10);
  canv->SetFrameFillStyle (0);
  canv->SetFrameLineStyle (0);
  canv->SetFrameBorderMode(0);
  canv->SetFrameBorderSize(10);      

  canv->cd();
  TPad *pad = 0;
  if (up){
    pad = new TPad("upper","pad",0, 0.26 ,1 ,1);
    pad->SetBottomMargin(0.05);
    pad->SetTopMargin(0.05);
    pad->Draw();
    pad->cd();
    return pad;
  }
  else {
    pad = new TPad("lower","pad",0, 0   ,1 ,0.26);  
    pad->SetTopMargin(0.05);
    pad->SetBottomMargin(0.24);
    pad->Draw();
    return pad;
  }

};

double getZpos(const unsigned iL){
  if (iL==0) return 322.103;//cm
  else if (iL==1) return 323.047;//cm
  else if (iL==2) return 325.073;//cm
  else if (iL==3) return 326.017;//cm
  else if (iL==4) return 328.043;//cm
  else if (iL==5) return 328.987;//cm
  else if (iL==6) return 331.013;//cm
  else if (iL==7) return 331.957;//cm
  else if (iL==8) return 333.983;//cm
  else if (iL==9) return 334.927;//cm
  else if (iL==10) return 336.953;//cm
  else if (iL==11) return 337.897;//cm
  else if (iL==12) return 339.923;//cm
  else if (iL==13) return 340.867;//cm
  else if (iL==14) return 342.893;//cm
  else if (iL==15) return 343.837;//cm
  else if (iL==16) return 345.863;//cm
  else if (iL==17) return 346.807;//cm
  else if (iL==18) return 348.833;//cm
  else if (iL==19) return 349.777;//cm
  else if (iL==20) return 351.803;//cm
  else if (iL==21) return 352.747;//cm
  else if (iL==22) return 354.773;//cm
  else if (iL==23) return 355.717;//cm
  else if (iL==24) return 357.743;//cm
  else if (iL==25) return 358.687;//cm
  else if (iL==26) return 360.713;//cm
  else if (iL==27) return 361.657;//cm
  else if (iL==28) return 367.699;//cm
  else if (iL==29) return 373.149;//cm
  else if (iL==30) return 378.599;//cm
  else if (iL==31) return 384.049;//cm
  else if (iL==32) return 389.499;//cm
  else if (iL==33) return 394.949;//cm
  else if (iL==34) return 400.399;//cm
  else if (iL==35) return 405.849;//cm
  else if (iL==36) return 411.299;//cm
  else if (iL==37) return 416.749;//cm
  else if (iL==38) return 422.199;//cm
  else if (iL==39) return 427.649;//cm
  else if (iL==40) return 436.199;//cm
  else if (iL==41) return 444.749;//cm
  else if (iL==42) return 453.299;//cm
  else if (iL==43) return 461.849;//cm
  else if (iL==44) return 470.399;//cm
  else if (iL==45) return 478.949;//cm
  else if (iL==46) return 487.499;//cm
  else if (iL==47) return 496.049;//cm
  else if (iL==48) return 504.599;//cm
  else return 0; 
}

double getDEta(const double & eta, const double & z, const double & dr){

  double r0 = z*tan(2*atan(exp(-1.*eta)));
  double theta1 = atan((r0+dr)/z); 
  double theta2 = atan((r0-dr)/z); 

  return -log(tan(theta2/2))+log(tan(theta1/2));
  
};
