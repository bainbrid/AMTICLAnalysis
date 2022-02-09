#include <algorithm>
#include <iomanip>
#include <stdlib.h>

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

int plotEreso(){

  const std::string pteta = "pt5_eta17";
  //std::string filePath = "D49_EM_trkFirst/";
  std::string filePath = "D49_Dummy_new/";

  std::string lPlotDir = filePath+"/"+pteta+"/LCPlots/";
  if (system(("mkdir -p "+lPlotDir).c_str())) return 1;

   
  TFile * fin = TFile::Open((filePath+"/CloseByPhotons/step3_"+pteta+"_FlatTracksters.root").c_str());
  if (!fin) return 1;
  else std::cout << " File " << fin->GetName() << " opened." << std::endl;
  fin->cd("pid");
  TTree *tree = (TTree*)gDirectory->Get("tree");
  if (!tree) return 1;
  SetTdrStyle();
  gStyle->SetOptStat(0);

  std::string lCutBase = "trackster>=0";
  unsigned nPhotons = tree->GetEntries(lCutBase.c_str());
  std::cout << " -- Number of photons: " << nPhotons << std::endl;
  
  
  const unsigned nV = 1;
  const unsigned nN = 4;//cut on rechits
  const unsigned nL = 28;//layers
  
  TCanvas *myc[nN];
  for (unsigned iN(0); iN<nN; ++iN){//loop on variables
    std::ostringstream label;
    label << "myc_" << iN;
    myc[iN] = new TCanvas(label.str().c_str(),label.str().c_str(),1);
    //myc[iN]->Divide(2,2);
  }

  std::string var[nV] = {"EtotOverECP"};
  //unsigned nBins[nV] = {20,50,50};
  //double binMin[nV] = {-0.5,0,0};
  //double binMax[nV] = {19.5,10,0.2};
  int color[8] = {1,2,3,4,6,7,8,9};



  TLatex lat;
  char buf[200];
  
  TH1F *hist[nV][nN];
  for (unsigned iV(0); iV<nV; ++iV){//loop on variables
    for (unsigned iN(0); iN<nN; ++iN){//loop on Nrechits
      std::ostringstream label;
      label << "h_" << var[iV] << "_" << iN+1 << "rechits";
      hist[iV][iN] = new TH1F(label.str().c_str(),";E^{sum}_{LC}/E_{CP} (GeV) ;per TS",100,0.5,1.5);
      hist[iV][iN]->Sumw2();
    }
  }

  for (unsigned iV(0); iV<nV; ++iV){//loop on variables
    for (unsigned iN(0); iN<nN; ++iN){//loop on Nrechits
      myc[iN]->cd();
      std::ostringstream lVar;

      for (unsigned iL(0); iL<nL; ++iL){//loop on layers

	if (iV==0) lVar << "lc_energy/cp_energy";
	lVar << ">>" << hist[iV][iN]->GetName();
	
	std::ostringstream lCut;
	lCut << "1./" << nPhotons << "*(" << lCutBase << " && lc_layer == " << iL+1 << " && lc_nrechits >= " << iN+1 << ")";

	tree->Draw(lVar.str().c_str(),lCut.str().c_str(),"colz");

	std::cout << lVar.str() << " " << var[iV] << " layer " << iL+1 << " nRecHits " << iN << " cut " << lCut.str() << " entries " << hist[iV][iL][iN]->GetEntries() << std::endl;
	sprintf(buf,"Layer %d, Nrechits = %d",iL+1,iN+1);
	lat.DrawLatexNDC(0.15,0.95,buf);
	
      }//loop on rechits
      
      /*      for (unsigned iN(0); iN<nN; ++iN){//loop on Nrechits
	hist[iV][iL][iN]->SetLineColor(color[iN]);
	hist[iV][iL][iN]->SetMarkerColor(color[iN]);
	hist[iV][iL][iN]->SetMarkerStyle(20+iN);
	hist[iV][iL][iN]->Draw(iN==0?"PE":"PEsame");
	sprintf(buf,"nRecHits #geq %d",iN+1);
	if (iV==0 && iL==0) leg->AddEntry(hist[iV][iL][iN],buf,"P");
      }
      leg->Draw("same");
      */

      hist[iR][iV]->Fit("gaus","","same");
      TF1 *fit = (TF1*)hist[iR][iV]->GetFunction("gaus");
      sprintf(buf,"#sigma/m=%.3f (%.3f/#sqrt{E})",fit->GetParameter(2)/fit->GetParameter(1),fit->GetParameter(2)/fit->GetParameter(1)*sqrt(E[iR]));
      lat.SetTextColor(color[iR]);
      lat.SetTextSize(0.03);
      lat.DrawLatexNDC(0.16,0.85-0.1*iR,buf);
      
      myc[iL]->Update();
      std::ostringstream label;
      label << lPlotDir << "/" << var[iV] << "_layer" ;
      if (iL<9) label << "0";
      label << iL+1 << ".pdf"; 
      myc[iL]->Print(label.str().c_str());
      
    }//loop on layers
  }//loop on variables
  

  return 0;

  
}
