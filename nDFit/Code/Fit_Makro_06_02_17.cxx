
#include <string>
#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TH1F.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TCanvas.h"
#include "TImage.h"

using namespace std;




using namespace std;

////////////////////////////////////////////////////////////////////////
//Fit Functions
////////////////////////////////////////////////////////////////////////


//Gauss
Double_t signal(Double_t *x, Double_t *par) {
  return par[0]*exp(-0.5*pow((x[0]-par[1])/par[2],2));
}



//Landau
Double_t landau(Double_t *x,  Double_t *par) { // 3 par
  return par[0]*(TMath::Landau(x[0], par[1], par[2])/par[2]);
}



//Landau_n
Double_t landau_n(Double_t *x,  Double_t *par) { // 3 par
  return par[0]*(TMath::Landau(2*par[1]-x[0], par[1], par[2])/par[2]);
}


//fit top masse

Double_t fit_mtop(Double_t *x,  Double_t *par) { 
return signal(x, par) + landau(x, &par[3]) + landau_n(x, &par[6]);
}

//fit W masse
Double_t fit_mw(Double_t *x,  Double_t *par) { 
return signal(x,par) + signal(x,&par[3]);
}
//fit Rbq
Double_t fit_rbq(Double_t *x,  Double_t *par) { 
return signal(x,par) + signal(x,&par[3]) + landau(x,&par[6]);
}















void run(){
	
	
	TFile *file1 = new TFile("/home/iwsatlas1/schulte/Test_Merged_files/Merged_ttbar.root");
	
	TH1F *h1 = (TH1F*)file1->Get("hist_klf_mtop_window");

	
	//Top 
	
	TF1 *ftop = new TF1("ftop",fit_mtop,125,210,9);
	
	ftop->SetParameters(10, 180, 0.5, 700, 145, 8, 500, 130 ,9);
	ftop->SetParLimits(3,0.0,1000000);
	ftop->SetParLimits(6,0.0,1000000);
    ftop->SetParLimits(7,100.0,1000000);
	  
		//c1->cd(1);

	TCanvas *c1 = new TCanvas();

	
	h1->Fit("ftop","RI","",130,210);
	h1->Fit("ftop","RI","",130,210);
	h1->Fit("ftop","RI","",130,210);
	h1->Fit("ftop","RI","",130,210);
	h1->Fit("ftop","RI","",130,210);
	h1->Draw();
	
	double p[9];
	
	for(int i = 0; i < 9; i++){
	 
		p[i] = ftop->GetParameter(i);
		};
		
		
	
	TFile *file2 = new TFile("/home/iwsatlas1/schulte/Test_Merged_files/Merged_ttbar_170.root");
	
	TH1F *h2 = (TH1F*)file2->Get("hist_klf_mtop_window");
	
	
	ftop->SetParameters(p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8]);
	ftop->SetParLimits(3,0.0,1000000);
	ftop->SetParLimits(6,0.0,1000000);
    ftop->SetParLimits(7,100.0,1000000);
    
  
	TCanvas *c2 = new TCanvas();

	h2->Fit("ftop","RI","",130,210);
	h2->Draw();
	
    
	  
	
	}
