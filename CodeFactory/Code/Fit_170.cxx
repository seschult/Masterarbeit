#include <string>
#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TH1F.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"


using namespace std;

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




void run(const char *File, const char *Save){
	
	
	//Top 
	TF1 *f00 = new TF1("f5",signal,125,210,3);
	TF1 *f01 = new TF1("f1",landau,125,210,3);
	TF1 *f02 = new TF1("f2",landau_n,125,210,3);
    TF1  *ftop = new TF1("ftop",fit_mtop,125,210,9);
	  
	 
	//MW 
    TF1  *f13 = new TF1("f3",signal,56,110,3);
    TF1 *f14 = new TF1("f4",signal,56,110,3);
	TF1 *fmw = new TF1("fmw",fit_mw,56,110,6);
	
	
	//RBQ 
	 TF1  *f25 = new TF1("f5",signal,0.3,3,3);
	 TF1  *f26 = new TF1("f6",signal,0.3,3,3);
	 TF1  *f27 = new TF1("f7",landau,0.3,3,3);
	 TF1  *frbq = new TF1("frbq",fit_rbq,0.3,3,9);
	
	
	
	////////////////////////////////////////////////////////////////////
	//Fit Top Masse
	/////////////////////////////////////////////////////////////////7//
	
	
	TFile *file1 = new TFile("~/Merged_ttbar.root","r");
  
	TH1F *h1 = (TH1F*)file1->Get("hist_klf_mtop_window");

	ftop->SetParameters(10, 170, 0.5, 400, 175, 8, 400, 120 ,9);
	ftop->SetParLimits(3,0.0,1000000);
	ftop->SetParLimits(6,0.0,1000000);
    ftop->SetParLimits(7,100.0,1000000);
	
		//c1->cd(1);
	double with = 600;
	double hight = 700;
	c1 = new TCanvas("c1","c1", with, hight);
	

	//Fit of the histogram   
	h1->Fit("ftop","RI","",130,210);
	h1->Draw();
	
	

	//Get fit parameters 
	double p1[9];
 
	for(int i = 0; i < 9; i++){
	 
		p1[i] = ftop->GetParameter(i);
		};
 
 
 
	TFile *file2 = new TFile(File,"r");
  
	TH1F *h2 = (TH1F*)file2->Get("hist_klf_mtop_window");

	ftop->SetParameters(p1[0],p1[1],p1[2],p1[3],p1[4],p1[5],p1[6],p1[7],p1[8]);
	ftop->SetParLimits(3,0.0,1000000);
	ftop->SetParLimits(6,0.0,1000000);
    ftop->SetParLimits(7,100.0,1000000);
	
		//c1->cd(1);
	double with2 = 600;
	double hight2 = 700;
	c2 = new TCanvas("c2","c2", with2, hight2);
	

	//Fit of the histogram   
	h2->Fit("ftop","RI","",130,210);
	h2->Draw();
	int npar_mtop = 9;
	double p2[9];
 
	for(int i = 0; i < 9; i++){
	 
		p2[i] = ftop->GetParameter(i);
		};
	
	
	
	
	
	Char_t Label[300];
 
	
	TFile *fend = new TFile(Save,"RECREATE");
	fend ->Write();
	TTree *t1 = new TTree("t1","Data");
	 /*

	//t1->Branch("Label",Label,"label/C");
	
	t1->Branch("fit_mtop_npar",npar_mtop,"fit_mtop_npar/I");
	t1->Branch("fit_mtoppar",p2,"fit_mtoppar[fit_mtop_npar]/D");
	//t1->Branch("fit_mtoperr",&p->err_mtop,"fit_mtoperr[fit_mtop_npar]/D");
	
	t1->Branch("fit_mW_npar",&p->npar_mw,"fit_mW_npar/I");
	t1->Branch("fit_mWpar",&p->par_mw,"fit_mWpar[fit_mW_npar]/D");
	t1->Branch("fit_mWerr",&p->err_mw,"fit_mWerr[fit_mW_npar]/D");
	
	t1->Branch("fit_rbq_npar",&p->npar_rbq,"fit_rbq_npar/I");
	t1->Branch("fit_rbq",&p->par_rbq,"fit_rbqpar[fit_rbq_npar]/D");
	t1->Branch("fit_rbqerr",&p->err_rbq,"fit_rbqerr[fit_rbq_npar]/D");
	
	//Nation[0]='';
	//sprintf(&(Nation[0]),"%s",A->Argv(1));
	t1->Fill();   
	fend->Write();
	fend->Close();
 
   */
 
 
 }
