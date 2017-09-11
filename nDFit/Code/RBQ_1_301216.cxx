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
#include "TLegend.h"
#include "TLatex.h"
#include "TAxis.h"
#include "TTree.h"


using namespace std;




/////////////////////////////////////////////////////////////////////////////////////////
//Definition of Functions
////////////////////////////////////////////////////////////////////////////////////////

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


//fit w_mass

Double_t fit_r(Double_t *x,  Double_t *par) { 
return signal(x,par) + signal(x,&par[3]) + landau(x,&par[6]);
}



void draw2(const char *File) {


 //Definition of functions

     TF1 *f0 = new TF1("f0",signal,0.3,3,3);
	 TF1 *f1 = new TF1("f1",signal,0.3,3,3);
	 TF1 *f2 = new TF1("f2",landau,0.3,3,3);
	 TF1 *fit = new TF1("fit",fit_r,0.3,3,9);
	 
	fit->SetParameters(10, 1, 2, 10, 2, 1, 10, 2 ,1);
	//fit->SetParameters(10, 1, 1, 1, 2, 1, 10, 2 ,15);
	//fit->SetParLimits(0,0.0,1000000);
  //  fit->SetParLimits(3,0.0,1000000);
   // fit->SetParLimits(7,0.0,1000000);
    
 // Open file and histogram   
	TFile *file = new TFile(File);
  
	TH1F *h1 = (TH1F*)file->Get("hist_klf_window_Rbq_reco");
	
	
 //Fit of the histogram 
    Double_t w = 600;
    Double_t h = 700;
	TCanvas *c1 = new TCanvas("c", "c", w, h);
	gStyle->SetOptStat(0);
	//c->SetWindowSize(w + (w - c->GetWw()), h + (h - c->GetWh()));
	TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
	pad1->SetBottomMargin(2); // Upper and lower plot are joined
    pad1->SetGridx();         // Vertical grid
    pad1->Draw();             // Draw the upper pad: pad1
    pad1->cd();               // pad1 becomes the current pad
	
	
	fit->SetLineColor(kRed);
	h1->Fit("fit","I","",0.3,3);
	//h1->Fit("fit","I","",0.301,2.99);
	h1->Draw();

	
	double chi2 = fit->GetChisquare();
	double NDF = fit->GetNDF();
	
	std::stringstream oss_Sep;
	oss_Sep   << setprecision(3) << chi2;

	TLatex l0;
	l0.SetTextAlign(9);
	l0.SetTextSize(0.038);
	l0.SetNDC();
	l0.DrawLatex(0.1, 0.92, ("chi^2: " + oss_Sep.str()).c_str());

	std::stringstream oss_NDF;
	oss_NDF << setprecision(3) << NDF;

	TLatex l1;
	l1.SetTextAlign(9);
	l1.SetTextSize(0.038);
	l1.SetNDC();
	l1.DrawLatex(0.1, 0.95, ("NDF: " + oss_NDF.str()).c_str());
 
 //Get fit parameters 
 double p[9];
 
 for(int i = 0; i < 9; i++){
	 
	 p[i] = fit->GetParameter(i);
 }
 
  
  
  
 
 
 //Draw functions
	
  
  
  f0->SetParameters(p[0],p[1],p[2]);
  f0->SetLineColor(3);
  f0->Draw("SAME");
 
  f1->SetParameters(p[3],p[4],p[5]);
  f1->SetLineColor(4); 
  f1->Draw("SAME");
  
 
  f2->SetParameters(p[6],p[7],p[8]);
  f2->SetLineColor(5);
  f2->Draw("SAME");
  
  
   TLegend *leg = new TLegend(0.6,0.7,0.78,0.9);
   //leg->SetHeader("The Legend Title","C"); // option "C" allows to center the header
   leg->AddEntry(h1,"Data","lep");
   leg->AddEntry(fit,"Fit","l");
   leg->AddEntry("f0","Gauss","l");
   leg->AddEntry("f1","Gauss","l");
   leg->AddEntry("f2","Landau","l");
   leg->Draw();
   
   
   c1->cd();          // Go back to the main canvas before defining pad2
   TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
   pad2->SetTopMargin(0);
   pad2->SetBottomMargin(0.2);
   pad2->SetGridx(); // vertical grid
   pad2->Draw();
   pad2->cd();   
       // pad2 becomes the current pad
    bin_tot = h1 -> GetSize();   
    double start = h1->GetXaxis()->GetBinCenter(0);
    double end = h1->GetXaxis()->GetBinCenter(bin_tot);
       
	//TH1F* h2 = new TH1F("h2", "test",);
	bin_tot = h1 -> GetSize();
	TH1F *h2 = new TH1F("h2", "dif", bin_tot, start, end);
	//bin_tot = h1 -> GetSize();
	//cout<< "TOTAL:" bin_tot << endl;
	
	
	 TF1 *func = h1->GetFunction("fit");
	 
 	for (int bin=1; bin<=bin_tot;bin++) {
    double x = h1->GetXaxis()->GetBinCenter(bin);
    double fval = func->Eval(x);
    double dif = h1->GetBinContent(bin);
    
    if(dif != 0){
    h2->SetBinContent(bin,h1->GetBinContent(bin)-fval);
	}else continue;
    //std::cout<< "func:"<<fval;
    //std::cout << "Bin:"<<dif;
    
    
};
	h2->Draw();
	
	
	
	/*
	for(int i = 0; i < bin_tot; i++ ){
		double bin = h1->GetBin(i);
		
		TAxis *xaxis = h1->GetXaxis(); 
        double binCenter = xaxis->GetBinCenter(bin);
		
		
		double func = fit->Eval(binCenter);
		
		double diff = func - bin;
		
		h2->Fill(binCenter);
	
		//printf ("Characters: %d \n", 'func', 65);
		
		};*/
	
		


/*
   // Define the ratio plot
   TH1F *h2 = (TH1F*)h1->Clone("h3");
   h2->SetLineColor(kBlack);
   h2->SetMinimum(0.8);  // Define Y ..
   h2->SetMaximum(1.35); // .. range
   //h2->Sumw2();
   h2->SetStats(0);      // No statistics on lower plot
   h2->Divide(fit);
   h2->SetMarkerStyle(21);
   h2->Draw("ep");       // Draw the ratio plot


*/
	
 TFile *fend = new TFile("Fit_Data.root","RECREATE");
 
 for(int i=0; i<9; i++){
  p[i]->Write();
};
   //TTree *t1 = new TTree("t1","Data");
   
	//t1->Fill();   
	//fend->Write();
	//fend->Close();

	

}
 


	


void draw(const char *File) {

{
 //Definition of functions

     TF1 *f0 = new TF1("f0",signal,0.3,3,3);
	 TF1 *f1 = new TF1("f1",signal,0.3,3,3);
	 TF1 *f2 = new TF1("f2",landau_n,0.3,3,3);
	 TF1 *fit = new TF1("fit",fit_r,0.3,3,9);
	 
	fit->SetParameters(10, 1, 2, 10, 0, 1, 10, 2 ,3);
	//fit->SetParameters(10, 1, 1, 1, 2, 1, 10, 2 ,15);
	//f0->SetParLimits(0,0.0,1000000);
    f0->SetParLimits(3,0.0,1000000);
    //f0->SetParLimits(6,0.0,1000000);
    
 // Open file and histogram   
	TFile *file = new TFile(File);
  
	TH1F *h1 = (TH1F*)file->Get("hist_klf_window_Rbq_reco");
	
	
 //Fit of the histogram   
	TCanvas *c1 = new TCanvas;
	h1->Fit("fit","","",0.3,3);
	h1->Draw();
 
 //Get fit parameters 
 double p[9];
 
 for(int i = 0; i < 9; i++){
	 
	 p[i] = fit->GetParameter(i);
 }
 
  
  
  
 
 
 //Draw functions

  
  
  f0->SetParameters(p[0],p[1],p[2]);
  f0->SetLineColor(3);
  f0->Draw("SAME");
 
  f1->SetParameters(p[3],p[4],p[5]);
  f1->SetLineColor(4);
  f1->Draw("SAME");
  
 
  f2->SetParameters(p[6],p[7],p[8]);
  f2->SetLineColor(5);
  f2->Draw("SAME");
  
  
  
   TFile *fend = new TFile("Fit_Data.root","RECREATE");
   TTree *t1 = new TTree("t1","Data");
   
	t1->Fill();   
	fend->Write();
	fend->Close();
	
   

	
	


	

}
 
  
	

  

  
}

