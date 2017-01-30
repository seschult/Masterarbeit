#include "ZROOT.h"


using namespace std;


////////////////////////////////////////////////////////////////////////
//Fit functions
////////////////////////////////////////////////////////////////////////
//Landau for fit of mtop
Double_t landau(Double_t *x,  Double_t *par) { // 3 par
  return par[0]*(TMath::Landau(x[0], par[1], par[2])/par[2]);
}
//Gauss
Double_t signal(Double_t *x, Double_t *par) {
  return par[0]*exp(-0.5*pow((x[0]-par[1])/par[2],2));
}
//fit Rbq and MW
Double_t fit(Double_t *x,  Double_t *par) { 
return signal(x,par) + signal(x,&par[3]);
}

//fit mtop
Double_t fit2(Double_t *x,  Double_t *par) { 
return signal(x,par) + landau(x,&par[3]);
}


////////////////////////////////////////////////////////////////////////
//BlaBlaBla




////////////////////////////////////////////////////////////////////////
//Class def
////////////////////////////////////////////////////////////////////////
class Background_fit{


public:	

Background_fit();
void top_Background();
void W_Background();
void rbq_Background();

	//MTOP
	TF1 *fmtop_Back;
	TF1 *fmtop_Landau;
	TF1 *fmtop_Gauss;
	
	TFile *file11;
	TFile *file12;
	TFile *file13;
	TFile *file14;
	TFile *file15;
	TH1F *h11; 
	TH1F *h12; 
	TH1F *h13; 
	TH1F *h14; 
	TH1F *h15; 
	
	
	//MW
	TF1 *fmW_Back;
	TF1 *fmW_Gauss1;
	TF1 *fmW_Gauss2;

	TH1F *h21; 
	TH1F *h22; 
	TH1F *h23; 
	TH1F *h24; 
	TH1F *h25; 
	
	
	//RBQ
	TF1 *frbq_Back;
	TF1 *frbq_Gauss1;
	TF1 *frbq_Gauss2;

	TH1F *h31; 
	TH1F *h32; 
	TH1F *h33; 
	TH1F *h34; 
	TH1F *h35; 
	
	
	TCanvas *c1;
	TCanvas *c2;
	TCanvas *c3;
	
	
	
	};
////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////
//Cunstructor
////////////////////////////////////////////////////////////////////////
Background_fit :: Background_fit(){
	
	
	//Mtop
	fmtop_Back = new TF1("fmtop_Back",fit2,125,220,6);
	fmtop_Back->SetParameters(10, 170, 0.5 ,400, 175, 300);  
	fmtop_Back->SetParLimits(0,0.0,1000000);
	fmtop_Back->SetParLimits(1,0.0,1000000);
	fmtop_Back->SetParLimits(3,0.0,1000000);
	
	fmtop_Gauss = new TF1("fmtop_Gauss",signal,125,220,3);
	fmtop_Landau = new TF1("fmtop_Landau",landau,125,220,3);
	
	 
	file11 = new TFile("/home/iwsatlas1/schulte/Merged_singleTop.root");
	file12 = new TFile("/home/iwsatlas1/schulte/Merged_ttV.root");
	file13 = new TFile("/home/iwsatlas1/schulte/Merged_ttbar.root");
	file14 = new TFile("/home/iwsatlas1/schulte/Merged_wjets.root");
	file15 = new TFile("/home/iwsatlas1/schulte/Merged_zjets.root");
  
	h11 = (TH1F*)file11->Get("hist_klf_mtop_window");

	h12 = (TH1F*)file12->Get("hist_klf_mtop_window");
	h13 = (TH1F*)file13->Get("hist_klf_mtop_window");
	h14 = (TH1F*)file14->Get("hist_klf_mtop_window");
	h15 = (TH1F*)file15->Get("hist_klf_mtop_window");
	
	
	//MW
	fmW_Back = new TF1("fmW_Back",fit,56,110,6);
	//fmW_Back->SetParameters(400, 190, 0.2);  
	fmW_Back->SetParameters(10, 80, 1, 10, 82, 10); 
	//fmW_Back->SetParLimits(0,0.0,1000000);
	 
	fmW_Gauss1 = new TF1("fmwGauss1",signal,56,110,3);
	fmW_Gauss2 = new TF1("fmwGauss2",signal,56,110,3);
;
	 
	h21 = (TH1F*)file11->Get("hist_klf_window_Whad_m");
	h22 = (TH1F*)file12->Get("hist_klf_window_Whad_m");
	h23 = (TH1F*)file13->Get("hist_klf_window_Whad_m");
	h24 = (TH1F*)file14->Get("hist_klf_window_Whad_m");
	h25 = (TH1F*)file15->Get("hist_klf_window_Whad_m");
	
	//Rbq
	frbq_Back = new TF1("frbq_Back",fit,0.3,3,6);
	frbq_Back->SetParameters(10, 0.9, 1, 10, 1.3, 10); 
	//frbq_Back->SetParameters(10, 0.9, 1, 10, 1.3, 10); 
	//frbq_Back->SetParameters(400, 190, 0.2);  
	frbq_Back->SetParLimits(0,0.0,1000000);
	frbq_Back->SetParLimits(3,0.0,1000000);
	
	frbq_Gauss1 = new TF1("frbqGauss1",signal,0.3,3,3);
	frbq_Gauss2 = new TF1("frbqGauss2",signal,0.3,3,3);
	 
	h31 = (TH1F*)file11->Get("hist_klf_window_Rbq_reco");
	h32 = (TH1F*)file12->Get("hist_klf_window_Rbq_reco");
	h33 = (TH1F*)file13->Get("hist_klf_window_Rbq_reco");
	h34 = (TH1F*)file14->Get("hist_klf_window_Rbq_reco");
	h35 = (TH1F*)file15->Get("hist_klf_window_Rbq_reco");
	
	
	
////////////////////////////////////////////////////////////////////////
//Add Histogramms
////////////////////////////////////////////////////////////////////////
	//mtop
	h11->Add(h12);
	//h11 = h12;
	h11->Add(h13);
	h11->Add(h14);
	h11->Add(h15);
	
	//mw
	h21->Add(h22);
	//h21 = h22;
	h21->Add(h23);
	h21->Add(h24);
	h21->Add(h25);
	
	//rbq
	h31->Add(h32);
	//h31 = h32;
	h31->Add(h33);
	h31->Add(h34);
	h31->Add(h35);
	
}
	
////////////////////////////////////////////////////////////////////////

void Background_fit :: top_Background(){
		
////////////////////////////////////////////////////////////////////////	
//Fit Histogramms mtop Background
////////////////////////////////////////////////////////////////////////
	double with = 600;
	double hight = 700;
	c1 = new TCanvas("c1","c1", with, hight);
	gStyle->SetOptStat(0);
	
  
	h11->Fit("fmtop_Back","RI","",130,210);
	h11->Draw();
	
	//double Chi_2 = f0->GetChisquare();
	
	
////////////////////////////////////////////////////////////////////////	
//Functions	
////////////////////////////////////////////////////////////////////////

	double p[6];
 
	for(int i = 0; i < 6; i++){
	 
		p[i] = fmtop_Back->GetParameter(i);
		};
 
	//Draw functions
	fmtop_Gauss->SetParameters(p[0],p[1],p[2]);
	fmtop_Gauss->SetLineColor(3);
	fmtop_Gauss->Draw("SAME");
 
	fmtop_Landau->SetParameters(p[3],p[4],p[5]);
	fmtop_Landau->SetLineColor(4);
	fmtop_Landau->Draw("SAME");	
	


 
////////////////////////////////////////////////////////////////////////
//chi2 , NDF and Prob
////////////////////////////////////////////////////////////////////////


	double chi2_1 = fmtop_Back->GetChisquare();
	double NDF_1 = fmtop_Back->GetNDF();
	double Prob = fmtop_Back->GetProb();
	
	std::stringstream oss_Sep1;
	oss_Sep1   << setprecision(3) << chi2_1;
	
	TLatex l10;
	l10.SetTextAlign(9);
	l10.SetTextSize(0.038);
	l10.SetNDC();
	l10.DrawLatex(0.1, 0.92, ("chi^2: " + oss_Sep1.str()).c_str());
	
	
	std::stringstream oss_NDF1;
	oss_NDF1 << setprecision(3) << NDF_1;
	
	TLatex l11;
	l11.SetTextAlign(9);
	l11.SetTextSize(0.038);
	l11.SetNDC();
	l11.DrawLatex(0.1, 0.95, ("NDF: " + oss_NDF1.str()).c_str());
	
	
	std::stringstream oss_Prob;
	oss_Prob << setprecision(3) << Prob;
	
	TLatex l12;
	l12.SetTextAlign(9);
	l12.SetTextSize(0.038);
	l12.SetNDC();
	l12.DrawLatex(0.1, 0.88, ("Prob: " + oss_Prob.str()).c_str());







////////////////////////////////////////////////////////////////////////	
	
	
	
	}


void Background_fit :: W_Background(){
		
////////////////////////////////////////////////////////////////////////	
//Fit Histogramms mtop Background
////////////////////////////////////////////////////////////////////////
	double with = 600;
	double hight = 700;
	c2 = new TCanvas("c2","c2", with, hight);
	gStyle->SetOptStat(0);
	
  
	h21->Fit("fmW_Back","RI","",56,110);
	h21->Draw();
	
	//double Chi_2 = f0->GetChisquare();



////////////////////////////////////////////////////////////////////////
//Show Gauss1 and Gauss2
////////////////////////////////////////////////////////////////////////

	
	double p[6];
 
	for(int i = 0; i < 6; i++){
	 
		p[i] = fmW_Back->GetParameter(i);
		};
 
	//Draw functions
	fmW_Gauss1->SetParameters(p[0],p[1],p[2]);
	fmW_Gauss1->SetLineColor(3);
	fmW_Gauss1->Draw("SAME");
 
	fmW_Gauss2->SetParameters(p[3],p[4],p[5]);
	fmW_Gauss2->SetLineColor(4);
	fmW_Gauss2->Draw("SAME");
 



 
////////////////////////////////////////////////////////////////////////
//chi2 , NDF and Prob
////////////////////////////////////////////////////////////////////////


	double chi2_1 = fmW_Back->GetChisquare();
	double NDF_1 = fmW_Back->GetNDF();
	double Prob = fmW_Back->GetProb();
	
	std::stringstream oss_Sep1;
	oss_Sep1   << setprecision(3) << chi2_1;
	
	TLatex l10;
	l10.SetTextAlign(9);
	l10.SetTextSize(0.038);
	l10.SetNDC();
	l10.DrawLatex(0.1, 0.92, ("chi^2: " + oss_Sep1.str()).c_str());
	
	
	std::stringstream oss_NDF1;
	oss_NDF1 << setprecision(3) << NDF_1;
	
	TLatex l11;
	l11.SetTextAlign(9);
	l11.SetTextSize(0.038);
	l11.SetNDC();
	l11.DrawLatex(0.1, 0.95, ("NDF: " + oss_NDF1.str()).c_str());
	
	
	std::stringstream oss_Prob;
	oss_Prob << setprecision(3) << Prob;
	
	TLatex l12;
	l12.SetTextAlign(9);
	l12.SetTextSize(0.038);
	l12.SetNDC();
	l12.DrawLatex(0.1, 0.88, ("Prob: " + oss_Prob.str()).c_str());








////////////////////////////////////////////////////////////////////////	
	
	
	
	}
	
	
void Background_fit :: rbq_Background(){
		
////////////////////////////////////////////////////////////////////////	
//Fit Histogramms mtop Background
////////////////////////////////////////////////////////////////////////
	double with = 600;
	double hight = 700;
	c3 = new TCanvas("c3","c3", with, hight);
	gStyle->SetOptStat(0);
	
  
	h31->Fit("frbq_Back","RI","",0.32,2.99);
	h31->Draw();
	
	//double Chi_2 = f0->GetChisquare();

////////////////////////////////////////////////////////////////////////
//Show Gauss1 and Gauss2
////////////////////////////////////////////////////////////////////////

	
	double p[6];
 
	for(int i = 0; i < 6; i++){
	 
		p[i] = fmW_Back->GetParameter(i);
		};
 
	//Draw functions
	frbq_Gauss1->SetParameters(p[0],p[1],p[2]);
	frbq_Gauss1->SetLineColor(3);
	frbq_Gauss1->Draw("SAME");
 
	frbq_Gauss2->SetParameters(p[3],p[4],p[5]);
	frbq_Gauss2->SetLineColor(4);
	frbq_Gauss2->Draw("SAME");
 
////////////////////////////////////////////////////////////////////////
//chi2 , NDF and Prob
////////////////////////////////////////////////////////////////////////


	double chi2_1 = frbq_Back->GetChisquare();
	double NDF_1 = frbq_Back->GetNDF();
	double Prob = frbq_Back->GetProb();
	
	std::stringstream oss_Sep1;
	oss_Sep1   << setprecision(3) << chi2_1;
	
	TLatex l10;
	l10.SetTextAlign(9);
	l10.SetTextSize(0.038);
	l10.SetNDC();
	l10.DrawLatex(0.1, 0.92, ("chi^2: " + oss_Sep1.str()).c_str());
	
	
	std::stringstream oss_NDF1;
	oss_NDF1 << setprecision(3) << NDF_1;
	
	TLatex l11;
	l11.SetTextAlign(9);
	l11.SetTextSize(0.038);
	l11.SetNDC();
	l11.DrawLatex(0.1, 0.95, ("NDF: " + oss_NDF1.str()).c_str());
	
	
	std::stringstream oss_Prob;
	oss_Prob << setprecision(3) << Prob;
	
	TLatex l12;
	l12.SetTextAlign(9);
	l12.SetTextSize(0.038);
	l12.SetNDC();
	l12.DrawLatex(0.1, 0.88, ("Prob: " + oss_Prob.str()).c_str());


////////////////////////////////////////////////////////////////////////	
	
	
	
	}	

////////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[]){
	TApplication* A =new TApplication("ssss",&argc, argv);

	Background_fit p;

	p.top_Background();
	p.W_Background();
	p.rbq_Background();

	
	A->Run();
	return 0;
	}






