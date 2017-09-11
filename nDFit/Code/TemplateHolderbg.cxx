#include "TTree.h"
#include "TF1.h"
#include "TH1F.h"
#include "TLine.h"
#include "TROOT.h"
#include "TPad.h"
#include "TMath.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TMinuit.h"
#include "TMinuitMinimizer.h"
#include "TAttMarker.h"
#include "Math/Functor.h"
#include "Math/Functor.h"

#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <math.h>
#include <stdio.h>
#include <string>
#include <sys/stat.h>
#include <unistd.h>
#include <iomanip>
#include <omp.h>

#include "Fill3D.h"
#include "FCN.h"
#include "TemplateHolder.h"

using namespace std;





Double_t TemplateHolder ::gauss(Double_t *x, Double_t *par) {return par[0]*exp(-0.5*pow((x[0]-par[1])/par[2],2));}
//Landau
Double_t TemplateHolder ::landau(Double_t *x,  Double_t *par) { return par[0]*(TMath::Landau(x[0], par[1], par[2])/par[2]);}


//topbg mass fit function
Double_t TemplateHolder ::mtop_functionbg(Double_t *x,  Double_t *par) { 
return gauss(x, par) + landau(x, &par[3]);}
//fit Wbg masse and Rbq
Double_t TemplateHolder ::fit_mwbg(Double_t *x,  Double_t *par) { 
return gauss(x,par) + gauss(x,&par[3]);}









TH1F* FitRatio(TF1* f, TH1F* h){
  TH1F* hratio=(TH1F*) h->Clone("hratio");
  for (int i = 0; i < hratio->GetNbinsX(); ++i)
    {
      double x=hratio->GetBinCenter(i+1);
      double valf=f->Eval(x);
      double valh=hratio->GetBinContent(i+1);
      double errh=hratio->GetBinError(i+1);
      double diff_sigma=(valh - valf)/errh;
      hratio->SetBinContent(i+1,diff_sigma);
      // double ratio=valf/valh;
      // double ratioerr=ratio*(errh/valh);
      // hratio->SetBinContent(i+1,ratio);
      // hratio->SetBinError(i+1,ratioerr);
    }
  //hratio->SetLineColor(f->GetLineColor());
  //hratio->SetLineWidth(1);
  return hratio;
}






TemplateHolderbg :: TemplateHolderbg (const char *File, double t, double j, double bg){
	
	fbg = bg;
	fbJSF=b;
	fJSF=j;
	
	
    fgaustop = new TF1("ftopgaus",this,&TemplateHolder::gauss,130,210,3);
    flandtop = new TF1("f1",this,&TemplateHolder::landau,130,210,3);
   
    ftopbg = new TF1("ftop",this,&TemplateHolder::mtop_functionbg,130,210,9);
    fgaussw = new TF1("f3",this,&TemplateHolder::gauss,56,110,3);
    fgaussw = new TF1("f4",this,&TemplateHolder::gauss,56,110,3);
	fw = new TF1("fmw",this,&TemplateHolder::fit_mwbg,56,110,6); 
	
	// Open file and histogram 
    puts(File);  
    file = new TFile(File,"r");
    htopbg = (TH1F*)file->Get("hist_klf_mtop_window");
    hwbg = (TH1F*)file->Get("hist_klf_window_Whad_m");
    hrbqbg = (TH1F*)file->Get("hist_klf_window_Rbq_reco");
    
    htopt->Add(htopbg);
    hwt->Add(hwbg);
    hrbqt->Add(rbqbg);
}


/*

void TemplateHolder::plot(std::string out){
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  double tsize = gStyle->GetLabelSize("X");

  TCanvas* c1 = new TCanvas("c1","c1", 800, 10, 800, 600);
  c1->cd();
  
  double padsize=0.3;
  TPad *pad1 = new TPad("pad1","pad1",0.,padsize,1.,1.);
  pad1->SetBottomMargin(0.03); 
  pad1->SetGridx();         // Vertical grid
  pad1->Draw();             // Draw the upper pad: pad1
  TPad *pad2 = new TPad("pad2","pad2",0.,0.,1.,padsize);
  pad2->SetTopMargin(0.045);
  pad2->SetBottomMargin(0.45);
  pad2->SetGridx(); // vertical grid
  pad2->SetGridy();
  pad2->Draw();
  
  pad1->cd();
  double fitmin = ftop->GetXmin();
  double fitmax = ftop->GetXmax(); 
  
  h1->GetXaxis()->SetRangeUser(fitmin,fitmax);
  h1->SetMaximum(1.2*h1->GetMaximum());
  //h_mlb->GetYaxis()->SetTitle(Form("Normalised events / %2.1f GeV", h_mlb->GetBinWidth(1))); 
  h1->GetYaxis()->SetTitle(Form("Events / %2.1f GeV", h1->GetBinWidth(1))); 
  h1->SetLabelOffset(0.3,"X");
  h1->SetLabelSize(0.06, "Y");
  h1->SetTitleSize(0.07, "Y");
  h1->SetTitleOffset(1.1,"Y");
  
  //Fit of the histogram   
  ftop->SetLineColor(2);
  ftop->SetLineWidth(3);
  h1->SetLineWidth(3);
	
  h1->Draw();
  ftop->Draw("Same");
 
	//Get fit parameters 
	////////////////////////////////////////////////////////////////////
	double p[9];
	
	double npar_mtop = 9;
 
	for(int i = 0; i < 9; i++){
	 
		p[i] = ftop->GetParameter(i);
		par_mtop[i] = p[i];
		err_mtop[i] = ftop->GetParError(i);
		};
 
	//Draw functions
	////////////////////////////////////////////////////////////////////
	f00->SetParameters(p[0],p[1],p[2]);
	f00->SetLineColor(6);
	f00->SetFillStyle(0); 
	f00->SetLineWidth(3);   
	f00->Draw("SAME"); 
	f01->SetParameters(p[3],p[4],p[5]);
	f01->SetLineColor(209);
	f01->SetFillStyle(0); 
	f01->SetLineWidth(3);  
	f01->Draw("SAME");  
	f02->SetParameters(p[6],p[7],p[8]);
	f02->SetLineColor(96);
	f02->SetFillStyle(0); 
	f02->SetLineWidth(3);  
	f02->Draw("SAME");
	
	//Get chi^2 
	////////////////////////////////////////////////////////////////////
	double chi2_0 = ftop->GetChisquare();
	double NDF_0 = ftop->GetNDF();
	double Prob = ftop->GetProb();
	double COM = chi2_0/NDF_0;
	
	std::stringstream oss_Sep;
	oss_Sep   << setprecision(3) << chi2_0;
	TLatex l00;
	l00.SetTextAlign(9);
	l00.SetTextSize(0.038);
	l00.SetLineWidth(2);
	l00.SetNDC();
	l00.DrawLatex(0.1358774,0.7424242, ("#chi^{2}: " + oss_Sep.str()).c_str());
	
	std::stringstream oss_NDF_0;
	oss_NDF_0 << setprecision(3) << NDF_0;
	TLatex l01;
	l01.SetTextAlign(9);
	l01.SetTextSize(0.038);
	l01.SetLineWidth(2);
	l01.SetNDC();
	l01.DrawLatex(.1358774,0.7811448, ("NDF: " + oss_NDF_0.str()).c_str());
	
	/*std::stringstream oss_Prob;
	oss_Prob << setprecision(3) << Prob;
	
	TLatex l03;
	l03.SetTextAlign(9);
	l03.SetTextSize(0.038);
	l03.SetLineWidth(2);
	l03.SetNDC();
	l03.DrawLatex(0.1369345,0.8299663, ("Prob: " + oss_Prob.str()).c_str());
	
	
	std::stringstream oss_COM;
	oss_COM   << setprecision(3) << COM;
	TLatex l04;
	l04.SetTextAlign(9);
	l04.SetTextSize(0.038);
	l04.SetLineWidth(2);
	l04.SetNDC();
	l04.DrawLatex(0.1358774,0.6952862, ("#chi^{2}/NDF: " + oss_COM.str()).c_str());
	
	TLatex lwork;
	lwork.SetTextAlign(9);
	lwork.SetTextSize(0.058);
	lwork.SetLineWidth(2);
	lwork.SetNDC();
	lwork.DrawLatex(0.1258774,0.8552862, ("ATLAS work-in-progress"));
	
	
	string Base = gSystem->BaseName(out.c_str());
	
	
    TLatex lfile;
	lfile.SetTextAlign(9);
	lfile.SetTextSize(0.058);
	lfile.SetLineWidth(2);
	lfile.SetNDC();
	lfile.DrawLatex(0.1258774,0.9552862, (Base.c_str()));
	////////////////////////////////////////////////////////////////////	
	//Legend
	////////////////////////////////////////////////////////////////////    
    TLegend *leg = new TLegend(0.7078059,0.6690754,0.8987342,0.8986966,NULL,"brNDC");
    leg->SetBorderSize(1);
    leg->SetLineColor(1);
    leg->SetLineStyle(1);
    leg->SetLineWidth(1);
    leg->SetFillColor(0);
    leg->SetFillStyle(1001);
    leg->AddEntry(h1,"Simulation","lep");
    leg->AddEntry(ftop,"Fit","l");
    leg->AddEntry(f00,"Gauss","l");
    leg->AddEntry(f01,"Landau","l");
    leg->AddEntry(f02,"Landau_n","l");
    leg->Draw();	
	////////////////////////////////////////////////////////////////////
   pad2->cd();   
       // pad2 becomes the current pad      
  TH1F* h_ratio = (TH1F*) FitRatio(ftop,h1);
  h_ratio->GetYaxis()->SetRangeUser(20,20);
  h_ratio->GetXaxis()->SetTitle("m_{top}^{reco}[GeV]");
  h_ratio->SetMaximum(100);    
  h_ratio->SetMinimum(-100);  
  h_ratio->SetLabelSize(0.14, "X");
  h_ratio->SetLabelOffset(0.05,"X");
  h_ratio->SetTitleOffset(1.2,"X");
  h_ratio->SetTitleSize(0.17, "X");
  h_ratio->SetLabelSize(0.10, "Y");
  h_ratio->SetNdivisions(305,"Y");  
  h_ratio->GetYaxis()->SetTitleOffset(1.0);
  h_ratio->GetYaxis()->SetTitleFont(63);
  h_ratio->GetYaxis()->SetTitleSize(20);
  h_ratio->GetYaxis()->SetTitle("(Sim.-fit)/#sigma");
  h_ratio->SetLineWidth(3);    
  h_ratio->SetMarkerStyle(20);	  
  h_ratio->SetMarkerSize(1);	  
  h_ratio->SetLineWidth(3); 
  h_ratio->Draw("P");

  TLine* tline_0 = new TLine(fitmin,0,fitmax,0);
  tline_0->SetLineColor(1);
  tline_0->SetLineStyle(1);
  tline_0->SetLineWidth(2);
  tline_0->Draw("SAME");


  c1->SaveAs((out+"mtop_reco.png").c_str());
 
//////////////////////////////////////////////////////////////////////// 
//PLOT M_W_reco  
////////////////////////////////////////////////////////////////////////  
  	


  TCanvas* cw = new TCanvas("cw","cw", 800, 10, 800, 600);
  cw->cd();
  

  TPad *padw1 = new TPad("padw1","padw1",0.,padsize,1.,1.);
  padw1->SetBottomMargin(0.03); 
  padw1->SetGridx();         // Vertical grid
  padw1->Draw();             // Draw the upper pad: pad1
  TPad *padw2 = new TPad("padw2","padw2",0.,0.,1.,padsize);
  padw2->SetTopMargin(0.045);
  padw2->SetBottomMargin(0.45);
  padw2->SetGridx(); // vertical grid
  padw2->SetGridy();
  padw2->Draw();
  
  padw1->cd();
  fitmin = fw->GetXmin();
  fitmax = fw->GetXmax(); 
  
  h2->GetXaxis()->SetRangeUser(fitmin,fitmax);
  h2->SetMaximum(1.2*h2->GetMaximum());
  //h_mlb->GetYaxis()->SetTitle(Form("Normalised events / %2.1f GeV", h_mlb->GetBinWidth(1))); 
  h2->GetYaxis()->SetTitle(Form("Events / %2.1f GeV", h2->GetBinWidth(1))); 
  h2->SetLabelOffset(0.3,"X");
  h2->SetLabelSize(0.06, "Y");
  h2->SetTitleSize(0.07, "Y");
  h2->SetTitleOffset(1.1,"Y");
  
  //Fit of the histogram   
  fw->SetLineColor(2);
  fw->SetLineWidth(3);
  h2->SetLineWidth(3);
	
  h2->Draw();
  fw->Draw("Same");
	//Get fit parameters 
	////////////////////////////////////////////////////////////////////
	double pw[6];
	
	double npar_mw = 6;
	
	for(int i = 0; i < 6; i++){
	 
		pw[i] = fw->GetParameter(i);
		par_mw[i] = pw[i];
		err_mw[i] = fw->GetParError(i);
		};
 
	//Draw functions
	////////////////////////////////////////////////////////////////////
	f13->SetParameters(pw[0],pw[1],pw[2]);
	f13->SetLineColor(6);
	f13->SetFillStyle(0); 
	f13->SetLineWidth(3);   
	f13->Draw("SAME"); 
	f14->SetParameters(pw[3],pw[4],pw[5]);
	f14->SetLineColor(209);
	f14->SetFillStyle(0); 
	f14->SetLineWidth(3);  
	f14->Draw("SAME");  
	
	
	//Get chi^2 
	////////////////////////////////////////////////////////////////////
	 chi2_0 = fw->GetChisquare();
	 NDF_0 = fw->GetNDF();
	 Prob = fw->GetProb();
	 COM = chi2_0/NDF_0;
	
	//std::stringstream oss_Sep;
	oss_Sep   << setprecision(3) << chi2_0;
	//TLatex l00;
	l00.SetTextAlign(9);
	l00.SetTextSize(0.038);
	l00.SetLineWidth(2);
	l00.SetNDC();
	l00.DrawLatex(0.1358774,0.7424242, ("#chi^{2}: " + oss_Sep.str()).c_str());
	
	//std::stringstream oss_NDF_0;
	oss_NDF_0 << setprecision(3) << NDF_0;
	//TLatex l01;
	l01.SetTextAlign(9);
	l01.SetTextSize(0.038);
	l01.SetLineWidth(2);
	l01.SetNDC();
	l01.DrawLatex(.1358774,0.7811448, ("NDF: " + oss_NDF_0.str()).c_str());
	
	/*std::stringstream oss_Prob;
	oss_Prob << setprecision(3) << Prob;
	
	TLatex l03;
	l03.SetTextAlign(9);
	l03.SetTextSize(0.038);
	l03.SetLineWidth(2);
	l03.SetNDC();
	l03.DrawLatex(0.1369345,0.8299663, ("Prob: " + oss_Prob.str()).c_str());
	
	
	//std::stringstream oss_COM;
	oss_COM   << setprecision(3) << COM;
	//TLatex l04;
	l04.SetTextAlign(9);
	l04.SetTextSize(0.038);
	l04.SetLineWidth(2);
	l04.SetNDC();
	l04.DrawLatex(0.1358774,0.6952862, ("#chi^{2}/NDF: " + oss_COM.str()).c_str());
	
	//TLatex lwork;
	lwork.SetTextAlign(9);
	lwork.SetTextSize(0.058);
	lwork.SetLineWidth(2);
	lwork.SetNDC();
	lwork.DrawLatex(0.1258774,0.8552862, ("ATLAS work-in-progress"));
	

    //TLatex lfile;
	lfile.SetTextAlign(9);
	lfile.SetTextSize(0.058);
	lfile.SetLineWidth(2);
	lfile.SetNDC();
	lfile.DrawLatex(0.1258774,0.9552862, (Base.c_str()));
	////////////////////////////////////////////////////////////////////	
	//Legend
	////////////////////////////////////////////////////////////////////    
    TLegend *leg1 = new TLegend(0.7078059,0.6690754,0.8987342,0.8986966,NULL,"brNDC");
    leg1->SetBorderSize(1);
    leg1->SetLineColor(1);
    leg1->SetLineStyle(1);
    leg1->SetLineWidth(1);
    leg1->SetFillColor(0);
    leg1->SetFillStyle(1001);
    leg1->AddEntry(h2,"Simulation","lep");
    leg1->AddEntry(fw,"Fit","l");
    leg1->AddEntry(f13,"Gauss","l");
    leg1->AddEntry(f14,"Gauss","l");
    leg1->Draw();	
	////////////////////////////////////////////////////////////////////
    padw2->cd();   
       // pad2 becomes the current pad      
  TH1F* h_ratiow = (TH1F*) FitRatio(fw,h2);
  h_ratiow->GetYaxis()->SetRangeUser(20,20);
  h_ratiow->GetXaxis()->SetTitle("m_{w}^{reco}[GeV]");
  h_ratiow->SetMaximum(100);    
  h_ratiow->SetMinimum(-100);  
  h_ratiow->SetLabelSize(0.14, "X");
  h_ratiow->SetLabelOffset(0.05,"X");
  h_ratiow->SetTitleOffset(1.2,"X");
  h_ratiow->SetTitleSize(0.17, "X");
  h_ratiow->SetLabelSize(0.10, "Y");
  h_ratiow->SetNdivisions(305,"Y");  
  h_ratiow->GetYaxis()->SetTitleOffset(1.0);
  h_ratiow->GetYaxis()->SetTitleFont(63);
  h_ratiow->GetYaxis()->SetTitleSize(20);
  h_ratiow->GetYaxis()->SetTitle("(Sim.-fit)/#sigma");
  h_ratiow->SetLineWidth(3);    
  h_ratiow->SetMarkerStyle(20);	  
  h_ratiow->SetMarkerSize(1);	  
  h_ratiow->SetLineWidth(3); 
  h_ratiow->Draw("P");

  tline_0 = new TLine(fitmin,0,fitmax,0);
  tline_0->SetLineColor(1);
  tline_0->SetLineStyle(1);
  tline_0->SetLineWidth(2);
  tline_0->Draw("SAME");

  

   cw->SaveAs((out+"mw_reco.png").c_str());
//////////////////////////////////////////////////////////////////////// 
//PLOT Rbq_reco  
////////////////////////////////////////////////////////////////////////  
  	   
  TCanvas* c3 = new TCanvas("c3","c3", 800, 10, 800, 600);
  c3->cd();
  
  
  TPad *padr1 = new TPad("padr1","padr1",0.,padsize,1.,1.);
  padr1->SetBottomMargin(0.03); 
  padr1->SetGridx();         // Vertical grid
  padr1->Draw();             // Draw the upper pad: pad1
  TPad *padr2 = new TPad("padr2","padr2",0.,0.,1.,padsize);
  padr2->SetTopMargin(0.045);
  padr2->SetBottomMargin(0.45);
  padr2->SetGridx(); // vertical grid
  padr2->SetGridy();
  padr2->Draw();
  
  padr1->cd();
  fitmin = frbq->GetXmin();
  fitmax = frbq->GetXmax(); 
  
  h3->GetXaxis()->SetRangeUser(fitmin,fitmax);
  h3->SetMaximum(1.2*h3->GetMaximum());
  //h_mlb->GetYaxis()->SetTitle(Form("Normalised events / %2.1f GeV", h_mlb->GetBinWidth(1))); 
  h3->GetYaxis()->SetTitle(Form("Events / %2.1f GeV", h1->GetBinWidth(1))); 
  h3->SetLabelOffset(0.3,"X");
  h3->SetLabelSize(0.06, "Y");
  h3->SetTitleSize(0.07, "Y");
  h3->SetTitleOffset(1.1,"Y");
  
  //Fit of the histogram   
  frbq->SetLineColor(2);
  frbq->SetLineWidth(3);
  h3->SetLineWidth(3);
	
  h3->Draw();
  frbq->Draw("Same");
  
	//Get fit parameters 
	////////////////////////////////////////////////////////////////////
	double pr[9];
	
	int npar_mrbq = 9;

 
	for(int i = 0; i < 9; i++){
	 
		pr[i] = frbq->GetParameter(i);
		par_mrbq[i] = pr[i];
		err_mrbq[i] = frbq->GetParError(i);
		};
 
	//Draw functions
	////////////////////////////////////////////////////////////////////
	f25->SetParameters(pr[0],pr[1],pr[2]);
	f25->SetLineColor(6);
	f25->SetFillStyle(0); 
	f25->SetLineWidth(3);   
	f25->Draw("SAME"); 
	f26->SetParameters(pr[3],pr[4],pr[5]);
	f26->SetLineColor(209);
	f26->SetFillStyle(0); 
	f26->SetLineWidth(3);  
	f26->Draw("SAME");  
	f27->SetParameters(pr[6],pr[7],pr[8]);
	f27->SetLineColor(96);
	f27->SetFillStyle(0); 
	f27->SetLineWidth(3);  
	f27->Draw("SAME");
	
	//Get chi^2 
	////////////////////////////////////////////////////////////////////
	chi2_0 = frbq->GetChisquare();
    NDF_0 = frbq->GetNDF();
    Prob = frbq->GetProb();
    COM = chi2_0/NDF_0;
	
	//std::stringstream oss_Sep;
	oss_Sep   << setprecision(3) << chi2_0;
	//TLatex l00;
	l00.SetTextAlign(9);
	l00.SetTextSize(0.038);
	l00.SetLineWidth(2);
	l00.SetNDC();
	l00.DrawLatex(0.1358774,0.7424242, ("#chi^{2}: " + oss_Sep.str()).c_str());
	
	//std::stringstream oss_NDF_0;
	oss_NDF_0 << setprecision(3) << NDF_0;
	//TLatex l01;
	l01.SetTextAlign(9);
	l01.SetTextSize(0.038);
	l01.SetLineWidth(2);
	l01.SetNDC();
	l01.DrawLatex(.1358774,0.7811448, ("NDF: " + oss_NDF_0.str()).c_str());
	
	/*std::stringstream oss_Prob;
	oss_Prob << setprecision(3) << Prob;
	
	TLatex l03;
	l03.SetTextAlign(9);
	l03.SetTextSize(0.038);
	l03.SetLineWidth(2);
	l03.SetNDC();
	l03.DrawLatex(0.1369345,0.8299663, ("Prob: " + oss_Prob.str()).c_str());
	
	
	//std::stringstream oss_COM;
	oss_COM   << setprecision(3) << COM;
	//TLatex l04;
	l04.SetTextAlign(9);
	l04.SetTextSize(0.038);
	l04.SetLineWidth(2);
	l04.SetNDC();
	l04.DrawLatex(0.1358774,0.6952862, ("#chi^{2}/NDF: " + oss_COM.str()).c_str());
	
	//TLatex lwork;
	lwork.SetTextAlign(9);
	lwork.SetTextSize(0.058);
	lwork.SetLineWidth(2);
	lwork.SetNDC();
	lwork.DrawLatex(0.1258774,0.8552862, ("ATLAS work-in-progress"));
	
    //TLatex lfile;
	lfile.SetTextAlign(9);
	lfile.SetTextSize(0.058);
	lfile.SetLineWidth(2);
	lfile.SetNDC();
	lfile.DrawLatex(0.1258774,0.9552862, (Base.c_str()));
	////////////////////////////////////////////////////////////////////	
	//Legend
	////////////////////////////////////////////////////////////////////    
    TLegend *leg2 = new TLegend(0.7078059,0.6690754,0.8987342,0.8986966,NULL,"brNDC");
    leg2->SetBorderSize(1);
    leg2->SetLineColor(1);
    leg2->SetLineStyle(1);
    leg2->SetLineWidth(1);
    leg2->SetFillColor(0);
    leg2->SetFillStyle(1001);
    leg2->AddEntry(h3,"Simulation","lep");
    leg2->AddEntry(frbq,"Fit","l");
    leg2->AddEntry(f25,"Gauss","l");
    leg2->AddEntry(f26,"Gauss","l");
    leg2->AddEntry(f27,"Landau_n","l");
    leg2->Draw();	
	////////////////////////////////////////////////////////////////////
   padr2->cd();   
       // pad2 becomes the current pad      
  TH1F* h_ratior = (TH1F*) FitRatio(frbq,h3);
  h_ratior->GetYaxis()->SetRangeUser(20,20);
  h_ratior->GetXaxis()->SetTitle("R_{bq}^{reco}[GeV]");
  h_ratior->SetMaximum(100);    
  h_ratior->SetMinimum(-100);  
  h_ratior->SetLabelSize(0.14, "X");
  h_ratior->SetLabelOffset(0.05,"X");
  h_ratior->SetTitleOffset(1.2,"X");
  h_ratior->SetTitleSize(0.17, "X");
  h_ratior->SetLabelSize(0.10, "Y");
  h_ratior->SetNdivisions(305,"Y");  
  h_ratior->GetYaxis()->SetTitleOffset(1.0);
  h_ratior->GetYaxis()->SetTitleFont(63);
  h_ratior->GetYaxis()->SetTitleSize(20);
  h_ratior->GetYaxis()->SetTitle("(Sim.-fit)/#sigma");
  h_ratior->SetLineWidth(3);    
  h_ratior->SetMarkerStyle(20);	  
  h_ratior->SetMarkerSize(1);	  
  h_ratior->SetLineWidth(3); 
  h_ratior->Draw("P");

  //TLine* tline_0 = new TLine(fitmin,0,fitmax,0);
  tline_0->SetLineColor(1);
  tline_0->SetLineStyle(1);
  tline_0->SetLineWidth(2);
  tline_0->Draw("SAME");
 
  c3->SaveAs((out+"rbq_reco.png").c_str());
 }
 
*/ 
 
