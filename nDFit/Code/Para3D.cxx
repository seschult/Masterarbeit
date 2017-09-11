#include "TSystem.h"
#include "TFile.h"
#include "TF1.h"
#include "TTree.h"
#include "TH1F.h"
#include "TROOT.h"
#include "TPad.h"
#include "TMath.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TMinuit.h"
#include "TAttMarker.h"

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
#include <omp.h>


using namespace std;
//Fit functions
//Gauss
Double_t gauss(Double_t *x, Double_t *par) {return par[0]*exp(-0.5*pow((x[0]-par[1])/par[2],2));}
//Landau
Double_t landau(Double_t *x,  Double_t *par) { return par[0]*(TMath::Landau(x[0], par[1], par[2])/par[2]);}
//Landau_n
Double_t landau_n(Double_t *x,  Double_t *par) {return par[0]*(TMath::Landau(2*par[1]-x[0], par[1], par[2])/par[2]);}
//top mass fit function
Double_t mtop_function(Double_t *x,  Double_t *par) { return gauss(x, par) + landau(x, &par[3]) + landau_n(x, &par[6]);}



inline bool check_file_exist (const std::string& name) {
	
  struct stat buffer;   
  return (stat (name.c_str(), &buffer) == 0); 
}


bool replace(std::string& str, const std::string& from, const std::string& to) {
	

    size_t start_pos = str.find(from);
    if(start_pos == std::string::npos)
        return false;
    str.replace(start_pos, from.length(), to);
    return true;
}




void big_ftop_fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
double chisq[200];
double chisqG=0;

#pragma omp parallel for 
for (int j=0;j< gALL.size();j++)
{
double mtoptemplate=gParameters[j][0];	

if (j==3) mtoptemplate=par[36];


chisq[j]=0;
double localpar[9];
for (int k=0;k<9;k++) localpar[k]=par[4*k+0]+par[4*k+1]*(mtoptemplate-172.5)+par[4*k+2]*(gParameters[j][1]-1.0)+par[4*k+3]*(gParameters[j][2]-1.0);
   gALL[j]->ftop->SetParameters(localpar);
   for (int i=1;i<gALL[j]->h1->GetNbinsX()+1;i++)
   {
   if (gALL[j]->h1->GetBinLowEdge(i)<130) continue;
   if (gALL[j]->h1->GetBinLowEdge(i)+gALL[j]->h1->GetBinWidth(i)>210) continue;
  double  v=std::abs(gALL[j]->ftop->Integral(gALL[j]->h1->GetBinLowEdge(i),gALL[j]->h1->GetBinLowEdge(i)+gALL[j]->h1->GetBinWidth(i))/(gALL[j]->h1->GetBinWidth(i)));
  // double v=gALL[j]->ftop->Eval(gALL[j]->h1->GetBinLowEdge(i)+0.5*gALL[j]->h1->GetBinWidth(i));
   if (gALL[j]->h1->GetBinContent(i)<0.001) continue;
   if (gALL[j]->h1->GetBinError(i)<0.000001) continue;
   chisq[j]+=pow((v-gALL[j]->h1->GetBinContent(i))/gALL[j]->h1->GetBinError(i),2);
   }
}
  


for (int k=0;k< gALL.size();k++) chisqG+=chisq[k];
   
   f = chisqG;
} 



std::vector<TemplateHolder*> gALL;
std::vector<std::vector<double> >  gParameters;
std::vector<double> gMTOP={170,171.5,172.5,173.5,175};
std::vector<double> gbJSF={0.96, 0.98, 1.0,   1.02, 1.04};
std::vector<double> gJSF={0.96, 0.98, 1.0,   1.02, 1.04};






//Constructor
TemplateHolder :: TemplateHolder (c){
    f00 = new TF1("f5",gauss,130,210,3);
    f01 = new TF1("f1",landau,130,210,3);
    f02 = new TF1("f2",landau_n,130,210,3);
    ftop = new TF1("ftop",mtop_function,130,210,9);    
    // Open file and histogram 
  
}


void TemplateHolder:: input(std::string InputFolder){
	
	for (std::vector<double>::iterator it1=gMTOP.begin(); it1!=gMTOP.end();it1++)
		for (std::vector<double>::iterator it2=gbJSF.begin(); it2!=gbJSF.end();it2++)
			for (std::vector<double>::iterator it3=gJSF.begin(); it3!=gJSF.end();it3++)
			{  
			  std::string file_name=Form("%s/Output_PlotFactory_2.4.24_JSF_%2.2f_bJSF_%2.2f/HistogramFolder/Merge_emujets_comb_2incl_1excl_nominal/Merged_ttbar_%4.1f.root",InputFolder.c_str(),*it3,*it2,*it1);  	
			  
			  replace(file_name,"172.5","");
			  replace(file_name,".5","p5");
		      replace(file_name,"1.00","1.0");
		      replace(file_name,"1.00","1.0");
		      replace(file_name,".0.root",".root");
				
			  if (!check_file_exist(file_name)) continue;
			  
			   puts(file_name.c_str());
		 
		  
			   std::vector<double> par;
			   par.push_back(*it1);	  
		       par.push_back(*it3);	  
		       par.push_back(*it2);	  
		  
		  
		       gParameters.push_back(par);
		       samples.push_back(file_name);
	
		     };
		     
	system("mkdir -p out");	     
	for(size_t i=0; i<samples.size(); i++)  gALL.push_back( new TemplateHolder(samples[i].c_str())); 
	  puts(File);  
    file = new TFile(File,"r");
    h1 = (TH1F*)file->Get("hist_klf_mtop_window");	     
	}

void TemplateHolder :: top_fit(std::vector<TemplateHolder*> gALL){

	TMinuit *gMinuit = new TMinuit(9*4+1);  //initialize TMinuit with a maximum of 5 params
	gMinuit->SetFCN(big_ftop_fcn);
	
	
	Double_t arglist[10];
	Int_t ierflg = 0;   
	arglist[0] = 1500000;
	arglist[1] = 0.0001;
    int i=0;

	/Set Fit Parameters
//1.10618e+04, 1.62284e+02, 1.28988e+01, 4.28790e+05, 1.87218e+02, 1.48783e+01, 4.46906e+05, 1.36078e+02, 1.02262e+01   
	gMinuit->mnparm(i, Form("a%i",i), 1.10618e+04, 1000.0, 1000,1000000,ierflg); i++;//zentralparameter
	gMinuit->mnparm(i, Form("a%i",i), 0.0, 0.1, 0,0,ierflg); i++;//tompmasse
	gMinuit->mnparm(i, Form("a%i",i), 0.0, 0.1, 0,0,ierflg); i++;//bJSf
	gMinuit->mnparm(i, Form("a%i",i), 0.0, 0.1, 0,0,ierflg); i++;//JSF
	gMinuit->mnparm(i, Form("a%i",i), 1.62284e+02, 10, 10,1000.0,ierflg); i++;
	gMinuit->mnparm(i, Form("a%i",i), 0.0, 0.1, 0,0,ierflg); i++;
	gMinuit->mnparm(i, Form("a%i",i), 0.0, 0.1, 0,0,ierflg); i++;
	gMinuit->mnparm(i, Form("a%i",i), 0.0, 0.1, 0,0,ierflg); i++;
	gMinuit->mnparm(i, Form("a%i",i), 1.28988e+01, 0.1, 0,0,ierflg); i++;
	gMinuit->mnparm(i, Form("a%i",i), 0.0, 0.1, 0,0,ierflg); i++;
	gMinuit->mnparm(i, Form("a%i",i), 0.0, 0.1, 0,0,ierflg); i++;
	gMinuit->mnparm(i, Form("a%i",i), 0.0, 0.1, 0,0,ierflg); i++;
	gMinuit->mnparm(i, Form("a%i",i), 4.28790e+05, 0.1, 1000,10000000,ierflg); i++;
	gMinuit->mnparm(i, Form("a%i",i), 0.0, 0.1, 0,0,ierflg); i++;
	gMinuit->mnparm(i, Form("a%i",i), 0.0, 0.1, 0,0,ierflg); i++;
	gMinuit->mnparm(i, Form("a%i",i), 0.0, 0.1, 0,0,ierflg); i++;
	gMinuit->mnparm(i, Form("a%i",i), 1.87218e+02, 0.1, 0,0,ierflg); i++;
	gMinuit->mnparm(i, Form("a%i",i), 0.0, 0.1, 0,0,ierflg); i++;
	gMinuit->mnparm(i, Form("a%i",i), 0.0, 0.1, 0,0,ierflg); i++;
	gMinuit->mnparm(i, Form("a%i",i), 0.0, 0.1, 0,0,ierflg); i++;
	gMinuit->mnparm(i, Form("a%i",i), 1.48783e+01, 0.1, 0,0,ierflg); i++;
	gMinuit->mnparm(i, Form("a%i",i), 0.0, 0.1, 0,0,ierflg); i++;
	gMinuit->mnparm(i, Form("a%i",i), 0.0, 0.1, 0,0,ierflg); i++;
	gMinuit->mnparm(i, Form("a%i",i), 0.0, 0.1, 0,0,ierflg); i++;
	gMinuit->mnparm(i, Form("a%i",i), 4.46906e+05, 0.1, 1000,1000000.0,ierflg); i++;
	gMinuit->mnparm(i, Form("a%i",i), 0.0, 0.1, 0,0,ierflg); i++;
	gMinuit->mnparm(i, Form("a%i",i), 0.0, 0.1, 0,0,ierflg); i++;
	gMinuit->mnparm(i, Form("a%i",i), 0.0, 0.1, 0,0,ierflg); i++;
	gMinuit->mnparm(i, Form("a%i",i),  1.36078e+02, 0.1, 0,0,ierflg); i++;
	gMinuit->mnparm(i, Form("a%i",i), 0.0, 0.1, 0,0,ierflg); i++;
	gMinuit->mnparm(i, Form("a%i",i), 0.0, 0.1, 0,0,ierflg); i++;
	gMinuit->mnparm(i, Form("a%i",i), 0.0, 0.1, 0,0,ierflg); i++;
	gMinuit->mnparm(i, Form("a%i",i), 1.02262e+01, 0.1, 0,0,ierflg); i++;
	gMinuit->mnparm(i, Form("a%i",i), 0.0, 0.1, 0,0,ierflg); i++;
	gMinuit->mnparm(i, Form("a%i",i), 0.0, 0.1, 0,0,ierflg); i++;
	gMinuit->mnparm(i, Form("a%i",i), 0.0, 0.1, 0,0,ierflg); i++;

	//TopMass
	gMinuit->mnparm(i, Form("a%i",i), 170.0, 0.1, 165,180,ierflg); i++;
	
	
	//Fit   
	gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);
	gMinuit->mnexcm("SIMPLEX", arglist ,2,ierflg);
	gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);
	gMinuit->mnexcm("SIMPLEX", arglist ,2,ierflg);
	gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);
	gMinuit->mnexcm("SIMPLEX", arglist ,2,ierflg);
	gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);
	gMinuit->mnrset(0);
	gMinuit->mnexcm("SIMPLEX", arglist ,2,ierflg);
	gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);
	gMinuit->mnrset(0);
	gMinuit->mnexcm("SIMPLEX", arglist ,2,ierflg);
	gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);
	gMinuit->mnrset(0);
	gMinuit->mnexcm("SIMPLEX", arglist ,2,ierflg);
	gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);

	gMinuit->mnexcm("HESSE", arglist ,2,ierflg);
	
}


void TemplateHolder :: mtop_plot(){ 
	

for(size_t i=0; i<samples.size(); i++){ 
	
	
	out = Form("./out_JSF_%2.2f_bJSF_%2.2f_mtop_%3.1f.png",gParameters[i][0],gParameters[i][1],gParameters[i][2]);   

	
	
  TCanvas* c1 = new TCanvas("c1","c1", 950, 900);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.32, 1, 1.0);
  pad1->SetBottomMargin(1); // Upper and lower plot are joined
  pad1->SetGridx();         // Vertical grid
  pad1->Draw();             // Draw the upper pad: pad1
  pad1->cd();               // pad1 becomes the current pad
  
  h1->GetXaxis()->SetLabelFont(63);
  h1->GetXaxis()->SetLabelSize(0);
  h1->GetXaxis()->SetTitleSize(0.0);
  h1->GetXaxis()->SetRangeUser(124,214);  
  h1->GetYaxis()->SetLabelFont(63);
  h1->GetYaxis()->SetLabelSize(28);
  h1->GetYaxis()->SetRangeUser(-73,1.1* h1->GetMaximum());
  h1->GetYaxis()->SetTitleOffset(1.45);
  h1->GetYaxis()->SetTitleFont(63);
  h1->GetYaxis()->SetTitleSize(28);
  h1->GetYaxis()->SetTitle("Events");
  h1->SetLineWidth(3);
  
  ftop->SetFillColor(19);
  ftop->SetFillStyle(0);
  ftop->SetLineColor(2);
  ftop->SetLineWidth(3);
 
  h1->Draw();
  ftop->Draw("SAME");
	
  double p[9];
  
  for(int i = 0; i < 9; i++) p[i] = ftop->GetParameter(i);
  
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
  f02->SetLineColor(12);
  f02->SetFillStyle(0); 
  f02->SetLineWidth(3);  
  f02->Draw("SAME");
  
  
  
   double chi2_0 = ftop->GetChisquare();
   double NDF_0 = ftop->GetNDF();
   double COM = chi2_0/NDF_0;
	
  
  
  TLatex l1;
  l1.SetTextAlign(9);
  l1.SetTextSize(0.058);
  l1.SetLineWidth(2);
  l1.SetNDC();
  l1.DrawLatex(0.1258774,0.8552862, ("ATLAS work-in-progress"));


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

  std::stringstream oss_COM;
  oss_COM   << setprecision(3) << COM;
  TLatex l04;
  l04.SetTextAlign(9);
  l04.SetTextSize(0.038);
  l04.SetLineWidth(2);
  l04.SetNDC();
  l04.DrawLatex(0.1358774,0.6952862, ("#chi^{2}/NDF: " + oss_COM.str()).c_str());








  
  //Legend
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
  leg->SetTextFont(63);
  leg->SetTextSize(25);
  leg->Draw();

  c1->cd(); 
  
  
  
  
//comparison plot between fit and MC
  TPad *pad2 = new TPad("pad2", "pad2", 0, 0.08, 1, 0.39);
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.2);
  pad2->SetGridx(); // vertical grid
  pad2->SetGridy();
  pad2->Draw();
  pad2->cd();   

  //Calculate the diffrence between fit and histogramm h1
  TH1F *htop = new TH1F("htop", "dif", h1 -> GetNbinsX(), 120, 220);  
  
  for (int bin=2; bin<=h1 -> GetNbinsX();bin++) 
	{    
    double x = h1->GetXaxis()->GetBinCenter(bin);
    double fval = ftop->Eval(x);
    double sub = h1->GetBinContent(bin)-fval;
    
    if(sub > -100 && x > 130 ){
		htop->SetBinContent(bin, sub);
		double   err = h1->GetBinError(bin);
		double   div = sub/err;
		htop->SetBinContent(bin, div);
		}else continue;
	};
  
  
  htop->GetYaxis()->SetRangeUser(-2.8,2.8);
  
  float lower_edge  = htop -> GetBinLowEdge(1);
  float number_bins = htop -> GetNbinsX();
  float upper_edge = htop -> GetBinLowEdge(number_bins) + htop->GetBinWidth(number_bins);
  
  htop->GetXaxis()->SetRangeUser(124,214);
  htop->GetXaxis()->SetLabelFont(63);
  htop->GetXaxis()->SetLabelSize(28);
  htop->GetXaxis()->SetTitle("m_{top}^{reco}[GeV]");
  htop->GetXaxis()->SetTitleSize(25);
  htop->GetXaxis()->SetTitleOffset(3.0);
  htop->GetXaxis()->SetTitleFont(63);
  htop->GetYaxis()->SetLabelFont(63);
  htop->GetYaxis()->SetLabelSize(28);
  htop->GetYaxis()->SetTitleSize(28);
  htop->GetYaxis()->SetTitleFont(63);
  htop->GetYaxis()->SetTitleOffset(1.4);
  htop->GetYaxis()->SetTitle("(Sim.-Fit)/Error[#sigma]");
  
  htop->SetLineWidth(3);   
  htop->SetMarkerStyle(15); 
  htop->Draw("P");
  
 
 //red line @ 0 in comparison plot
  TF1 *norm1 = new TF1("fa1","0", lower_edge, upper_edge);
  norm1 -> SetLineColor(kRed);
  norm1 -> SetLineStyle(1);
  norm1 -> SetLineWidth(3);
  norm1->Draw("SAME");

  c1->SaveAs(out.c_str());
  std::string Name = file->GetName();
  
  };
}
	

