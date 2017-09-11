



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


	

class TemplateHolder{
	public: 	
	TemplateHolder(const char *File); //constructor		
	void mtop_plot(std::string out); //constructor
	void top_fit();
	
	TFile *file;//
	TH1F *h1; //
	TF1 *ftop;//
	TF1 *f00;//
	TF1 *f01;//
	TF1 *f02;//		
};






//Constructor
TemplateHolder :: TemplateHolder (const char *File){
    f00 = new TF1("f5",gauss,130,210,3);
    f01 = new TF1("f1",landau,130,210,3);
    f02 = new TF1("f2",landau_n,130,210,3);
    ftop = new TF1("ftop",mtop_function,130,210,9);    
    // Open file and histogram 
    puts(File);  
    file = new TFile(File,"r");
    h1 = (TH1F*)file->Get("hist_klf_mtop_window");
}


std::vector<TemplateHolder*> gALL;
std::vector<std::vector<double> >  gParameters;
std::vector<std::vector<std::string> >  gParametersout;
std::vector<double> gMTOP={170,171.5,172.5,173.5,175};
std::vector<double> gbJSF={0.96, 0.98, 1.0,   1.02, 1.04};
std::vector<double> gJSF={0.96, 0.98, 1.0,   1.02, 1.04};
std::vector<std::string> parout;


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

/*

  4.19640e+06   4.59254e+01   2.00100e+00   4.19640e+06
   2  a1           2.40930e+04   2.67981e-01   1.14884e-02   2.40930e+04
   3  a2          -6.89471e+03   6.77375e+02   1.10979e+01  -6.89471e+03
   4  a3          -3.12643e+02   2.65583e-01   1.36819e-03  -3.12643e+02
   5  a4           5.23514e+03   8.48202e-02   2.49631e-03   5.23514e+03
   6  a5           3.03757e+01   4.93148e-04   1.44843e-05   3.03757e+01
   7  a6           6.31281e+01   1.81862e+00   2.03043e-02   6.31281e+01
   8  a7           8.12061e-01   4.94936e-04   2.51007e-06   8.12061e-01
   9  a8           3.64765e+03   6.50285e-02   1.73933e-03   3.64765e+03
  10  a9           2.14348e+01   3.81078e-04   1.02209e-05   2.14348e+01
  11  a10          6.26267e+00   1.42928e+00   1.56238e-02   6.26267e+00
  12  a11          2.57372e-01   3.78362e-04   1.92776e-06   2.57372e-01
  13  a12          2.18391e+08   1.03226e+04   1.04137e+02   2.18391e+08
  14  a13          1.69712e+06   5.99055e+01   8.09249e-01   1.69712e+06
  15  a14          4.29849e+06   3.34363e+05   2.23200e+02   4.29849e+06
  16  a15          4.36293e+05   6.10954e+01   2.97774e-01   4.36293e+05
  17  a16          2.24191e+04   8.01251e-01   1.06903e-02   2.24191e+04
  18  a17          1.56974e+02   4.59926e-03   7.48512e-05   1.56974e+02
  19  a18          2.81399e+02   2.20199e+01   1.87716e-01   2.81399e+02
  20  a19          2.79414e+01   4.69495e-03   2.33064e-05   2.79414e+01
  21  a20         -1.37058e+04   5.13763e-01   6.53543e-03  -1.37058e+04
  22  a21         -6.31535e+01   3.08940e-03   3.01139e-05  -6.31535e+01
  23  a22          1.72407e+02   1.41924e+01   1.20620e-01   1.72407e+02
  24  a23          1.72915e+01   3.07220e-03   1.54220e-05   1.72915e+01
  25  a24         -1.16842e+07   3.02985e+03   1.51272e+01  -1.16842e+07
  26  a25         -8.80409e+04   1.76586e+01   8.84065e-02  -8.80409e+04
  27  a26          1.83437e+05   8.23490e+04   2.14030e+02   1.83437e+05
  28  a27         -1.81404e+04   1.74533e+01   8.84655e-02  -1.81404e+04
  29  a28         -5.01715e+02   2.02879e-01   1.03405e-03  -5.01715e+02
  30  a29         -3.18121e+00   1.18364e-03   6.02946e-06  -3.18121e+00
  31  a30          3.62313e+01   3.90329e+00   4.92543e-02   3.62313e+01
  32  a31          5.29738e-01   1.18823e-03   6.05097e-06   5.29738e-01
  33  a32          2.32415e+02   1.00732e-01   5.16715e-04   2.32415e+02
  34  a33          1.22134e+00   5.93424e-04   3.01290e-06   1.22134e+00
  35  a34          1.43310e+00   2.06276e+00   2.45469e-02   1.43310e+00
  36  a35         -7.39338e-02   5.91137e-04   3.02254e-06  -7.39338e-02
*/




void TemplateHolder :: top_fit(){
    TMinuit *gMinuit = new TMinuit(9*4+1);  //initialize TMinuit with a maximum of 5 params
	gMinuit->SetFCN(big_ftop_fcn);
	Double_t arglist[10];
	Int_t ierflg = 0;   
	arglist[0] = 1500000;
	arglist[1] = 0.0001;
    int i=0;
   //Set Fit Parameters
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
  /////////////////////////////////////////////////////////////////////////////////////////////
}	






void TemplateHolder :: mtop_plot(std::string out){ 
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
}

	












int main(int argc, char* argv[])
{

	
	ROOT::EnableThreadSafety();
	
	
if (argc>2) exit(1);
  std::string InputFolder    = argv[1];
  std::vector<std::string> samples;  
  for (std::vector<double>::iterator it1=gMTOP.begin(); it1!=gMTOP.end();it1++)
	for (std::vector<double>::iterator it2=gbJSF.begin(); it2!=gbJSF.end();it2++)
		for (std::vector<double>::iterator it3=gJSF.begin(); it3!=gJSF.end();it3++)
		{  
		  std::string file_name=Form("%s/Output_PlotFactory_2.4.24_JSF_%2.2f_bJSF_%2.2f/HistogramFolder/Merge_emujets_comb_2incl_1excl_nominal/Merged_ttbar_%4.1f.root",InputFolder.c_str(),*it3,*it2,*it1);  
		  //puts(file_name.c_str());
		  replace(file_name,"172.5","");
		  replace(file_name,".5","p5");
		  replace(file_name,"1.00","1.0");
		  replace(file_name,"1.00","1.0");
		  replace(file_name,".0.root",".root");
		  
		  puts(file_name.c_str());
		  
		  if (!check_file_exist(file_name)) continue;
		  
		  puts(file_name.c_str());
		 
		  
		  std::vector<double> par;
		  
		  par.push_back(*it1);	  
		  par.push_back(*it3);	  
		  par.push_back(*it2);	  
		  
		  
		  gParameters.push_back(par);
		  samples.push_back(file_name);
		
		
		
		
		
		  std::ostringstream strs0;
		  std::ostringstream strs1;
		  std::ostringstream strs2;
		  strs0 << *it3;
		  strs1 << *it2;
		  strs2 << *it1;
		  parout.push_back(strs0.str());
		  parout.push_back(strs1.str());
		  parout.push_back(strs2.str());
		  gParametersout.push_back(parout);
		 
				
		}	  



system("mkdir -p out");
for(size_t i=0; i<samples.size(); i++)  gALL.push_back( new TemplateHolder(samples[i].c_str()));
///////////////////////////////////////////////////////////////////////////////////////////////



//gALL[0]->top_fit("h1");
gALL[0]->top_fit();

for(size_t i=0; i<samples.size(); i++) { 
	
	
	
	
	
//	std::string Nameout =  gParametersout[i][1]+"_"+  gParametersout[i][2] +"_"+  gParametersout[i][3];
	//gALL[i]->mtop_plot(Form("./out_%i.png",i));
	gALL[i]->mtop_plot(Form("./out_JSF_%2.2f_bJSF_%2.2f_mtop_%3.1f.png",gParameters[i][0],gParameters[i][1],gParameters[i][2]));
	//gALL[i]->mtop_plot(Nameout);

};



for(size_t i=0; i<samples.size(); i++) { 
	
	std::string Nameout =  gParametersout[i][2];
	
	std::cout<<Nameout<<std::endl;
	};

return 0;
}

