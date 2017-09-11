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
	void mtop_fit(); //constructor
	TFile *file;//
	TH1F *h1; //
	TF1 *ftop;//
	TF1 *f00;//
	TF1 *f01;//
	TF1 *f02;//	
};

//Fit Functions
//Gauss
Double_t gauss(Double_t *x, Double_t *par) {return par[0]*exp(-0.5*pow((x[0]-par[1])/par[2],2));}
//Landau
Double_t landau(Double_t *x,  Double_t *par) { return par[0]*(TMath::Landau(x[0], par[1], par[2])/par[2]);}
//Landau_n
Double_t landau_n(Double_t *x,  Double_t *par) {return par[0]*(TMath::Landau(2*par[1]-x[0], par[1], par[2])/par[2]);}
//top mass fit function
Double_t mtop_function(Double_t *x,  Double_t *par) { return gauss(x, par) + landau(x, &par[3]) + landau_n(x, &par[6]);}

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
  
  //Fit of the histogram   
  ftop->SetFillColor(19);
  ftop->SetFillStyle(0);
  ftop->SetLineColor(2);
  ftop->SetLineWidth(3);
  h1->SetLineWidth(3);
  
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
  
  TLatex l05;
  l05.SetTextAlign(9);
  l05.SetTextSize(0.058);
  l05.SetLineWidth(2);
  l05.SetNDC();
  l05.DrawLatex(0.1258774,0.8552862, ("ATLAS work-in-progress"));
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
  TPad *pad2 = new TPad("pad2", "pad2", 0, 0.08, 1, 0.39);
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.2);
  pad2->SetGridx(); // vertical grid
  pad2->SetGridy();
  pad2->Draw();
  pad2->cd();   

  //Calculate the diffrence between Fitfunction and Histogramm h1
  TH1F *htop = new TH1F("htop", "dif", h1 -> GetNbinsX(), 120, 220);  
  for (int bin=2; bin<=h1 -> GetNbinsX();bin++) {    
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
  htop->Draw();
  
  TF1 *norm1 = new TF1("fa1","0", lower_edge, upper_edge);
  norm1 -> SetLineColor(kRed);
  norm1 -> SetLineStyle(1);
  norm1 -> SetLineWidth(3);
  norm1->Draw("SAME");

  c1->SaveAs(out.c_str());
  std::string Name = file->GetName();
}

	
std::vector<TemplateHolder*> gALL;
std::vector<std::vector<double> >  gParameters;
std::vector<double> gMTOP={170,171.5,172.5,173.5,175};
std::vector<double> gbJSF={0.96, 0.98, 1.0,   1.02, 1.04};
std::vector<double> gJSF={0.96, 0.98, 1.0,   1.02, 1.04};
std::vector<std::vector<std::string> >  gParametersout;
std::vector<std::string> parout;

void big_ftop_fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
double chisq[200];
double chisqG=0;

#pragma omp parallel for 
for (int j=0;j< gALL.size();j++)
{
chisq[j]=0;
double localpar[9];
for (int k=0;k<9;k++) localpar[k]=par[4*k+0]+par[4*k+1]*(gParameters[j][0]-172.5)+par[4*k+2]*(gParameters[j][1]-1.0)+par[4*k+3]*(gParameters[j][2]-1.0);
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


void TemplateHolder:: mtop_fit(){
	   TMinuit *gMinuit = new TMinuit(9*4);  //initialize TMinuit with a maximum of 5 params
   gMinuit->SetFCN(big_ftop_fcn);
   Double_t arglist[10];
   Int_t ierflg = 0;   
   arglist[0] = 15000;
   arglist[1] = 0.001;
   int i=0;
//1.10618e+04, 1.62284e+02, 1.28988e+01, 4.28790e+05, 1.87218e+02, 1.48783e+01, 4.46906e+05, 1.36078e+02, 1.02262e+01   
gMinuit->mnparm(i, Form("a%i",i), 1.10618e+04, 0.1, 0,0,ierflg); i++;
gMinuit->mnparm(i, Form("a%i",i), 0.0, 0.1, 0,0,ierflg); i++;
gMinuit->mnparm(i, Form("a%i",i), 0.0, 0.1, 0,0,ierflg); i++;
gMinuit->mnparm(i, Form("a%i",i), 0.0, 0.1, 0,0,ierflg); i++;
gMinuit->mnparm(i, Form("a%i",i), 1.62284e+02, 0.1, 0,0,ierflg); i++;
gMinuit->mnparm(i, Form("a%i",i), 0.0, 0.1, 0,0,ierflg); i++;
gMinuit->mnparm(i, Form("a%i",i), 0.0, 0.1, 0,0,ierflg); i++;
gMinuit->mnparm(i, Form("a%i",i), 0.0, 0.1, 0,0,ierflg); i++;
gMinuit->mnparm(i, Form("a%i",i), 1.28988e+01, 0.1, 0,0,ierflg); i++;
gMinuit->mnparm(i, Form("a%i",i), 0.0, 0.1, 0,0,ierflg); i++;
gMinuit->mnparm(i, Form("a%i",i), 0.0, 0.1, 0,0,ierflg); i++;
gMinuit->mnparm(i, Form("a%i",i), 0.0, 0.1, 0,0,ierflg); i++;
gMinuit->mnparm(i, Form("a%i",i), 4.28790e+05, 0.1, 0,0,ierflg); i++;
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
gMinuit->mnparm(i, Form("a%i",i), 4.46906e+05, 0.1, 0,0,ierflg); i++;
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
   
   gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);
   gMinuit->mnexcm("SIMPLEX", arglist ,2,ierflg);
   gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);
   gMinuit->mnexcm("SIMPLEX", arglist ,2,ierflg);
   gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);
   gMinuit->mnexcm("HESSE", arglist ,2,ierflg);
}













int main(int argc, char* argv[])
{
if (argc>2) exit(1);


ROOT::EnableThreadSafety();

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
  par.push_back(*it3);	  
  par.push_back(*it2);	  
  par.push_back(*it1);	  
  gParameters.push_back(par);
  samples.push_back(file_name);
}	  

system("mkdir -p out");
  for(size_t i=0; i<samples.size(); i++)  gALL.push_back( new TemplateHolder(samples[i].c_str()));
///////////////////////////////////////////////////////////////////////////////////////////////




  /////////////////////////////////////////////////////////////////////////////////////////////
for(size_t i=0; i<samples.size(); i++) {
gALL[i]->mtop_fit();
//std::string Output = "~/" +  From("out_%i.png", i) ;

 gALL[i]->mtop_plot(Form("./out_%i.png",i));
 
};
return 0;



}
