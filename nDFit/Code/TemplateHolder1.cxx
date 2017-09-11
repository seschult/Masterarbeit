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

using namespace std;

////////////////////////////////////////////////////////////////////////////////////////////////////
class TemplateHolder{


public:
void plot(std::string out);
double fmtop;
double fJSF;
double fbJSF;
TFile *file;//
TH1F *h1; //
TH1F *h2; //
TH1F *h3; //
TF1 *ftop;//
TF1 *fw;//
TF1 *frbq;
TF1 *f13;//
TF1 *f14;//
TF1 *f00;//
TF1 *f01;//
TF1 *f02;//		
TF1 *f25;		
TF1 *f26;		
TF1 *f27;

Double_t gauss(Double_t *x, Double_t *par);// {return par[0]*exp(-0.5*pow((x[0]-par[1])/par[2],2));}
//Landau
Double_t landau(Double_t *x,  Double_t *par);// { return par[0]*(TMath::Landau(x[0], par[1], par[2])/par[2]);}
//Landau_n
Double_t landau_n(Double_t *x,  Double_t *par);// {return par[0]*(TMath::Landau(2*par[1]-x[0], par[1], par[2])/par[2]);}
//top mass fit function
Double_t mtop_function(Double_t *x,  Double_t *par) ;//{ return gauss(x, par) + landau(x, &par[3]) + landau_n(x, &par[6]);}
//fit W masse
Double_t fit_mw(Double_t *x,  Double_t *par);// { return gauss(x,par) + gauss(x,&par[3]);}
//fit Rbq
Double_t fit_rbq(Double_t *x,  Double_t *par) ;//{ return gauss(x,par) + gauss(x,&par[3]) + landau(x,&par[6]);}
TemplateHolder (const char *File, double t, double j, double b);
};

class FCN{
	

public:
std::vector<TemplateHolder*> gALL;

FCN(std::vector<TemplateHolder*> susi);
//Fit functions
//Gauss
double big_ftop_fcn(const Double_t *par) const
{
double f;	
double chisq[200];
double chisqG=0;
#pragma omp parallel for 
for (int j=0;j< gALL.size();j++)
{
	double mtoptemplate=gALL[j]->fmtop;
//if (j==3) mtoptemplate=par[36];
	chisq[j]=0;
	double localpar[9];
    for (int k=0;k<9;k++) localpar[k]=par[4*k+0]+par[4*k+1]*(mtoptemplate-172.5)+par[4*k+2]*(gALL[j]->fJSF-1.0)+par[4*k+3]*(gALL[j]->fbJSF-1.0);
    gALL[j]->ftop->SetParameters(localpar);
    for (int i=1;i<gALL[j]->h1->GetNbinsX()+1;i++)
	{
		if (gALL[j]->h1->GetBinLowEdge(i)<130) continue;
		if (gALL[j]->h1->GetBinLowEdge(i)+gALL[j]->h1->GetBinWidth(i)>210) continue;
		//double  v=std::abs(gALL[j]->ftop->Integral(gALL[j]->h1->GetBinLowEdge(i),gALL[j]->h1->GetBinLowEdge(i)+gALL[j]->h1->GetBinWidth(i))/(gALL[j]->h1->GetBinWidth(i)));
		double v=gALL[j]->ftop->Eval(gALL[j]->h1->GetBinCenter(i));
		if (gALL[j]->h1->GetBinContent(i)<0.001) continue;
		if (gALL[j]->h1->GetBinError(i)<0.000001) continue;
		chisq[j]+=pow((v-gALL[j]->h1->GetBinContent(i))/gALL[j]->h1->GetBinError(i),2);
   }
}
  
for (int k=0;k< gALL.size();k++) chisqG+=chisq[k];
   f = chisqG;
return f;
}  
double big_fw_fcn(const Double_t *par) const
{
double f;	
double chisq[200];
double chisqG=0;

#pragma omp parallel for 
for (int j=0;j< gALL.size();j++)
{
	



chisq[j]=0;
double localpar[6];

for (int k=0;k<6;k++) localpar[k]=par[4*k+0]+par[4*k+1]*(gALL[j]->fmtop-172.5)+par[4*k+2]*(gALL[j]->fJSF-1.0)+par[4*k+3]*(gALL[j]->fbJSF-1.0);

   
//printf("GG=%f\n",GG);
   
   gALL[j]->fw->SetParameters(localpar);
   
   for (int i=1;i<gALL[j]->h2->GetNbinsX()+1;i++)
   {
   if (gALL[j]->h2->GetBinLowEdge(i)<58) continue;
   if (gALL[j]->h2
   ->GetBinLowEdge(i)+gALL[j]->h2->GetBinWidth(i)>108) continue;
   
   //double  v=std::abs(gALL[j]->ftop->Integral(gALL[j]->h1->GetBinLowEdge(i),gALL[j]->h1->GetBinLowEdge(i)+gALL[j]->h1->GetBinWidth(i))/(gALL[j]->h1->GetBinWidth(i)));
double v=gALL[j]->fw->Eval(gALL[j]->h2->GetBinCenter(i));
   if (gALL[j]->h2->GetBinContent(i)<0.001) continue;
   if (gALL[j]->h2->GetBinError(i)<0.000001) continue;
   chisq[j]+=pow((v-gALL[j]->h2->GetBinContent(i))/gALL[j]->h2->GetBinError(i),2);
   }
}
  


for (int k=0;k< gALL.size();k++) chisqG+=chisq[k];
   
   f = chisqG;
   return f;
}   
double big_frbq_fcn(const Double_t *par) const
{



double f;
double chisq[200];
double chisqG=0;

#pragma omp parallel for 
for (int j=0;j< gALL.size();j++)
{
double rbqtemplate=gALL[j]->fmtop;	

//if (j==3) mtoptemplate=par[36];



chisq[j]=0;
double localpar[9];

for (int k=0;k<9;k++) localpar[k]=par[4*k+0]+par[4*k+1]*(rbqtemplate-172.5)+par[4*k+2]*(gALL[j]->fJSF-1.0)+par[4*k+3]*(gALL[j]->fbJSF-1.0);

   
   
   gALL[j]->frbq->SetParameters(localpar);
   
   for (int i=1;i<gALL[j]->h3->GetNbinsX()+1;i++)
   {
   if (gALL[j]->h3->GetBinLowEdge(i)<0.3) continue;
   if (gALL[j]->h3->GetBinLowEdge(i)+gALL[j]->h3->GetBinWidth(i)>3.9) continue;
   
   //double  v=std::abs(gALL[j]->frbq->Integral(gALL[j]->h3->GetBinLowEdge(i),gALL[j]->h3->GetBinLowEdge(i)+gALL[j]->h3->GetBinWidth(i))/(gALL[j]->h3->GetBinWidth(i)));
double v=gALL[j]->frbq->Eval(gALL[j]->h3->GetBinCenter(i));
   if (gALL[j]->h3->GetBinContent(i)<0.001) continue;
   if (gALL[j]->h3->GetBinError(i)<0.000001) continue;
   chisq[j]+=pow((v-gALL[j]->h3->GetBinContent(i))/gALL[j]->h3->GetBinError(i),2);
   }
}
  


for (int k=0;k< gALL.size();k++) chisqG+=chisq[k];
   
   f = chisqG;
   return f;
}   


//void big_all_fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag){
double big_all_fcn(const Double_t *par) const{
	double f;
	Int_t iflagq=0;

	int ftop_npar=36;
	double chi2_top= big_ftop_fcn(par);
	int fw_npar=24;
	double chi2_w=big_fw_fcn(&par[36]);
	int frbq_npar=36;
	double chi2_rbq=big_frbq_fcn(&par[36+24]);


	f=chi2_rbq+chi2_top+chi2_w;
	return f;
}	

};

FCN::FCN(std::vector<TemplateHolder*> susi){

gALL = susi;
 	
}



//////////////////////////////////////////////////////////////////////////////////////////////////////
class Fill3D{

public:
Fill3D(std::string InputFolder,std::vector<double> gMTOP,std::vector<double> gbJSF,std::vector<double> gJSF );
void fillup();
void top_fit();
std::vector<TemplateHolder*> gALL;
//std::vector<double> gMTOP={170,171.5,172.5,173.5,175};
//std::vector<double> gbJSF={0.96, 0.98, 1.0, 1.02, 1.04};
//std::vector<double> gJSF={0.96, 0.98, 1.0, 1.02, 1.04};

bool replace(std::string& str, const std::string& from, const std::string& to) {
	

    size_t start_pos = str.find(from);
    if(start_pos == std::string::npos)
        return false;
    str.replace(start_pos, from.length(), to);
    return true;
}

bool check_file_exist (const std::string& name) {
	
  struct stat buffer;   
  return (stat (name.c_str(), &buffer) == 0); 
}


};	

Fill3D::Fill3D(std::string InputFolder,std::vector<double> gMTOP,std::vector<double> gbJSF,std::vector<double> gJSF  ){

	for (std::vector<double>::iterator it1=gMTOP.begin(); it1!=gMTOP.end();it1++)
		for (std::vector<double>::iterator it2=gbJSF.begin(); it2!=gbJSF.end();it2++)
			for (std::vector<double>::iterator it3=gJSF.begin(); it3!=gJSF.end();it3++)
			{  
			  std::string file_name=Form("%s/Output_PlotFactory_2.4.24_JSF_%2.2f_bJSF_%2.2f/HistogramFolder/Merge_emujets_comb_2incl_1excl_nominal/Merged_ttbar_%4.1f.root",InputFolder.c_str(),*it3,*it2,*it1);  
			  //puts(file_name.c_str());
			  replace(file_name,"_172.5","");
			  replace(file_name,".5","p5");
			  replace(file_name,"1.00","1.0");
			  replace(file_name,"1.00","1.0");
			  replace(file_name,".0.root",".root");
			  
			  puts(file_name.c_str());
			  
			  if (!check_file_exist(file_name)) continue;
			  
			  puts(file_name.c_str());
			 gALL.push_back( new TemplateHolder(file_name.c_str(),*it1,*it3,*it2));
			  
			};
}


void Fill3D:: top_fit(){


	
	TMinuitMinimizer *aMinuit = new TMinuitMinimizer("wow",6*4+9*4+9*4); 
	 FCN* Mf3=new FCN(gALL);
    
    ROOT::Math::Functor f(Mf3,&FCN::big_all_fcn,96);
    
  
  
	aMinuit->SetFunction(f);
	Double_t arglist[10];
	Double_t arglistw[10];
	Double_t arglistr[10];
	Int_t ierflg = 0;   
	Int_t ierflgw = 0;   
	Int_t ierflgr = 0;   
	arglist[0] = 1500000;
	arglist[1] = 0.0001;
	arglistw[0] = 1500000;
	arglistw[1] = 0.0001;
	arglistr[0] = 1500000;
	arglistr[1] = 0.0001;
	int q=0;
  
  	aMinuit->SetLimitedVariable(q, Form("a%i",q), 1.10618e+04, 1000.0, 1000,1000000); q++;//zentralparameter
	aMinuit->SetLimitedVariable(q, Form("a%i",q), 0.0, 0.1, -10000,10000); q++;//tompmasse
	aMinuit->SetLimitedVariable(q, Form("a%i",q), 0.0, 0.1, -10000,10000); q++;//bJSf
	aMinuit->SetLimitedVariable(q, Form("a%i",q), 0.0, 0.1, -10000,10000); q++;//JSF
	aMinuit->SetLimitedVariable(q, Form("a%i",q), 1.62284e+02, 10, 10,1000.0); q++;
	aMinuit->SetLimitedVariable(q, Form("a%i",q), 0.0, 0.1, -10000,10000); q++;
	aMinuit->SetLimitedVariable(q, Form("a%i",q), 0.0, 0.1, -10000,10000); q++;
	aMinuit->SetLimitedVariable(q, Form("a%i",q), 0.0, 0.1, -10000,10000); q++;
	aMinuit->SetLimitedVariable(q, Form("a%i",q), 1.28988e+01, 0.1, 0,1000000); q++;
	aMinuit->SetLimitedVariable(q, Form("a%i",q), 0.0, 0.1, -10000,10000); q++;
	aMinuit->SetLimitedVariable(q, Form("a%i",q), 0.0, 0.1, -10000,10000); q++;
	aMinuit->SetLimitedVariable(q, Form("a%i",q), 0.0, 0.1, -10000,10000); q++;
	aMinuit->SetLimitedVariable(q, Form("a%i",q), 4.28790e+05, 0.1, 1000,10000000); q++;
	aMinuit->SetLimitedVariable(q, Form("a%i",q), 0.0, 0.1, -10000,10000); q++;
	aMinuit->SetLimitedVariable(q, Form("a%i",q), 0.0, 0.1, -10000,10000); q++;
	aMinuit->SetLimitedVariable(q, Form("a%i",q), 0.0, 0.1, -10000,10000); q++;
	aMinuit->SetLimitedVariable(q, Form("a%i",q), 1.87218e+02, 0.1, 0,1000000); q++;
	aMinuit->SetLimitedVariable(q, Form("a%i",q), 0.0, 0.1, -10000,10000); q++;
	aMinuit->SetLimitedVariable(q, Form("a%i",q), 0.0, 0.1, -10000,10000); q++;
	aMinuit->SetLimitedVariable(q, Form("a%i",q), 0.0, 0.1, -10000,10000); q++;
	aMinuit->SetLimitedVariable(q, Form("a%i",q), 1.48783e+01, 0.1, 0,1000000); q++;
	aMinuit->SetLimitedVariable(q, Form("a%i",q), 0.0, 0.1, -10000,10000); q++;
	aMinuit->SetLimitedVariable(q, Form("a%i",q), 0.0, 0.1, -10000,10000); q++;
	aMinuit->SetLimitedVariable(q, Form("a%i",q), 0.0, 0.1, -10000,10000); q++;
	aMinuit->SetLimitedVariable(q, Form("a%i",q), 4.46906e+05, 0.1, 1000,1000000.0); q++;
	aMinuit->SetLimitedVariable(q, Form("a%i",q), 0.0, 0.1, -10000,10000); q++;
	aMinuit->SetLimitedVariable(q, Form("a%i",q), 0.0, 0.1, -10000,10000); q++;
	aMinuit->SetLimitedVariable(q, Form("a%i",q), 0.0, 0.1, -10000,10000); q++;
	aMinuit->SetLimitedVariable(q, Form("a%i",q),  1.36078e+02, 0.1, 0,1000000); q++;
	aMinuit->SetLimitedVariable(q, Form("a%i",q), 0.0, 0.1, -10000,10000); q++;
	aMinuit->SetLimitedVariable(q, Form("a%i",q), 0.0, 0.1, -10000,10000); q++;
	aMinuit->SetLimitedVariable(q, Form("a%i",q), 0.0, 0.1, -10000,10000); q++;
	aMinuit->SetLimitedVariable(q, Form("a%i",q), 1.02262e+01, 0.1, 0,1000000); q++;
	aMinuit->SetLimitedVariable(q, Form("a%i",q), 0.0, 0.1, -10000,10000); q++;
	aMinuit->SetLimitedVariable(q, Form("a%i",q), 0.0, 0.1, -10000,10000); q++;
	aMinuit->SetLimitedVariable(q, Form("a%i",q), 0.0, 0.1, -10000,10000); q++;

  
  
    aMinuit->SetLimitedVariable(q, Form("a%i",q),  2.48333e+04  , 10., 0.0,10000000); q++;//zentralparameter
	aMinuit->SetLimitedVariable(q, Form("a%i",q), -5.36697e+02, 0.1, -100000,10000000); q++;//aMasse
	aMinuit->SetLimitedVariable(q, Form("a%i",q),0, 0.1, -100000,10000000); q++;//bJSf
	aMinuit->SetLimitedVariable(q, Form("a%i",q), 0, 0.1, -1000,10000000); q++;//JSF
	aMinuit->SetLimitedVariable(q, Form("a%i",q), 8.10386e+01   , 10, 40,100); q++;
	aMinuit->SetLimitedVariable(q, Form("a%i",q),   2.01724e-02 , 0.1, 0,1000000); q++;
	aMinuit->SetLimitedVariable(q, Form("a%i",q), 0, 0.1, 0,1000000); q++;
	aMinuit->SetLimitedVariable(q, Form("a%i",q), 0.0, 0.1, -10000,10000); q++;
	aMinuit->SetLimitedVariable(q, Form("a%i",q), 6.79767e+00 , 0.1, 1.0,15); q++;
	aMinuit->SetLimitedVariable(q, Form("a%i",q), 0.0, 0.1, -10000,10000); q++;
	aMinuit->SetLimitedVariable(q, Form("a%i",q), 0.0, 0.1, -10000,10000); q++;
	aMinuit->SetLimitedVariable(q, Form("a%i",q), 0.0, 0.1, -10000,10000); q++;
	
	aMinuit->SetLimitedVariable(q, Form("a%i",q), 1.77091e+04, 1.0, 100.0,10000000); q++;//zentralparameter
	aMinuit->SetLimitedVariable(q, Form("a%i",q), -3.28105e+02 , 0.1, -1000,1000); q++;//aMasse
	aMinuit->SetLimitedVariable(q, Form("a%i",q),0, 0.1, -100000,100000); q++;//bJSf
	aMinuit->SetLimitedVariable(q, Form("a%i",q), 0, 0.1, -100000,100000); q++;//JSF
	aMinuit->SetLimitedVariable(q, Form("a%i",q), 8.39734e+01, 5, 70,120); q++;
	aMinuit->SetLimitedVariable(q, Form("a%i",q), -2.22161e-02, 0.1, 0,1000000); q++;
	aMinuit->SetLimitedVariable(q, Form("a%i",q), 0., 0.1, 0,1000000); q++;
	aMinuit->SetLimitedVariable(q, Form("a%i",q), 0.0, 0.1, -10000,10000); q++;
	aMinuit->SetLimitedVariable(q, Form("a%i",q), 1.70114e+01, 0.01, 10.0,50); q++;
	aMinuit->SetLimitedVariable(q, Form("a%i",q), -4.48616e-02, 0.1, 0,1000000); q++;
	aMinuit->SetLimitedVariable(q, Form("a%i",q), 0.0, 0.1, -10000,10000); q++;
	aMinuit->SetLimitedVariable(q, Form("a%i",q), 0.0, 0.1, -10000,10000); q++;


  
  	aMinuit->SetLimitedVariable(q, Form("a%i",q), 4.51153e+03, 1000.0, 1000,1000000); q++;//zentralparameter
	aMinuit->SetLimitedVariable(q, Form("a%i",q), 0.0, 0.01, -1000000,10000000); q++;//tompmasse
	aMinuit->SetLimitedVariable(q, Form("a%i",q), 0.0, 0.01, -1000000,10000000); q++;//bJSf
	aMinuit->SetLimitedVariable(q, Form("a%i",q), 0.0, 0.01, -1000000,10000000); q++;//JSF
	aMinuit->SetLimitedVariable(q, Form("a%i",q), 1.92480e+00, 10, 0,1000.0); q++;
	aMinuit->SetLimitedVariable(q, Form("a%i",q), 0.0, 0.1, -10000,10000); q++;
	aMinuit->SetLimitedVariable(q, Form("a%i",q), 0.0, 0.1, -10000,10000); q++;
	aMinuit->SetLimitedVariable(q, Form("a%i",q), 0.0, 0.1, -10000,10000); q++;
	aMinuit->SetLimitedVariable(q, Form("a%i",q), 9.04483e-01, 0.1, 0,1000000); q++;
	aMinuit->SetLimitedVariable(q, Form("a%i",q), 0.0, 0.1, -10000,10000); q++;
	aMinuit->SetLimitedVariable(q, Form("a%i",q), 0.0, 0.1, -10000,10000); q++;
	aMinuit->SetLimitedVariable(q, Form("a%i",q), 0.0, 0.1, -10000,10000); q++;
	aMinuit->SetLimitedVariable(q, Form("a%i",q), 6.16592e+03 , 0.1, 1000,10000000); q++;
	aMinuit->SetLimitedVariable(q, Form("a%i",q), 0.0, 0.1, -10000,10000); q++;
	aMinuit->SetLimitedVariable(q, Form("a%i",q), 0.0, 0.1, -10000,10000); q++;
	aMinuit->SetLimitedVariable(q, Form("a%i",q), 0.0, 0.1, -10000,10000); q++;
	aMinuit->SetLimitedVariable(q, Form("a%i",q), 1.17703e+00, 0.1, 0,1000000); q++;
	aMinuit->SetLimitedVariable(q, Form("a%i",q), 0.0, 0.1, -10000,10000); q++;
	aMinuit->SetLimitedVariable(q, Form("a%i",q), 0.0, 0.1, -10000,10000); q++;
	aMinuit->SetLimitedVariable(q, Form("a%i",q), 0.0, 0.1, -10000,10000); q++;
	aMinuit->SetLimitedVariable(q, Form("a%i",q), 5.52227e-01, 0.1, 0,1000000); q++;
	aMinuit->SetLimitedVariable(q, Form("a%i",q), 0.0, 0.1, -10000,10000); q++;
	aMinuit->SetLimitedVariable(q, Form("a%i",q), 0.0, 0.1, -10000,10000); q++;
	aMinuit->SetLimitedVariable(q, Form("a%i",q), 0.0, 0.1, -10000,10000); q++;
	aMinuit->SetLimitedVariable(q, Form("a%i",q), 3.10282e+04, 0.1, 1000,1000000.0); q++;
	aMinuit->SetLimitedVariable(q, Form("a%i",q), 0.0, 0.1, -10000,10000); q++;
	aMinuit->SetLimitedVariable(q, Form("a%i",q), 0.0, 0.1, -10000,10000); q++;
	aMinuit->SetLimitedVariable(q, Form("a%i",q), 0.0, 0.1, -10000,10000); q++;
	aMinuit->SetLimitedVariable(q, Form("a%i",q), 7.36829e-01, 0.1, 0,1000000); q++;
	aMinuit->SetLimitedVariable(q, Form("a%i",q), 0.0, 0.1, -10000,10000); q++;
	aMinuit->SetLimitedVariable(q, Form("a%i",q), 0.0, 0.1, -10000,10000); q++;
	aMinuit->SetLimitedVariable(q, Form("a%i",q), 0.0, 0.1, -10000,10000); q++;
	aMinuit->SetLimitedVariable(q, Form("a%i",q), 2.44623e-01, 0.1, 0,1000000); q++;
	aMinuit->SetLimitedVariable(q, Form("a%i",q), 0.0, 0.1, -10000,10000); q++;
	aMinuit->SetLimitedVariable(q, Form("a%i",q), 0.0, 0.1, -10000,10000); q++;
	aMinuit->SetLimitedVariable(q, Form("a%i",q), 0.0, 0.1, -10000,10000); q++;

  
  
  
  
  	aMinuit->Minimize();
	
  /*
  	aMinuit->mnexcm("MIGRAD", arglistr ,2);
	aMinuit->mnexcm("SIMPLEX", arglistr ,2);
	aMinuit->mnexcm("MIGRAD", arglistr ,2);
	aMinuit->mnexcm("SIMPLEX", arglistr ,2);
	aMinuit->mnexcm("MIGRAD", arglistr ,2);
	aMinuit->mnexcm("SIMPLEX", arglistr ,2);
	aMinuit->mnexcm("MIGRAD", arglistr ,2);
	aMinuit->mnrset(0);
	aMinuit->mnexcm("SIMPLEX", arglistr ,2);
	aMinuit->mnexcm("MIGRAD", arglistr ,2);
	aMinuit->mnrset(0);
	aMinuit->mnexcm("SIMPLEX", arglistr ,2);
	aMinuit->mnexcm("MIGRAD", arglistr ,2);
	aMinuit->mnrset(0);
	aMinuit->mnexcm("SIMPLEX", arglistr ,2);
	aMinuit->mnexcm("MIGRAD", arglistr ,2);

	aMinuit->mnexcm("HESSE", arglistr ,2);
*/




	
	
}

/////////////////////////////////////////////////////////////////////////////////////////////////////



Double_t TemplateHolder ::gauss(Double_t *x, Double_t *par) {return par[0]*exp(-0.5*pow((x[0]-par[1])/par[2],2));}
//Landau
Double_t TemplateHolder ::landau(Double_t *x,  Double_t *par) { return par[0]*(TMath::Landau(x[0], par[1], par[2])/par[2]);}
//Landau_n
Double_t TemplateHolder ::landau_n(Double_t *x,  Double_t *par) {return par[0]*(TMath::Landau(2*par[1]-x[0], par[1], par[2])/par[2]);}
//top mass fit function
Double_t TemplateHolder ::mtop_function(Double_t *x,  Double_t *par) { return gauss(x, par) + landau(x, &par[3]) + landau_n(x, &par[6]);}
//fit W masse
Double_t TemplateHolder ::fit_mw(Double_t *x,  Double_t *par) { 
return gauss(x,par) + gauss(x,&par[3]);
}
//fit Rbq
Double_t TemplateHolder ::fit_rbq(Double_t *x,  Double_t *par) { 
return gauss(x,par) + gauss(x,&par[3]) + landau(x,&par[6]);
}

TemplateHolder :: TemplateHolder (const char *File, double t, double j, double b){
	
	fmtop=t;
	fbJSF=b;
	fJSF=j;
    f00 = new TF1("f5",this,&TemplateHolder::gauss,130,210,3);
    f01 = new TF1("f1",this,&TemplateHolder::landau,130,210,3);
    f02 = new TF1("f2",this,&TemplateHolder::landau_n,130,210,3);
    ftop = new TF1("ftop",this,&TemplateHolder::mtop_function,130,210,9);
    f13 = new TF1("f3",this,&TemplateHolder::gauss,56,110,3);
    f14 = new TF1("f4",this,&TemplateHolder::gauss,56,110,3);
	fw = new TF1("fmw",this,&TemplateHolder::fit_mw,56,110,6); 
	f25 = new TF1("f5",this,&TemplateHolder::gauss,0.3,3,3);
	f26 = new TF1("f6",this,&TemplateHolder::gauss,0.3,3,3);
	f27 = new TF1("f7",this,&TemplateHolder::landau,0.3,3,3);
    frbq = new TF1("frbq",this,&TemplateHolder::fit_rbq,0.3,3,9);
	// Open file and histogram 
    puts(File);  
    file = new TFile(File,"r");
    h1 = (TH1F*)file->Get("hist_klf_mtop_window");
    h2 = (TH1F*)file->Get("hist_klf_window_Whad_m");
    h3 = (TH1F*)file->Get("hist_klf_window_Rbq_reco");
}

void TemplateHolder::plot(std::string out){
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
  
  //Legend
  TLegend *leg1 = new TLegend(0.7078059,0.6690754,0.8987342,0.8986966,NULL,"brNDC");
  leg1->SetBorderSize(1);
  leg1->SetLineColor(1);
  leg1->SetLineStyle(1);
  leg1->SetLineWidth(1);
  leg1->SetFillColor(0);
  leg1->SetFillStyle(1001);
  leg1->AddEntry(h1,"Simulation","lep");
  leg1->AddEntry(ftop,"Fit","l");
  leg1->AddEntry(f00,"Gauss","l");
  leg1->AddEntry(f01,"Landau","l");
  leg1->AddEntry(f02,"Landau_n","l");
  leg1->SetTextFont(63);
  leg1->SetTextSize(25);
  leg1->Draw();
  
  c1->cd();
  //comparison plot between fit and MC
  TPad *pad2 = new TPad("pad2", "pad2", 0, 0.08, 1, 0.39);
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.2);
  pad2->SetGridx(); // vertical grid
  pad2->SetGridy();
  pad2->Draw();
  pad2->cd();
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
  c1->SaveAs((out+"mtop.png").c_str());
  c1->SaveAs((out+"mtop.root").c_str());
  std::string Name = file->GetName();
  
  //plot Mw 
	TCanvas* c2 = new TCanvas("c2","c2", 950, 900);
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);
	double x = h2->GetMaximum();
    double limit = x + x*0.1;
	TPad *pad3 = new TPad("pad1", "pad1", 0, 0.32, 1, 1.0);
	pad3->SetBottomMargin(2); // Upper and lower plot are joined
    pad3->SetGridx();         // Vertical grid
    pad3->Draw();             // Draw the upper pad: pad1
    pad3->cd();               // pad1 becomes the current pad
	h2->GetXaxis()->SetLabelFont(63);
	h2->GetXaxis()->SetLabelSize(0); // labels will be 14 pixels
	h2->GetXaxis()->SetTitleSize(0.0);
    h2->GetXaxis()->SetTitleOffset(1.11);
    h2->GetXaxis()->SetTitleFont(26);
    h2->GetXaxis()->SetTitle("m_{w}^{reco}[GeV]"); // labels will be 14 pixels
	h2->GetYaxis()->SetLabelFont(63);
	h2->GetYaxis()->SetLabelSize(26);
	h2->GetYaxis()->SetRangeUser(-120,limit);
	h2->GetYaxis()->SetTitleOffset(1.40);
	h2->GetYaxis()->SetTitleFont(63);
	h2->GetYaxis()->SetTitleSize(26);
	h2->GetYaxis()->SetTitle("Events");
	fw->SetFillColor(19);
	fw->SetFillStyle(0);
	fw->SetLineColor(2);
    fw->SetLineWidth(3);   
    h2->SetLineWidth(3);
	h2->Draw();
	fw->Draw("SAME");
	//Get fit parameters 
	double pw[6];
	int npar_mw = 6;
	for(int i = 0; i < 6; i++){
		pw[i] = fw->GetParameter(i);
	double par_mw[100];
	double err_mw[100];
	 par_mw[i] = pw[i];
	 err_mw[i] = fw->GetParError(i);
		};
//Draw functions
	f13->SetParameters(pw[0],pw[1],pw[2]);
	f13->SetLineColor(6);
	f13->SetFillStyle(0); 
	f13->SetLineWidth(3);
	f13->Draw("SAME");
	f14->SetParameters(pw[3],pw[4],pw[5]);
	f14->SetLineColor(8);
	f14->SetFillStyle(0); 
	f14->SetLineWidth(3);
	f14->Draw("SAME");
/////////////////////////////////////////////////////////////////////	
   TLegend *leg = new TLegend(0.7078059,0.6690754,0.8987342,0.8986966,NULL,"brNDC");
   leg->SetBorderSize(1);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(1001);
   leg->AddEntry(h2,"Data","lep");
   leg->AddEntry(fw,"Fit","l");
   leg->AddEntry(f13,"Gauss","l");
   leg->AddEntry(f14,"Gauss","l");
   leg->Draw();
   c2->cd(); 
   TPad *pad21 = new TPad("pad2", "pad2",0, 0.08, 1, 0.39);
   pad21->SetTopMargin(0);
   pad21->SetBottomMargin(0.2);
   pad21->SetGridx(); // vertical grid
   pad21->SetGridy(); // vertical grid
   pad21->Draw();
   pad21->cd();   
   double bin_tot = h2 -> GetSize();  
   double start = h2->GetXaxis()->GetBinCenter(0);
   double end = h2->GetXaxis()->GetBinCenter(bin_tot);
   double div;
   double err;  
   TH1F *hw = new TH1F("hw", "dif", bin_tot, 56, 115);
   hw->GetXaxis()->SetLabelFont(63);
   hw->GetXaxis()->SetLabelSize(26); // labels will be 14 pixels
   hw->GetYaxis()->SetLabelFont(63);
   hw->GetYaxis()->SetLabelSize(26);
   for (int bin=2; bin<=bin_tot;bin++) {
		double x = h2->GetXaxis()->GetBinCenter(bin);
		double fval = fw->Eval(x);
		double dif = h2->GetBinContent(bin);
		double sub = h2->GetBinContent(bin)-fval;
    
		if(sub > -100 && x > 56){
			hw->SetBinContent(bin, sub);
			err = h2->GetBinError(bin);
			div = sub/err;
			hw->SetBinContent(bin, div);
			}else continue;
		};
  	int nbins = hw -> GetNbinsX();
	float lower_edge1  = hw -> GetBinLowEdge(1);
	float bin_width1  = hw -> GetBinWidth(1);
	float number_bins1 = hw -> GetNbinsX();
	float upper_edge1 = hw -> GetBinLowEdge(number_bins) + hw->GetBinWidth(number_bins);      
	hw->GetXaxis()->SetRangeUser(56,115);
    hw->GetXaxis()->SetLabelFont(63);
    hw->GetXaxis()->SetLabelSize(26);
    hw->GetXaxis()->SetTitle("m_{w}^{reco}[GeV]");  
    hw->GetXaxis()->SetTitleSize(26);
    hw->GetXaxis()->SetTitleOffset(1.0);
    hw->GetXaxis()->SetTitleFont(63);
    hw->GetYaxis()->SetLabelFont(63);
    hw->GetYaxis()->SetLabelSize(26);
    hw->GetYaxis()->SetTitleSize(26);
    hw->GetYaxis()->SetTitleFont(63);
    hw->GetYaxis()->SetRangeUser(-3.4,3.4);	      
    hw->GetYaxis()->SetTitleOffset(1.4);
	hw->GetYaxis()->SetTitle("(Sim.-Fit)/Error[#sigma]");
    hw->SetLineWidth(3);  
    hw->SetMarkerStyle(15); 
	hw->Draw("P");
	TF1 *norm11 = new TF1("fa1","0", lower_edge, upper_edge);
	norm11-> SetLineColor(kRed);
	norm11-> SetLineStyle(1);
	norm11 -> SetLineWidth(3);
	norm11->Draw("SAME");
	c2->SaveAs((out+"mw.png").c_str());
	
//rbq

	TCanvas *c3 = new TCanvas("c3","c3", 950, 900);
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);
	TPad *padr = new TPad("pad1", "pad1", 0, 0.32, 1, 1.0);
	padr->SetBottomMargin(1); // Upper and lower plot are joined
    padr->SetGridx();         // Vertical grid
    padr->Draw();             // Draw the upper pad: pad1
    padr->cd();               // pad1 becomes the current pad
	h3->GetXaxis()->SetLabelFont(63);
	h3->GetXaxis()->SetLabelSize(0); // labels will be 14 pixels
	h3->GetXaxis()->SetTitleSize(0.0);
    h3->GetXaxis()->SetTitleOffset(1.11);
    h3->GetXaxis()->SetTitleFont(63);
	h3->GetXaxis()->SetTitle("R_{bq}^{reco}"); // labels will be 14 pixels
	h3->GetYaxis()->SetLabelFont(63);
	h3->GetYaxis()->SetLabelSize(26);
	double xr = h3->GetMaximum();
    double limitr = xr + xr*0.1;
	h3->GetYaxis()->SetRangeUser(-120,limit);
	h3->GetYaxis()->SetTitleFont(63);
	h3->GetYaxis()->SetTitleSize(26);
	h3->GetYaxis()->SetTitleOffset(1.40);
	h3->GetYaxis()->SetTitle("Events");
	//Fit of the histogram
	frbq->SetFillColor(19);
	frbq->SetFillStyle(0);
	frbq->SetLineColor(2);
    frbq->SetLineWidth(3);   
    h3->SetLineWidth(3);
	h3->GetXaxis()->SetRangeUser(0.1,4);  
	h3->Draw();
	//Get fit parameters 
	double pr[9];
	double npar_rbq = 9;
	double par_rbq[100];
	double err_rbq[100];
	for(int i = 0; i < 9; i++){
	 	p[i] = frbq->GetParameter(i);
		par_rbq[i] = pr[i];
		err_rbq[i] = frbq->GetParError(i);
		};
	//Draw functions
	f25->SetParameters(p[0],p[1],p[2]);
	f25->SetLineColor(209);
	f25->SetLineWidth(3);
	f25->Draw("SAME");
	f26->SetParameters(p[3],p[4],p[5]);
	f26->SetLineColor(6);
	f26->SetLineWidth(3);
	f26->Draw("SAME");
	f27->SetParameters(p[6],p[7],p[8]);
	f27->SetLineColor(7);
	f27->SetLineWidth(3);
	f27->Draw("SAME");
	frbq->Draw("SAME");
   TLegend *legr = new TLegend(0.7078059,0.6690754,0.8987342,0.8986966,NULL,"brNDC");
   legr->SetBorderSize(1);
   legr->SetLineColor(1);
   legr->SetLineStyle(1);
   legr->SetLineWidth(1);
   legr->SetFillColor(0);
   legr->SetFillStyle(1001);
   legr->AddEntry(h3,"Data","lep");
   legr->AddEntry(frbq,"Fit","l");
   legr->AddEntry(f25,"Gauss","l");
   legr->AddEntry(f26,"Gauss","l");
   legr->AddEntry(f27,"Landau","l");
   legr->Draw();
   c3->cd(); 
   TPad *padr2 = new TPad("pad2", "pad2", 0, 0.08, 1, 0.39);
   padr2->SetTopMargin(0);
   padr2->SetBottomMargin(0.2);
   padr2->SetGridx(); // vertical grid
   padr2->SetGridy(); // vertical grid
   padr2->Draw();
   padr2->cd();      
   double bin_totr = h3 -> GetSize();  
   double startr = h3->GetXaxis()->GetBinCenter(0);
   double endr = h3->GetXaxis()->GetBinCenter(bin_tot);
   double errr;
   double divr;
   TH1F *hrbq = new TH1F("hrbq", "dif", bin_tot, 0, 3);
   for (int bin=2; bin<=bin_tot;bin++) {
		double x = h3->GetXaxis()->GetBinCenter(bin);
		double fval = frbq->Eval(x);
		double dif = h3->GetBinContent(bin);
		double sub = h3->GetBinContent(bin)-fval;
        if(sub > -100 && x > 0.3){
			hrbq->SetBinContent(bin, sub);
			err = h3->GetBinError(bin);
			div = sub/err;
			hrbq->SetBinContent(bin, div);
			}else continue;
	   };
	int nbinsr = hrbq -> GetNbinsX();
	float lower_edger  = hrbq -> GetBinLowEdge(1);
	float bin_widthr   = hrbq -> GetBinWidth(1);
	float number_binsr = hrbq -> GetNbinsX();
	float upper_edger = hrbq -> GetBinLowEdge(number_bins) + hrbq->GetBinWidth(number_bins);
    hrbq->GetXaxis()->SetRangeUser(0.3,4);
    hrbq->GetXaxis()->SetLabelFont(63);
    hrbq->GetXaxis()->SetLabelSize(26);
    hrbq->GetXaxis()->SetTitle("R_{bq}^{reco}");
    hrbq->GetXaxis()->SetTitleSize(26);
    hrbq->GetXaxis()->SetTitleOffset(3.2);
    hrbq->GetXaxis()->SetTitleFont(63);
    hrbq->GetYaxis()->SetLabelFont(63);
    hrbq->GetYaxis()->SetLabelSize(26);
    hrbq->GetYaxis()->SetTitleSize(26);
    hrbq->GetYaxis()->SetTitleFont(63);
    hrbq->GetYaxis()->SetRangeUser(-2.4,2.4);	
    hrbq->GetYaxis()->SetTitleSize(26);
    hrbq->GetYaxis()->SetTitleOffset(1.4);
    hrbq->GetYaxis()->SetTitleFont(63);
    hrbq->GetYaxis()->SetTitle("(Sim.-Fit)/Error[#sigma]");
    hrbq->SetLineWidth(3);   
    hrbq->SetMarkerStyle(15); 
    hrbq->Draw("P");
    TF1 *normr = new TF1("fa1","0", lower_edge, upper_edge);
	normr -> SetLineColor(kRed);
	normr -> SetLineStyle(1);
	normr -> SetLineWidth(3);
	normr->Draw("SAME");
    c3->SaveAs((out+"rbq.png").c_str());	
 }
///////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////



int main(int argc, char* argv[])
{
	
	// uncomment me on 6 ROOT::EnableThreadSafety();
if (argc>2) exit(1);	
std::vector<double> gMTOP={170,171.5,172.5,173.5,175};
std::vector<double> gbJSF={0.96, 0.98, 1.0, 1.02, 1.04};
std::vector<double> gJSF={0.96, 0.98, 1.0, 1.02, 1.04};
Fill3D* topquark = new Fill3D(argv[1],gMTOP,gbJSF,gJSF);
topquark->top_fit();

for(int i = 0; i < topquark->gALL.size(); i++){
   topquark->gALL[i]->plot(Form("./out_mass_%2.2f_JSF_%2.2f_bJSF_%3.1f",topquark->gALL[i]->fmtop,topquark->gALL[i]->fJSF,topquark->gALL[i]->fbJSF));
};
return 0;
}
