#include "Fill3D.h"
#include "FCN.h"
#include "TemplateHolder.h"

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



FCN::FCN(std::vector<TemplateHolder*> susi){

gALL = susi;
 	
}


double FCN::big_ftop_fcn(const Double_t *par) const
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



double FCN::big_fw_fcn(const Double_t *par) const
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



double FCN::big_frbq_fcn(const Double_t *par) const
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



double FCN::big_all_fcn(const Double_t *par) const
{
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






double FCN::big_ftop_fcnbg(const Double_t *par) const
{
double f;	
double chisq[200];
double chisqG=0;
#pragma omp parallel for 
for (int j=0;j< gbackground.size();j++)
{
	double mtoptemplate=gbackground[j]->ftopbg;
//if (j==3) mtoptemplate=par[36];
	chisq[j]=0;
	double localpar[9];
    for (int k=0;k<9;k++) localpar[k]=par[4*k+0]+par[4*k+1]*(mtoptemplate-172.5)+par[4*k+2]*(gALL[j]->fJSF-1.0)+par[4*k+3]*(gALL[j]->fbJSF-1.0);
    gALL[j]->ftop->SetParameters(localpar);
    for (int i=1;i<gbackground[j]->htopt->GetNbinsX()+1;i++)
	{
		if (gbackground[j]->htopt->GetBinLowEdge(i)<130) continue;
		if (gbackground[j]->htopt->GetBinLowEdge(i)+gbackground[j]->htopt->GetBinWidth(i)>210) continue;
		//double  v=std::abs(gALL[j]->ftop->Integral(gALL[j]->h1->GetBinLowEdge(i),gALL[j]->h1->GetBinLowEdge(i)+gALL[j]->h1->GetBinWidth(i))/(gALL[j]->h1->GetBinWidth(i)));
		double v=gbackground[j]->ftopbg->Eval(gbackground[j]->htopt->GetBinCenter(i));
		if (gbackground[j]->htopt->GetBinContent(i)<0.001) continue;
		if (gbackground[j]->htopt->GetBinError(i)<0.000001) continue;
		chisq[j]+=pow((v-gparameter[j]->htopt->GetBinContent(i))/gparameter[j]->htopt->GetBinError(i),2);
   }
}
  
for (int k=0;k< gbackground.size();k++) chisqG+=chisq[k];
   f = chisqG;
return f;
} 









double FCN::big_all_fcnbg(const Double_t *par) const
{
	double f;
	Int_t iflagq=0;

	int ftop_npar=24;
	double chi2_topbg= big_ftop_fcnbg(par);
	int fw_npar=24;
	double chi2_w=big_fw_fcnbg(&par[24]);
	int frbq_npar=24;
	double chi2_rbqbg=big_frbq_fcnbg(&par[24+24]);


	f=chi2_rbqbg+chi2_topbg+chi2_wbg;
	return f;
}	
