#ifndef TemplateHolder_H_
#define TemplateHolder_H_

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


class TemplateHolder{


public:
void plot(std::string out);
void save(std::string output);
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

double par_mrbq[9];
double err_mrbq[9];
double par_mw[6];
double err_mw[6];
double err_mtop[9];
double par_mtop[9];


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


#endif
