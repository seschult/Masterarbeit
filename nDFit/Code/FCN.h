#ifndef FCN_H_
#define FCN_H_


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

#include "TemplateHolder.h"
using namespace std;


class FCN{
	

public:
std::vector<TemplateHolder*> gALL;

FCN(std::vector<TemplateHolder*> susi);
//Fit functions
//Gauss

double big_ftop_fcn(const Double_t *par) const;
 
double big_fw_fcn(const Double_t *par) const;
  
double big_frbq_fcn(const Double_t *par) const;

//void big_all_fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag){
double big_all_fcn(const Double_t *par) const;


};


#endif
