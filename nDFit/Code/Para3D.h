#ifndef Para3D_H_
#define Para3D_H_

#include "TF1.h"
#include "TH1F.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TFile.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TChain.h"
#include "TGraphErrors.h"
#include "TSystemDirectory.h"
#include "TSystemFile.h"
#include "TSystem.h"
#include "TList.h"

#include "TMinuit.h"
#include "Para3D.h"
#include <string>






using namespace std;




class TemplateHolder{
public: 	
TemplateHolder(); //constructor	
void top_fit();	
void mtop_plot(); //constructor
void input (std::string InputFolder);
TFile *file;//
TH1F *h1; //
TF1 *ftop;//
TF1 *f00;//
TF1 *f01;//
TF1 *f02;//	
std::vector<TemplateHolder*> gALL;
std::vector<std::string> samples; 	
std::string out;





	
};



#endif
