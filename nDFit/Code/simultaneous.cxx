//#include "PlotFactory_13TeV/TemplatePara.h"
//#include "PlotFactory_13TeV/AtlasStyle.h"
//#include "PlotFactory_13TeV/ConfigClass.h"

#include "TSystem.h"
#include "TFile.h"
#include "TF1.h"
#include "TTree.h"
#include "TH1.h"
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
#include <iomanip>
#include <fstream>
#include <math.h>
#include <stdio.h>


#ifndef TemplatePara_H_
#define TemplatePara_H_

#include "TF1.h"
#include "TH1F.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TMinuit.h"

#include <string>

using namespace std;

class TemplatePara{
		    //class TemplatePara: public TObject	
	//members of class
	//{

	public: 
	
        //constructor
	TemplatePara(); //constructor
	TemplatePara(const char *File, const string plot_dir, const string AnalysisType); //constructor
	
	////////////////////////////////////////////////////////////////////
	//For Save
	double par_mtop[20];
	double err_mtop[20];
	int npar_mtop;
	
	double par_mw[20];
	double err_mw[20];
	int npar_mw;
	
	double par_rbq[20];
	double err_rbq[20];
	int npar_rbq;

	double par_mlb[20];
	double err_mlb[20];
	int npar_mlb;
	
	double par_rbq_dilep[20];
	double err_rbq_dilep[20];
	int npar_rbq_dilep;

	double mass;
	////////////////////////////////////////////////////////////////////
	std::string ma;
	
	void mtop_fit(); //constructor
	void mw_fit(); //constructor
	void rbq_fit(); //constructor
	void mlb_fit(); //constructor
	void rbq_dilep_fit(); //constructor
	
	
	string out_dir;//

	TFile *file;//
	
	TH1F *h1; //
	TH1F *h2;//
	TH1F *h3;//
	TH1F *h_mlb;//
	TH1F *h_rbq_dilep;//

	TF1 *ftop;//
	TF1 *fmw;//
	TF1 *frbq;//
	TF1 *f_mlb;//
	TF1 *f_rbq_dilep;//

	TF1 *f00;//
	TF1 *f01;//
	TF1 *f02;//
	TF1 *f13;//
	TF1 *f14;//
	TF1 *f25;//
	TF1 *f26;//
	TF1 *f27;//
	TF1 *f1_mlb;//
	TF1 *f2_mlb;//
	TF1 *f3_mlb;//
	TF1 *f1_rbq_dilep;//
	TF1 *f2_rbq_dilep;//
	TF1 *f3_rbq_dilep;//
		
	TCanvas *c1;//
	TCanvas *c2;//
	TCanvas *c3;//
	TCanvas *c_mlb;//
	TCanvas *c_rbq_dilep;//
	
	
	
	//ClassDef(TemplatePara,0)
	
	
};

#endif


using namespace std;

////////////////////////////////////////////////////////////////////////
//Fit Functions
////////////////////////////////////////////////////////////////////////

//Gauss
Double_t gauss(Double_t *x, Double_t *par) {
  return par[0]*exp(-0.5*pow((x[0]-par[1])/par[2],2));
}
//Landau
Double_t landau(Double_t *x,  Double_t *par) { // 3 par
  return par[0]*(TMath::Landau(x[0], par[1], par[2])/par[2]);
}
//Landau_n
Double_t landau_n(Double_t *x,  Double_t *par) { // 3 par
  return par[0]*(TMath::Landau(2*par[1]-x[0], par[1], par[2])/par[2]);
}


//top mass fit function
Double_t mtop_function(Double_t *x,  Double_t *par) { 
return gauss(x, par) + landau(x, &par[3]) + landau_n(x, &par[6]);
}
//W mass fit function
Double_t mw_function(Double_t *x,  Double_t *par) { 
return gauss(x,par) + gauss(x,&par[3]);
}
//Rbq fit function
Double_t rbq_function(Double_t *x,  Double_t *par) { 
return gauss(x,par) + gauss(x,&par[3]) + landau(x,&par[6]);
}
//mlb fit function
Double_t mlb_function(Double_t* x, Double_t* par) {
  return gauss(x,par) + landau(x,&par[3]);
}
//Dilepton Rbq fit function
Double_t rbq_dilep_function(Double_t* x, Double_t* par) {
  return gauss(x,par) + gauss(x,&par[3]) + landau(x,&par[6]);
}
Double_t linpar_function(Double_t* x, Double_t* par) {
  return par[0]*(x[0] - 172.5) + par[1];
}


////////////////////////////////////////////////////////////////////////
//Fit

TH1F* gTEMP_HISTO;
TF1*  gTEMP_FUNCTION;
void ftop_fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);







////////////////////////////////////////////////////////////////////////
//Ratio plot
////////////////////////////////////////////////////////////////////////

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
  hratio->SetLineColor(f->GetLineColor());
  hratio->SetLineWidth(2);
  return hratio;
}



////////////////////////////////////////////////////////////////////////
//Class def
////////////////////////////////////////////////////////////////////////


//ClassImp(TemplatePara)

////////////////////////////////////////////////////////////////////////
//Default constructor
////////////////////////////////////////////////////////////////////////


TemplatePara :: TemplatePara (){


  out_dir="";

    h1=0; //
    h2=0;//
    h3=0;//
    file=0;//
    
    f00=0;//
    f01=0;//
    f02=0;//
    f13=0;//
    f14=0;//
    f25=0;//
    f26=0;//
    f27=0;//
        
    ftop=0;//
    fmw=0;//
    frbq=0;//
    
    //map<string,TF1*>    fFunctions=0;
    c1=0;//
    c2=0;//
    c3=0;//


    //Dilepton
    h_mlb=0; //
    h_rbq_dilep=0; //
    
    f_mlb=0;//
    f_rbq_dilep=0;//
    
    f1_mlb=0;//
    f2_mlb=0;//
    f3_mlb=0;//
    f1_rbq_dilep=0;//
    f2_rbq_dilep=0;//
    f3_rbq_dilep=0;//
    
    c_mlb=0;//
    c_rbq_dilep=0;//
   
}




/////////////////////////////////////////////////////////////////////////
//Constructor
////////////////////////////////////////////////////////////////////////

TemplatePara :: TemplatePara (const char *File, const string plot_dir, const string AnalysisType){
  string fname(gSystem->BaseName(File));



  if (fname=="Merged_ttbar_170.root"){ mass = 170;
	  ma = "170";
	};
  if (fname=="Merged_ttbar_171p5.root"){ mass = 171.5;
	ma = "171p5";
    };
  if (fname=="Merged_ttbar_173p5.root"){ mass = 173.5;
	ma =  "173p5";
	};
  if (fname=="Merged_ttbar_175.root"){ mass = 175;
	ma = "175";
	};
  if (fname=="Merged_ttbar.root"){ mass = 172.5;
	ma = "172p5";
	};
  
  out_dir = plot_dir;

  if(AnalysisType == "lepjets"){
    //Top 
    f00 = new TF1("f5",gauss,130,210,3);
    f01 = new TF1("f1",landau,130,210,3);
    f02 = new TF1("f2",landau_n,130,210,3);
    ftop = new TF1("ftop",mtop_function,130,210,9);    
    
    //MW 
    f13 = new TF1("f3",gauss,56,105,3);
    f14 = new TF1("f4",gauss,56,105,3);
    fmw = new TF1("fmw",mw_function,56,105,6);    
    
    //RBQ 
    f25 = new TF1("f5",gauss,0.3,3,3);
    f26 = new TF1("f6",gauss,0.3,3,3);
    f27 = new TF1("f7",landau,0.3,3,3);
    frbq = new TF1("frbq",rbq_function,0.3,3,9);
        
    
    // Set initial parameters and limits 
    /* ftop->SetParameters(1318, 167, -1.05, 4340, 182, 1.27, 54200, 142 ,9.56);
    ftop->SetParLimits(3,0.0,1000000);
    ftop->SetParLimits(6,0.0,1000000);
    ftop->SetParLimits(7,100.0,1000000);
    
    fmw->SetParameters(1742, 78, 6.8, 1616, 79, 19); 
    
    frbq->SetParameters(600, 1.5, 0.7, 1000, 1, 0.5, 4525, 0.95 ,0.3);
    frbq->SetParLimits(0,0.0,1000000);
    frbq->SetParLimits(3,0.0,1000000);
    frbq->SetParLimits(6,0.0,1000000);
    
    */   
    // Open file and histogram 
    puts(File);  
    file = new TFile(File,"r");
    
    h1 = (TH1F*)file->Get("hist_klf_mtop_window");
    h2 = (TH1F*)file->Get("hist_klf_window_Whad_m");
    h3 = (TH1F*)file->Get("hist_klf_window_Rbq_reco");
  }
  else if(AnalysisType == "dilepton"){
    //mlb 
    f_mlb = new TF1("f_mlb",mlb_function,30,170,6);    
    f1_mlb = new TF1("f1_mlb",gauss,30,170,3);
    f2_mlb = new TF1("f2_mlb",landau,30,170,3);
    //f3_mlb = new TF1("f3_mlb",landau_n,30,170,3);
    
    //rbq_dilepton
    f_rbq_dilep = new TF1("f_rbq_dilep",rbq_dilep_function,0.3,3,9);
    f1_rbq_dilep = new TF1("f1_rbq_dilep",gauss,0.3,3,3);
    f2_rbq_dilep = new TF1("f2_rbq_dilep",gauss,0.3,3,3);
    f3_rbq_dilep = new TF1("f3_rbq_dilep",landau,0.3,3,3);
        
    
    // Set initial parameters and limits 
    f_mlb->SetParameters(1000,92,10,1000,70,27);//gaus+landau
    //f_mlb->SetParLimits(3,0.0,1000000);
    //f_mlb_func->SetParLimits(0,0.01,0.04);
    //f_mlb_func->SetParLimits(1,70,110);
    //f_mlb_func->SetParLimits(3,0.05,0.6);
    
    f_rbq_dilep->SetParameters(600, 1.5, 0.7, 1000, 1, 0.5, 4525, 0.95 ,0.3);
    f_rbq_dilep->SetParLimits(0,0.0,1000000);
    f_rbq_dilep->SetParLimits(3,0.0,1000000);
    f_rbq_dilep->SetParLimits(6,0.0,1000000);
    
        
    // Open file and histogram 
    puts(File);  
    file = new TFile(File,"r");
    
    h_mlb = (TH1F*)file->Get("hist_mlb_logscale"); 
    h_rbq_dilep = (TH1F*)file->Get("hist_Rbq"); 
  }
  
}





//130,210

void ftop_fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
   
  // TH1F* fH=(TH1F*)(gTEMP_HISTO->Clone("Adddddddd"));
   gTEMP_FUNCTION->SetParameters(par);
   double chisq=0;
   //loop over   all mqsses beginfor (mtop=170 ... 175) {
   
   for (int i=1;i<gTEMP_HISTO->GetNbinsX()+1;i++)
   {
   if (gTEMP_HISTO->GetBinLowEdge(i)<130) continue;
   if (gTEMP_HISTO->GetBinLowEdge(i)+gTEMP_HISTO->GetBinWidth(i)>210) continue;
   
   double v=std::abs(gTEMP_FUNCTION->Integral(gTEMP_HISTO->GetBinLowEdge(i),gTEMP_HISTO->GetBinLowEdge(i)+gTEMP_HISTO->GetBinWidth(i))/(gTEMP_HISTO->GetBinWidth(i)));
   if (gTEMP_HISTO->GetBinContent(i)<0.001) continue;
   if (gTEMP_HISTO->GetBinError(i)<0.000001) continue;
   chisq+=pow((v-gTEMP_HISTO->GetBinContent(i))/gTEMP_HISTO->GetBinError(i),2);
   
   }
   
   //fH->Eval(gTEMP_FUNCTION,"R");
    //chisq=fH->Chi2Test(gTEMP_HISTO,"CHI2");    
   f = chisq;
}

////////////////////////////////////////////////////////////////////////
//Fit top Mass
////////////////////////////////////////////////////////////////////////

void TemplatePara :: mtop_fit(){ 
  
 


  double bin = h1->GetMaximum();

  

  double mean = h1->GetMean();
  double RMS  = h1->GetRMS();
  double skew = h1->GetSkewness();
  double kurtosis = h1->GetKurtosis();
  
  cout<<"mean"<< mean<< endl;
  cout<<"RMS"<< RMS<< endl;
  cout<<"skew"<< skew<< endl;

/*

  ftop->SetParameters(2*bin, mean , RMS/2, bin*0.75, mean + RMS , RMS/3, bin*0.75,130 ,RMS/2);



  ftop->SetParLimits(0,0,10000000);
  ftop->SetParLimits(1,130,180);
  ftop->SetParLimits(2,0,1000*RMS);
  ftop->SetParLimits(3,0,100000000);
  ftop->SetParLimits(4,130,240);
  ftop->SetParLimits(5,0,10000);
  //ftop->SetParLimits(6,0,1000000);
  ftop->SetParLimits(7,90,150);
  ftop->SetParLimits(8,0,100
  );
  */
    ftop->SetParLimits(0,1.0,100000000);
	ftop->SetParLimits(3,0.0,1000000);
	ftop->SetParLimits(6,0.0,1000000);
    ftop->SetParLimits(7,100.0,1000000);

if ( std::abs(mass - 170)<0.001)  
ftop->SetParameters(5.74282e+03, 1.65763e+02, 9.86866e+00, 7.50921e+05, 1.64739e+02, 1.55647e+01, 5.33607e+05, 1.31487e+02, 1.32490e+01);

if ( std::abs(mass - 171.5)<0.001)  
ftop->SetParameters(1.15047e+04, 1.62162e+02, 1.22716e+01, 4.47290e+05, 1.87020e+02, 1.52078e+01, 4.43965e+05, 1.37312e+02, 9.43957e+00);


if ( std::abs(mass - 172.5)<0.001)  
ftop->SetParameters(1.10618e+04, 1.62284e+02, 1.28988e+01, 4.28790e+05, 1.87218e+02, 1.48783e+01, 4.46906e+05, 1.36078e+02, 1.02262e+01);


if ( std::abs(mass - 173.5)<0.001)  
ftop->SetParameters(6.10094e+03, 1.67810e+02, 1.09465e+01, 6.79509e+05, 1.69965e+02, 1.75019e+01, 4.66939e+05, 1.34979e+02, 1.21092e+01);


if ( std::abs(mass - 175)<0.001)  
{
ftop->SetParameters(
              6.12031e+03,
              1.68487e+02,
              9.18321e+00,
              4.64806e+05,
              1.79979e+02,
              1.20146e+01,
              5.44530e+05,
             1.44945e+02,
   9.53774e+00);


}




 cout << mass <<endl;
  //SetAtlasStyle();

  double with = 950;
  double hight = 900;
  c1 = new TCanvas("c1","c1", with, hight);
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
  
  double x = h1->GetMaximum();
  double limit = x + x*0.1;
  h1->GetYaxis()->SetRangeUser(-73,limit);
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


	cout<< mass<< endl;
	cout<< mass<< endl;
	cout<< mass<< endl;
	cout<< mass<< endl;
	cout<< mass<< endl;
	cout<< mass<< endl;
	cout<< mass<< endl;
	cout<< mass<< endl;
	cout<< mass<< endl;
	cout<< mass<< endl;
	cout<< mass<< endl;
	cout<< mass<< endl;
	cout<< mass<< endl;
	cout<< mass<< endl;
	cout<< mass<< endl;
	cout<< "top"<< endl;
	
  h1->Fit("ftop","RI","",130,210);
  h1->Fit("ftop","RI","",130,210);
  h1->Fit("ftop","RI","",130,210);
  h1->Fit("ftop","RI","",130,210);
  h1->Fit("ftop","RI","",130,210);


  
  h1->Draw();
  
	
  
  
  //Get fit parameters 
  ////////////////////////////////////////////////////////////////////
  double p[9];
  
  npar_mtop = 9;
  
  for(int i = 0; i < 9; i++){
    
    p[i] = ftop->GetParameter(i);
    par_mtop[i] = p[i];
    err_mtop[i] = ftop->GetParError(i);
  };
  
    
    
///////////////////////////////////////////////////////////////////////////////////////////////
   TMinuit *gMinuit = new TMinuit(9);  //initialize TMinuit with a maximum of 5 params
   
   gTEMP_HISTO=h1;
   gTEMP_FUNCTION=ftop;
   
   gMinuit->SetFCN(ftop_fcn);
   Double_t arglist[10];
   Int_t ierflg = 0;
   
   arglist[0] = 15000;
   arglist[1] = 0.001;
   
   for (int i=0;i<9;i++)    gMinuit->mnparm(i, Form("a%i",i), par_mtop[i], err_mtop[i]+1.0, 0,0,ierflg);
   gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);
   gMinuit->mnexcm("SIMPLEX", arglist ,2,ierflg);
   gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);
   gMinuit->mnexcm("SIMPLEX", arglist ,2,ierflg);
   gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);
   gMinuit->mnexcm("HESSE", arglist ,2,ierflg);
  /////////////////////////////////////////////////////////////////////////////////////////////
    
   
   for (int i=0;i<9;i++) 
   {
	 double currentValue,  currentError;
   	gMinuit->GetParameter(i, currentValue, currentError);
   ftop->SetParameter(i,currentValue);
}
     
   h1->Fit("ftop","RI","",130,210);


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
  f02->SetLineColor(12);
  f02->SetFillStyle(0); 
  f02->SetLineWidth(3);  
  f02->Draw("SAME");
  
  
  
  //Get chi^2 
  ////////////////////////////////////////////////////////////////////
  double chi2_0 = ftop->GetChisquare();
  double NDF_0 = ftop->GetNDF();
  double Prob = ftop->GetProb();
  double COM = chi2_0/NDF_0;
  
  stringstream oss_Sep;
  oss_Sep   << setprecision(3) << chi2_0;
  TLatex l00;
  l00.SetTextAlign(9);
  l00.SetTextSize(0.048);
  l00.SetLineWidth(2);
  l00.SetNDC();
  //l00.DrawLatex(0.1358774,0.6824242, ("#chi^{2}: " + oss_Sep.str()).c_str());
  
  
  stringstream oss_NDF_0;
  oss_NDF_0 << setprecision(3) << NDF_0;
  
  TLatex l01;
  l01.SetTextAlign(9);
  l01.SetTextSize(0.048);
  l01.SetLineWidth(2);
  l01.SetNDC();
  //l01.DrawLatex(.1358774,0.7211448, ("NDF: " + oss_NDF_0.str()).c_str());
  
  stringstream oss_Prob;
  oss_Prob << setprecision(3) << Prob;
  
  TLatex l03;
  l03.SetTextAlign(9);
  l03.SetTextSize(0.048);
  l03.SetLineWidth(2);
  l03.SetNDC();
  //l03.DrawLatex(0.1369345,0.7699663, ("Prob: " + oss_Prob.str()).c_str());
  
  
  
  stringstream oss_COM;
  oss_COM   << setprecision(3) << COM;
  TLatex l04;
  l04.SetTextAlign(9);
  l04.SetTextSize(0.048);
  l04.SetLineWidth(2);
  l04.SetNDC();
  //l04.DrawLatex(0.1358774,0.6352862, ("#chi^{2}/NDF: " + oss_COM.str()).c_str());
  //l04.DrawLatex(0.1,0.77, ("#chi^{2}/NDF = "+ oss_Sep.str() + "/" +oss_NDF_0.str() + " = " + oss_COM.str()).c_str());
  l04.DrawLatex(0.12,0.8, ("#chi^{2}/NDF = "+ oss_Sep.str() + "/" +oss_NDF_0.str() + " = " + oss_COM.str()).c_str());
  TLatex l05;
  l05.SetTextAlign(9);
  l05.SetTextSize(0.058);
  l05.SetLineWidth(2);
  l05.SetNDC();
  l05.DrawLatex(0.1258774,0.8552862, ("ATLAS work-in-progress"));
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
  leg->SetTextFont(63);
  leg->SetTextSize(25);
  leg->Draw();
  
  
  
  ////////////////////////////////////////////////////////////////////
  c1->cd(); 
  TPad *pad2 = new TPad("pad2", "pad2", 0, 0.08, 1, 0.39);
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.2);
  pad2->SetGridx(); // vertical grid
  pad2->SetGridy();
  pad2->Draw();
  pad2->cd();   
  // pad2 becomes the current pad
  
  
  //Calculate the diffrence between Fitfunction and Histogramm h1
  ////////////////////////////////////////////////////////////////////
  
  double bin_tot = h1 -> GetSize();  
  // double start = h1->GetXaxis()->GetBinCenter(0);
  // double end = h1->GetXaxis()->GetBinCenter(bin_tot);
  double err;
  double div;
  
  
  
  TH1F *htop = new TH1F("htop", "dif", bin_tot, 120, 220);
  
  
  
  
  
  
  
  
  
  for (int bin=2; bin<=bin_tot;bin++) {
    
    double x = h1->GetXaxis()->GetBinCenter(bin);
    double fval = ftop->Eval(x);
    //double dif = h1->GetBinContent(bin);
    double sub = h1->GetBinContent(bin)-fval;
    if(sub > -100 && x > 130 ){
      htop->SetBinContent(bin, sub);
      err = h1->GetBinError(bin);
      div = sub/err;
      htop->SetBinContent(bin, div);
    }else continue;
    
  };
  
  // double x1 = htop->GetMaximum();
  // double limit1 = x1 + x1*0.12;
  htop->GetYaxis()->SetRangeUser(-2.8,2.8);
  
  
  // int nbins = htop -> GetNbinsX();
  
  float lower_edge  = htop -> GetBinLowEdge(1);
  //float bin_width   = htop -> GetBinWidth(1);
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
  
  
  //htop->GetYaxis()->SetRangeUser(-3.4,3.4);	      
  
  htop->GetYaxis()->SetTitleOffset(1.4);
  
  htop->GetYaxis()->SetTitle("(Sim.-Fit)/Error[#sigma]");
  htop->SetLineWidth(3);  
  
  
  
  htop->Draw();
  
  TF1 *norm1 = new TF1("fa1","0", lower_edge, upper_edge);
  norm1 -> SetLineColor(kRed);
  norm1 -> SetLineStyle(1);
  norm1 -> SetLineWidth(3);
  norm1->Draw("SAME");
  
  string Name = file->GetName();
  //string Base = gSystem->BaseName(Name.c_str());
  
  //c1 -> Print((out_dir + Base + "mtop.pdf").c_str());
  c1 -> Print((out_dir + "/" + ma + "_temp_mtop.png").c_str());
  
  
  
  
  
}




////////////////////////////////////////////////////////////////////////
//Fit mw
////////////////////////////////////////////////////////////////////////

void TemplatePara :: mw_fit(){



  double bin = h2->GetMaximum();

  

  double mean = h2->GetMean();
  double RMS  = h2->GetRMS();
  double skew = h2->GetSkewness();
  double kurtosis = h2->GetKurtosis();
  
  cout<<"mean"<< mean<< endl;
  cout<<"RMS"<< RMS<< endl;
  cout<<"skew"<< skew<< endl;



  fmw->SetParameters(2*bin, mean-RMS  , RMS/2, bin*0.75, mean+RMS, RMS*2);
  fmw->SetParLimits(2,0,1000*RMS);
  fmw->SetParLimits(5,0,1000*RMS);
  

  fmw->SetParLimits(4,mean-RMS/1000, 100.0);
  fmw->SetParLimits(1,0,mean-RMS/1000);





























  //SetAtlasStyle();

  double with = 950;
  double hight = 900;
  c2 = new TCanvas("c2","c2", with, hight);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  
  
  
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.32, 1, 1.0);
  pad1->SetBottomMargin(2); // Upper and lower plot are joined
  pad1->SetGridx();         // Vertical grid
  pad1->Draw();             // Draw the upper pad: pad1
  pad1->cd();               // pad1 becomes the current pad
  
  
  h2->GetXaxis()->SetLabelFont(63);
  h2->GetXaxis()->SetLabelSize(0); // labels will be 14 pixels
  h2->GetXaxis()->SetTitleSize(0.0);
  h2->GetXaxis()->SetTitleOffset(1.45);
  h2->GetXaxis()->SetTitleFont(28);
  h2->GetXaxis()->SetTitle("m_{w}^{reco}[GeV]"); // labels will be 14 pixels
  
  h2->GetXaxis()->SetRangeUser(50,115);
  
  h2->GetYaxis()->SetLabelFont(63);
  h2->GetYaxis()->SetLabelSize(28);
  
  double x = h2->GetMaximum();
  double limit = x + x*0.1;
  h2->GetYaxis()->SetRangeUser(-120,limit);
  h2->GetYaxis()->SetTitleOffset(1.45);
  
  h2->GetYaxis()->SetTitleFont(63);
  h2->GetYaxis()->SetTitleSize(28);
  h2->GetYaxis()->SetTitle("Events");
  
  //Fit of the histogram   
  
  
  fmw->SetFillColor(19);
  fmw->SetFillStyle(0);
  fmw->SetLineColor(2);
  fmw->SetLineWidth(3);   
  h2->SetLineWidth(3);
  
  
  

	cout<< mass<< endl;
	cout<< mass<< endl;
	cout<< mass<< endl;
	cout<< mass<< endl;
	cout<< mass<< endl;
	cout<< mass<< endl;
	cout<< mass<< endl;
	cout<< mass<< endl;
	cout<< mass<< endl;
	cout<< mass<< endl;
	cout<< mass<< endl;
	cout<< mass<< endl;
	cout<< mass<< endl;
	cout<< mass<< endl;
	cout<< mass<< endl;
	cout<< "W"<< endl;
  
  h2->Fit("fmw","","",56,105);
  h2->Fit("fmw","","",56,105);
  h2->Fit("fmw","","",56,105);
  h2->Fit("fmw","","",56,105);
  h2->Fit("fmw","","",56,105);




  h2->Draw();
  
  
  //Get fit parameters 
  double p[6];
  npar_mw = 6;
  
  for(int i = 0; i < 6; i++){
    
    p[i] = fmw->GetParameter(i);
    par_mw[i] = p[i];
    err_mw[i] = fmw->GetParError(i);
  };
  
  
  //Draw functions
  f13->SetParameters(p[0],p[1],p[2]);
  f13->SetLineColor(6);
  f13->SetFillStyle(0); 
  f13->SetLineWidth(3);
  f13->Draw("SAME");
  
  f14->SetParameters(p[3],p[4],p[5]);
  f14->SetLineColor(209);
  f14->SetFillStyle(0); 
  f14->SetLineWidth(3);
  f14->Draw("SAME");
  
  
  
  
  //Get chi^2
  ////////////////////////////////////////////////////////////////////////
  
  
  
  double chi2_0 = fmw->GetChisquare();
  double NDF_0 = fmw->GetNDF();
  double Prob = fmw->GetProb();
  double COM = chi2_0/NDF_0;
  
  
  
  stringstream oss_Sep;
  oss_Sep   << setprecision(3) << chi2_0;
  TLatex l00;
  l00.SetTextAlign(9);
  l00.SetTextSize(0.048);
  l00.SetLineWidth(2);
  l00.SetNDC();
  //l00.DrawLatex(0.1358774,0.6824242, ("#chi^{2}: " + oss_Sep.str()).c_str());
  
  
  stringstream oss_NDF_0;
  oss_NDF_0 << setprecision(3) << NDF_0;
  
  TLatex l01;
  l01.SetTextAlign(9);
  l01.SetTextSize(0.048);
  l01.SetLineWidth(2);
  l01.SetNDC();
  //l01.DrawLatex(.1358774,0.7211448, ("NDF: " + oss_NDF_0.str()).c_str());
  
  stringstream oss_Prob;
  oss_Prob << setprecision(3) << Prob;
  
  TLatex l03;
  l03.SetTextAlign(9);
  l03.SetTextSize(0.048);
  l03.SetLineWidth(2);
  l03.SetNDC();
  //l03.DrawLatex(0.1369345,0.7699663, ("Prob: " + oss_Prob.str()).c_str());
  
  
  
  stringstream oss_COM;
  oss_COM   << setprecision(3) << COM;
  TLatex l04;
  l04.SetTextAlign(9);
  l04.SetTextSize(0.048);
  l04.SetLineWidth(2);
  l04.SetNDC();
  //l04.DrawLatex(0.1358774,0.6352862, ("#chi^{2}/NDF: " + oss_COM.str()).c_str());
  l04.DrawLatex(0.12,0.8, ("#chi^{2}/NDF = "+ oss_Sep.str() + "/" +oss_NDF_0.str() + " = " + oss_COM.str()).c_str());
  
  
  TLatex l05;
  l05.SetTextAlign(9);
  l05.SetTextSize(0.058);
  l05.SetLineWidth(2);
  l05.SetNDC();
  l05.DrawLatex(0.1258774,0.8552862, ("ATLAS work-in-progress"));
  ////////////////////////////////////////////////////////////////////
  
  
  /////////////////////////////////////////////////////////////////////	
  TLegend *leg = new TLegend(0.7078059,0.6690754,0.8987342,0.8986966,NULL,"brNDC");
  leg->SetBorderSize(1);
  leg->SetLineColor(1);
  leg->SetLineStyle(1);
  leg->SetLineWidth(1);
  leg->SetFillColor(0);
  leg->SetFillStyle(1001);
  leg->AddEntry(h2,"Simulation","lep");
  leg->AddEntry(fmw,"Fit","l");
  leg->AddEntry(f13,"Gauss","l");
  leg->AddEntry(f14,"Gauss","l");
  leg->SetTextFont(63);
  leg->SetTextSize(25);
  
  leg->Draw();
  
  
  
  
  c2->cd(); 
  TPad *pad2 = new TPad("pad2", "pad2",0, 0.08, 1, 0.39);
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.2);
  pad2->SetGridx(); // vertical grid
  pad2->SetGridy(); // vertical grid
  pad2->Draw();
  pad2->cd();   
  // pad2 becomes the current pad
  
  
  //Calculate the diffrence between Fitfunction and Histogramm h1
  
  
  //Calculate the diffrence between Fitfunction and Histogramm h1
  
  double bin_tot = h2 -> GetSize();  
  //double start = h2->GetXaxis()->GetBinCenter(0);
  //double end = h2->GetXaxis()->GetBinCenter(bin_tot);
  double div;
  double err;
  
  TH1F *hw = new TH1F("hw", "dif", bin_tot, 45, 120);
  hw->GetXaxis()->SetLabelFont(63);
  hw->GetXaxis()->SetLabelSize(26); // labels will be 14 pixels
  hw->GetYaxis()->SetLabelFont(63);
  hw->GetYaxis()->SetLabelSize(26);
  //TF1 *func = h1->GetFunction("ftop");
  
  for (int bin=2; bin<=bin_tot;bin++) {
    
    double x = h2->GetXaxis()->GetBinCenter(bin);
    double fval = fmw->Eval(x);
    // double dif = h2->GetBinContent(bin);
    double sub = h2->GetBinContent(bin)-fval;
    
    if(sub > -100 && x > 54 && x < 105){
      hw->SetBinContent(bin, sub);
      err = h2->GetBinError(bin);
      div = sub/err;
      hw->SetBinContent(bin, div);
    }else continue;
    
  };
  
  
  // int nbins = hw -> GetNbinsX();
  
  float lower_edge  = hw -> GetBinLowEdge(1);
  // float bin_width   = hw -> GetBinWidth(1);
  float number_bins = hw -> GetNbinsX();
  float upper_edge = hw -> GetBinLowEdge(number_bins) + hw->GetBinWidth(number_bins);
  
  hw->GetXaxis()->SetRangeUser(50,115);
  
  hw->GetXaxis()->SetLabelFont(63);
  hw->GetXaxis()->SetLabelSize(28);
  hw->GetXaxis()->SetTitle("m_{w}^{reco}[GeV]");
  hw->GetXaxis()->SetTitleSize(28);
  hw->GetXaxis()->SetTitleOffset(3.0);
  hw->GetXaxis()->SetTitleFont(63);
  hw->GetYaxis()->SetLabelFont(63);
  hw->GetYaxis()->SetLabelSize(28);
  hw->GetYaxis()->SetTitleSize(28);
  hw->GetYaxis()->SetTitleFont(63);
  
  
  hw->GetYaxis()->SetRangeUser(-2.8,2.8);	      
  
  hw->GetYaxis()->SetTitleOffset(1.45);
  
  hw->GetYaxis()->SetTitle("(Sim.-Fit)/Error[#sigma]");
  hw->SetLineWidth(3);  
  
  
  
  hw->Draw();
  
  TF1 *norm1 = new TF1("fa1","0", lower_edge, upper_edge);
  norm1 -> SetLineColor(kRed);
  norm1 -> SetLineStyle(1);
  norm1 -> SetLineWidth(3);
  norm1->Draw("SAME");
  
  
  string Name = file->GetName();
  //string Base = gSystem->BaseName(Name.c_str());
  
  //c2 -> Print((out_dir + Base + "mw.pdf").c_str());
  c2 -> Print((out_dir + "/" + ma + "_temp_mw.png").c_str());
}



////////////////////////////////////////////////////////////////////////
//Fit rbq
////////////////////////////////////////////////////////////////////////
	
void TemplatePara :: rbq_fit(){ 





  double bin = h3->GetMaximum();

  

  double mean = h3->GetMean();
  double RMS  = h3->GetRMS();
  double skew = h3->GetSkewness();
  
  
  cout<<"mean"<< mean<< endl;
  cout<<"RMS"<< RMS<< endl;
  cout<<"skew"<< skew<< endl;



  frbq->SetParameters(2*bin, mean , RMS/2, bin*0.75, mean + RMS , RMS/3, bin*0.75, mean+RMS ,RMS/2);
  
  
  

















 // SetAtlasStyle();

  double with = 950;
  double hight = 900;
  c3 = new TCanvas("c3","c3", with, hight);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.32, 1, 1.0);
  pad1->SetBottomMargin(1); // Upper and lower plot are joined
  pad1->SetGridx();         // Vertical grid
  pad1->Draw();             // Draw the upper pad: pad1
  pad1->cd();               // pad1 becomes the current pad
  
  h3->GetXaxis()->SetLabelFont(63);
  h3->GetXaxis()->SetLabelSize(0); // labels will be 14 pixels
  h3->GetXaxis()->SetTitleSize(0.0);
  h3->GetXaxis()->SetTitleOffset(1.45);
  h3->GetXaxis()->SetTitleFont(63);
  h3->GetXaxis()->SetTitle("R_{bq}^{reco}"); // labels will be 14 pixels
  
  
  
  h3->GetYaxis()->SetLabelFont(63);
  
  h3->GetYaxis()->SetLabelSize(28);
  
  double x = h3->GetMaximum();
  double limit = x + x*0.1;
  h3->GetYaxis()->SetRangeUser(-120,limit);
  
  
  h3->GetYaxis()->SetTitleFont(63);
  h3->GetYaxis()->SetTitleSize(28);
  h3->GetYaxis()->SetTitleOffset(1.45);
  h3->GetYaxis()->SetTitle("Events");
  
  
  //Fit of the histogram
  
  frbq->SetFillColor(19);
  frbq->SetFillStyle(0);
  frbq->SetLineColor(2);
  frbq->SetLineWidth(3);   
  h3->SetLineWidth(3);
  
  h3->GetXaxis()->SetRangeUser(0.1,4);  
  
  //Fit of the histogram   




  h3->Fit("frbq","I","",0.3,3);
  h3->Fit("frbq","I","",0.3,3);
  h3->Fit("frbq","I","",0.3,3);
  h3->Fit("frbq","I","",0.3,3);
  h3->Fit("frbq","I","",0.3,3);






  h3->Draw();
  
  
  
  
  
  //Get fit parameters 
  double p[9];
  
  npar_rbq = 9;
  
  for(int i = 0; i < 9; i++){
    
    p[i] = frbq->GetParameter(i);
    par_rbq[i] = p[i];
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
  f27->SetLineColor(12);
  f27->SetLineWidth(3);
  f27->Draw("SAME");
  
  
  //////////////////////////////////////////////////////////////////////
  double chi2_0 = frbq->GetChisquare();
  double NDF_0 = frbq->GetNDF();
  double Prob = frbq->GetProb();
  double COM = chi2_0/NDF_0;
  
  stringstream oss_Sep;
  oss_Sep   << setprecision(3) << chi2_0;
  TLatex l00;
  l00.SetTextAlign(9);
  l00.SetTextSize(0.048);
  l00.SetLineWidth(2);
  l00.SetNDC();
  //l00.DrawLatex(0.5058774,0.6824242, ("#chi^{2}: " + oss_Sep.str()).c_str());
  
  
  stringstream oss_NDF_0;
  oss_NDF_0 << setprecision(3) << NDF_0;
  
  TLatex l01;
  l01.SetTextAlign(9);
  l01.SetTextSize(0.048);
  l01.SetLineWidth(2);
  l01.SetNDC();
  //l01.DrawLatex(.5058774,0.7211448, ("NDF: " + oss_NDF_0.str()).c_str());
  
  stringstream oss_Prob;
  oss_Prob << setprecision(3) << Prob;
  
  TLatex l03;
  l03.SetTextAlign(9);
  l03.SetTextSize(0.048);
  l03.SetLineWidth(2);
  l03.SetNDC();
  //l03.DrawLatex(0.5069345,0.7699663, ("Prob: " + oss_Prob.str()).c_str());
  
  
  
  stringstream oss_COM;
  oss_COM   << setprecision(3) << COM;
  TLatex l04;
  l04.SetTextAlign(9);
  l04.SetTextSize(0.048);
  l04.SetLineWidth(2);
  l04.SetNDC();
  l04.DrawLatex(0.5558774,0.5652862, ("#chi^{2}/NDF = "+ oss_Sep.str() + "/" +oss_NDF_0.str() + " = " + oss_COM.str()).c_str());
  
  
  TLatex l05;
  l05.SetTextAlign(9);
  l05.SetTextSize(0.058);
  l05.SetLineWidth(2);
  l05.SetNDC();
  l05.DrawLatex(0.1258774,0.8552862, ("ATLAS work-in-progress"));
  ////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////
  TLegend *leg = new TLegend(0.7078059,0.6690754,0.8987342,0.8986966,NULL,"brNDC");
  leg->SetBorderSize(1);
  leg->SetLineColor(1);
  leg->SetLineStyle(1);
  leg->SetLineWidth(1);
  leg->SetFillColor(0);
  leg->SetFillStyle(1001);
  leg->AddEntry(h3,"Simulation","lep");
  leg->AddEntry(frbq,"Fit","l");
  leg->AddEntry(f25,"Gauss","l");
  leg->AddEntry(f26,"Gauss","l");
  leg->AddEntry(f27,"Landau","l");
  leg->SetTextFont(63);
  leg->SetTextSize(25);
  
  leg->Draw();
  
  
  
  c3->cd(); 
  TPad *pad2 = new TPad("pad2", "pad2", 0, 0.08, 1, 0.39);
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.2);
  pad2->SetGridx(); // vertical grid
  pad2->SetGridy(); // vertical grid
  pad2->Draw();
  pad2->cd();   
  // pad2 becomes the current pad
  
  
  //Calculate the diffrence between Fitfunction and Histogramm h1
  
  double bin_tot = h3 -> GetSize();  
  // double start = h3->GetXaxis()->GetBinCenter(0);
  // double end = h3->GetXaxis()->GetBinCenter(bin_tot);
  double err;
  double div;
  
  
  TH1F *hrbq = new TH1F("hrbq", "dif", bin_tot, 0, 3);
  
  
  
  
  for (int bin=2; bin<=bin_tot;bin++) {
    
    double x = h3->GetXaxis()->GetBinCenter(bin);
    double fval = frbq->Eval(x);
    // double dif = h3->GetBinContent(bin);
    double sub = h3->GetBinContent(bin)-fval;
    
    if(sub > -100 && x > 0.3){
      hrbq->SetBinContent(bin, sub);
      err = h3->GetBinError(bin);
      div = sub/err;
      hrbq->SetBinContent(bin, div);
    }else continue;
    
  };
  
  
  
  
	
  // int nbins = hrbq -> GetNbinsX();
  
  float lower_edge  = hrbq -> GetBinLowEdge(1);
  // float bin_width   = hrbq -> GetBinWidth(1);
  float number_bins = hrbq -> GetNbinsX();
  float upper_edge = hrbq -> GetBinLowEdge(number_bins) + hrbq->GetBinWidth(number_bins);
  
  
  
  
  
  
  
  hrbq->GetXaxis()->SetRangeUser(0.1,4);
  
  hrbq->GetXaxis()->SetLabelFont(63);
  hrbq->GetXaxis()->SetLabelSize(28);
  hrbq->GetXaxis()->SetTitle("R_{bq}^{reco}");
  hrbq->GetXaxis()->SetTitleSize(28);
  hrbq->GetXaxis()->SetTitleOffset(3.0);
  hrbq->GetXaxis()->SetTitleFont(63);
  
  
  
  hrbq->GetYaxis()->SetLabelFont(63);
  hrbq->GetYaxis()->SetLabelSize(28);
  hrbq->GetYaxis()->SetTitleSize(28);
  hrbq->GetYaxis()->SetTitleFont(63);
  
  
  hrbq->GetYaxis()->SetRangeUser(-3.4,3.4);	
  
  hrbq->GetYaxis()->SetTitleSize(28);
  hrbq->GetYaxis()->SetTitleOffset(1.45);
  hrbq->GetYaxis()->SetTitleFont(63);
  
  hrbq->GetYaxis()->SetTitle("(Sim.-Fit)/Error[#sigma]");
  hrbq->SetLineWidth(3);  
  
  hrbq->Draw();
  
  
  TF1 *norm1 = new TF1("fa1","0", lower_edge, upper_edge);
  norm1 -> SetLineColor(kRed);
  norm1 -> SetLineStyle(1);
  norm1 -> SetLineWidth(3);
  norm1->Draw("SAME");
  
  
  string Name = file->GetName();
 // string Base = gSystem->BaseName(Name.c_str());
  
  //c3 -> Print((out_dir + Base + "rbq.pdf").c_str());
  c3 -> Print((out_dir + "/" + ma + "_temp_rbq.png").c_str());

}	

	
	
	
	
////////////////////////////////////////////////////////////////////////
//Fit mlb
////////////////////////////////////////////////////////////////////////

void TemplatePara :: mlb_fit(){ 
  
  
 // SetAtlasStyle();
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  double tsize = gStyle->GetLabelSize("X");

  c_mlb = new TCanvas("c_mlb","c_mlb", 800, 10, 600, 500);
  c_mlb->cd();
  
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

  pad1->cd();               // pad1 becomes the current pad
  
  double fitmin = f_mlb->GetXmin();
  double fitmax = f_mlb->GetXmax();

  h_mlb->GetXaxis()->SetRangeUser(fitmin,fitmax);
  h_mlb->SetMaximum(1.2*h_mlb->GetMaximum());
  //h_mlb->GetYaxis()->SetTitle(Form("Normalised events / %2.1f GeV", h_mlb->GetBinWidth(1))); 
  h_mlb->GetYaxis()->SetTitle(Form("Events / %2.1f GeV", h_mlb->GetBinWidth(1))); 
  h_mlb->SetLabelOffset(0.3,"X");
  h_mlb->SetLabelSize(0.06, "Y");
  h_mlb->SetTitleSize(0.07, "Y");
  h_mlb->SetTitleOffset(1.1,"Y");
  
  //Fit of the histogram   
  f_mlb->SetLineColor(2);
  f_mlb->SetLineWidth(3);
  h_mlb->SetLineWidth(3);
	
  h_mlb->Draw();
  h_mlb->Fit("f_mlb","RI0","",fitmin,fitmax);
  
  
  //Get fit parameters 
  ////////////////////////////////////////////////////////////////////
  npar_mlb = 6;
  
  double p[npar_mlb];
  
  for(int i = 0; i < npar_mlb; i++){
    
    p[i] = f_mlb->GetParameter(i);
    par_mlb[i] = p[i];
    err_mlb[i] = f_mlb->GetParError(i);
  };
  
  //Draw functions
  ////////////////////////////////////////////////////////////////////
  f1_mlb->SetParameters(p[0],p[1],p[2]);
  f1_mlb->SetLineColor(6);
  f1_mlb->SetFillStyle(0); 
  f1_mlb->SetLineWidth(3);   
  f1_mlb->Draw("SAME");
  
  f2_mlb->SetParameters(p[3],p[4],p[5]);
  f2_mlb->SetLineColor(209);
  f2_mlb->SetFillStyle(0); 
  f2_mlb->SetLineWidth(3);  
  f2_mlb->Draw("SAME");
  
  //f3_mlb->SetParameters(p[6],p[7],p[8]);
  //f3_mlb->SetLineColor(12);
  //f3_mlb->SetFillStyle(0); 
  //f3_mlb->SetLineWidth(3);  
  //f3_mlb->Draw("SAME");

  h_mlb->Draw("same");
  f_mlb->Draw("same");
  
  //Get chi^2 
  ////////////////////////////////////////////////////////////////////
  double chi2_0 = f_mlb->GetChisquare();
  double NDF_0 = f_mlb->GetNDF();
  double Prob = f_mlb->GetProb();
  double COM = chi2_0/NDF_0;
  
  stringstream oss_chi2_0; oss_chi2_0 << setprecision(3) << chi2_0;
  stringstream oss_NDF_0; oss_NDF_0 << NDF_0;
  stringstream oss_Prob; oss_Prob << setprecision(3) << Prob;
  stringstream oss_COM; oss_COM << setprecision(3) << COM;
  
  TLatex l04;
  l04.SetTextAlign(9);
  l04.SetTextSize(0.048);
  l04.SetLineWidth(2);
  l04.SetNDC();
  //l04.DrawLatex(0.1358774,0.6352862, ("#chi^{2}/NDF: " + oss_COM.str()).c_str());
  //l04.DrawLatex(0.1,0.77, ("#chi^{2}/NDF = "+ oss_chi2_0.str() + "/" +oss_NDF_0.str() + " = " + oss_COM.str()).c_str());
  l04.DrawLatex(0.19,0.8, ("#chi^{2}/NDF = "+ oss_chi2_0.str() + "/" +oss_NDF_0.str() + " = " + oss_COM.str()).c_str());
  TLatex l05;
  l05.SetTextAlign(9);
  l05.SetTextSize(0.058);
  l05.SetLineWidth(2);
  l05.SetNDC();
  l05.DrawLatex(0.19,0.855, ("ATLAS work-in-progress"));
  ////////////////////////////////////////////////////////////////////
  
  //Legend
  ////////////////////////////////////////////////////////////////////    
  TLegend *leg = new TLegend(0.7078059,0.6690754,0.8987342,0.8986966,NULL,"brNDC");
  leg->SetLineColor(kNone);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(tsize);
  leg->AddEntry(h_mlb,"Simulation","lep");
  leg->AddEntry(f_mlb,"Fit","l");
  leg->AddEntry(f1_mlb,"Gauss","l");
  leg->AddEntry(f2_mlb,"Landau","l");
  //leg->AddEntry(f3_mlb,"Landau_n","l");
  //leg->SetTextFont(63);
  leg->Draw("same");
  
  
  
  ////////////////////////////////////////////////////////////////////
  //c_mlb->cd(); 
  pad2->cd();    // pad2 becomes the current pad
  
  
  //Calculate the diffrence between Fitfunction and Histogramm
  ////////////////////////////////////////////////////////////////////
  
  TH1F* h_ratio = (TH1F*) FitRatio(f_mlb,h_mlb);

  
  h_ratio->GetXaxis()->SetTitle("m_{lb}^{reco}[GeV]");
  h_ratio->GetYaxis()->SetTitle("(Sim.-fit)/#sigma");

  h_ratio->SetMaximum(4);    
  h_ratio->SetMinimum(-4);  
  
  h_ratio->SetLabelSize(0.14, "X");
  h_ratio->SetLabelOffset(0.05,"X");
  h_ratio->SetTitleOffset(1.2,"X");
  h_ratio->SetTitleSize(0.17, "X");

  h_ratio->SetLabelSize(0.10, "Y");
  h_ratio->SetTitleOffset(0.5,"Y");
  h_ratio->SetTitleSize(0.15, "Y");
  h_ratio->SetNdivisions(305,"Y");

  h_ratio->SetLineWidth(3);    
    
  h_ratio->Draw("hist");

  TLine* tline_0 = new TLine(fitmin,0,fitmax,0);
  tline_0->SetLineColor(1);
  tline_0->SetLineStyle(1);
  tline_0->SetLineWidth(2);
  tline_0->Draw("SAME");

  h_ratio->Draw("hist same");
  
  //string Name = file->GetName();
  //string Base = gSystem->BaseName(Name.c_str());
  
  stringstream oss_mass; oss_mass << mass;

  c_mlb->Print((out_dir + oss_mass.str() + "_mlb.pdf").c_str());
  //c_mlb -> Print((out_dir + ma + "_temp_mlb.png").c_str());
  
  
  cout<<"#############################################################################################################################################################"<<endl;
  cout<<"#############################################################################################################################################################"<<endl;
  cout<<"#############################################################################################################################################################"<<endl;
  cout<<"#############################################################################################################################################################"<<endl;
  cout<<"#############################################################################################################################################################"<<endl;
  cout<<"#############################################################################################################################################################"<<endl;
  cout<<"#############################################################################################################################################################"<<endl;
  cout<<"#############################################################################################################################################################"<<endl;
  cout<<"#############################################################################################################################################################"<<endl;
  cout<<"#############################################################################################################################################################"<<endl;
  cout<<"#############################################################################################################################################################"<<endl;
  cout<<"#############################################################################################################################################################"<<endl;
  cout<<"#############################################################################################################################################################"<<endl;
  
}

int in_out(const char* in, const char* out, const string plot_dir, const string AnalysisType){  
  
  
  if(AnalysisType == "lepjets"){
    ////////////////////////////////////////////////////////////////////
    TemplatePara* p= new TemplatePara(in, plot_dir, AnalysisType);
    
    p->mtop_fit();
    
    p->mw_fit();
    
    p->rbq_fit();
    
    ////////////////////////////////////////////////////////////////////
    TFile *fend = new TFile(out,"RECREATE");
    TTree *t1 = new TTree("t1","Data");
    
    
    t1->Branch("mass",&p->mass,"mass/D");
    
    t1->Branch("fit_mtop_npar",&p->npar_mtop,"fit_mtop_npar/I");
    t1->Branch("fit_mtop_par",&p->par_mtop,"fit_mtop_par[fit_mtop_npar]/D");
    t1->Branch("fit_mtop_err",&p->err_mtop,"fit_mtop_err[fit_mtop_npar]/D");
    
    t1->Branch("fit_mW_npar",&p->npar_mw,"fit_mW_npar/I");
    t1->Branch("fit_mW_par",&p->par_mw,"fit_mW_par[fit_mW_npar]/D");
    t1->Branch("fit_mW_err",&p->err_mw,"fit_mW_err[fit_mW_npar]/D");
    
    t1->Branch("fit_rbq_npar",&p->npar_rbq,"fit_rbq_npar/I");
    t1->Branch("fit_rbq_par",&p->par_rbq,"fit_rbq_par[fit_rbq_npar]/D");
    t1->Branch("fit_rbq_err",&p->err_rbq,"fit_rbq_err[fit_rbq_npar]/D");
    
    t1->Fill();   
    fend->Write();
    fend->Close();
    ////////////////////////////////////////////////////////////////////	
  }
  else if(AnalysisType == "dilepton"){
    ////////////////////////////////////////////////////////////////////
    TemplatePara* p= new TemplatePara(in, plot_dir, AnalysisType);
    
    p->mlb_fit();
    
    //p->rbq_dilep_fit();
    
    ////////////////////////////////////////////////////////////////////
    TFile *fend = new TFile(out,"RECREATE");
    TTree *t1 = new TTree("t1","Data");
    
    
    t1->Branch("mass",&p->mass,"mass/D");
    
    t1->Branch("fit_mlb_npar",&p->npar_mlb,"fit_mlb_npar/I");
    t1->Branch("fit_mlb_par",&p->par_mlb,"fit_mlb_par[fit_mlb_npar]/D");
    t1->Branch("fit_mlb_err",&p->err_mlb,"fit_mlb_err[fit_mlb_npar]/D");
    
    //t1->Branch("fit_rbq_dilep_npar",&p->npar_rbq_dilep,"fit_rbq_dilep_npar/I");
    //t1->Branch("fit_rbq_dilep_par",&p->par_rbq_dilep,"fit_rbq_dilep_par[fit_rbq_dilep_npar]/D");
    //t1->Branch("fit_rbq_dilep_err",&p->err_rbq_dilep,"fit_rbq_dilep_err[fit_rbq_dilep_npar]/D");
    
    t1->Fill();   
    fend->Write();
    fend->Close();
    ////////////////////////////////////////////////////////////////////	
  }
  else
    cout << endl << AnalysisType << " is an unknown analysis type!" << endl << endl;

  return 0;
}	

std::vector<TemplatePara*> gALL;
std::vector<std::vector<double> >  gParameters;

std::vector<double> gMTOP={170,171.5,172.5,173.5,175};
std::vector<double> gbJSF={0.96, 0.98, 1.0,   1.02, 1.04};
std::vector<double> gJSF={0.96, 0.98, 1.0,   1.02, 1.04};

/*
std::vector<double> gMTOP={170,171.5,172.5,173.5,175};
std::vector<double> gbJSF={1.0};
std::vector<double> gJSF={1.0};
*/



void big_ftop_fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
   
//par = A0, B0 .....
double chisq=0;
int j;
//double mtop[5]={170,171,172.5,173,175};

///gParameters[gALL.size()-1][0]=par[37];

for (j=0;j< gALL.size();j++)
{

double localpar[9];
for (int k=0;k<9;k++)
localpar[k]=par[4*k+0]+par[4*k+1]*(gParameters[j][0]-172.5)+par[4*k+2]*(gParameters[j][1]-1.0)+par[4*k+3]*(gParameters[j][2]-1.0);
/*
localpar[1]=par[2]+par[3]*(gParameters[j][0]-172.5);
localpar[2]=par[4]+par[5]*(gParameters[j][0]-172.5);
localpar[3]=par[6]+par[7]*(gParameters[j][0]-172.5);
localpar[4]=par[8]+par[9]*(gParameters[j][0]-172.5);
localpar[5]=par[10]+par[11]*(gParameters[j][0]-172.5);
localpar[6]=par[12]+par[13]*(gParameters[j][0]-172.5);
localpar[7]=par[14]+par[15]*(gParameters[j][0]-172.5);
localpar[8]=par[16]+par[17]*(gParameters[j][0]-172.5);
*/


   gALL[j]->ftop->SetParameters(localpar);
   
   for (int i=1;i<gALL[j]->h1->GetNbinsX()+1;i++)
   {
   if (gALL[j]->h1->GetBinLowEdge(i)<130) continue;
   if (gALL[j]->h1->GetBinLowEdge(i)+gALL[j]->h1->GetBinWidth(i)>210) continue;
   
   double v=std::abs(gALL[j]->ftop->Integral(gALL[j]->h1->GetBinLowEdge(i),gALL[j]->h1->GetBinLowEdge(i)+gALL[j]->h1->GetBinWidth(i))/(gALL[j]->h1->GetBinWidth(i)));
   if (gALL[j]->h1->GetBinContent(i)<0.001) continue;
   if (gALL[j]->h1->GetBinError(i)<0.000001) continue;
   chisq+=pow((v-gALL[j]->h1->GetBinContent(i))/gALL[j]->h1->GetBinError(i),2);
   
   }
}   
   f = chisq;
}
#include <sys/stat.h>
#include <unistd.h>
#include <string>
inline bool exists_test3 (const std::string& name) {
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

int main(int argc, char* argv[]){
  
  string InputFolder    = argv[1];
  //string OutputFolder   = argv[2];
  //string ConfigFileName = argv[3];
  
  
  //Output_PlotFactory_2.4.24_JSF_*_bJSF_*/HistogramFolder/Merge_emujets_comb_2incl_1excl_nominal/Merged_ttbar*
  //JSF 0.96 0.98 1.0   1.02 1.04
  //bJSF 0.96 0.98 1.0   1.02 1.04
  
  
  vector<string> samples;
  
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
  
  
  
  
  if (!exists_test3(file_name)) continue;
	  puts(file_name.c_str());
std::vector<double> par;
par.push_back(*it3);	  
par.push_back(*it2);	  
par.push_back(*it1);	  
gParameters.push_back(par);

samples.push_back(file_name);


  }	  
/*
std::vector<double> pardata;
par.push_back(0.0);	  
par.push_back(bjetscf_data);	  
par.push_back(jetscf_data);	  
gParameters.push_back(pardata);

samples.push_back(file_name);

*/
  //exit(1);
  
  //samples.push_back(InputFolder+"/Merged_ttbar_170.root");
  //samples.push_back(InputFolder+"/Merged_ttbar_171p5.root");
  //samples.push_back(InputFolder+"/Merged_ttbar.root");
  //samples.push_back(InputFolder+"/Merged_ttbar_173p5.root");
  //samples.push_back(InputFolder+"/Merged_ttbar_175.root");


  for(size_t i=0; i<samples.size(); i++)
  {
    string IN =samples[i];
   // string OUT = OutputFolder + "/SepFitPars_" + samples[i] + ".root";
  //  in_out(IN.c_str(), OUT.c_str(), OutputFolder, AnalysisType);
  system("mkdir -p out");
  gALL.push_back( new TemplatePara(IN.c_str(), std::string("./out"), std::string("lepjets")));
  
  }
  
   
    
///////////////////////////////////////////////////////////////////////////////////////////////
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
  /////////////////////////////////////////////////////////////////////////////////////////////

for(size_t i=0; i<samples.size(); i++)
{
 TCanvas* c1 = new TCanvas("c1","c1", 1024, 768);
 c1->cd();
 gALL[i]->h1->Draw();
 gALL[i]->ftop->Draw("PLSAME");
 c1->SaveAs(Form("./out_%i.png",i));
    
}  
  
  
  return 0;
}













