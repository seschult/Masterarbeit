#include "ZROOT.h"
#include "mtop_fit.h"





using namespace std;

////////////////////////////////////////////////////////////////////////
//Fit Functions
////////////////////////////////////////////////////////////////////////


//Gauss
Double_t signal(Double_t *x, Double_t *par) {
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


//fit top masse

Double_t fit_mtop(Double_t *x,  Double_t *par) { 
return signal(x, par) + landau(x, &par[3]) + landau_n(x, &par[6]);
}

//fit W masse
Double_t fit_mw(Double_t *x,  Double_t *par) { 
return signal(x,par) + signal(x,&par[3]);
}
//fit Rbq
Double_t fit_rbq(Double_t *x,  Double_t *par) { 
return signal(x,par) + signal(x,&par[3]) + landau(x,&par[6]);
}

////////////////////////////////////////////////////////////////////////
//Class def
////////////////////////////////////////////////////////////////////////



ClassImp(mtop_fit)

////////////////////////////////////////////////////////////////////////
//Default constructor
////////////////////////////////////////////////////////////////////////


mtop_fit :: mtop_fit (){


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

//std::map<std::string,TF1*>    fFunctions=0;
	c1=0;//
	c2=0;//
	c3=0;//


}

/////////////////////////////////////////////////////////////////////////
//Constructor
////////////////////////////////////////////////////////////////////////

mtop_fit :: mtop_fit (const char *File){
	
	//Top 
	f00 = new TF1("f5",signal,125,210,3);
	f01 = new TF1("f1",landau,125,210,3);
	f02 = new TF1("f2",landau_n,125,210,3);
    ftop = new TF1("ftop",fit_mtop,125,210,9);
	  
	 
	//MW 
    f13 = new TF1("f3",signal,56,110,3);
    f14 = new TF1("f4",signal,56,110,3);
	fmw = new TF1("fmw",fit_mw,56,110,6);
	
	
	//RBQ 
	f25 = new TF1("f5",signal,0.3,3,3);
	f26 = new TF1("f6",signal,0.3,3,3);
	f27 = new TF1("f7",landau,0.3,3,3);
    frbq = new TF1("frbq",fit_rbq,0.3,3,9);
	
	
	
	 
	ftop->SetParameters(1318, 167, -1.05, 4340, 182, 1.27, 54200, 142 ,9.56);
	ftop->SetParLimits(3,0.0,1000000);
	ftop->SetParLimits(6,0.0,1000000);
    ftop->SetParLimits(7,100.0,1000000);
	  
	fmw->SetParameters(1742, 78, 6.8, 1616, 79, 19); 
	
	frbq->SetParameters(600, 1.5, 0.7, 1000, 1, 0.5, 4525, 0.95 ,0.3);
	frbq->SetParLimits(0,0.0,1000000);
    frbq->SetParLimits(3,0.0,1000000);
    frbq->SetParLimits(6,0.0,1000000);
	
    
    
	// Open file and histogram 
    puts(File);  
	file = new TFile(File,"r");
  
	h1 = (TH1F*)file->Get("hist_klf_mtop_window");
	h2 = (TH1F*)file->Get("hist_klf_window_Whad_m");
	h3 = (TH1F*)file->Get("hist_klf_window_Rbq_reco");
	
	
}
	
////////////////////////////////////////////////////////////////////////
//Fit top Masse
////////////////////////////////////////////////////////////////////////

void mtop_fit :: top_fit(){ 
	
	
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
	h1->GetXaxis()->SetLabelSize(20);
	h1->GetXaxis()->SetTitleSize(0.035);
    h1->GetXaxis()->SetTitleOffset(1.11);
    h1->GetXaxis()->SetTitleFont(42);
	h1->GetXaxis()->SetTitle("m_{top}^{reco}[GeV]"); 
	h1->GetXaxis()->SetRangeUser(124,214);
	h1->GetYaxis()->SetLabelFont(63);
	h1->GetYaxis()->SetLabelSize(20);
	h1->GetYaxis()->SetRangeUser(-120,19300);
	h1->GetYaxis()->SetTitleOffset(1.40);
	h1->GetYaxis()->SetTitle("Entries");
	
	//Fit of the histogram   
	ftop->SetFillColor(19);
	ftop->SetFillStyle(0);
	ftop->SetLineColor(2);
    ftop->SetLineWidth(3);
     h1->SetLineWidth(3);
	
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
 
	//Draw functions
	////////////////////////////////////////////////////////////////////
	f00->SetParameters(p[0],p[1],p[2]);
	f00->SetLineColor(6);
	f00->SetFillStyle(0); 
	f00->SetLineWidth(3);   
	f00->Draw("SAME");
 
	f01->SetParameters(p[3],p[4],p[5]);
	f01->SetLineColor(7);
	f01->SetFillStyle(0); 
	f01->SetLineWidth(3);  
	f01->Draw("SAME");
  
	f02->SetParameters(p[6],p[7],p[8]);
	f02->SetLineColor(8);
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
	
	std::stringstream oss_Prob;
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
   c1->cd(); 
   TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
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
    double start = h1->GetXaxis()->GetBinCenter(0);
    double end = h1->GetXaxis()->GetBinCenter(bin_tot);
    double err;
    double div;
    
    
   
	TH1F *htop = new TH1F("htop", "dif", bin_tot, 120, 220);

	  htop->GetXaxis()->SetRangeUser(124,214);
      htop->GetXaxis()->SetLabelFont(42);
      htop->GetXaxis()->SetLabelSize(0.08);
      htop->GetXaxis()->SetTitle("m_{top}^{reco}[GeV]");
      htop->GetXaxis()->SetTitleSize(0.08);
      htop->GetXaxis()->SetTitleOffset(1.18);
      htop->GetXaxis()->SetTitleFont(42);
      htop->GetYaxis()->SetLabelFont(42);
      htop->GetYaxis()->SetLabelSize(0.08);
      htop->GetYaxis()->SetTitleSize(0.035);
      htop->GetYaxis()->SetTitleFont(42);
      htop->GetYaxis()->SetRangeUser(-2.9,2.9);	      
      htop->GetYaxis()->SetTitleSize(0.08);
	  htop->GetYaxis()->SetTitleOffset(0.29);
	  htop->GetYaxis()->SetTitleFont(42);
      htop->GetYaxis()->SetTitle("(Sim.-Fit)/Error[#sigma]");
      htop->SetLineWidth(3);  
	 
 	 
 	for (int bin=2; bin<=bin_tot;bin++) {
		
		double x = h1->GetXaxis()->GetBinCenter(bin);
		double fval = ftop->Eval(x);
		double dif = h1->GetBinContent(bin);
		double sub = h1->GetBinContent(bin)-fval;
    
		if(sub > -100){
			htop->SetBinContent(bin, sub);
			err = h1->GetBinError(bin);
			div = sub/err;
			htop->SetBinContent(bin, div);
			
			
			}else continue;
		
		};

	//Red line
	////////////////////////////////////////////////////////////////////
	int nbins = htop -> GetNbinsX();

	float lower_edge  = htop -> GetBinLowEdge(1);
	float bin_width   = htop -> GetBinWidth(1);
	float number_bins = htop -> GetNbinsX();
	float upper_edge = htop -> GetBinLowEdge(number_bins) + htop->GetBinWidth(number_bins);

	htop->Draw();

	TF1 *norm1 = new TF1("fa1","0", lower_edge, upper_edge);
	norm1 -> SetLineColor(kRed);
	norm1 -> SetLineStyle(1);
	norm1 -> SetLineWidth(3);
	norm1->Draw("SAME");
	
	string Name = file->GetName();
	string Base = gSystem->BaseName(Name.c_str());
	
	 c1 -> Print(("/home/iwsatlas1/schulte/plots_15217/"+ Base +"mtop.png").c_str());

 
 
 
 
}




////////////////////////////////////////////////////////////////////////
//Fit mw
////////////////////////////////////////////////////////////////////////

void mtop_fit :: mw_fit(){
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
	h2->GetXaxis()->SetLabelSize(20); // labels will be 14 pixels
	h2->GetXaxis()->SetTitleSize(0.035);
    h2->GetXaxis()->SetTitleOffset(1.11);
    h2->GetXaxis()->SetTitleFont(42);
	h2->GetXaxis()->SetTitle("m_{w}^{reco}[GeV]"); // labels will be 14 pixels
	h2->GetYaxis()->SetLabelFont(63);
	h2->GetYaxis()->SetLabelSize(20);
	h2->GetYaxis()->SetRangeUser(-120,34030);
	h2->GetYaxis()->SetTitleOffset(1.40);
	
	h2->GetYaxis()->SetTitle("Entries");
	
	//Fit of the histogram   
	
	
	fmw->SetFillColor(19);
	fmw->SetFillStyle(0);
	fmw->SetLineColor(3);
    fmw->SetLineWidth(3);   
    h2->SetLineWidth(3);
	
	   
	h2->Fit("fmw","","",56,110);
	

	
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
	f13->SetLineColor(3);
	f13->SetFillStyle(0); 
	f13->SetLineWidth(3);
	f13->Draw("SAME");
 
	f14->SetParameters(p[3],p[4],p[5]);
	f14->SetLineColor(7);
	f14->SetFillStyle(0); 
	f14->SetLineWidth(3);
	f14->Draw("SAME");


//Get chi^2
////////////////////////////////////////////////////////////////////////
	
	

	double chi2_0 = fmw->GetChisquare();
	double NDF_0 = fmw->GetNDF();
	double Prob = fmw->GetProb();
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
	
	std::stringstream oss_Prob;
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
	////////////////////////////////////////////////////////////////////
	
	
   /////////////////////////////////////////////////////////////////////	
   TLegend *leg = new TLegend(0.7078059,0.6690754,0.8987342,0.8986966,NULL,"brNDC");
   leg->SetBorderSize(1);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(1001);
   leg->AddEntry(h2,"Data","lep");
   leg->AddEntry(fmw,"Fit","l");
   leg->AddEntry(f13,"Gauss","l");
   leg->AddEntry(f14,"Gauss","l");
  
   leg->Draw();
	
 
	
	
   c2->cd(); 
   TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
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
    double start = h2->GetXaxis()->GetBinCenter(0);
    double end = h2->GetXaxis()->GetBinCenter(bin_tot);
    double div;
    double err;
    
	TH1F *hw = new TH1F("hw", "dif", bin_tot, 40, 120);
	hw->GetXaxis()->SetLabelFont(63);
	hw->GetXaxis()->SetLabelSize(14); // labels will be 14 pixels
	hw->GetYaxis()->SetLabelFont(63);
	hw->GetYaxis()->SetLabelSize(14);
	//TF1 *func = h1->GetFunction("ftop");
	 
 	for (int bin=2; bin<=bin_tot;bin++) {
		
		double x = h2->GetXaxis()->GetBinCenter(bin);
		double fval = fmw->Eval(x);
		double dif = h2->GetBinContent(bin);
		double sub = h2->GetBinContent(bin)-fval;
    
		if(sub > -100){
			hw->SetBinContent(bin, sub);
			err = h2->GetBinError(bin);
			div = sub/err;
			hw->SetBinContent(bin, div);
			}else continue;
		
		};
	
	
 	int nbins = hw -> GetNbinsX();

	float lower_edge  = hw -> GetBinLowEdge(1);
	float bin_width   = hw -> GetBinWidth(1);
	float number_bins = hw -> GetNbinsX();
	float upper_edge = hw -> GetBinLowEdge(number_bins) + hw->GetBinWidth(number_bins);
       
	  hw->GetXaxis()->SetRangeUser(54,115);
      hw->GetXaxis()->SetLabelFont(42);
      hw->GetXaxis()->SetLabelSize(0.08);
      hw->GetXaxis()->SetTitle("m_{w}^{reco}[GeV]");
      hw->GetXaxis()->SetTitleSize(0.08);
      hw->GetXaxis()->SetTitleOffset(1.18);
      hw->GetXaxis()->SetTitleFont(42);
      hw->GetYaxis()->SetLabelFont(42);
      hw->GetYaxis()->SetLabelSize(0.08);
      hw->GetYaxis()->SetTitleSize(0.035);
      hw->GetYaxis()->SetTitleFont(42);
      hw->GetYaxis()->SetRangeUser(-3.4,3.4);	      
      hw->GetYaxis()->SetTitleSize(0.08);
	  hw->GetYaxis()->SetTitleOffset(0.29);
	  hw->GetYaxis()->SetTitleFont(42);
      hw->GetYaxis()->SetTitle("(Sim.-Fit)/Error[#sigma]");
      hw->SetLineWidth(3);  
	
	
	
	hw->Draw();
	
	TF1 *norm1 = new TF1("fa1","0", lower_edge, upper_edge);
	norm1 -> SetLineColor(kRed);
	norm1 -> SetLineStyle(1);
	norm1 -> SetLineWidth(3);
	norm1->Draw("SAME");
	

	string Name = file->GetName();
	string Base = gSystem->BaseName(Name.c_str());
	
	 //c2 -> Print(("~/plots15217/" + Base +".png").c_str());
 	 c2 -> Print(("/home/iwsatlas1/schulte/plots_15217/"+ Base +"mw.png").c_str());
 
}
 


////////////////////////////////////////////////////////////////////////
//Fit rbq
////////////////////////////////////////////////////////////////////////
	
void mtop_fit :: rbq_fit(){ 
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
	h3->GetXaxis()->SetLabelSize(20); // labels will be 14 pixels
	h3->GetXaxis()->SetTitleSize(0.035);
    h3->GetXaxis()->SetTitleOffset(1.11);
    h3->GetXaxis()->SetTitleFont(42);
	h3->GetXaxis()->SetTitle("R_{bq}^{reco}"); // labels will be 14 pixels
	h3->GetYaxis()->SetLabelFont(63);
	h3->GetYaxis()->SetLabelSize(20);
	h3->GetYaxis()->SetRangeUser(-120,20030);
	
	h3->GetYaxis()->SetTitleOffset(1.40);
	h3->GetYaxis()->SetTitle("Entries");
	
	
	//Fit of the histogram
	
	frbq->SetFillColor(19);
	frbq->SetFillStyle(0);
	frbq->SetLineColor(2);
    frbq->SetLineWidth(3);   
    h3->SetLineWidth(3);
	
	h3->GetXaxis()->SetRangeUser(0.1,4);  
	
	//Fit of the histogram   
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
	f25->SetLineColor(3);
	f25->SetLineWidth(3);
	f25->Draw("SAME");
 
	f26->SetParameters(p[3],p[4],p[5]);
	f26->SetLineColor(7);
	f26->SetLineWidth(3);
	f26->Draw("SAME");
  

	f27->SetParameters(p[6],p[7],p[8]);
	f27->SetLineColor(9);
	f27->SetLineWidth(3);
	f27->Draw("SAME");
	
	//////////////////////////////////////////////////////////////////////
	double chi2_0 = fmw->GetChisquare();
	double NDF_0 = fmw->GetNDF();
	double Prob = fmw->GetProb();
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
	
	std::stringstream oss_Prob;
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
	
	//////////////////////////////////////////////////////////////////////
   TLegend *leg = new TLegend(0.7078059,0.6690754,0.8987342,0.8986966,NULL,"brNDC");
   leg->SetBorderSize(1);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(1001);
   leg->AddEntry(h3,"Data","lep");
   leg->AddEntry(frbq,"Fit","l");
   leg->AddEntry(f25,"Gauss","l");
   leg->AddEntry(f26,"Gauss","l");
   leg->AddEntry(f27,"Landau","l");
   leg->Draw();

 
 
   c3->cd(); 
   TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
   pad2->SetTopMargin(0);
   pad2->SetBottomMargin(0.2);
   pad2->SetGridx(); // vertical grid
   pad2->SetGridy(); // vertical grid
   pad2->Draw();
   pad2->cd();   
       // pad2 becomes the current pad
       
       
	//Calculate the diffrence between Fitfunction and Histogramm h1
	
	double bin_tot = h3 -> GetSize();  
    double start = h3->GetXaxis()->GetBinCenter(0);
    double end = h3->GetXaxis()->GetBinCenter(bin_tot);
    double err;
    double div;
    
	
	TH1F *hrbq = new TH1F("hrbq", "dif", bin_tot, 0, 3);
	hrbq->GetXaxis()->SetLabelFont(63);
	hrbq->GetXaxis()->SetLabelSize(14); // labels will be 14 pixels
	hrbq->GetYaxis()->SetLabelFont(63);
	hrbq->GetYaxis()->SetLabelSize(14);
	
	
 
 	for (int bin=2; bin<=bin_tot;bin++) {
		
		double x = h3->GetXaxis()->GetBinCenter(bin);
		double fval = frbq->Eval(x);
		double dif = h3->GetBinContent(bin);
		double sub = h3->GetBinContent(bin)-fval;
    
		if(sub > -100 && dif !=0){
			hrbq->SetBinContent(bin, sub);
			err = h3->GetBinError(bin);
			div = sub/err;
			hrbq->SetBinContent(bin, div);
			}else continue;
		
		};
		
		
		
		
	
 	int nbins = hrbq -> GetNbinsX();

	float lower_edge  = hrbq -> GetBinLowEdge(1);
	float bin_width   = hrbq -> GetBinWidth(1);
	float number_bins = hrbq -> GetNbinsX();
	float upper_edge = hrbq -> GetBinLowEdge(number_bins) + hrbq->GetBinWidth(number_bins);
       
       
       
	
	
	
	
	  hrbq->GetXaxis()->SetRangeUser(0.1,4);
      hrbq->GetXaxis()->SetLabelFont(42);
      hrbq->GetXaxis()->SetLabelSize(0.08);
      hrbq->GetXaxis()->SetTitle("R_{bq}^{reco}");
      hrbq->GetXaxis()->SetTitleSize(0.08);
      hrbq->GetXaxis()->SetTitleOffset(1.18);
      hrbq->GetXaxis()->SetTitleFont(42);
 
 

      hrbq->GetYaxis()->SetLabelFont(42);
      hrbq->GetYaxis()->SetLabelSize(0.08);
      hrbq->GetYaxis()->SetTitleSize(0.035);
      hrbq->GetYaxis()->SetTitleFont(42);
      hrbq->GetYaxis()->SetRangeUser(-5.4,5.4);	
      
      hrbq->GetYaxis()->SetTitleSize(0.08);
	  hrbq->GetYaxis()->SetTitleOffset(0.29);
	  hrbq->GetYaxis()->SetTitleFont(42);
      
      hrbq->GetYaxis()->SetTitle("(Sim.-Fit)/Error[#sigma]");
      hrbq->SetLineWidth(3);  
		
	  hrbq->Draw();
 
 
 	TF1 *norm1 = new TF1("fa1","0", lower_edge, upper_edge);
	norm1 -> SetLineColor(kRed);
	norm1 -> SetLineStyle(1);
	norm1 -> SetLineWidth(3);
	norm1->Draw("SAME");
		
		
	string Name = file->GetName();
	string Base = gSystem->BaseName(Name.c_str());
	
	// c3 -> Print(("~/plots15217/" + Base +".png").c_str());
 	 c3 -> Print(("/home/iwsatlas1/schulte/plots_15217/"+ Base +"rbq.png").c_str());
}	
	
	
	
	

	
	
	
	
	
