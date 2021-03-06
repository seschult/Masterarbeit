#include "ZROOT.h"






using namespace std;

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





class mtop_fit{
	
	//members of class
	
	public:
	
	mtop_fit(); //constructor
	void top_fit(const char* File); //constructor
	void mw_fit(); //constructor
	void rbq_fit(); //constructor
	
	
	
	mtop_fit(const char *File); //constructor
	 //constructor
	
	//Gauss

	TH1F *h1;
	TH1F *h2;
	TH1F *h3;
	TFile *file;
	
	
	TF1 *f00;
	TF1 *f01;
	TF1 *f02;
	TF1 *f13;
	TF1 *f14;
	TF1 *f25;
	TF1 *f26;
	TF1 *f27;
	
	
	TF1 *ftop;
	TF1 *fmw;
	TF1 *frbq;

//std::map<std::string,TF1*>    fFunctions;
	TCanvas *c1;
	TCanvas *c2;
	TCanvas *c3;
	
	

	
	
	
};







mtop_fit :: mtop_fit (const char *File){
	std::cout<<"constructor ok1";
	//Top 
	f00 = new TF1("f5",signal,125,210,3);
	f01 = new TF1("f1",landau,125,210,3);
	f02 = new TF1("f2",landau_n,125,210,3);
    ftop = new TF1("ftop",fit_mtop,125,210,9);
	  
	 std::cout<<"constructor 2";
	//MW 
    f13 = new TF1("f3",signal,56,110,3);
    f14 = new TF1("f4",signal,56,110,3);
	fmw = new TF1("fmw",fit_mw,56,110,6);
	
	
	//RBQ 
	f25 = new TF1("f5",signal,0.3,3,3);
	f26 = new TF1("f6",signal,0.3,3,3);
	f27 = new TF1("f7",landau,0.3,3,3);
    frbq = new TF1("frbq",fit_rbq,0.3,3,9);
	
	
	
	 
	//ftop->SetParameters(3260, 176, -45, 40000, 189, 4.4, 700000, 160 ,10);
	//ftop->SetParameters(1318, 167, -1.05, 4340, 182, 1.27, 54200, 142 ,9.56);
	ftop->SetParameters(1318, 172, 1.05, 4340, 182, 1.27, 54200, 132 ,9.56);
	
	ftop->SetParLimits(3,0.0,1000000);
	ftop->SetParLimits(6,0.0,1000000);
    ftop->SetParLimits(7,100.0,1000000);
	  
	 
	fmw->SetParameters(1742, 78, 6.8, 1616, 79, 19); 
	
	
	
	//frbq->SetParameters(10, 1, 2, 10, 2, 1, 10, 2 ,1);
	frbq->SetParameters(600, 1.5, 0.7, 1000, 1, 0.5, 4525, 0.95 ,0.3);
	
	//fit->SetParameters(10, 1, 1, 1, 2, 1, 10, 2 ,15);
	frbq->SetParLimits(0,0.0,1000000);
    frbq->SetParLimits(3,0.0,1000000);
    frbq->SetParLimits(6,0.0,1000000);
	
	
    
    
    puts(File);
	// Open file and histogram   
	file = new TFile(File,"r");
	cout << File;
	h1 = (TH1F*)file->Get("hist_klf_mtop_window");
	h2 = (TH1F*)file->Get("hist_klf_window_Whad_m");
	h3 = (TH1F*)file->Get("hist_klf_window_Rbq_reco");
			
}
	


void mtop_fit :: top_fit(const char* File){ 
	
	
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
	h1->GetXaxis()->SetLabelSize(0); // labels will be 14 pixels
	h1->GetXaxis()->SetTitleSize(0.0);
    h1->GetXaxis()->SetTitleOffset(1.0);
    h1->GetXaxis()->SetTitleFont(42);
	h1->GetXaxis()->SetTitle("m_{top}[GeV]"); // labels will be 14 pixels
	h1->GetYaxis()->SetLabelFont(63);
	
	
	
	h1->GetYaxis()->SetLabelSize(26);
	
   double x = h1->GetMaximum();
   double limit = x + x*0.1;
	
	h1->GetYaxis()->SetRangeUser(-430,limit);
	
	h1->GetYaxis()->SetTitleSize(26);
	h1->GetYaxis()->SetTitle("Events");
	h1->GetYaxis()->SetTitleOffset(0.5);
	//Fit of the histogram
	ftop->SetFillColor(19);
	ftop->SetFillStyle(0);
	ftop->SetLineColor(2);
    ftop->SetLineWidth(3);   
    h1->SetLineWidth(3);   
	h1->Fit("ftop","RI","",130,210);
	
	h1->GetXaxis()->SetRangeUser(124,214);
	h1->Draw();
	
	
		
	
	

	//Get fit parameters 
	double p[9];
 
	for(int i = 0; i < 9; i++){
	 
		p[i] = ftop->GetParameter(i);
		};
 
 
 
	//Draw functions
	f00->SetParameters(p[0],p[1],p[2]);
	f00->SetLineColor(3);
	f00->SetFillStyle(0); 
	f00->SetLineWidth(3);   
	f00->Draw("SAME");
 
	f01->SetParameters(p[3],p[4],p[5]);
	f01->SetLineColor(7);
	f01->SetFillStyle(0); 
	f01->SetLineWidth(3);  
	f01->Draw("SAME");
  
	f02->SetParameters(p[6],p[7],p[8]);
	f02->SetLineColor(9);
	f02->SetFillStyle(0); 
	f02->SetLineWidth(3);  
	f02->Draw("SAME");
	
	double chi2_0 = ftop->GetChisquare();
	double NDF_0 = ftop->GetNDF();
	double Prob = ftop->GetProb();
	double COM = chi2_0/NDF_0;
	
	
	
	//Output chi2/NDF/Prob
	////////////////////////////////////////////////////////////////////
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
	l01.DrawLatex(.1458774,0.7811448, ("NDF: " + oss_NDF_0.str()).c_str());
	
	std::stringstream oss_Prob;
	oss_Prob << setprecision(3) << Prob;
	
	TLatex l03;
	l03.SetTextAlign(9);
	l03.SetTextSize(0.038);
	l03.SetLineWidth(2);
	l03.SetNDC();
	l03.DrawLatex(0.1469345,0.8299663, ("Prob: " + oss_Prob.str()).c_str());
	
	
	
	std::stringstream oss_COM;
	oss_COM   << setprecision(3) << COM;
	TLatex l04;
	l04.SetTextAlign(9);
	l04.SetTextSize(0.038);
	l04.SetLineWidth(2);
	l04.SetNDC();
	l04.DrawLatex(0.1458774,0.6952862, ("#chi^{2}/NDF: " + oss_COM.str()).c_str());
	////////////////////////////////////////////////////////////////////
	
   //TLegend *leg = new TLegend(0.6,0.7,0.78,0.9);
   TLegend *leg = new TLegend(0.7078059,0.6690754,0.8987342,0.8986966,NULL,"brNDC");
   leg->SetBorderSize(1);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(1001);
   leg->AddEntry(h1,"Data","lep");
   leg->AddEntry(ftop,"Fit","l");
   leg->AddEntry(f00,"Gauss","l");
   leg->AddEntry(f01,"Landau","l");
   leg->AddEntry(f02,"Landau_n","l");
   leg->Draw();
   

   c1->cd(); 
   TPad *pad2 = new TPad("pad2", "pad2", 0, 0.08, 1, 0.39);
   pad2->SetTopMargin(0);
   pad2->SetBottomMargin(0.2);
   pad2->Range(0,0,0,0.0);
   pad2->SetGridx(); // vertical grid
   pad2->Draw();
   pad2->cd(); 
 
 
     
       
	//Calculate the diffrence between Fitfunction and Histogramm h1
	std::cout<<"constructor 37";
	double bin_tot = h1 -> GetSize();  
    double start = h1->GetXaxis()->GetBinCenter(0);
    double end = h1->GetXaxis()->GetBinCenter(bin_tot);
    double err;
    double div;
    
    
    std::cout<<"constructor 38";
	TH1F *htop = new TH1F("htop", "dif", bin_tot, 120, 220);
	

	  htop->GetXaxis()->SetRangeUser(124,214);
      htop->GetXaxis()->SetLabelFont(63);
      htop->GetXaxis()->SetLabelSize(0.11);
      
      
     
      htop->GetXaxis()->SetTitle("m_{top}[GeV]");
      htop->GetXaxis()->SetTitleSize(0.10);
      htop->GetXaxis()->SetTitleOffset(1.0);
      htop->GetXaxis()->SetTitleFont(42);
 
 

      htop->GetYaxis()->SetLabelFont(63);
      htop->GetYaxis()->SetLabelSize(26);
      htop->GetYaxis()->SetTitleSize(0.11);
      htop->GetYaxis()->SetTitleFont(42);
      	
      
      htop->GetYaxis()->SetTitleSize(0.11);
	  htop->GetYaxis()->SetTitleOffset(0.4);
	  htop->GetYaxis()->SetTitleFont(63);
      
      htop->GetYaxis()->SetTitle("(Data-Fit)/Error[#sigma]");
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
		
		
		double x1 = htop->GetMaximum();
		double limit1 = x1+ x1*0.1;
	
	htop->GetYaxis()->SetRangeUser(-limit1,limit1);	
	
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


	string Name = "massetop";
	std::cout<<Name;
	c1->Print(Form("%s%s%s","~/",Name.c_str(), ".png" ));
	
		
	
  TFile *mtop_par = new TFile("~/mtop_par.root","update");
  
    double label_top; 

  TTree *t1 = (TTree*)(mtop_par->Get("t1"));//
  if (!t1) 
  {
	  
  t1=new TTree("t1","Data");

	
  t1->Branch("mtop",&label_top,"label_top/D");		
  t1->Branch("p0",&(p[0]),"par0/D");	
  t1->Branch("p1",&(p[1]),"par1/D");	
  t1->Branch("p2",&(p[2]),"par2/D");	
  t1->Branch("p3",&(p[3]),"par3/D");	
  t1->Branch("p4",&(p[4]),"par4/D");	
  t1->Branch("p5",&(p[5]),"par5/D");	
  t1->Branch("p6",&p[6],"par6/D");	
  t1->Branch("p7",&p[7],"par7/D");	
  t1->Branch("p8",&p[8],"par8/D");	
 	
  
  
  
  
}
else
{


	
  t1->SetBranchAddress("mtop",&label_top);		
  t1->SetBranchAddress("p0",&(p[0]));	
  t1->SetBranchAddress("p1",&(p[1]));	
  t1->SetBranchAddress("p2",&(p[2]));	
  t1->SetBranchAddress("p3",&(p[3]));	
  t1->SetBranchAddress("p4",&(p[4]));	
  t1->SetBranchAddress("p5",&(p[5]));	
  t1->SetBranchAddress("p6",&p[6]);	
  t1->SetBranchAddress("p7",&p[7]);	
  t1->SetBranchAddress("p8",&p[8]);	
 	
t1->GetEntry(0);

}

  //Get Labels
  

  
     //puts(File);
      std::string fname(gSystem->BaseName(File));
     
     if (fname=="Merged_ttbar_170.root") label_top = 170;
     if (fname=="Merged_ttbar_171p5.root") label_top = 171.5;
     if (fname=="Merged_ttbar_173p5.root") label_top = 173.5;
     if (fname=="Merged_ttbar_175.root") label_top = 175;
     if (fname=="Merged_ttbar.root") label_top = 172.5;
     
     

       
 
		





  t1->Fill();   
  t1->Write("",TObject::kOverwrite);
  mtop_par->Write();
  mtop_par->Close();	
	
	
	
}









void mtop_fit :: mw_fit(){
	double with = 950;
	double hight = 900;
	c2 = new TCanvas("c2","c2", with, hight);
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);
	

	TPad *pad1 = new TPad("pad1", "pad1", 0, 0.32, 1, 1.0);
	pad1->SetBottomMargin(1); // Upper and lower plot are joined
    pad1->SetGridx();         // Vertical grid
    pad1->Draw();             // Draw the upper pad: pad1
    pad1->cd();               // pad1 becomes the current pad
	
		
	h2->GetXaxis()->SetLabelFont(63);
	h2->GetXaxis()->SetLabelSize(0); // labels will be 14 pixels
	h2->GetXaxis()->SetTitleSize(0.035);
    h2->GetXaxis()->SetTitleOffset(1.11);
    h2->GetXaxis()->SetTitleFont(42);
	h2->GetXaxis()->SetTitle("m_{w}[GeV]"); // labels will be 14 pixels
	h2->GetYaxis()->SetLabelFont(63);
	h2->GetYaxis()->SetLabelSize(25);
	
	double x = h2->GetMaximum();
   double limit = x + x*0.1;
	
	h2->GetYaxis()->SetRangeUser(-120,limit);
	h2->GetYaxis()->SetTitle("Entries");
	
	
	//Fit of the histogram
	
	fmw->SetFillColor(19);
	fmw->SetFillStyle(0);
	fmw->SetLineColor(2);
    fmw->SetLineWidth(3);   
    h2->SetLineWidth(3);
	
	   
	h2->Fit("fmw","","",56,110);
	
	h2->GetXaxis()->SetRangeUser(54,115);
	h2->Draw();
	
	//Get fit parameters 
	double p[6];
 
	for(int i = 0; i < 6; i++){
	 
		p[i] = fmw->GetParameter(i);
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


	////////////////////////////////////////////////////////////////////
	
	

	double chi2_0 = fmw->GetChisquare();
	double NDF_0 = fmw->GetNDF();
	double Prob = fmw->GetProb();
	double COM = chi2_0/NDF_0;
	
	
	
	//Output chi2/NDF/Prob
	////////////////////////////////////////////////////////////////////
	std::stringstream oss_Sep;
	oss_Sep   << setprecision(3) << chi2_0;
	TLatex l00;
	l00.SetTextAlign(9);
	l00.SetTextSize(0.038);
	l00.SetLineWidth(2);
	l00.SetNDC();
	l00.DrawLatex(0.1458774,0.7424242, ("#chi^{2}: " + oss_Sep.str()).c_str());
	
	
	std::stringstream oss_NDF_0;
	oss_NDF_0 << setprecision(3) << NDF_0;
	
	TLatex l01;
	l01.SetTextAlign(9);
	l01.SetTextSize(0.038);
	l01.SetLineWidth(2);
	l01.SetNDC();
	l01.DrawLatex(.1458774,0.7811448, ("NDF: " + oss_NDF_0.str()).c_str());
	
	std::stringstream oss_Prob;
	oss_Prob << setprecision(3) << Prob;
	
	TLatex l03;
	l03.SetTextAlign(9);
	l03.SetTextSize(0.038);
	l03.SetLineWidth(2);
	l03.SetNDC();
	l03.DrawLatex(0.1469345,0.8299663, ("Prob: " + oss_Prob.str()).c_str());
	
	
	
	std::stringstream oss_COM;
	oss_COM   << setprecision(3) << COM;
	TLatex l04;
	l04.SetTextAlign(9);
	l04.SetTextSize(0.038);
	l04.SetLineWidth(2);
	l04.SetNDC();
	l04.DrawLatex(0.1458774,0.6952862, ("#chi^{2}/NDF: " + oss_COM.str()).c_str());
	////////////////////////////////////////////////////////////////////
	
	
	////////////////////////////////////////////////////////////////////
	
	
	
 
	//TLegend* leg = new TLegend(0.6,0.7,0.78,0.9);
   //leg->SetHeader("The Legend Title","C"); // option "C" allows to center the header
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
   TPad *pad2 = new TPad("pad2", "pad2", 0, 0.08, 1, 0.38);
   pad2->SetTopMargin(0);
   pad2->SetBottomMargin(0.2);
   pad2->SetGridx(); // vertical grid
   pad2->Draw();
   pad2->cd();   
       // pad2 becomes the current pad
       
       
	//Calculate the diffrence between Fitfunction and Histogramm h1
	
	double bin_tot = h2 -> GetSize();  
    double start = h2->GetXaxis()->GetBinCenter(0);
    double end = h2->GetXaxis()->GetBinCenter(bin_tot);
    double div;
    double err;
    
	TH1F *hw = new TH1F("hw", "dif", bin_tot, 40, 120);
	hw->GetXaxis()->SetLabelFont(63);
	hw->GetXaxis()->SetLabelSize(0); // labels will be 14 pixels
	hw->GetYaxis()->SetLabelFont(63);
	hw->GetYaxis()->SetLabelSize(25);
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
       
       
       
	   double x1 = hw->GetMaximum();
       double limit1 = x1 + x1*0.1;
	
	
	
	  hw->GetXaxis()->SetRangeUser(54,115);
      hw->GetXaxis()->SetLabelFont(42);
      hw->GetXaxis()->SetLabelSize(0.11);
      
      
     
      hw->GetXaxis()->SetTitle("m_{w}[GeV]");
      hw->GetXaxis()->SetTitleSize(0.08);
      hw->GetXaxis()->SetTitleOffset(1.18);
      hw->GetXaxis()->SetTitleFont(42);
 
 

      hw->GetYaxis()->SetLabelFont(42);
      hw->GetYaxis()->SetLabelSize(0.11);
      hw->GetYaxis()->SetTitleSize(0.035);
      hw->GetYaxis()->SetTitleFont(42);
      hw->GetYaxis()->SetRangeUser(-limit1,limit1);	
      
      hw->GetYaxis()->SetTitleSize(0.08);
	  hw->GetYaxis()->SetTitleOffset(0.29);
	  hw->GetYaxis()->SetTitleFont(42);
      
      hw->GetYaxis()->SetTitle("(Data-Fit)/Error[#sigma]");
      hw->SetLineWidth(3);  
	
	
	
	hw->Draw();
	
	TF1 *norm1 = new TF1("fa1","0", lower_edge, upper_edge);
	norm1 -> SetLineColor(kRed);
	norm1 -> SetLineStyle(1);
	norm1 -> SetLineWidth(3);
	norm1->Draw("SAME");
	
 
  string Name = "massew";
	std::cout<<Name;
	c2->Print(Form("%s%s%s","~/",Name.c_str(), ".png" ));
 
 
}








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
	h3->GetXaxis()->SetLabelSize(0); // labels will be 14 pixels
	h3->GetXaxis()->SetTitleSize(0.035);
    h3->GetXaxis()->SetTitleOffset(1.11);
    h3->GetXaxis()->SetTitleFont(42);
	h3->GetXaxis()->SetTitle("R_{bq}"); // labels will be 14 pixels
	h3->GetYaxis()->SetLabelFont(63);
	h3->GetYaxis()->SetLabelSize(25);
	
	
	double x = h3->GetMaximum();
    double limit = x + x*0.1;
	h3->GetYaxis()->SetRangeUser(-120,limit);
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
	
	//npar_rbq = 9;
 
	for(int i = 0; i < 9; i++){
	 
		p[i] = frbq->GetParameter(i);
		//par_rbq[i] = p[i];
		//err_rbq[i] = frbq->GetParError(i);
		
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
	
	
	
	//Output chi2/NDF/Prob
	////////////////////////////////////////////////////////////////////
	std::stringstream oss_Sep;
	oss_Sep   << setprecision(3) << chi2_0;
	TLatex l00;
	l00.SetTextAlign(9);
	l00.SetTextSize(0.038);
	l00.SetLineWidth(2);
	l00.SetNDC();
	l00.DrawLatex(0.1458774,0.7424242, ("#chi^{2}: " + oss_Sep.str()).c_str());
	
	
	std::stringstream oss_NDF_0;
	oss_NDF_0 << setprecision(3) << NDF_0;
	
	TLatex l01;
	l01.SetTextAlign(9);
	l01.SetTextSize(0.038);
	l01.SetLineWidth(2);
	l01.SetNDC();
	l01.DrawLatex(.1458774,0.7811448, ("NDF: " + oss_NDF_0.str()).c_str());
	
	std::stringstream oss_Prob;
	oss_Prob << setprecision(3) << Prob;
	
	TLatex l03;
	l03.SetTextAlign(9);
	l03.SetTextSize(0.038);
	l03.SetLineWidth(2);
	l03.SetNDC();
	l03.DrawLatex(0.1469345,0.8299663, ("Prob: " + oss_Prob.str()).c_str());
	
	
	
	std::stringstream oss_COM;
	oss_COM   << setprecision(3) << COM;
	TLatex l04;
	l04.SetTextAlign(9);
	l04.SetTextSize(0.038);
	l04.SetLineWidth(2);
	l04.SetNDC();
	l04.DrawLatex(0.1458774,0.6952862, ("#chi^{2}/NDF: " + oss_COM.str()).c_str());
	////////////////////////////////////////////////////////////////////
	
	//////////////////////////////////////////////////////////////////////
 
    //leg->SetHeader("The Legend Title","C"); // option "C" allows to center the header
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
   TPad *pad2 = new TPad("pad2", "pad2", 0, 0.08, 1, 0.38);
   pad2->SetTopMargin(0);
   pad2->SetBottomMargin(0.2);
   pad2->SetGridx(); // vertical grid
   pad2->Draw();
   pad2->cd();   
       // pad2 becomes the current pad
       
       
	//Calculate the diffrence between Fitfunction and Histogramm h1
	
	double bin_tot = h3 -> GetSize();  
    double start = h3->GetXaxis()->GetBinCenter(0);
    double end = h3->GetXaxis()->GetBinCenter(bin_tot);
    double err;
    double div;
    
	
	TH1F *hrbq = new TH1F("hw", "dif", bin_tot, 0, 3);
	hrbq->GetXaxis()->SetLabelFont(63);
	hrbq->GetXaxis()->SetLabelSize(0.11); // labels will be 14 pixels
	hrbq->GetYaxis()->SetLabelFont(63);
	hrbq->GetYaxis()->SetLabelSize(0.11);
	
	
 
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
      hrbq->GetXaxis()->SetLabelSize(0.11);
      
      
     
      hrbq->GetXaxis()->SetTitle("Rbq");
      hrbq->GetXaxis()->SetTitleSize(0.08);
      hrbq->GetXaxis()->SetTitleOffset(1.18);
      hrbq->GetXaxis()->SetTitleFont(42);
 
 

      hrbq->GetYaxis()->SetLabelFont(42);
      hrbq->GetYaxis()->SetLabelSize(0.11);
      hrbq->GetYaxis()->SetTitleSize(0.035);
      hrbq->GetYaxis()->SetTitleFont(42);
      
      double x1 = hrbq->GetMaximum();
      double limit1 = x1 + x1*0.1;
      hrbq->GetYaxis()->SetRangeUser(-limit1, limit1);	
      
      hrbq->GetYaxis()->SetTitleSize(0.08);
	  hrbq->GetYaxis()->SetTitleOffset(0.29);
	  hrbq->GetYaxis()->SetTitleFont(42);
      
      hrbq->GetYaxis()->SetTitle("(Data-Fit)/Error[#sigma]");
      hrbq->SetLineWidth(3);  
		
	  hrbq->Draw();
 
 
 	TF1 *norm1 = new TF1("fa1","0", lower_edge, upper_edge);
	norm1 -> SetLineColor(kRed);
	norm1 -> SetLineStyle(1);
	norm1 -> SetLineWidth(3);
	norm1->Draw("SAME");
	
 
	string Name = file->GetName();
	//std::cout<<Name;
	c3->Print(Form("%s%s%s","~/",Name.c_str(), ".png" ));
 
}	




	
	
	
//top_mass
	
	




	

	
	
	
	

	
	



int main(int argc, char* argv[]){
	TApplication* A =new TApplication("ssss",&argc, argv);
	
	//TFile* F= new TFile("1.root","recreate");
	
	
	
	std::cout<<"Main start";
	mtop_fit p(A->Argv(1));
	std::cout<<"constructor ok";
	p.top_fit(A->Argv(1));
	p.mw_fit();
	p.rbq_fit();

	
	
	
	A->Run();
	return 0;
}
