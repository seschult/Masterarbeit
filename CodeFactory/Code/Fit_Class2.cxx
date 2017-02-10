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
	void top_fit(); //constructor
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
	
	
	
	 
	//ftop->SetParameters(10, 170, 0.5, 400, 175, 8, 400, 120 ,9);
	ftop->SetParameters(1318, 167, -1.05, 4340, 182, 1.27, 54200, 142 ,9.56);
	ftop->SetParLimits(3,0.0,1000000);
	ftop->SetParLimits(6,0.0,1000000);
    ftop->SetParLimits(7,100.0,1000000);
	  
	fmw->SetParameters(10, 80, 1, 10, 82, 10); 
	fmw->SetParameters(1742, 78, 6.8, 1616, 79, 19); 
	
	
	
	//frbq->SetParameters(10, 1, 2, 10, 2, 1, 10, 2 ,1);
	frbq->SetParameters(600, 1.5, 0.7, 1000, 1, 0.5, 4525, 0.95 ,0.3);
	
	//fit->SetParameters(10, 1, 1, 1, 2, 1, 10, 2 ,15);
	frbq->SetParLimits(0,0.0,1000000);
    frbq->SetParLimits(3,0.0,1000000);
    frbq->SetParLimits(6,0.0,1000000);
	
	
	//f0->SetParameters(10, 150, 2, 400, 170, 8, 100, 50 ,100);
    
    
    puts(File);
	// Open file and histogram   
	file = new TFile(File,"r");
  
	h1 = (TH1F*)file->Get("hist_klf_mtop_window");
	h2 = (TH1F*)file->Get("hist_klf_window_Whad_m");
	h3 = (TH1F*)file->Get("hist_klf_window_Rbq_reco");
	
	/*double with = 2000;
	double hight = 500;
	c1 = new TCanvas("c","c", with, hight);
	c1->Divide(3,2);
	gStyle->SetOptStat(0); //New
	*/
	
	//fFunctions["mtop_f0"]=new TF1("mtop_f0",signal,0,400,3);
	//fFunctions["mtop_f0fsfssfs"]=new TF1("mtop_f0ff",signal,0,400,3);
		
	
}
	


void mtop_fit :: top_fit(){ 
	
	//c1->cd(1);
	double with = 600;
	double hight = 700;
	c1 = new TCanvas("c1","c1", with, hight);
	gStyle->SetOptStat(0);
	
	TPad *pad1 = new TPad("pad1", "pad1", 0, 0.32, 1, 1.0);
	pad1->SetBottomMargin(5); // Upper and lower plot are joined
    pad1->SetGridx();         // Vertical grid
    pad1->Draw();             // Draw the upper pad: pad1
    pad1->cd();               // pad1 becomes the current pad
	
	h1->GetXaxis()->SetLabelFont(63);
	h1->GetXaxis()->SetLabelSize(14); // labels will be 14 pixels
	h1->GetYaxis()->SetLabelFont(63);
	h1->GetYaxis()->SetLabelSize(14);
	
	//Fit of the histogram   
	h1->Fit("ftop","RI","",130,210);
	h1->Draw();
	
	

	//Get fit parameters 
	double p[9];
 
	for(int i = 0; i < 9; i++){
	 
		p[i] = ftop->GetParameter(i);
		};
 
	//Draw functions
	f00->SetParameters(p[0],p[1],p[2]);
	f00->SetLineColor(3);
	f00->Draw("SAME");
 
	f01->SetParameters(p[3],p[4],p[5]);
	f01->SetLineColor(7);
	f01->Draw("SAME");
  

	f02->SetParameters(p[6],p[7],p[8]);
	f02->SetLineColor(9);
	f02->Draw("SAME");
	
	double chi2_0 = ftop->GetChisquare();
	double NDF_0 = ftop->GetNDF();
	double Prob = ftop->GetProb();

	//Output chi2/NDF/Prob
	////////////////////////////////////////////////////////////////////
	std::stringstream oss_Sep;
	oss_Sep   << setprecision(3) << chi2_0;
	
	TLatex l00;
	l00.SetTextAlign(9);
	l00.SetTextSize(0.038);
	l00.SetNDC();
	l00.DrawLatex(0.1, 0.92, ("chi^2: " + oss_Sep.str()).c_str());
	
	
	std::stringstream oss_NDF_0;
	oss_NDF_0 << setprecision(3) << NDF_0;
	
	TLatex l01;
	l01.SetTextAlign(9);
	l01.SetTextSize(0.038);
	l01.SetNDC();
	l01.DrawLatex(0.1, 0.95, ("NDF: " + oss_NDF_0.str()).c_str());
	
	std::stringstream oss_Prob;
	oss_Prob << setprecision(3) << Prob;
	
	TLatex l03;
	l03.SetTextAlign(9);
	l03.SetTextSize(0.038);
	l03.SetNDC();
	l03.DrawLatex(0.1, 0.88, ("Prob: " + oss_Prob.str()).c_str());
	////////////////////////////////////////////////////////////////////
	
	    
   TLegend *leg = new TLegend(0.6,0.7,0.78,0.9);
   //leg->SetHeader("The Legend Title","C"); // option "C" allows to center the header
   leg->AddEntry(h1,"Data","lep");
   leg->AddEntry(ftop,"Fit","l");
   leg->AddEntry(f00,"Gauss","l");
   leg->AddEntry(f01,"Landau","l");
   leg->AddEntry(f02,"Landau_n","l");
   leg->Draw();
   
	
   c1->cd(); 
   TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
   pad2->SetTopMargin(0);
   pad2->SetBottomMargin(0.2);
   pad2->SetGridx(); // vertical grid
   pad2->Draw();
   pad2->cd();   
       // pad2 becomes the current pad
       
 
     
       
	//Calculate the diffrence between Fitfunction and Histogramm h1
	
	double bin_tot = h1 -> GetSize();  
    double start = h1->GetXaxis()->GetBinCenter(0);
    double end = h1->GetXaxis()->GetBinCenter(bin_tot);
    double err;
    double div;
    
    
	TH1F *htop = new TH1F("htop", "dif", bin_tot, 120, 220);
	htop->GetXaxis()->SetLabelFont(63);
	htop->GetXaxis()->SetLabelSize(14); // labels will be 14 pixels
	htop->GetYaxis()->SetLabelFont(63);
	htop->GetYaxis()->SetLabelSize(14);
	
	//TF1 *func = h1->GetFunction("ftop");
	 
 	for (int bin=2; bin<=bin_tot;bin++) {
		
		double x = h1->GetXaxis()->GetBinCenter(bin);
		double fval = ftop->Eval(x);
		double dif = h1->GetBinContent(bin);
		double sub = h1->GetBinContent(bin)-fval;
    
		if(sub > -100){
			htop->SetBinContent(bin, sub);
			err = htop->GetBinError(bin);
			div = sub/err;
			htop->SetBinContent(bin, div);
			
			
			}else continue;
		
		};
		
	
	int nbins = htop -> GetNbinsX();

	float lower_edge  = htop -> GetBinLowEdge(1);
	float bin_width   = htop -> GetBinWidth(1);
	float number_bins = htop -> GetNbinsX();
	float upper_edge = htop -> GetBinLowEdge(number_bins) + htop->GetBinWidth(number_bins);
       
       
       
	TF1 *norm1 = new TF1("fa1","1", lower_edge, upper_edge);
	norm1 -> SetLineColor(kRed);
	norm1 -> SetLineStyle(1);
	norm1 -> SetLineWidth(2);
	norm1->Draw();
	
	htop->Draw("Same");
	
	
	
 
}


void mtop_fit :: mw_fit(){
	double with = 600;
	double hight = 700;
	c2 = new TCanvas("c2","c2", with, hight);
	gStyle->SetOptStat(0);
	//c1->cd(2);
	TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
	pad1->SetBottomMargin(2); // Upper and lower plot are joined
    pad1->SetGridx();         // Vertical grid
    pad1->Draw();             // Draw the upper pad: pad1
    pad1->cd();               // pad1 becomes the current pad
	
	
	h2->GetXaxis()->SetLabelFont(63);
	h2->GetXaxis()->SetLabelSize(14); // labels will be 14 pixels
	h2->GetYaxis()->SetLabelFont(63);
	h2->GetYaxis()->SetLabelSize(14);
	
	//Fit of the histogram   
	h2->Fit("fmw","","",56,110);
	h2->Draw();
	
	//Get fit parameters 
	double p[6];
 
	for(int i = 0; i < 6; i++){
	 
		p[i] = fmw->GetParameter(i);
		};
 
	//Draw functions
	f13->SetParameters(p[0],p[1],p[2]);
	f13->SetLineColor(3);
	f13->Draw("SAME");
 
	f14->SetParameters(p[3],p[4],p[5]);
	f14->SetLineColor(7);
	f14->Draw("SAME");


	////////////////////////////////////////////////////////////////////
	double chi2_1 = fmw->GetChisquare();
	double NDF_1 = fmw->GetNDF();
	double Prob = fmw->GetProb();
	
	std::stringstream oss_Sep1;
	oss_Sep1   << setprecision(3) << chi2_1;
	
	TLatex l10;
	l10.SetTextAlign(9);
	l10.SetTextSize(0.038);
	l10.SetNDC();
	l10.DrawLatex(0.1, 0.92, ("chi^2: " + oss_Sep1.str()).c_str());
	
	
	std::stringstream oss_NDF1;
	oss_NDF1 << setprecision(3) << NDF_1;
	
	TLatex l11;
	l11.SetTextAlign(9);
	l11.SetTextSize(0.038);
	l11.SetNDC();
	l11.DrawLatex(0.1, 0.95, ("NDF: " + oss_NDF1.str()).c_str());
	
	
	std::stringstream oss_Prob;
	oss_Prob << setprecision(3) << Prob;
	
	TLatex l12;
	l12.SetTextAlign(9);
	l12.SetTextSize(0.038);
	l12.SetNDC();
	l12.DrawLatex(0.1, 0.88, ("Prob: " + oss_Prob.str()).c_str());
	
	
	////////////////////////////////////////////////////////////////////
	
	
	
 
	TLegend* leg = new TLegend(0.6,0.7,0.78,0.9);
   //leg->SetHeader("The Legend Title","C"); // option "C" allows to center the header
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
			err = hw->GetBinError(bin);
			div = sub/err;
			hw->SetBinContent(bin, div);
			}else continue;
		
		};
	
	
 	int nbins = hw -> GetNbinsX();

	float lower_edge  = hw -> GetBinLowEdge(1);
	float bin_width   = hw -> GetBinWidth(1);
	float number_bins = hw -> GetNbinsX();
	float upper_edge = hw -> GetBinLowEdge(number_bins) + hw->GetBinWidth(number_bins);
       
       
       
	TF1 *norm1 = new TF1("fa1","1", lower_edge, upper_edge);
	norm1 -> SetLineColor(kRed);
	norm1 -> SetLineStyle(1);
	norm1 -> SetLineWidth(2);
	norm1->Draw();
	
	hw->Draw("Same");
	
	
	
 
}
	
void mtop_fit :: rbq_fit(){ 
	double with = 600;
	double hight = 700;
	c3 = new TCanvas("c3","c3", with, hight);
	gStyle->SetOptStat(0);
	//c1->cd(3);
	TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
	pad1->SetBottomMargin(2); // Upper and lower plot are joined
    pad1->SetGridx();         // Vertical grid
    pad1->Draw();             // Draw the upper pad: pad1
    pad1->cd();               // pad1 becomes the current pad
	
	
    h3->GetXaxis()->SetLabelFont(63);
	h3->GetXaxis()->SetLabelSize(14); // labels will be 14 pixels
	h3->GetYaxis()->SetLabelFont(63);
	h3->GetYaxis()->SetLabelSize(14);
	
	//Fit of the histogram   
	h3->Fit("frbq","I","",0.3,3);
	h3->Fit("frbq","I","",0.3,3);
	h3->Fit("frbq","I","",0.3,3);
	h3->Fit("frbq","I","",0.3,3);
	h3->Fit("frbq","I","",0.3,3);
	
	h3->Draw();
	


	//Get fit parameters 
	double p[9];
 
	for(int i = 0; i < 9; i++){
	 
		p[i] = frbq->GetParameter(i);
		};
 
	//Draw functions
	f25->SetParameters(p[0],p[1],p[2]);
	f25->SetLineColor(3);
	f25->Draw("SAME");
 
	f26->SetParameters(p[3],p[4],p[5]);
	f26->SetLineColor(7);
	f26->Draw("SAME");
  

	f27->SetParameters(p[6],p[7],p[8]);
	f27->SetLineColor(9);
	f27->Draw("SAME");
	
	//////////////////////////////////////////////////////////////////////
	double chi2_2 = frbq->GetChisquare();
	double NDF_2 = frbq->GetNDF();
	double Prob = frbq->GetProb();
	
	std::stringstream oss_Sep2;
	oss_Sep2   << setprecision(3) << chi2_2;
	
	TLatex l20;
	l20.SetTextAlign(9);
	l20.SetTextSize(0.038);
	l20.SetNDC();
	l20.DrawLatex(0.1, 0.92, ("chi^2: " + oss_Sep2.str()).c_str());
	
	
	std::stringstream oss_NDF2;
	oss_NDF2 << setprecision(3) << NDF_2;
	
	TLatex l21;
	l21.SetTextAlign(9);
	l21.SetTextSize(0.038);
	l21.SetNDC();
	l21.DrawLatex(0.1, 0.95, ("NDF: " + oss_NDF2.str()).c_str());
	
	
	std::stringstream oss_Prob;
	oss_Prob << setprecision(3) << Prob;
	
	TLatex l22;
	l22.SetTextAlign(9);
	l22.SetTextSize(0.038);
	l22.SetNDC();
	l22.DrawLatex(0.1, 0.88, ("Prob: " + oss_Prob.str()).c_str());
	
	//////////////////////////////////////////////////////////////////////
 
	TLegend* leg = new TLegend(0.6,0.7,0.78,0.9);
   //leg->SetHeader("The Legend Title","C"); // option "C" allows to center the header
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
	hrbq->GetXaxis()->SetLabelSize(14); // labels will be 14 pixels
	hrbq->GetYaxis()->SetLabelFont(63);
	
	  hrbq->SetTitle("Dif Data-Fit");
 
 	for (int bin=2; bin<=bin_tot;bin++) {
		
		double x = h3->GetXaxis()->GetBinCenter(bin);
		double fval = frbq->Eval(x);
		double dif = h3->GetBinContent(bin);
		double sub = h3->GetBinContent(bin)-fval;
    
		if(sub > -100 && dif !=0){
			hrbq->SetBinContent(bin, sub);
			err = hrbq->GetBinError(bin);
			div = sub/err;
			hrbq->SetBinContent(bin, div);
			
			}else continue;
		
		};

		
		
 	int nbins = hrbq -> GetNbinsX();

	float lower_edge  = hrbq -> GetBinLowEdge(1);
	float bin_width   = hrbq -> GetBinWidth(1);
	float number_bins = hrbq -> GetNbinsX();
	float upper_edge = hrbq -> GetBinLowEdge(number_bins) + hrbq->GetBinWidth(number_bins);
       
       
       
	TF1 *norm1 = new TF1("fa1","1", lower_edge, upper_edge);
	norm1 -> SetLineColor(kRed);
	norm1 -> SetLineStyle(1);
	norm1 -> SetLineWidth(2);
	norm1->Draw();
	
	hrbq->Draw("Same");
	
	
	
}	
	
	
	
	

	
	
	
	
	
	/*
	 * n order to avoid that, a class can include a special function called its constructor, which is automatically called whenever a new object of this class is created, allowing the class to initialize member variables or allocate storage.

This constructor function is declared just like a regular member function, but with a name that matches the class name and without any return type; not even void.

The Rectangle class above can easily be improved by implementing a constructor:
	 * 
	 * */
	
	
	
	
	 //obj of class like int a
	/*After declaraiton of the class, any of the public members of object rect can be accessed s if they were a normal function 
	 * .between object ame and member name: rect.set_values();
myarea = rect.area(); */




int main(int argc, char* argv[]){
	TApplication* A =new TApplication("ssss",&argc, argv);
	
	//TFile* F= new TFile("1.root","recreate");
	
	mtop_fit p(A->Argv(1));

	p.top_fit();
	p.mw_fit();
	p.rbq_fit();

	
	
	
	
	/*TFile *fend = new TFile(A->Argv(2),"RECREATE");
	TTree *t1 = new TTree("t1","Data");
	

	t1->Branch("Label",Label,"label/C");
	
	t1->Branch("fit_mtop_npar",&p->npar_mtop,"fit_mtop_npar/I");
	t1->Branch("fit_mtoppar",&p->par_mtop,"fit_mtoppar[fit_mtop_npar]/D");
	t1->Branch("fit_mtoperr",&p->err_mtop,"fit_mtoperr[fit_mtop_npar]/D");
	
	t1->Branch("fit_mW_npar",&p->npar_mw,"fit_mW_npar/I");
	t1->Branch("fit_mWpar",&p->par_mw,"fit_mWpar[fit_mW_npar]/D");
	t1->Branch("fit_mWerr",&p->err_mw,"fit_mWerr[fit_mW_npar]/D");
	
	t1->Branch("fit_rbq_npar",&p->npar_rbq,"fit_rbq_npar/I");
	t1->Branch("fit_rbq",&p->par_rbq,"fit_rbqpar[fit_rbq_npar]/D");
	t1->Branch("fit_rbqerr",&p->err_rbq,"fit_rbqerr[fit_rbq_npar]/D");
	
	//Nation[0]='';
	//sprintf(&(Nation[0]),"%s",A->Argv(1));
	t1->Fill();   
	fend->Write();
	fend->Close();




	
	
	//F->cd();
	//p.Write("FFF");
	//F->Close();
	*/
	
	
	
	A->Run();
	return 0;
	}
