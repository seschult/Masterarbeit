#include "ZROOT.h"
#include "mtop_fit.h"





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





ClassImp(mtop_fit)


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
	
	
	
	 
	ftop->SetParameters(10, 170, 0.5, 400, 175, 8, 400, 120 ,9);
	ftop->SetParLimits(3,0.0,1000000);
	ftop->SetParLimits(6,0.0,1000000);
    ftop->SetParLimits(7,100.0,1000000);
	  
	fmw->SetParameters(10, 80, 1, 10, 82, 10); 
	
	
	frbq->SetParameters(10, 1, 2, 10, 2, 1, 10, 2 ,1);
	//fit->SetParameters(10, 1, 1, 1, 2, 1, 10, 2 ,15);
	//frbq->SetParLimits(0,0.0,1000000);
    //frbq->SetParLimits(3,0.0,1000000);
    //frbq->SetParLimits(6,0.0,1000000);
	
	
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
			}else continue;
		
		};
		htop->Draw();
	
	
	
	
	
 
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
			}else continue;
		
		};
		hw->Draw();
	
 
 
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
			}else continue;
		
		};
		hrbq->Draw();
 
}	
	
	
	
	

	
	
	
	
	
