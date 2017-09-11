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







class mtop_fit{
	
	//members of class
	
	public:
	
	
	mtop_fit(const char *Merge){
		fit_mtop(Merge);
		fit_mw(Merge);
		fit_Rbq(Merge);
		
		
		}; //constructor
	
	 //constructor
	
	//Gauss


	

};





void mtop_fit:: fit_mtop (const char *File){
	
	
	
	
	 TF1 *f5 = new TF1("f5",signal,0,400,3);
	 TF1 *f1 = new TF1("f1",landau,0,400,3);
	 TF1 *f2 = new TF1("f2",landau_n,-100,400,3);
	 TF1 *f0 = new TF1("f0",fit_mtop,0,400,9);
	 
	f0->SetParameters(10, 150, 0.5, 400, 170, 8, 100, 50 ,100);  
	//f0->SetParameters(10, 150, 2, 400, 170, 8, 100, 50 ,100);
    
    
    puts(File);
 // Open file and histogram   
	TFile *file = new TFile(File,"r");
  
	TH1F *h1 = (TH1F*)file->Get("hist_klf_mtop_param");
	
	
 //Fit of the histogram   
	TCanvas *c1 = new TCanvas;
	h1->Fit("f0","","",0,500);
	h1->Draw();
	double Chi_2 = f0->GetChisquare();


 //Get fit parameters 
 double p[9];
 
 for(int i = 0; i < 8; i++){
	 
	 p[i] = f0->GetParameter(i);
 }
 
  
  
  
 
 
 //Draw functions
 
 //
  
  
  f5->SetParameters(p[0],p[1],p[2]);
  f5->SetLineColor(3);
  f5->Draw("SAME");
 
  f1->SetParameters(p[3],p[4],p[5]);
  f1->SetLineColor(7);
  f1->Draw("SAME");
  

  f2->SetParameters(1000000*p[6],p[7],p[8]);
  f2->SetLineColor(9);
  

  f2->Draw("SAME");
  
  
 TLegend* leg = new TLegend(0.4,0.7,0.78,0.9);
	
	leg->AddEntry("f0","Fit","l");
	leg->AddEntry("f5","Gauss","l");
	leg->AddEntry("f1","Landau","l");
	leg->AddEntry("f2","Inv_Landau","l");
	leg->AddEntry("(TObject*)0", "Chi_2 = 5100.06","");

	leg->Draw();
	
	
			 
		}

	
void mtop_fit:: fit_mw(const char *File){
	
	
	   TF1 *f0 = new TF1("f0",signal,0,400,3);
	 TF1 *f1 = new TF1("f1",signal,0,400,3);
	// TF1 *f2 = new TF1("f2",landau_n,-100,400,3);
	 TF1 *fit = new TF1("fit",fit_mw,0,400,9);
	 
	fit->SetParameters(10, 80, 1, 10, 82, 1);
	//fit->SetParameters(10, 79, 2, 10, 82, 2);
    
 // Open file and histogram   
	TFile *file = new TFile(File);
  
	TH1F *h1 = (TH1F*)file->Get("hist_klf_original_Whad_m");
	
	
 //Fit of the histogram   
	TCanvas *c1 = new TCanvas;
	h1->Fit("fit","","",0,450);
	h1->Draw();
 
 //Get fit parameters 
 double p[6];
 
 for(int i = 0; i < 6; i++){
	 
 p[i] = fit->GetParameter(i);
 
 }
 
  
  
  
 
 
 //Draw functions
 
 //TCanvas *c2 = new TCanvas;
  
  
  f0->SetParameters(p[0],p[1],p[2]);
  f0->SetLineColor(3);
  f0->Draw("SAME");
 
  f1->SetParameters(p[3],p[4],p[5]);
  f1->SetLineColor(4);
  f1->Draw("SAME");
  
 
  //f2->SetParameters(p[6],p[7],p[8]);
  //f2->SetLineColor(5);
  //f2->Draw("SAME");
  
  
 leg = new TLegend(0.4,0.7,0.78,0.9);
	
	leg->AddEntry("fit","Fit","l");
	leg->AddEntry("f0","Gauss","l");
	leg->AddEntry("f1","Gauss_2","l");
	

	leg->AddEntry("(TObject*)0", "Chi_2 = 1","");

	leg->Draw("SAME");
  
  
  
  //gSystem->ProcessEvents();
  //TImage *img = TImage::Create();
  //img->FromPad(c);
  //img->WriteImage("canvas.png");*/
 
}
	
		
void mtop_fit:: fit_Rbq(const char *File){
	
 //Definition of functions

     TF1 *f0 = new TF1("f0",signal,0,5,3);
	 TF1 *f1 = new TF1("f1",signal,0,5,3);
	 TF1 *f2 = new TF1("f2",landau_n,-1,5,3);
	 TF1 *fit = new TF1("fit",fit_r,0,5,9);
	 
	//fit->SetParameters(10, 0, 2, 10, 0, 1, 10, 2 ,3);
	fit->SetParameters(10, 0, 2, 10, -1, 1, 10, 2 ,3);
	
    
 // Open file and histogram   
	TFile *file = new TFile(File);
  
	TH1F *h1 = (TH1F*)file->Get("hist_klf_original_Rbq_reco");
	
	
 //Fit of the histogram   
	TCanvas *c1 = new TCanvas;
	h1->Fit("fit","","",0,5);
	h1->Draw();
 
 //Get fit parameters 
 double p[9];
 
 for(int i = 0; i < 8; i++){
	 
	 p[i] = fit->GetParameter(i);
 }
 
  
  
  
 
 
 //Draw functions
 
 //TCanvas *c2 = new TCanvas;
  
  
  f0->SetParameters(p[0],p[1],p[2]);
  f0->SetLineColor(3);
  f0->Draw("SAME");
 
  f1->SetParameters(p[3],p[4],p[5]);
  f1->SetLineColor(4);
  f1->Draw("SAME");
  
 
  f2->SetParameters(p[6],p[7],p[8]);
  f2->SetLineColor(5);
  f2->Draw("SAME");
  
  
  
   
   
     
 leg = new TLegend(0.4,0.7,0.78,0.9);
	
	leg->AddEntry("fit","Fit","l");
	
	


	
	leg->AddEntry("f0","Gauss","l");
	leg->AddEntry("f1","Gauss2","l");
	leg->AddEntry("f2","Landau","l");
leg->AddEntry("(TObject*)0", "Chi_2 = 9596.44","");
	leg->Draw("SAME");
  
  
  
  //gSystem->ProcessEvents();
  //TImage *img = TImage::Create();
  //img->FromPad(c);
  //img->WriteImage("canvas.png");*/
 
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

	mtop_fit top(A->Argv(1));
	
	A->Run();
	return 0;
	}
