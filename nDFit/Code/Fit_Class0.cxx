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
	
	c1 = new TCanvas();
	c1->Divide(3,1);
	
	
	//fFunctions["mtop_f0"]=new TF1("mtop_f0",signal,0,400,3);
	//fFunctions["mtop_f0fsfssfs"]=new TF1("mtop_f0ff",signal,0,400,3);
		
	
}
	


void mtop_fit :: top_fit(){ 
	
	c1->cd(1);
	
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
 
}


void mtop_fit :: mw_fit(){
	
	c1->cd(2);
	
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
 
}
	
void mtop_fit :: rbq_fit(){ 
	
	c1->cd(3);
	
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

	mtop_fit p(A->Argv(1));

	p.top_fit();
	p.mw_fit();
	p.rbq_fit();

	
	A->Run();
	return 0;
	}
