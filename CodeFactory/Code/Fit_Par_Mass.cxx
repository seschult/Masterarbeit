  #include "ZROOT.h"
  
  
  //Tapplication wie FitClass3
  
int main(int argc, char* argv[]){
	TApplication* A =new TApplication("ssss",&argc, argv);
////////////////////////////////////////////////////////////////////////
//Get Entries
////////////////////////////////////////////////////////////////////////  
  // Declaration of leaf types
   Char_t          Label[7];
   Int_t           fit_mtop_npar;
   Double_t        fit_mtoppar[9];   //[fit_mtop_npar]
   Double_t        fit_mtoperr[9];   //[fit_mtop_npar]
   Int_t           fit_mW_npar;
   Double_t        fit_mWpar[6];   //[fit_mW_npar]
   Double_t        fit_mWerr[6];   //[fit_mW_npar]
   Int_t           fit_rbq_npar;
   Double_t        fit_rbq[9];   //[fit_rbq_npar]
   Double_t        fit_rbqerr[9];   //[fit_rbq_npar]

   // List of branches
   TBranch        *b_label;   //!
   TBranch        *b_fit_mtop_npar;   //!
   TBranch        *b_fit_mtoppar;   //!
   TBranch        *b_fit_mtoperr;   //!
   TBranch        *b_fit_mW_npar;   //!
   TBranch        *b_fit_mWpar;   //!
   TBranch        *b_fit_mWerr;   //!
   TBranch        *b_fit_rbq_npar;   //!
   TBranch        *b_fit_rbq;   //!
   TBranch        *b_fit_rbqerr;   //!
   
   
	TChain * fChain = new TChain("t1","");
    fChain->Add(A->Argv(1));
      
   

   fChain->SetBranchAddress("Label", Label, &b_label);
   fChain->SetBranchAddress("fit_mtop_npar", &fit_mtop_npar, &b_fit_mtop_npar);
   fChain->SetBranchAddress("fit_mtoppar", fit_mtoppar, &b_fit_mtoppar);
   fChain->SetBranchAddress("fit_mtoperr", fit_mtoperr, &b_fit_mtoperr);
   fChain->SetBranchAddress("fit_mW_npar", &fit_mW_npar, &b_fit_mW_npar);
   fChain->SetBranchAddress("fit_mWpar", fit_mWpar, &b_fit_mWpar);
   fChain->SetBranchAddress("fit_mWerr", fit_mWerr, &b_fit_mWerr);
   fChain->SetBranchAddress("fit_rbq_npar", &fit_rbq_npar, &b_fit_rbq_npar);
   fChain->SetBranchAddress("fit_rbq", fit_rbq, &b_fit_rbq);
   fChain->SetBranchAddress("fit_rbqerr", fit_rbqerr, &b_fit_rbqerr);
   
   
 
 
 
   fChain->GetEntry(0);

////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////

int n = fChain->GetEntriesFast();



//n =nummber of entries

//mtop
Double_t x[n+10], mtopp0[n+10];

Double_t mtop[20][20];
Double_t mtoperr[20][20];






Double_t mtoperr1[n+10], massy[n+10];  // x=mass y= Parameter

x[0] = 170;//Fix me
x[1] = 171.5;
x[2] = 173.5;
x[3] = 175;
x[4] = 177;

    



//apl option draw
for(int i = 0; i < fChain->GetEntriesFast(); ++i) {
	
		fChain->GetEntry(i);
		
		
		for(int j = 0; j < 9; ++j){
			
			mtop[i][j] = fit_mtoppar[j];
			mtoperr[i][j] = fit_mtoperr[j];
			
			};	
		
		
		
    };
    
    
    
    
   TGraphErrors *gr[20];
 
   
   for(int j = 0; j < 9; ++j){
	   
   gr[j] = new TGraphErrors(n,x,mtop[j],0,mtoperr[j]);
   
	};

	TCanvas *c[20];
   for(int j = 0; j < 9; ++j){
   c[j] = new TCanvas(Form("%i",j),Form("%i",j),200,10,700,500);
   gr[j]->Draw("AC*");
   };



A->Run();


return 0;
}
