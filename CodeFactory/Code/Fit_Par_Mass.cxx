  #include "ZROOT.h"
  
  
  //Tapplication wie FitClass3
  
int main(){
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
    fChain->Add("~/Errors.root/t1");
      
   

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
   
   
 
 
 
   fChain->GetEntry(0)

////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////

TGraph ABC
TGraph errors
//apl option draw
for(int n = 0; n < fChain->GetEntriesFast(); ++n) {
    fChain->GetEntry(n);
    
    ABC -> SetPoint(N, masse, Parameterxy)
    }
   // then do something with various vector[i]'s


return 0;
}
