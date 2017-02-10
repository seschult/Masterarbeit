#include "ZROOT.h"
#include "mtop_fit.h"





int in_out(const char*   in, const char* out)
{
	
	
	Char_t Label[300];
   	
	////////////////////////////////////////////////////////////////////
	//mtop_fit p(A->Argv(1));
	mtop_fit* p= new mtop_fit(in);

	p->top_fit();
	
	p->mw_fit();
	
	p->rbq_fit();



	////////////////////////////////////////////////////////////////////
	TFile *fend = new TFile(out,"RECREATE");
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

	////////////////////////////////////////////////////////////////////	
	
	
	
	
	
	
}	






int main(int argc, char* argv[]){
	TApplication* A =new TApplication("ssss",&argc, argv);



	std::vector<std::string> OUT;
	std::vector<std::string> IN;
	const char *ext=".root";

   TSystemDirectory dir(A->Argv(1), A->Argv(1));
   TList *files = dir.GetListOfFiles();
   
   if (files) {
      TSystemFile *file;
      TString fname;
      TIter next(files);
	  TH1F *h1;
   
      while ((file=(TSystemFile*)next())) {
         fname = file->GetName();
       
         
     
         if (!file->IsDirectory() && fname.EndsWith(ext)) {
           IN.push_back(std::string(A->Argv(1))+"/"+gSystem->BaseName(fname));
           OUT.push_back(std::string(A->Argv(2))+"/out"+gSystem->BaseName(fname));
            

         }
	 }
}








for (int i=0;i<IN.size();i++) in_out(IN[i].c_str(),OUT[i].c_str());
   
   
	

	A->Run();
	return 0;
	}

