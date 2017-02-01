#include "ZROOT.h"
#include "mtop_fit.h"
#include "ZBestNumber.h"






int main(int argc, char* argv[]){
	TApplication* A =new TApplication("ssss",&argc, argv);




   
   
	Char_t Label[300];
   	
	////////////////////////////////////////////////////////////////////
	//mtop_fit p(A->Argv(1));
	mtop_fit* p= new mtop_fit(A->Argv(1));

	p->top_fit();
	
	p->mw_fit();
	
	p->rbq_fit();



	////////////////////////////////////////////////////////////////////
	TFile *fend = new TFile(A->Argv(2),"RECREATE");
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
	

	A->Run();
	return 0;
	}
