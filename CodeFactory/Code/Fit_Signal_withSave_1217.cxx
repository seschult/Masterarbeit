#include "ZROOT.h"
#include "mtop_fit.h"
#include "ZBestNumber.h"






int main(int argc, char* argv[]){
	TApplication* A =new TApplication("ssss",&argc, argv);



const Int_t kMaxTrack = 500;
   Int_t ntrack;
   Int_t stat[kMaxTrack];
   Int_t sign[kMaxTrack];
   Float_t px[kMaxTrack];
   Float_t py[kMaxTrack];
   Float_t pz[kMaxTrack];
   Float_t pt[kMaxTrack];
   Float_t zv[kMaxTrack];
   Float_t chi2[kMaxTrack];
   Double_t sumstat;
   
   
   Char_t Nation[300];
   
	
	
	
	ntrack = 5;
	
	
	//mtop_fit p(A->Argv(1));
	mtop_fit* p= new mtop_fit(A->Argv(1));

	p->top_fit();
	
	p->mw_fit();
	
	p->rbq_fit();

	TFile *fend = new TFile(A->Argv(2),"RECREATE");
	TTree *t1 = new TTree("t1","Data");
	
	

	
	
	t1->Branch("ntrack",&ntrack,"ntrack/I");
	t1->Branch("stat",stat,"stat[ntrack]/I");
	t1->Branch("sign",sign,"sign[ntrack]/I");
	
	t1->Branch("px",px,"px[ntrack]/F");
	t1->Branch("py",py,"py[ntrack]/F");
	t1->Branch("pz",pz,"pz[ntrack]/F");
	t1->Branch("Nation",Nation,"Nation/C");
	
	
	
	//
	t1->Branch("fit_mW_npar",&p->npar_mw,"fit_mW_npar/I");
	t1->Branch("fit_mWpar",p->par_mw,"fit_mWpar[fit_mW_npar]/D");
	
	//Nation[0]='';
	sprintf(&(Nation[0]),"%s",A->Argv(1));
	t1->Fill();   
	fend->Write();
	fend->Close();

		
	

	A->Run();
	return 0;
	}
