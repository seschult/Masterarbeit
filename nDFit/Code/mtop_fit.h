#include "ZROOT.h"




class mtop_fit: public TObject
	
	//members of class
	
	    
	{
	public: 
	
	
	////////////////////////////////////////////////////////////////////
	//For Save
	double par_mtop[20];
	double err_mtop[20];
	int npar_mtop;
	
	double par_mw[20];
	double err_mw[20];
	int npar_mw;
	
	double par_rbq[20];
	double err_rbq[20];
	int npar_rbq;
	////////////////////////////////////////////////////////////////////
	
	
	mtop_fit(); //constructor
	void top_fit(); //constructor
	void mw_fit(); //constructor
	void rbq_fit(); //constructor
	void save(const char* out);
	
	
	mtop_fit(const char *File); //constructor
	 //constructor
	

	int fA;//
	TH1F *h1; //
	TH1F *h2;//
	TH1F *h3;//
	TFile *file;//
	
	TF1 *f00;//
	TF1 *f01;//
	TF1 *f02;//
	TF1 *f13;//
	TF1 *f14;//
	TF1 *f25;//
	TF1 *f26;//
	TF1 *f27;//
	
	
	TF1 *ftop;//
	TF1 *fw;//
	TF1 *frbq;//


	TCanvas *c1;//
	TCanvas *cw;//
	TCanvas *c3;//
	
	
	
	ClassDef(mtop_fit,0)
	
	
};
