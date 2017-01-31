#include "ZROOT.h"




class mtop_fit: public TObject
	
	//members of class
	
	    
	{
	public: 
	mtop_fit(); //constructor
	void top_fit(); //constructor
	void mw_fit(); //constructor
	void rbq_fit(); //constructor
	
	mtop_fit(const char *File); //constructor
	 //constructor
	
	//Gauss
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
	TF1 *fmw;//
	TF1 *frbq;//

//std::map<std::string,TF1*>    fFunctions;
	TCanvas *c1;//
	TCanvas *c2;//
	TCanvas *c3;//
	    ClassDef(mtop_fit,0)
	
	
};
