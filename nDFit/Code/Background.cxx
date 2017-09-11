#include "TString.h"
#include "TMinuitMinimizer.h"
#include "TMatrixD.h"


#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <math.h>
#include <stdio.h>
#include <string>
#include <sys/stat.h>
#include <unistd.h>
#include <iomanip>
#include <omp.h>

#include "Fill3D.h"
#include "TemplateHolder.h"
#include "FCN.h"

using namespace std;


bool Fill3D:: check_file_exist (const std::string& name) {
	
  struct stat buffer;   
  return (stat (name.c_str(), &buffer) == 0); 
}

bool Fill3D::replace(std::string& str, const std::string& from, const std::string& to) {
	size_t start_pos = str.find(from);
    if(start_pos == std::string::npos)
        return false;
    str.replace(start_pos, from.length(), to);
    return true;
}







/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////

#ifndef LJBackground_H_
#define LJBackground_H_



#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <math.h>
#include <stdio.h>
#include <string>
#include <sys/stat.h>
#include <unistd.h>
#include <iomanip>
#include <omp.h>
#include "TMinuitMinimizer.h"


#include "TemplateHolder.h"

using namespace std;


	
class LJBackground{

public:
Fill3D(std::string InputFolder,std::string Channel ,std::vector<double> gMTOP,std::vector<double> gbJSF,std::vector<double> gJSF );
void fillup();
void top_fit();
void save_par(std::string output);
std::vector<TemplateHolder*> gALL;
TMinuitMinimizer *aMinuit; 

bool replace(std::string& str, const std::string& from, const std::string& to); 

bool check_file_exist (const std::string& name); 


};


#endif 






void::Fill3D fillbg(std::string InputFolder, std::string Channel, std::vector<double> gMTOP,std::vector<double> gbJSF,std::vector<double> gJSF  ){
std::vector<std::string> gbackg = {diboson, singleTop, ttV, wjets, zjets};

	for (std::vector<double>::iterator it1=gbackg.begin(); it1!=gbackg.end();it1++)
		for (std::vector<double>::iterator it2=gbJSF.begin(); it2!=gbJSF.end();it2++)
			for (std::vector<double>::iterator it3=gJSF.begin(); it3!=gJSF.end();it3++)
				

				 std::string file_namebg=Form("%s/Output_PlotFactory_2.4.24_JSF_%2.2f_bJSF_%2.2f/HistogramFolder/%s/Merged_ttbar_%4.1f.root",InputFolder.c_str(),*it3,*it2,*it1,Channel.c_str()); 
				 puts(file_namebg.c_str());
			  
			     if (!check_file_exist(file_name)) continue;
				
				gBackground.push.back(new TemplateHolderbg(file_namebg.c_str()*it1,*it3,*it2);
				
				}
			    
			
		

}

/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////


void Fill3D:: bg_fit(){


	
	bgMinuit = new TMinuitMinimizer("bgwow",6*4+6*4+6*4); 
	FCN* Mfbg=new FCN(gBackground);
    
    ROOT::Math::Functor f(Mfbg,&FCN::big_all_fcnbg,72);
    
  
  
	bgMinuit->SetFunction(f);
	Double_t arglist[10];
	Double_t arglistw[10];
	Double_t arglistr[10];
	Int_t ierflg = 0;   
	Int_t ierflgw = 0;   
	Int_t ierflgr = 0;   
	arglist[0] = 1500000;
	arglist[1] = 0.0001;
	arglistw[0] = 1500000;
	arglistw[1] = 0.0001;
	arglistr[0] = 1500000;
	arglistr[1] = 0.0001;
	int q=0;
  
  	bgMinuit->SetLimitedVariable(q, Form("mtop%i",q), 1.10618e+04, 1000.0, 1000,1000000); q++;//zentralparameter
	bgMinuit->SetLimitedVariable(q, Form("mtop%i",q), 0.0, 0.1, -10000,10000); q++;//tompmasse
	bgMinuit->SetLimitedVariable(q, Form("mtop%i",q), 0.0, 0.1, -10000,10000); q++;//bJSf
	bgMinuit->SetLimitedVariable(q, Form("mtop%i",q), 0.0, 0.1, -10000,10000); q++;//JSF
	bgMinuit->SetLimitedVariable(q, Form("mtop%i",q), 1.62284e+02, 10, 10,1000.0); q++;
	bgMinuit->SetLimitedVariable(q, Form("mtop%i",q), 0.0, 0.1, -10000,10000); q++;
	bgMinuit->SetLimitedVariable(q, Form("mtop%i",q), 0.0, 0.1, -10000,10000); q++;
	bgMinuit->SetLimitedVariable(q, Form("mtop%i",q), 0.0, 0.1, -10000,10000); q++;
	bgMinuit->SetLimitedVariable(q, Form("mtop%i",q), 1.28988e+01, 0.1, 0,1000000); q++;
	bgMinuit->SetLimitedVariable(q, Form("mtop%i",q), 0.0, 0.1, -10000,10000); q++;
	bgMinuit->SetLimitedVariable(q, Form("mtop%i",q), 0.0, 0.1, -10000,10000); q++;
	bgMinuit->SetLimitedVariable(q, Form("mtop%i",q), 0.0, 0.1, -10000,10000); q++;
	bgMinuit->SetLimitedVariable(q, Form("mtop%i",q), 4.28790e+05, 0.1, 1000,10000000); q++;
	bgMinuit->SetLimitedVariable(q, Form("mtop%i",q), 0.0, 0.1, -10000,10000); q++;
	bgMinuit->SetLimitedVariable(q, Form("mtop%i",q), 0.0, 0.1, -10000,10000); q++;
	bgMinuit->SetLimitedVariable(q, Form("mtop%i",q), 0.0, 0.1, -10000,10000); q++;
	bgMinuit->SetLimitedVariable(q, Form("mtop%i",q), 1.87218e+02, 0.1, 0,1000000); q++;
	bgMinuit->SetLimitedVariable(q, Form("mtop%i",q), 0.0, 0.1, -10000,10000); q++;
	bgMinuit->SetLimitedVariable(q, Form("mtop%i",q), 0.0, 0.1, -10000,10000); q++;
	bgMinuit->SetLimitedVariable(q, Form("mtop%i",q), 0.0, 0.1, -10000,10000); q++;
	bgMinuit->SetLimitedVariable(q, Form("mtop%i",q), 1.48783e+01, 0.1, 0,1000000); q++;
	bgMinuit->SetLimitedVariable(q, Form("mtop%i",q), 0.0, 0.1, -10000,10000); q++;
	bgMinuit->SetLimitedVariable(q, Form("mtop%i",q), 0.0, 0.1, -10000,10000); q++;
	bgMinuit->SetLimitedVariable(q, Form("mtop%i",q), 0.0, 0.1, -10000,10000); q++;


  
  
    bgMinuit->SetLimitedVariable(q, Form("mw%i",q),  2.48333e+04  , 10., 0.0,10000000); q++;//zentralparameter
	bgMinuit->SetLimitedVariable(q, Form("mw%i",q), -5.36697e+02, 0.1, -100000,10000000); q++;//aMasse
	bgMinuit->SetLimitedVariable(q, Form("mw%i",q),0, 0.1, -100000,10000000); q++;//bJSf
	bgMinuit->SetLimitedVariable(q, Form("mw%i",q), 0, 0.1, -1000,10000000); q++;//JSF
	bgMinuit->SetLimitedVariable(q, Form("mw%i",q), 8.10386e+01   , 10, 40,100); q++;
	bgMinuit->SetLimitedVariable(q, Form("mw%i",q),   2.01724e-02 , 0.1, 0,1000000); q++;
	bgMinuit->SetLimitedVariable(q, Form("mw%i",q), 0, 0.1, 0,1000000); q++;
	bgMinuit->SetLimitedVariable(q, Form("mw%i",q), 0.0, 0.1, -10000,10000); q++;
	bgMinuit->SetLimitedVariable(q, Form("mw%i",q), 6.79767e+00 , 0.1, 1.0,15); q++;
	bgMinuit->SetLimitedVariable(q, Form("mw%i",q), 0.0, 0.1, -10000,10000); q++;
	bgMinuit->SetLimitedVariable(q, Form("mw%i",q), 0.0, 0.1, -10000,10000); q++;
	bgMinuit->SetLimitedVariable(q, Form("mw%i",q), 0.0, 0.1, -10000,10000); q++;
	bgMinuit->SetLimitedVariable(q, Form("mw%i",q), 1.77091e+04, 1.0, 100.0,10000000); q++;//zentralparameter
	bgMinuit->SetLimitedVariable(q, Form("mw%i",q), -3.28105e+02 , 0.1, -1000,1000); q++;//aMasse
	bgMinuit->SetLimitedVariable(q, Form("mw%i",q),0, 0.1, -100000,100000); q++;//bJSf
	bgMinuit->SetLimitedVariable(q, Form("mw%i",q), 0, 0.1, -100000,100000); q++;//JSF
	bgMinuit->SetLimitedVariable(q, Form("mw%i",q), 8.39734e+01, 5, 70,120); q++;
	bgMinuit->SetLimitedVariable(q, Form("mw%i",q), -2.22161e-02, 0.1, 0,1000000); q++;
	bgMinuit->SetLimitedVariable(q, Form("mw%i",q), 0., 0.1, 0,1000000); q++;
	bgMinuit->SetLimitedVariable(q, Form("mw%i",q), 0.0, 0.1, -10000,10000); q++;
	bgMinuit->SetLimitedVariable(q, Form("mw%i",q), 1.70114e+01, 0.01, 10.0,50); q++;
	bgMinuit->SetLimitedVariable(q, Form("mw%i",q), -4.48616e-02, 0.1, 0,1000000); q++;
	bgMinuit->SetLimitedVariable(q, Form("mw%i",q), 0.0, 0.1, -10000,10000); q++;
	bgMinuit->SetLimitedVariable(q, Form("mw%i",q), 0.0, 0.1, -10000,10000); q++;


  
  	bgMinuit->SetLimitedVariable(q, Form("rbq%i",q), 4.51153e+03, 1000.0, 1000,1000000); q++;//zentralparameter
	bgMinuit->SetLimitedVariable(q, Form("rbq%i",q), 0.0, 0.01, -1000000,10000000); q++;//tompmasse
	bgMinuit->SetLimitedVariable(q, Form("rbq%i",q), 0.0, 0.01, -1000000,10000000); q++;//bJSf
	bgMinuit->SetLimitedVariable(q, Form("rbq%i",q), 0.0, 0.01, -1000000,10000000); q++;//JSF
	bgMinuit->SetLimitedVariable(q, Form("rbq%i",q), 1.92480e+00, 10, 0,1000.0); q++;
	bgMinuit->SetLimitedVariable(q, Form("rbq%i",q), 0.0, 0.1, -10000,10000); q++;
	bgMinuit->SetLimitedVariable(q, Form("rbq%i",q), 0.0, 0.1, -10000,10000); q++;
	bgMinuit->SetLimitedVariable(q, Form("rbq%i",q), 0.0, 0.1, -10000,10000); q++;
	bgMinuit->SetLimitedVariable(q, Form("rbq%i",q), 9.04483e-01, 0.1, 0,1000000); q++;
	bgMinuit->SetLimitedVariable(q, Form("rbq%i",q), 0.0, 0.1, -10000,10000); q++;
	bgMinuit->SetLimitedVariable(q, Form("rbq%i",q), 0.0, 0.1, -10000,10000); q++;
	bgMinuit->SetLimitedVariable(q, Form("rbq%i",q), 0.0, 0.1, -10000,10000); q++;
	bgMinuit->SetLimitedVariable(q, Form("rbq%i",q), 6.16592e+03 , 0.1, 1000,10000000); q++;
	bgMinuit->SetLimitedVariable(q, Form("rbq%i",q), 0.0, 0.1, -10000,10000); q++;
	bgMinuit->SetLimitedVariable(q, Form("rbq%i",q), 0.0, 0.1, -10000,10000); q++;
	bgMinuit->SetLimitedVariable(q, Form("rbq%i",q), 0.0, 0.1, -10000,10000); q++;
	bgMinuit->SetLimitedVariable(q, Form("rbq%i",q), 1.17703e+00, 0.1, 0,1000000); q++;
	bgMinuit->SetLimitedVariable(q, Form("rbq%i",q), 0.0, 0.1, -10000,10000); q++;
	bgMinuit->SetLimitedVariable(q, Form("rbq%i",q), 0.0, 0.1, -10000,10000); q++;
	bgMinuit->SetLimitedVariable(q, Form("rbq%i",q), 0.0, 0.1, -10000,10000); q++;
	bgMinuit->SetLimitedVariable(q, Form("rbq%i",q), 5.52227e-01, 0.1, 0,1000000); q++;
	bgMinuit->SetLimitedVariable(q, Form("rbq%i",q), 0.0, 0.1, -10000,10000); q++;
	bgMinuit->SetLimitedVariable(q, Form("rbq%i",q), 0.0, 0.1, -10000,10000); q++;
	bgMinuit->SetLimitedVariable(q, Form("rbq%i",q), 0.0, 0.1, -10000,10000); q++;


  
  
  
  
  	bgMinuit->Minimize();
  	bgMinuit->Hesse();

	

	
}

void Fill3D::save_par(std::string output){
	
	std::string outpath = output + "/parameter.root";
	
	TFile *fend = new TFile(outpath.c_str(),"RECREATE");
	

	std::vector<double> Parameter;
	std::vector<double> Err;
	std::vector<std::string> Name;
	
	const double *P = aMinuit->X();
	const double *E = aMinuit->Errors();
	double cov[200000];
	
	
	aMinuit->GetCovMatrix(cov);
	TMatrixD* CMatrix = new TMatrixD(aMinuit->NDim(), aMinuit->NDim(), cov);
	if (P) 
	
	{
	for(int i =0; i< aMinuit->NDim(); i++){
	
		Err.push_back(E[i]);   //(aMinuit->X())[i]
		Name.push_back(aMinuit->VariableName(i));   //(aMinuit->X())[i]
		};
		
	}
	
	else puts("No parameters");	

		
	
		TTree *t1 = new TTree("t1","Data");
		TBranch *bParameter = 0;
		TBranch *bErrors = 0;
		TBranch *bNames = 0;
		
		t1->Branch("Names",&Name);
		t1->Branch("Parameter",&Parameter);
		t1->Branch("Errors",&Err);
	
	
		t1->Fill();
	

	t1->Write();
	CMatrix->Write();


	
	
	fend->Write();
	fend->Close();
	
	
	
	}
