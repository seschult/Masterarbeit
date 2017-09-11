#include "ZROOT.h"
#include "mtop_fit.h"
#include "ZBestNumber.h"





	
	
	
	
	
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
	//TApplication* A =new TApplication("ssss",&argc, argv);

	//
	
	TFile* F= new TFile("1.root","recreate");
	
	//mtop_fit p(A->Argv(1));
	mtop_fit* p= new mtop_fit(argv[1]);

	p->top_fit();
	p->mw_fit();
	p->rbq_fit();
p->fA=10;
	
	F->cd();
	p->Write("FFF");
	ZBestNumber *Q= new ZBestNumber(2.0);
	Q->Write("QQQ");
	F->Close();


TFile* F2= new TFile("1.root","read");
ZBestNumber *Z=(ZBestNumber*)(F2->Get("QQQ"));
Z->Print();
printf("%f\n",Z->fV);
mtop_fit* q=(mtop_fit*)(F2->Get("FFF"));
q->Print();
if (q) if (q->c1) { puts("OK"); printf("%f\n",q->fA);}


	//A->Run();
	return 0;
	}
