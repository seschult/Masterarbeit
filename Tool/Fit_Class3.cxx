#include "ZROOT.h"
#include "mtop_fit.h"





	
	
	
	
	
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
	mtop_fit p(argv[1]);

	p.top_fit();
	p.mw_fit();
	p.rbq_fit();

	
	F->cd();
	p.Write("FFF");
	F->Close();
	//A->Run();
	return 0;
	}
