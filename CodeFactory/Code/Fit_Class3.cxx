#include "ZROOT.h"
#include "mtop_fit.h"
#include "ZBestNumber.h"






int main(int argc, char* argv[]){
	TApplication* A =new TApplication("ssss",&argc, argv);

	//
	
	
	
	//mtop_fit p(A->Argv(1));
	mtop_fit* p= new mtop_fit(A->Argv(1));

	p->top_fit();
	p->mw_fit();
	p->rbq_fit();

	
	

	A->Run();
	return 0;
	}
