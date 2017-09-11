#ifdef __CINT__
#else
#include "ZBestNumber.h"
#include "TFile.h"
int main()
#endif
{

TFile* q= new TFile("1111.root","recreate");


    ZBestNumber l=ZBestNumber(0.0016);

l.Write();

q->Close();


}
//0.095±0.006+0.007−0.006
