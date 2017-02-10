
#include <string>
#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TH1F.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TCanvas.h"
#include "TImage.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TAxis.h"
#include "TTree.h"
	#include "TSystem.h"
	
using namespace std;

void run(const char* File){
	


  //const char* inDir = "$ROOTSYS/tutorials";

  char* dir = gSystem->ExpandPathName(File);
  void* dirp = gSystem->OpenDirectory(dir);

  const char* entry;
  const char* filename[100];
  Int_t n = 0;
  TString str;

  while((entry = (char*)gSystem->GetDirEntry(dirp))) {
    str = entry;
    if(str.EndsWith(ext))
      filename[n++] = gSystem->ConcatFileName(dir, entry);
  }

  for (Int_t i = 0; i < n; i++)
    Printf("file -> %s", filename[i]);
}
	
	
	
	
	
	
	
	}
