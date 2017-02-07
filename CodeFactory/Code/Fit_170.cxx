#include <string>
#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TH1F.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TSystem.h"

using namespace std;


void list_files(const char *dirname, const char *ext=".root")
{
   TSystemDirectory dir(dirname, dirname);
   TList *files = dir.GetListOfFiles();
   
   if (files) {
      TSystemFile *file;
      TString fname;
      TIter next(files);
	  TH1F *h1;
   
      while ((file=(TSystemFile*)next())) {
         fname = file->GetName();
         
         if (!file->IsDirectory() && fname.EndsWith(ext)) {
            cout << fname.Data() << endl;
            
            
            TFile *file1 = new TFile(fname,"r");
            h1 = (TH1F*)file1->Get("hist_klf_mtop_window");
            h1->Draw();
            
         }
	 }
}
}
