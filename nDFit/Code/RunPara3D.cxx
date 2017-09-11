#include "TSystem.h"
#include "TFile.h"
#include "TF1.h"
#include "TTree.h"
#include "TH1F.h"
#include "TROOT.h"
#include "TPad.h"
#include "TMath.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TMinuit.h"
#include "TAttMarker.h"

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
#include <omp.h>
#include "Para3D.h"











int main(int argc, char* argv[])
{
	
	ROOT::EnableThreadSafety();
	if (argc>2) exit(1);
	
	//std::string InputFolder = argv[1];
	//std::string OutputFolder = argv[1];

	TemplateHolder* top_quark = new TemplateHolder();
	
	top_quark->input(argv[1]);
		  
	top_quark->top_fit();  
	
	top_quark->mtop_plot();
	
	
	
return 0;	
}
