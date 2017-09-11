#ifndef Fill3D_H_
#define Fill3D_H_



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


	
class Fill3D{

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
