/*
 * ZROOT.h
 *
 * Copyright 2014,2015 Andrii Verbytskyi <andriish@mppmu.mpg.de>
 * Max-Planck Institut für Physik
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 */
#ifndef ZROOT_H
#define ZROOT_H
typedef long double zdouble;
#define DEFAULT_FFORMAT "%6.6Lf"
#define _USE_MATH_DEFINES
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <iostream>
#include <iterator>
#include <algorithm>
#include <string.h>
#include <sstream>
#include <iostream>
#include <limits>
#include <cfloat>
#include <cmath>
#include <math.h>
#include <set>
#include <map>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <cassert>
#include <stdlib.h>

#include "TROOT.h"
#include "TMatrixT.h"
#include "TMatrixTSym.h"
#include "TVectorT.h"
#include "TTree.h"
#include "TBranch.h"
#include "TH1F.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TCanvas.h"
#include "TImage.h"
#include "TF1.h"
#include "TApplication.h"
#include "TMath.h"
#include "TLegend.h"
#include "TFile.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TChain.h"
#include "TGraphErrors.h"
#include "TSystemDirectory.h"
#include "TSystemFile.h"
#include "TSystem.h"
#include "TAttMarker.h"
#include "TMinuit.h"


typedef TMatrixTSym<double> ZBestMatrixTSym;
typedef TMatrixT<double>    ZBestMatrixT;
typedef TVectorT<double>    ZBestVectorT;
#endif
