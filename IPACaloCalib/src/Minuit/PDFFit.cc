//Author: S Middleton
//Date: August 2019
//Purpose: PDF Functions for minuit fitting

//ROOT:
#include "Math/VectorUtil.h"
#include "TMath.h"
#include "Math/Math.h"
#include <TSystem.h>
#include <TROOT.h>
#include <TObjString.h>

//Utilities:
#include "PDFFit.hh"

//Minuit
#include <Minuit2/FCNBase.h>
#include <Minuit2/MnMigrad.h>
#include <Minuit2/MnMinos.h>
#include <Minuit2/MnStrategy.h>
#include <Minuit2/MnUserParameters.h>
#include <Minuit2/FunctionMinimum.h>

using namespace mu2e;

double PDFFit::operator() (const std::vector<double> &x) const
{
  	long double F = 0;
  	for (size_t i=0;i<this->seeds.size();i++){
      double const& seed = seeds[i];
  		F+=(seeds[i]*this->caloE[i]-this->trackerE)/this->errors[i]
    }
    return F;
  }

