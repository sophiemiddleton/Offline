#include "Math/VectorUtil.h"
#include "TMath.h"
#include "Math/Math.h"
#include <TSystem.h>
#include <TROOT.h>
#include <TObjString.h>
#include "RecoDataProducts/inc/ComboHit.hh"
#include "IPACaloCalib/inc/Minimizer.hh"

#include <Minuit2/FCNBase.h>
#include <Minuit2/MnMigrad.h>
#include <Minuit2/MnMinos.h>
#include <Minuit2/MnStrategy.h>
#include <Minuit2/MnUserParameters.h>
#include <Minuit2/FunctionMinimum.h>


using namespace mu2e;

double Minimizer::operator() (const std::vector<double> &x) const
{
  
  double totF = 0;
  for(size_t e=0;e<this->Nevents;e++){
	double par_i = 0;
	for (size_t c=0;c<this->Nclusters;c++){
		//int cluster_number = this->cluster_list[c];
		par_i += x[c]*this->EoP - this->avEoP;
	}
	totF += par_i;
	  
  }
  return totF;
}
