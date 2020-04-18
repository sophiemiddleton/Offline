#ifndef FITTER_HH
#define FITTER_HH

//ROOT
#include "TMath.h"
#include "TF1.h"
#include "TH1F.h"

//Minuit
#include <Minuit2/FCNBase.h>

using namespace mu2e;
  struct CaloFitResult{
    	public:
    		std::vector<std::string> names;
    		std::vector<double> bestfit; //the constants
    		std::vector<double> bestfiterrors; //errors on those constants
    		std::vector<double> bestfitcov; //covariences on the constants
    		double NLL;
	};

namespace Fitter {
	CaloFitResult DoFit(std::vector<double> seeds, std::vector<double> errors);
}


#endif

