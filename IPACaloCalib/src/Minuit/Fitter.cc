//Author : S Middleton
//Date : August 2019
//Purpose: calls  minuit fitting to cosmic track seed. Input is CosmicTrackSeed, can then derive parameters from CosmicTrack stored there.

//ROOT:
#include "Math/VectorUtil.h"
#include "TMath.h"
#include "Math/Math.h"
#include "Fitter.hh"
#include "PDFFit.hh"

//For Drift:
#include <TSystem.h>
#include <TROOT.h>
#include <TObjString.h>

//Minuit
#include <Minuit2/FCNBase.h>
#include <Minuit2/MnMigrad.h>
#include <Minuit2/MnMinos.h>
#include <Minuit2/MnStrategy.h>
#include <Minuit2/MnPrint.h>
#include <Minuit2/MnUserParameters.h>
#include <Minuit2/FunctionMinimum.h>

using namespace mu2e;

namespace MinuitDriftFitter{
	 CaloFitResult DoFit(std::vector<double> seeds, std::vector<double> errors,
		 std::vector<double> caloE, std::vector<double> trackerE){

	std::vector<double> errors(674,0);
	std::vector<double> seed(674,0);
	CaloFitResult FitResult;

	//Constrain to mean = 0 for 4 parameters (T0 might not be so...)
	std::vector<double> constraint_means(674,1);
	std::vector<double> constraints(674,1);

	PDFFit fit(td::vector<double> seeds, std::vector<double> errors,
		std::vector<double> caloE, std::vector<double> trackerE,
		constraint_means, constraints);

	//Initiate Minuit Fit:
	ROOT::Minuit2::MnStrategy mnStrategy(2);
	ROOT::Minuit2::MnUserParameters params(seed,errors);
	ROOT::Minuit2::MnMigrad migrad(fit,params,mnStrategy);

	//Set Limits as tracker dimensions:
	for(size_t i=0;i<seeds.size(); i++){
		migrad.SetLimits((signed) i, 0, 1);
	}
	int maxfcn = 100;
	double tolerance = 1000;

	//Define Minimization method as "MIGRAD" (see minuit documentation)
	ROOT::Minuit2::FunctionMinimum min = migrad(maxfcn, tolerance);
	if(_diag > 1){
		ROOT::Minuit2::MnPrint::SetLevel(3);
		ROOT::Minuit2::operator<<(cout, min);
	}

	//Will be the results of the fit routine:
	ROOT::Minuit2::MnUserParameters results = min.UserParameters();
	double minval = min.Fval();

	//Define name for parameters
	FitResult.bestfit = results.Params();
	FitResult.bestfiterrors = results.Errors();

	//Store Minuit Covarience if exists:
	if(min.HasValidCovariance()) FitResult.bestfitcov = min.UserCovariance().Data();

	//Name Parameters:
	for(size_t i=0;i<seeds.size(); i++){
		FitResult.names.push_back("i");
	}
	FitResult.NLL = minval;
	return FitResult;

}

