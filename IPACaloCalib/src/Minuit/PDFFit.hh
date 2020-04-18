#ifndef PDFFit_HH
#define PDFFit_HH

#include "TMath.h"
#include "TF1.h"
#include "TH1F.h"

#include <Minuit2/FCNBase.h>


using namespace mu2e;

    class PDFFit : public ROOT::Minuit2::FCNBase {
    	public:
        std::vector<double> seeds;
        std::vector<double> errors;
     		std::vector<double> caloE;
        std::vector<double> trackerE;
    		std::vector<double> constraint_means;
    		std::vector<double> constraints;

    		int nparams = 674;

        PDFFit(std::vector<double> _seeds, std::vector<double> _errors,
        std::vector<double> _caloE, std::vector<double> _trackerE,
        std::vector<double> _constraint_means,
        std::vector<double> _constraints) :
        seeds(_seeds), errors(_errors), caloE(_caloE),
        _trackerE(trackerE), _constraint_means(constraint_means),
        _constraints(constraints){};
        double Up() const { return 0.5; };
        double operator() (const std::vector<double> &x) const;

    };


};

#endif

