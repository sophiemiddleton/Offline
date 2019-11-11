#ifndef _IPACALOCALIB_MINIMIZER_HH
#define _IPACALOCALIB_MINIMIZER_HH
#include "DataProducts/inc/XYZVec.hh"
//ROOT
#include "TMath.h"
#include "TF1.h"
#include "TH1F.h"
//Minuit
#include <Minuit2/FCNBase.h>


using namespace mu2e;

class Minimizer : public ROOT::Minuit2::FCNBase {
  public:
    
    double Nevents;
    double Nclusters;
    double EoP;
    double avEoP;
    std::vector<double> cluster_list;
    std::vector<double> constraint_means;
    std::vector<double> constraints;
    int nparams =674; 
    
     Minimizer(double &_Nevents, double &_Nclusters, double &_EoP, double &_avEoP, std::vector<double> &_cluster_list, std::vector<double> &_constraint_means, std::vector<double> &_constraints) : Nevents(_Nevents), Nclusters(_Nclusters), EoP(_EoP), avEoP(_avEoP), cluster_list(_cluster_list),constraint_means(_constraint_means), constraints(_constraints){};
   
    double Up() const { return 0.5; };
    double operator() (const std::vector<double> &x) const;
    
    
};


#endif

  
