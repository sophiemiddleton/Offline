//Author: S Middleton
//Date: Dec 2019
//Purpose: Parameter Translate the Cosmic Track Fit to BTRk infrastructure

/*
Provisional parameterization of the KF linear trajectory class. The (geometric) components will be:
 	-theta and phi of the direction vector WRT to the tracker coordinate system -->Track direction (called theta and phi here)
	-signed distance of closest approach (DOCA) to the tracker system origin ---->Signed distance (called d0 here)
	- azimuth of the point of closest approach (POCA) vector --> Measure of Distance (called phi0 here)

*/



#ifndef COSMICLINEPARAMS_HH
#define COSMICLINEPARAMS_HH
#include "BTrk/TrkBase/TrkParams.hh"

// Class interface //
class CosmicLineParams : public TrkParams {
public:

  enum ParIndex {d0Index=0, phi0Index, thetaIndex, phiIndex,  NHLXPRM};

  CosmicLineParams(const CLHEP::HepVector&, const CLHEP::HepSymMatrix&);
  ~CosmicLineParams();
  //position: --> DOCA, POCA
  double d0() const                              {return parvec[d0Index];}
  double phi0() const                            {return parvec[phi0Index];}
  //direction: ---> angles
  double theta() const                              {return parvec[thetaIndex];}
  double phi() const                          {return parvec[phiIndex];}

  const CLHEP::HepVector& params() const                {return parvec;}
  CLHEP::HepVector& params()                            {return parvec;}
  const CLHEP::HepSymMatrix& covariance() const            {return parcov;}
  CLHEP::HepSymMatrix& covariance()                        {return parcov;}
//Distance Parameters:
  void setD0(double in)                          {parvec[d0Index] = in;} 
  void setPhi0(double in)                        {parvec[phi0Index] = in;} 
//Direction Parameters:
  void setTheta(double in)                          {parvec[ThetaIndex] = in;} 
  void setPhi(double in)                      {parvec[PhiIndex] = in;} 

  void setError(const CLHEP::HepSymMatrix& in)          {parcov = in;}

  void print(std::ostream& o) const;		// Print parameters on one line
  void printAll(std::ostream& o) const;	// Print parameters and error matrix

private:	
};

// Output operator, useful for debugging
std::ostream& operator<<(std::ostream& o, const CosmicLineParams& cosmic);

#endif






