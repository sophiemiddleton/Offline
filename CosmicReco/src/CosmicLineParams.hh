//Author: S Middleton
//Date: Nov 2019
//Purpose: Parameter Translate the Cosmic Track Fit to BTRk infrastructure
#ifndef COSMICLINEPARAMS_HH
#define COSMICLINEPARAMS_HH
#include "BTrk/TrkBase/TrkParams.hh"

// Class interface //
class CosmicLineParams : public TrkParams {
public:
  enum ParIndex {A0Index ==0, A1Index=1, B0Index, B1Index, T0Index};

  CosmicLineParams(const CLHEP::HepVector&, const CLHEP::HepSymMatrix&);
  ~CosmicLineParams();
  
  double A0() const                              {return parvec[A0Index];}
  double A1() const                            {return parvec[A1Index];}
  double B0() const                           {return parvec[B0Index];}
  double B1() const                              {return parvec[B1Index];}
  //double T0() const                          {return parvec[T0Index];}

  const CLHEP::HepVector& params() const                {return parvec;}
  CLHEP::HepVector& params()                            {return parvec;}
  const CLHEP::HepSymMatrix& covariance() const            {return parcov;}
  CLHEP::HepSymMatrix& covariance()                        {return parcov;}

  void setA0(double in)                          {parvec[A0Index] = in;} 
  void setA1(double in)                        {parvec[A1Index] = in;} 
  void setB0(double in)                       {parvec[B0Index] = in;} 
  void setB1(double in)                          {parvec[B1Index] = in;} 
  //void setT0(double in)                      {parvec[T0Index] = in;} 
  void setError(const CLHEP::HepSymMatrix& in)          {parcov = in;}

  void print(std::ostream& o) const;		// Print parameters on one line
  void printAll(std::ostream& o) const;	// Print parameters and error matrix

private:	
};

// Output operator, useful for debugging
std::ostream& operator<<(std::ostream& o, const CosmicLineParams& cosmic_line);

#endif
