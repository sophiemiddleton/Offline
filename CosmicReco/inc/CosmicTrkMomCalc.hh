
#ifndef COSMICTRKMOMCALC_HH
#define COSMICTRKMOMCALC_HH

class TrkSimpTraj;
class CosmicLineTraj;
class BField;
class HelixTraj;

#include "CLHEP/Vector/ThreeVector.h"

class BbrVectorErr;
class DifNumber;

// Class interface //
class CosmicTrkMomCalc{

public:

  CosmicTrkMomCalc();
  virtual ~CosmicTrkMomCalc();

  //********************************
  // The calculator functions:
  //********************************

  static double       ptMom( const TrkSimpTraj&, const BField&, double fltlen);
  static CLHEP::Hep3Vector   vecMom(const TrkSimpTraj&, const BField&, double fltlen);
  static BbrVectorErr errMom(const TrkSimpTraj&, const BField&, double fltlen);
  static int          charge(const TrkSimpTraj&, const BField&, double fltlen);

// calculate the df form of the momentum magnitude.  This can be done for generic
// trajectories, but is hideously inefficient
  static DifNumber momMag(const HelixTraj& theTraj, const BField&);
  static DifNumber momMagCosmic(const CosmicLineTraj& theTraj, const BField&);
  // Interface to vertexing algorithms (M.Bondioli 7/17/98)
  // covariance matrices of the track at fixed flight length 
  static CLHEP::HepMatrix    posmomCov(const TrkSimpTraj&,const BField&,
				double fltlen);
  static void         getAllCovs(const TrkSimpTraj&,const BField&,
				 double fltlen,
				 CLHEP::HepSymMatrix& xxCov,
				 CLHEP::HepSymMatrix& ppCov,
				 CLHEP::HepMatrix&    xpCov);
  // invert 6x6
  static bool weightToCov(const CLHEP::HepSymMatrix& inXX,const CLHEP::HepSymMatrix& inPP,const CLHEP::HepMatrix& inXP,
			  CLHEP::HepSymMatrix& outXX,     CLHEP::HepSymMatrix& outPP,     CLHEP::HepMatrix& outXP);

  // accessors to 2nd derivative of chi2 wrt x and p.
  // x0 and p0 are filled with the pos and 3mom around which expansion 
  // takes place, whilst Weights are filled with the 2nd deriv of chi2.
  // so that:
  //  dx=x-x0   dp=p-p0
  //
  //  chi2(x,p)=0.5 dx^t*Wxx*dx + dx^t*Wxp*dp + 0.5 dp^t*Wpp*dp
  // where:
  //  pos and mom are 3-dim Vectors,
  //  xxWeight ppWeight and xpWeight are 3 by 3 matrices
  //
  // takes place, whilst Weights are filled with the 2nd deriv of chi2.
  // so that:
  //  dx=x-x0   dp=p-p0
  //
  //  chi2(x,p)=0.5 dx^t*Wxx*dx + dx^t*Wxp*dp + 0.5 dp^t*Wpp*dp
  // where:
  //  pos and mom are 3-dim Vectors,
  //  xxWeight ppWeight and xpWeight are 3 by 3 matrices
  //
  static void         getAllWeights(const TrkSimpTraj&, const BField&,
				    double fltlen,
				    CLHEP::HepVector& pos,
				    CLHEP::HepVector& mom,
				    CLHEP::HepSymMatrix& xxWeight,
				    CLHEP::HepSymMatrix& ppWeight,
				    CLHEP::HepMatrix&    xpWeight);

private:

  static double       calcCurvPtMom(const CLHEP::Hep3Vector&,  
                                    double curvature, 
				    const BField&);
  static double       calcCosmicPtMom(const CLHEP::Hep3Vector&,  
                                    double curvature, 
				    const BField&);
  static CLHEP::Hep3Vector   calcCurvVecMom(const CLHEP::Hep3Vector&, 
                                     double curvature, 
				     const BField&);
    static CLHEP::Hep3Vector   calcCosmicVecMom(const CLHEP::Hep3Vector&, 
                                     double curvature, 
				     const BField&);
  static BbrVectorErr calcCurvErrMom(const TrkSimpTraj&, 
				     const BField&,  
                                     double flt); 
  static BbrVectorErr calcNeutErrMom(const TrkSimpTraj&, const BField&,  
                                     double flt);
  static BbrVectorErr calcCosmicErrMom(const TrkSimpTraj&, const BField&,  
                                     double flt);
  static int          calcCurvCharge(const CLHEP::Hep3Vector&, 
                                     double curvature, 
				     const BField&);        

  static CLHEP::HepMatrix    calcCurvPosmomCov(const TrkSimpTraj&,const BField&,
					double fltlen);
  static CLHEP::HepMatrix    calcNeutPosmomCov(const TrkSimpTraj&,const BField&,
					double fltlen);
  static CLHEP::HepMatrix    calcCosmicPosmomCov(const TrkSimpTraj&,const BField&,
					double fltlen);
  static void         calcCurvAllCovs(const TrkSimpTraj&,const BField&,
				      double fltlen,
				      CLHEP::HepSymMatrix& xxCov,
				      CLHEP::HepSymMatrix& ppCov,
				      CLHEP::HepMatrix&    xpCov);

  static void         calcCurvAllCovsOLD(const TrkSimpTraj&,const BField&,
				      double fltlen,
				      CLHEP::HepSymMatrix& xxCov,
				      CLHEP::HepSymMatrix& ppCov,
				      CLHEP::HepMatrix&    xpCov);

  static void         calcNeutAllCovs(const TrkSimpTraj&,const BField&,
				      double fltlen,
				      CLHEP::HepSymMatrix& xxCov,
				      CLHEP::HepSymMatrix& ppCov,
				      CLHEP::HepMatrix&    xpCov);

 
  static void         calcNeutAllCovs(const TrkSimpTraj&,const BField&,
				      double fltlen,
				      CLHEP::HepVector& x0,CLHEP::HepVector& p0,
				      CLHEP::HepSymMatrix& xxCov,
				      CLHEP::HepSymMatrix& ppCov,
				      CLHEP::HepMatrix&    xpCov);

 static void         calcCosmicAllCovs(const TrkSimpTraj&,const BField&,
				      double fltlen,
				      CLHEP::HepSymMatrix& xxCov,
				      CLHEP::HepSymMatrix& ppCov,
				      CLHEP::HepMatrix&    xpCov);


  static void         calcCurvAllWeights(const TrkSimpTraj&,const BField&,
					 double fltlen,
					 CLHEP::HepVector& pos,
					 CLHEP::HepVector& mom,
					 CLHEP::HepSymMatrix& xxWeight,
					 CLHEP::HepSymMatrix& ppWeight,
					 CLHEP::HepMatrix&    xpWeight);

  static void         calcCurvAllWeightsOLD(const TrkSimpTraj&,const BField&,
					 double fltlen,
					 CLHEP::HepVector& pos,
					 CLHEP::HepVector& mom,
					 CLHEP::HepSymMatrix& xxWeight,
					 CLHEP::HepSymMatrix& ppWeight,
					 CLHEP::HepMatrix&    xpWeight);

  static void         calcNeutAllWeights(const TrkSimpTraj&,const BField&,
					 double fltlen,
					 CLHEP::HepVector& pos,
					 CLHEP::HepVector& mom,
					 CLHEP::HepSymMatrix& xxWeight,
					 CLHEP::HepSymMatrix& ppWeight,
					 CLHEP::HepMatrix&    xpWeight);
  static void         calcCosmicAllWeights(const TrkSimpTraj&,const BField&,
					 double fltlen,
					 CLHEP::HepVector& pos,
					 CLHEP::HepVector& mom,
					 CLHEP::HepSymMatrix& xxWeight,
					 CLHEP::HepSymMatrix& ppWeight,
					 CLHEP::HepMatrix&    xpWeight);

  static int          nearestInt(double);

};

#endif