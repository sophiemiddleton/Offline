//Author: S Middleton
//Date: Nov 2019
//Purpose: To allow cosmic fitting to be compatible with Kalman/BTrk

/*What we need: 
need to store parameters, positon, direction, 6 element convarience matrix



*/
#ifndef COSMICLINETRAJ_HH
#define COSMICLINETRAJ_HH

#include "BTrk/TrkBase/TrkDifTraj.hh"
#include "BTrk/TrkBase/TrkKalDeriv.hh"
#include "BTrk/TrkBase/TrkParams.hh"
#include "BTrk/BbrGeom/HepPoint.h"
#include <vector>

class TrkVisitor;
#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Matrix/SymMatrix.h"
#include <iosfwd>


class CosmicLineTraj : public TrkDifTraj, public TrkKalDeriv {

public:

  static HepPoint _theOrigin; // define the origin as a HepPoint


public:

  //-----------------------
  // Constructors and such
  //-----------------------
  CosmicLineTraj(const CLHEP::HepVector& params, const CLHEP::HepSymMatrix& cov,
              const double startRange = -99999.,const double endRange =99999.,
              const HepPoint& refpoint = _theOrigin);
  CosmicLineTraj(const TrkParams& params,
              const double startRange = -99999.,const double endRange =99999.,
              const HepPoint& refpoint = _theOrigin);
  virtual ~CosmicLineTraj();

  virtual CosmicLineTraj* clone() const = 0;

  //--------------------------------------------------
  //  Access to parameters, errors and reference point
  //--------------------------------------------------
  TrkParams*                 parameters()                   {return &_dtparams;}
  const TrkParams*           parameters() const             {return &_dtparams;}
  virtual const CosmicLineTraj* localTrajectory(double fltLen, double& localFlt)
    const;
  const HepPoint&            referencePoint() const         {return _refpoint;}
  virtual void               print(std::ostream& os) const;
  virtual                    void printAll(std::ostream& os) const;

  //--------------------------------------------------
  // Change contents
  //--------------------------------------------------
  // Change the reference point and the parameters
  void           changePoint(const HepPoint& newpoint,double& fltlen);
  // Set the ref point and don't change the params
  void           setPoint(const HepPoint& newpoint)    {_refpoint = newpoint;}
  // inversion function: this inverts both the flight range and the parameters
  // so that the same traj is described but going in the opposite direction.
  CosmicLineTraj& invert();
  // invert the track parameters passed in newparams.
  // Returns true in flags if the inversion requires a sign change in the
  // covariance matrix as well.
  virtual void invertParams(TrkParams* newparams, std::vector<bool>& flags) const = 0;
  // Provide function which translates the reference point of parameters
  virtual TranslateParams paramFunction() const = 0;

  //--------------------------------------------------
  // Visitor access for momentum functions
  //--------------------------------------------------

  virtual void visitAccept(TrkVisitor* vis) const = 0;

  bool operator==(const CosmicLineTraj&) const; // return equivalence, not identy


protected:
  TrkParams _dtparams;
  HepPoint _refpoint; // reference point for parameters
private:
  // Preempt 
  CosmicLineTraj(const CosmicLineTraj &);
  CosmicLineTraj&   operator= (const CosmicLineTraj&);
};
#endif
