#include "BTrk/BaBar/BaBar.hh"
#include <assert.h>
#include <math.h>
#include <limits.h>

#include "BTrk/BaBar/Constants.hh"
#include "BTrk/BbrGeom/HepPoint.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Matrix/SymMatrix.h"
#include "BTrk/TrkBase/HelixTraj.hh"
#include "CosmicReco/inc/CosmicTrkVisitor.hh"
#include "BTrk/difAlgebra/DifNumber.hh"
#include "BTrk/difAlgebra/DifPoint.hh"
#include "BTrk/difAlgebra/DifVector.hh"
#include "BTrk/BbrGeom/BbrAngle.hh"
#include "CosmicReco/inc/CosmicLineParams.hh"
#include "CosmicReco/inc/CosmicLineTraj.hh"
#include "CosmicReco/inc/CosmicTrkMomVisitor.hh"
#include "BTrk/BaBar/ErrLog.hh"
using std::endl;
using std::ostream;
using namespace CLHEP;
#ifndef M_2PI
#define M_2PI 2*M_PI
#endif

//reference point --> is this the center of the tracker system?

CosmicLineTraj::CosmicLineTraj(const HepVector& pvec, const HepSymMatrix& pcov,
                     double lowlim, double hilim, const HepPoint& refpoint) :
  TrkSimpTraj(pvec, pcov, lowlim,hilim,refpoint)
{
  //  Make sure the dimensions of the input matrix and vector are correct
  if( pvec.num_row() != NHLXPRM || pcov.num_row() != NHLXPRM ){
    ErrMsg(fatal) 
      << "CosmicLineTraj: incorrect constructor vector/matrix dimension" << endmsg;
  }

}


CosmicLineTraj::CosmicLineTraj(const CosmicLineParams& inpar,
                     double lowlim, double hilim, const HepPoint& refpoint) :
  TrkSimpTraj(inpar.params(), inpar.covariance(), lowlim,hilim,refpoint)
{
  
}

CosmicLineTraj::CosmicLineTraj(const TrkParams& inpar,
                     double lowlim, double hilim, const HepPoint& refpoint) :
  TrkSimpTraj(inpar, lowlim,hilim,refpoint)
{
  assert(inpar.parameter().num_row()==NHLXPRM);
}

CosmicLineTraj::CosmicLineTraj( const CosmicLineTraj& h )
  : TrkSimpTraj(h.parameters()->parameter(), h.parameters()->covariance(), 
                h.lowRange(),h.hiRange(),h.referencePoint())
{
}

CosmicLineTraj*
CosmicLineTraj::clone() const //TODO-causing compilation errors due to non-complete type
{
  return new CosmicLineTraj(*this);
}

CosmicLineTraj&
CosmicLineTraj::operator=(const CosmicLineTraj& h)
{
  if( &h != this ){
    Trajectory::operator=(h);
    _dtparams = *(h.parameters());
    _refpoint = h._refpoint;
  }
  return *this;
}

CosmicLineTraj::~CosmicLineTraj()
{
}

//Z: returns the z position i.e. position along DS axis.
double
CosmicLineTraj::z(const double& f) const  
{
  return  (d0()*tan(theta())*cos(phi() - phi0())+referencePoint().z()f*cos(theta()); 
}

//Z0: returns the POCA z0 value.
CosmicLineTraj::z0() const {
 return (d0()*tan(theta())*sin(phi() - phi0());

}

//zflight: returns the projection of the flight along the z axis.
double
CosmicLineTraj::zFlight(double zpos, double z0) const { 
  return (zpos - z0)/cos(theta());
}

//Position returns a 3 element position (x,y,z) 'point' at flight length 'f"
HepPoint
CosmicLineTraj::position(double f) const  
{
  double sphi0 = sin(phi0());
  double cphi0 = cos(phi0());
  //keep definition the same as for Helix for consistancy:
  double x_pos = -1*d0()*sphi0+referencePoint().x() + cos(phi())*sin(theta())*f;
  double y_pos = d0()*cphi0+referencePoint().y()+  sin(theta())*sin(phi())*f;
  double z_pos = d0()*tan(phi()-phi0())+referencePoint().z() + cos(theta())*f;
  return HepPoint(x_pos, y_pos, z_pos);
}

//Direction returns the unit vector for track direction
Hep3Vector
CosmicLineTraj::direction(double f) const 
{

double x_dir = cos(phi())*sin(theta());
double y_dir = sin(theta())*sin(phi());
double z_dir = cos(theta());
return Hep3Vector (x_dir, y_dir, z_dir);
}


//DelDirect is the 2nd deriv wrt flight lengths
Hep3Vector
CosmicLineTraj::delDirect( double fltLen ) const
{  
  double delX = 0;
  double delY = 0; 
  return Hep3Vector(delX, delY, 0.0);
}

double
CosmicLineTraj::distTo1stError(double s, double tol, int pathDir) const 
{
  return 9999;
}

double
CosmicLineTraj::distTo2ndError(double s, double tol, int pathDir) const 
{
  return 9999;
}

void
CosmicLineTraj::getInfo(double fltLen, HepPoint& pos, Hep3Vector& dir) const
{
  pos = position(fltLen);
  dir = direction(fltLen);
}

void
CosmicLineTraj::getInfo(double fltLen, HepPoint& pos, Hep3Vector& dir, 
                   Hep3Vector& delDir) const 
{
  pos = position(fltLen);
  dir = direction(fltLen);
  delDir = delDirect(fltLen);
}

//derivDeflect: considers scattering in 2 orthogonal directions.
HepMatrix
CosmicLineTraj::derivDeflect(double fltlen, deflectDirection idirect) const //TODO
{
//  This function computes the column matrix of derrivatives for the change
//  in parameters for a change in the direction of a track at a point along
//  its flight, holding the momentum and position constant.  The effects for
//  changes in 2 perpendicular directions (theta1 = theta and
//  theta2 = phi()*sin(theta) can sometimes be added, as scattering in these
//  are uncorrelated. These axes are track specific. as cosmics are not always coming along the same track direction it is necessary to have difference parameterization than that used for the helixal case.
//
  HepMatrix ddflct(NHLXPRM,1);

  switch (idirect) {
  case theta1: 
    
    ddflct(thetaIndex+1,1) = 1;
    ddflct(d0Index+1,1) = 0;
    ddflct(phi0Index+1,1) = 0;
    ddflct(phiIndex+1,1) = 0;
    break;

  case theta2: 
    
    ddflct(thetaIndex+1,1) = 0;
    ddflct(d0Index+1,1) = 0;
    ddflct(phi0Index+1,1) = 1/sin(theta());
    ddflct(phiIndex+1,1) = 1/sin(theta());
    break;
  }

  return ddflct;
}

//DerivDisplace: the effect of a change in position: this is gnored as it is a second order effect.
HepMatrix
CosmicLineTraj::derivDisplace(double fltlen, deflectDirection idirect) const
{
 HepMatrix ddflct(NHLXPRM,1);

 switch (idirect) {
  case theta1:
   
    ddflct(thetaIndex+1,1) = 0;
    ddflct(d0Index+1,1) = 0;
    ddflct(phi0Index+1,1) = 0;
    ddflct(phiIndex+1,1) = 0;
    break;
  case theta2:
    
    ddflct(thetaIndex+1,1) = 0;
    ddflct(d0Index+1,1) = 0;
    ddflct(phi0Index+1,1) = 0;
    ddflct(phiIndex+1,1) = 0;
    break;
  }

  return ddflct;
}

//derivPFract: Function computes column marix of derivatives for parameters from a fractional change in the track momentium. Holds postion and direction constant. The momentum change could from energy loss.
HepMatrix CosmicLineTraj::derivPFract(double fltLen) const {

     HepMatrix dmomfrac(NHLXPRM,1);
    
     dmomfrac(phiIndex+1,1) = 0;
     dmomfrac(theta0Index+1,1) = 0;
     dmomfrac(d0Index+1,1) = 0;
     dmomfrac(phi0Index+1,1) = 0;
     return dmomfrac;
     

}

void
CosmicLineTraj::getDFInfo(double flt, DifPoint& pos, DifVector& dir,
                     DifVector& delDir) const
{
  //Provides difNum version of information for calculation of derivatives.
  //  All arithmetic operations have been replaced by +=, etc. versions 
  //  for speed.

  DifNumber phi0Df(phi0(), phi0Index+1, NHLXPRM);
  phi0Df.setIndepPar( parameters() );
  DifNumber d0Df(d0(), d0Index+1, NHLXPRM);
  d0Df.setIndepPar( parameters() );
  DifNumber thetaDf(theta(), thetaIndex+1, NHLXPRM);
  thetaDf.setIndepPar( parameters() );
  DifNumber phiDf(phi(), phiIndex+1, NHLXPRM);
  phiDf.setIndepPar( parameters() );

  static DifNumber sTheta; 
  thetaDf.cosAndSin(cTheta, sTheta);
  static DifNumber sinPhi0, cosPhi0; 
  phi0Df.cosAndSin(cosPhi0, sinPhi0);
  static DifNumber sinPhi, cosPhi; 
  phi0Df.cosAndSin(cosPhi, sinPhi);

  bool lref = (referencePoint().x() != 0. || referencePoint().y() != 0. ||
               referencePoint().z() != 0.);


  DifNumber x =  sTheta*flt*cosPhi - d0Df*sinPhi0;
  DifNumber y =  sTheta*flt*sinPhi + d0Df*cosPhi0;
  DifNumber z =  d0Df*tan(phiDf-phi0Df) + flt*cTheta;
  pos =  DifPoint(x, y, z);
  dir = DifVector( cosPhi*sTheta, sinPhi*sTheta, cTheta);
  delDir = DifVector(0., 0., 0.);

  if (lref) {
    DifNumber px(referencePoint().x());
    DifNumber py(referencePoint().y());
    DifNumber pz(referencePoint().z());
    pos.x += px;
    pos.y += py;
    pos.z += pz;
  }

}



void CosmicLineTraj::getDFInfo2(double flt, DifPoint& pos, DifVector& dir) const //TODO
{
  //Provides difNum version of information for calculation of derivatives.
  //  All arithmetic operations have been replaced by +=, etc. versions 
  //  for speed.


  DifNumber phi0Df(phi0(), phi0Index+1, NHLXPRM);
  phi0Df.setIndepPar( parameters() );
  DifNumber d0Df(d0(), d0Index+1, NHLXPRM);
  d0Df.setIndepPar( parameters() );
  DifNumber thetaDf(theta(), thetaIndex+1, NHLXPRM);
  thetaDf.setIndepPar( parameters() );
  DifNumber phiDf(phi(), phiIndex+1, NHLXPRM);
  phiDf.setIndepPar( parameters() );

  static DifNumber sTheta; 
  thetaDf.cosAndSin(cTheta, sTheta);
  static DifNumber sinPhi0, cosPhi0; 
  phi0Df.cosAndSin(cosPhi0, sinPhi0);
  static DifNumber sinPhi, cosPhi; 
  phi0Df.cosAndSin(cosPhi, sinPhi);

  bool lref = (referencePoint().x() != 0. || referencePoint().y() != 0. ||
               referencePoint().z() != 0.);


  DifNumber x =  sTheta*flt*cosPhi - d0Df*sinPhi0;
  DifNumber y =  sTheta*flt*sinPhi + d0Df*cosPhi0;
  DifNumber z =  d0Df*tan(phiDf-phi0Df) + flt*cTheta;
  pos =  DifPoint(x, y, z);
  dir = DifVector( cosPhi*sTheta, sinPhi*sTheta, cTheta);
  delDir = DifVector(0., 0., 0.);

  if (lref) {
    DifNumber px(referencePoint().x());
    DifNumber py(referencePoint().y());
    DifNumber pz(referencePoint().z());
    pos.x += px;
    pos.y += py;
    pos.z += pz;
  }


}


double
CosmicLineTraj::curvature() const   //TODO 
{

  return 1;
}


void CosmicLineTraj::visitAccept(TrkVisitor* vis) const
{ 
   vis->trkVisitCosmicLineTraj(this); 
}

double
CosmicLineTraj::angle(const double& f) const 
{
  return BbrAngle(phi0() + f*phi0());
}

void
CosmicLineTraj::printAll(ostream& os) const
{
  os  << "CosmicLineTraj with range "
      << lowRange() <<" to " << hiRange() << " and parameters " << endl
      << "d0= " << d0() << " phi0= "
      << phi0() << " theta = "
      << theta() << " phi =  "
      << phi()  << endl;
}

void
CosmicLineTraj::print(ostream& os) const
{
  Trajectory::print(os << "CosmicLineTraj" );
}
