#include "BTrk/BaBar/BaBar.hh"
#include <assert.h>
#include <math.h>
#include <limits.h>

#include "BTrk/BaBar/Constants.hh"
#include "BTrk/BbrGeom/HepPoint.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Matrix/SymMatrix.h"
#include "BTrk/TrkBase/CosmicLineTraj.hh"
#include "BTrk/TrkBase/TrkVisitor.hh"
#include "BTrk/difAlgebra/DifNumber.hh"
#include "BTrk/difAlgebra/DifPoint.hh"
#include "BTrk/difAlgebra/DifVector.hh"
#include "BTrk/BbrGeom/BbrAngle.hh"
#include "BTrk/TrkBase/CosmicLineParams.hh"
#include "BTrk/BaBar/ErrLog.hh"
using std::endl;
using std::ostream;
using namespace CLHEP;


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
CosmicLineTraj::clone() const
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

double
CosmicLineTraj::z(const double& f) const //TODO 
{
  return  f + referencePoint().z();
}

double
CosmicLineTraj::zFlight(double zpos) const { //TODO 
  return (zpos);
}

HepPoint
CosmicLineTraj::position( double f) const //TODO check this maths makes sense
{
  double sphi0 = sin(phi0());
  double cphi0 = cos(phi0());
 
  double x_pos = d0()*cphi0+referencePoint().x();
  double y_pos = sphi0+referencePoint().y();
  double z_pos = f+referencePoint().z();
  return HepPoint(x_pos, y_pos, z_pos);
}


Hep3Vector
CosmicLineTraj::direction( double f) const //TODO - check maths
{
 double x_dir = cos(theta());
 double y_dir = cos(theta())/tan(phi());
 double z_dir = sqrt(pow(sin(theta()),2) - pow(cos(theta()),2)/pow(tan(phi()),2));
 return Hep3Vector (x_dir, y_dir, z_dir);
}

Hep3Vector
CosmicLineTraj::delDirect( double fltLen ) const 
{
 return Hep3Vector(0,0, 0.0);
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
CosmicLineTraj::getInfo(double fltLen, HepPoint& pos, Hep3Vector& dir, //TODO
                   Hep3Vector& delDir) const
{
  pos = position(fltLen);
  dir = direction(fltLen);
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
  
  static DifNumber cDip;
  dipDf.cosAndSin(cDip, dir.z);
  static DifNumber sinPhi0, cosPhi0;
  phi0Df.cosAndSin(cosPhi0, sinPhi0);

  bool lref = (referencePoint().x() != 0. || referencePoint().y() != 0. ||
               referencePoint().z() != 0.);
/*TODO---> what is alpha?
  DifNumber alphaDf = cDip;
  alphaDf *= omegaDf;
  alphaDf *= flt;
  alphaDf += phi0Df;

  alphaDf.mod(-Constants::pi, Constants::pi);
  alphaDf.cosAndSin(dir.x, dir.y);

  pos.x = dir.y;
  pos.x -= sinPhi0;
  pos.x /= omegaDf;
  DifNumber temp = d0Df;
  temp *= sinPhi0;
  pos.x -= temp;

  pos.y = cosPhi0;
  pos.y -= dir.x;
  pos.y /= omegaDf;
  temp = d0Df;
  temp *= cosPhi0;
  pos.y += temp;

  pos.z = flt;
  pos.z *= dir.z;
  pos.z += z0Df;

  if (lref) {
    DifNumber px(referencePoint().x());
    DifNumber py(referencePoint().y());
    DifNumber pz(referencePoint().z());
    pos.x += px;
    pos.y += py;
    pos.z += pz;
  }

  dir.x *= cDip;
  dir.y *= cDip;
*/
}

void
CosmicLineTraj::getDFInfo(double flt, DifPoint& pos, DifVector& dir, //TODO
                     DifVector& delDir) const
{
  //Provides difNum version of information for calculation of derivatives.
  //  All arithmetic operations have been replaced by +=, etc. versions 
  //  for speed.

  // Create difNumber versions of parameters
  DifNumber phi0Df(phi0(), phi0Index+1, NHLXPRM);
  DifNumber d0Df(d0(), d0Index+1, NHLXPRM);
  DifNumber thetaDf(theta(), thetaIndex+1, NHLXPRM);
  DifNumber phiDf(phi(), phiIndex+1, NHLXPRM);
  
  phi0Df.setIndepPar( parameters() );
  d0Df.setIndepPar( parameters() );
  thetaDf.setIndepPar( parameters() );
  phiDf.setIndepPar( parameters() );
  
/*
  static DifNumber cDip;
  dipDf.cosAndSin(cDip, dir.z);
  static DifNumber sinPhi0, cosPhi0;
  phi0Df.cosAndSin(cosPhi0, sinPhi0);

  bool lref = (referencePoint().x() != 0. || referencePoint().y() != 0. ||
	       referencePoint().z() != 0.);

  DifNumber alphaDf = cDip;
  alphaDf *= omegaDf;
  alphaDf *= flt;
  alphaDf += phi0Df;
alphaDf.mod(-Constants::pi, Constants::pi);
  alphaDf.cosAndSin(dir.x, dir.y);

  pos.x = dir.y;
  pos.x -= sinPhi0;
  pos.x /= omegaDf;
  DifNumber temp = d0Df;
  temp *= sinPhi0;
  pos.x -= temp;

  pos.y = cosPhi0;
  pos.y -= dir.x;
  pos.y /= omegaDf;
  temp = d0Df;
  temp *= cosPhi0;
  pos.y += temp;

  pos.z = flt;
  pos.z *= dir.z;
  pos.z += z0Df;

  if (lref) {
    DifNumber px(referencePoint().x());
    DifNumber py(referencePoint().y());
    DifNumber pz(referencePoint().z());
    pos.x += px;
    pos.y += py;
    pos.z += pz;
  }

  delDir.x = -omegaDf;
  delDir.x *= cDip;
  delDir.x *= cDip;
  delDir.x *= dir.y;

  delDir.y =  omegaDf;
  delDir.y *= cDip;
  delDir.y *= cDip;
  delDir.y *= dir.x;

  delDir.z = 0.;

  dir.x *= cDip;
  dir.y *= cDip;
*/ //TODO
}

double
CosmicLineTraj::curvature(double ) const   //TODO
{
return 1;
}

double
CosmicLineTraj::phi0() const 
{
  return BbrAngle(parameters()->parameter()[phi0Index]).rad();
}

void
CosmicLineTraj::paramFunc(const HepPoint& oldpoint,const HepPoint& newpoint,
                     const HepVector& oldvec,const HepSymMatrix& oldcov,
                     HepVector& newvec,HepSymMatrix& newcov,
                     double fltlen) //TODO
{

  HepVector parvec(oldvec);
  newvec = parvec;

  double delx = newpoint.x()-oldpoint.x();
  double dely = newpoint.y()-oldpoint.y();
  double delz = newpoint.z()-oldpoint.z();
/*

  double rad = 1./parvec[omegaIndex];//TODO
  double rad2 = rad*rad;
  double delta = rad + parvec[d0Index];
  double cos0 = cos(parvec[phi0Index]);
  double sin0 = sin(parvec[phi0Index]);

  double perp = delx*sin0-dely*cos0;
  double para = delx*cos0+dely*sin0;
  double tand = parvec[tanDipIndex];
  double oldphi  = parvec[phi0Index] +
    fltlen*parvec[omegaIndex]/sqrt(1.+tand*tand);
// delta
  double newdelta2 = delta*delta + delx*delx + dely*dely +
    2.0*delta*perp;
// assume delta, newdelta have the same sign
  double newdelta = delta>0 ? sqrt(newdelta2) : -sqrt(newdelta2);
  double invdelta = 1.0/newdelta;
  double invdelta2 = 1.0/newdelta2;
// d0
  newvec[d0Index] = newdelta - rad;
// phi0; check that we don't get the wrong wrapping. Atan2 has 2Pi ambiguity, not pi
  double newphi = atan2(sin0+delx/delta,cos0-dely/delta);
  while(fabs(newphi - oldphi)>M_PI)
    if(newphi > oldphi)
      newphi -= M_2PI;
    else
      newphi += M_2PI;
  newvec[phi0Index] = newphi;
  double delphi = newphi-parvec[phi0Index];
//z0
  newvec[z0Index] += tand*rad*(delphi) - delz;
// now covariance: first, compute the rotation matrix
// start with 0: lots of terms are zero
  static HepMatrix covrot(NHLXPRM,NHLXPRM,0);

//theta:
  covrot(thetaIndex+1,thetaIndex+1) = 1.0;
// d0
  covrot(d0Index+1,d0Index+1) = 1;
  covrot(d0Index+1,phi0Index+1) = ?;
// phi0
  covrot(phi0Index+1,d0Index+1) = ?;
  covrot(phi0Index+1,phi0Index+1) = 1.0;
// z0
  covrot(phiIndex+1,phiIndex+1) = 1.0;

//  Apply the rotation
  newcov = oldcov.similarity(covrot);
*/ //TODO
}

void CosmicLineTraj::invertParams(TrkParams* params, std::vector<bool>& flags) const
{
  assert(1==0);
}

double
CosmicLineTraj::angle(const double& f) const //TODO
{
  return BbrAngle(phi0() + arc(f));
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
