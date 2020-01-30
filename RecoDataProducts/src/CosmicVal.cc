//Author: S Middleton
//Date: Dec 2019
//Purpose: Store parameters for BTrk parameterizaiton of cosmic track

#include "RecoDataProducts/inc/CosmicVal.hh"
using CLHEP::HepVector;
using CLHEP::HepSymMatrix;

namespace mu2e {
  CosmicVal::CosmicVal(HepVector const& pvec) {
    for(int irow=0; irow < CosmicParams::NHLXPRM; ++irow)
      _pars[irow] = pvec[irow];
  }

  CosmicVal& CosmicVal::operator = (HepVector const& pvec) {
    for(int irow=0; irow < CosmicParams::NHLXPRM; ++irow)
      _pars[irow] = pvec[irow];
    return *this;
  }

  void CosmicVal::hepVector(HepVector& pvec) const {
    pvec = HepVector(CosmicParams::NHLXPRM,0);
    for(int irow=0; irow < CosmicParams::NHLXPRM; ++irow)
      pvec[irow] = _pars[irow];
  }

  CosmicCov::CosmicCov()  {
    for(size_t ipar=0;ipar<12;++ipar) //TODO - this number --> is it correct?
      _cov[ipar]=0;
  }

  CosmicCov::CosmicCov(HepSymMatrix const& pcov) {
    size_t index(0);
    for(int irow=0; irow < CosmicParams::NHLXPRM; ++irow){
      for(int icol=0; icol <= irow; ++icol){
	_cov[index] = pcov(irow+1,icol+1); // CLHEP convention!
	++index;
      }
    }
  }

  CosmicCov& CosmicCov::operator = (HepSymMatrix const& pcov) {
    size_t index(0);
    for(int irow=0; irow < CosmicParams::NHLXPRM; ++irow){
      for(int icol=0; icol <= irow; ++icol){
	_cov[index] = pcov(irow+1,icol+1); // CLHEP convention!
	++index;
      }
    }
    return *this;
  }

  void CosmicCov::symMatrix(HepSymMatrix& pcov) const {
  // dimension the matrix appropriately
    if(pcov.num_row() != 4)pcov = HepSymMatrix(4,0);
    size_t index(0);
    for(int irow=0; irow < CosmicParams::NHLXPRM; ++irow){
      for(int icol=0; icol <= irow; ++icol){
	pcov(irow+1,icol+1) =  _cov[index]; // CLHEP convention!
	++index;
      }
    }
  }


  void CosmicVal::direction(float fltlen, XYZVec& dir) const { 
   double x_dir = cos(phi0())*sin(theta());
   double y_dir = sin(theta())*sin(phi0());
   double z_dir = cos(theta());

   dir(x_dir, y_dir, z_dir);
  }

  void CosmicVal::position(float f, XYZVec& pos) const { 
   double sphi0 = sin(phi0());
  double cphi0 = cos(phi0());
  //keep definition the same as for Helix for consistancy:
  double x_pos = -1*d0()*sphi0+ cos(phi0())*sin(theta())*f  ;
  double y_pos = d0()*cphi0+  sin(theta())*sin(phi0())*f;
  double z_pos = z0() + cos(theta())*f;

   dir(x_pos, y_pos, z_pos);
  }

}
