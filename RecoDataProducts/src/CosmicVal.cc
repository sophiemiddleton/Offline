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
   double x_dir = cos(phi());
   double y_dir = cos(theta());
   double z_dir = cos(phi())*tan(phi());
 
   dir(x_dir, y_dir, z_dir);
  }

  void CosmicVal::position(float fltlen, XYZVec& pos) const { //TODO
   	pos(0,0,0);
  }

}
