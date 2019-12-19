#ifndef RecoDataProducts_CosmicVal_hh
#define RecoDataProducts_CosmicVal_hh

// Mu2e
#include "DataProducts/inc/Helicity.hh"
#include "DataProducts/inc/XYZVec.hh"
// BTrk
#include "CosmicReco/inc/CosmicLineParams.hh"
// CLHEP
#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Matrix/SymMatrix.h"
// C includes 
#include <Rtypes.h>
#include <math.h>

namespace mu2e {

  struct CosmicVal {
    CosmicVal() { for(size_t ipar=0;ipar<4;++ipar)_pars[ipar] = 0;} //Note 4 params now, not 5!
    CosmicVal(CLHEP::HepVector const& pvec);
    CosmicVal& operator = (CLHEP::HepVector const& pvec);

    Float_t d0() const { return _pars[CosmicLineParams::d0Index]; }
    Float_t phi0() const { return _pars[CosmicLineParams::phi0Index]; }
    Float_t theta() const { return _pars[CosmicLineParams::thetaIndex]; }
    Float_t phi() const { return _pars[CosmicLineParams::phiIndex]; }
   
    void position(float fltlen,XYZVec& pos) const; 
    void position(const XYZVec& pos, float& fltlen) const; // to go from XYZVec to fltlen
    void direction(float fltlen,XYZVec& pos) const; 
  
    Float_t& d0() { return _pars[CosmicLineParams::d0Index]; }
    Float_t& phi0() { return _pars[CosmicLineParams::phi0Index]; }
    Float_t& theta() { return _pars[CosmicLineParams::thetaIndex]; }
    Float_t& phi() { return _pars[CosmicLineParams::phiIndex]; }
    
    void hepVector(CLHEP::HepVector& pvec) const;
    Float_t _pars[4]; 
  };

  struct CosmicCov {
// lower-diagonal array.  Unfortunately genreflex can't handle std::Array, so I must
// hard-code the limit
    Float_t _cov[12]; //TODO --> I just scaled proportionally, not sure if 12 is enough/too much here?
    CosmicCov();
    CosmicCov(CLHEP::HepSymMatrix const& pcov);
    CosmicCov& operator = (CLHEP::HepSymMatrix const& pcov);
    void symMatrix(CLHEP::HepSymMatrix& pcov) const;
  };

} // namespace mu2e

#endif /* RecoDataProducts_CosmicVal_hh */
