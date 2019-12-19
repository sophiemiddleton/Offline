#include "BTrk/BaBar/BaBar.hh"
#include "BTrk/BbrGeom/BbrAngle.hh"
#include "CosmicReco/inc/CosmicLineParams.hh"
#include <assert.h>
#include <iostream>
using std::endl;
using std::ostream;
using namespace CLHEP;

CosmicLineParams::CosmicLineParams(const HepVector& params, const HepSymMatrix& pcov) 
  :  TrkParams(params,pcov) { 

    assert(parvec.num_row() == NHLXPRM); //TODO
    assert(parcov.num_row() == NHLXPRM);
    parvec[phi0Index] = BbrAngle(parvec[phi0Index]).rad();//TODO--> is this ok still?
}

CosmicLineParams::~CosmicLineParams(){}

void CosmicLineParams::print(ostream& o) const {
  o << "theta="  << theta() << " phi="   << phi() 
    << " d0 =" << d0() << " phi0 = " << phi0();
}

void CosmicLineParams::printAll(ostream& o) const {
  print(o);
  o << endl << covariance() << endl;
}

ostream& operator<<(ostream& o, const CosmicLineParams& cosmic) {
  cosmic.print(o);
  return o;
}
