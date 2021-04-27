// Mu2e includes
#include "ConfigTools/inc/SimpleConfig.hh"
#include "ConfigTools/inc/checkForStale.hh"
#include "DetectorSolenoidGeom/inc/DetectorSolenoid.hh"
#include "BeamlineGeom/inc/DSPA.hh"
#include "GeometryService/inc/DSPAMaker.hh"

// C++ includes
#include <algorithm>
#include <iostream>
#include <vector>

// CLHEP includes
#include "CLHEP/Vector/ThreeVector.h"

// Other includes
#include "cetlib_except/exception.h"

namespace mu2e {

    std::unique_ptr<DSPA> DSPAMaker::make(const SimpleConfig& c, const DetectorSolenoid& ds ) {

    std::unique_ptr<DSPA> dspa ( new DSPA() );

    dspa->_r4          = c.getDouble("dspa.r4");

    dspa->_halfLength4 = c.getDouble("dspa.halfLength4");

    dspa->_position    = CLHEP::Hep3Vector( ds.position().x(),0,c.getDouble("dspa.z0"));

    dspa->_mat4        = c.getString("dspa.materialName");


    return dspa;

  } // make()

} // namespace mu2e
