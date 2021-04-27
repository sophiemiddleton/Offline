

#ifndef BEAMLINEGEOM_DSPA_HH
#define BEAMLINEGEOM_DSPA_HH

#include "CLHEP/Vector/ThreeVector.h"

#include "Mu2eInterfaces/inc/Detector.hh"

namespace mu2e {

  class DSPAMaker;

  class DSPA : virtual public Detector {
  public:

    double r4() const { return _r4; }

    double halfLength4() const { return _halfLength4; }

    const CLHEP::Hep3Vector& position() const { return _position; }

    std::string material4() const { return _mat4; }

    int version() const { return _version; }

    //----------------------------------------------------------------
  private:
    friend class DSPAMaker;

    DSPA();

    double _r4;

    double _halfLength4;
    
    std::string _mat4;

    CLHEP::Hep3Vector _position;

    int _version;

    // Needed for persistency
    //    template<class T> friend class art::Wrapper;
    //    DSPA() {}
  };
}

#endif/*BEAMLINEGEOM_DSPA_HH*/
