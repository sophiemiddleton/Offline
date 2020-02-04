//
// Track definition object
//
// Original author David Brown, LBNL
//
#include "BTrk/BaBar/BaBar.hh"
#include "CosmicReco/inc/CosmicTrkDef.hh"


using CLHEP::Hep3Vector;
using CLHEP::HepSymMatrix;
using CLHEP::HepVector;
namespace mu2e 
{

    CosmicTrkDef::CosmicTrkDef(TimeCluster const& tclust, const CosmicLineTraj& cosmic,
      TrkParticle const& tpart, TrkFitDirection const& fdir) :
    _timeCluster(tclust),_c0(cosmic),_tpart(tpart),_fdir(fdir)
  {}

  CosmicTrkDef::CosmicTrkDef(const CosmicTrkDef& other ) : 
    _timeCluster(other._timeCluster),
    _c0(other._c0),  _tpart(other._tpart),
    _fdir(other._fdir)
  {}
  
  CosmicTrkDef&
  CosmicTrkDef::operator = (const CosmicTrkDef& other) {
    if(this != &other){
      _timeCluster = other._timeCluster;
      _c0 = other._c0;
      _tpart = other._tpart;
      _fdir = other._fdir;
    }
    return *this;
  }
    
  CosmicTrkDef::~CosmicTrkDef(){}
}
