#include "BTrk/BaBar/BbrCollectionUtils.hh"
#include "CosmicReco/inc/CosmicKalFitData.hh"

namespace mu2e {

//-----------------------------------------------------------------------------
  CosmicKalFitData::CosmicKalFitData() {
    cosmicTraj   = NULL;
    //    tpart       = TrkParticle::muons;        // TODO --> need both types...
    fdir        = TrkFitDirection::downstream; // TODO ---> do we need this?
    //seedcol       = 0;
    cosmicKalSeed = 0;
    cosmicSeed     = 0;

  }

//-----------------------------------------------------------------------------
  CosmicKalFitData::~CosmicKalFitData() {}

//-----------------------------------------------------------------------------
  void CosmicKalFitData::deleteTrack() {
      krep = NULL; 
  } 

//-----------------------------------------------------------------------------
  void CosmicKalFitData::init() {
    deleteTrack();
    cosmicTraj    = NULL;

}


}
