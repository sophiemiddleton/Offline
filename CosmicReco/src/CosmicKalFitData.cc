
#include "BTrk/BaBar/BbrCollectionUtils.hh"
#include "TrkReco/inc/CosmicKalFitData.hh"

namespace mu2e {

//-----------------------------------------------------------------------------
  CosmicKalFitData::CosmicKalFitData() {
    cosmicTraj   = NULL;
    //    tpart       = TrkParticle::muons;        // TODO --> need both types...
    fdir        = TrkFitDirection::downstream; // TODO ---> do we need this?
    cosmictrackseedcol       = 0;
    cosmicKalSeed 	  = 0;
    cosmicSeed     = 0;
   
  }

//-----------------------------------------------------------------------------
  CosmicKalFitData::~CosmicKalFitData() {}

//-----------------------------------------------------------------------------
  void CosmicKalFitData::deleteTrack() {
      cosmicrep = NULL; 
  } 

//-----------------------------------------------------------------------------
  void CosmicKalFitData::init() {
    deleteTrack();
    cosmicTraj    = NULL;
     
}


}
