#ifndef CosmicKalFitData_HH
#define CosmicKalFitData_HH

#include "RecoDataProducts/inc/StrawHitPositionCollection.hh"
#include "RecoDataProducts/inc/StrawHitFlagCollection.hh"
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/TrkFitDirection.hh"
#include "RecoDataProducts/inc/CosmicTrackSeed.hh"
#include "RecoDataProducts/inc/CosmicTrackSeedCollection.hh"
#include "RecoDataProducts/inc/CosmicKalSeed.hh"
#include "BTrkData/inc/Doublet.hh"

#include "BTrk/TrkBase/TrkT0.hh"
#include "BTrk/BaBar/BaBar.hh"

#include "BTrkData/inc/TrkStrawHit.hh"
#include "BTrk/KalmanTrack/KalRep.hh"
#include "BTrk/TrkBase/TrkParticle.hh"
#include "CosmicReco/inc/CosmicLineTraj.hh"

namespace art {
  class Event;
}

namespace mu2e {
  
  struct CosmicKalFitData {
    struct Diag_t {
     double    doca;   
    };
    
    const art::Event*                 event;
    
    TrkFitDirection                   fdir; //TODO --> do we need this?
    CosmicTrackSeedCollection	      cosmictrackseedcol;
    const CosmicSeedTrack*            cosmicSeed;      
    const CosmicKalSeed*              cosmicKalSeed;       
    CosmicLineTraj*                   cosmicTraj;      // initial parameterization of the track
    
    int                               fitType;        // 0:seed 1:final TODO ---> do we need this?
    
    Diag_t                            diag;

    CosmicKalFitData();
    ~CosmicKalFitData();

    void    deleteTrack ();
    void    init        ();
  };

} 

#endif
