#ifndef CosmicKalFitData_HH
#define CosmicKalFitData_HH

#include "RecoDataProducts/inc/StrawHitPositionCollection.hh"
#include "RecoDataProducts/inc/StrawHitFlagCollection.hh"
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/TrkFitDirection.hh"
#include "RecoDataProducts/inc/CosmicTrackSeed.hh"
#include "RecoDataProducts/inc/CosmicKalSeed.hh"
#include "BTrkData/inc/Doublet.hh"

#include "BTrk/TrkBase/TrkT0.hh"
#include "BTrk/BaBar/BaBar.hh"

#include "BTrkData/inc/TrkStrawHit.hh"
#include "BTrk/KalmanTrack/KalRep.hh"
#include "BTrk/TrkBase/TrkParticle.hh"
#include "BTrk/TrkBase/CosmicLineTraj.hh"

namespace art {
  class Event;
}

namespace mu2e {

  struct CosmicKalFitData {
    struct Diag_t {
     double    doca;   
    };

    const art::Event*                 event;
    KalRep*                           krep;           // Kalman rep, owned by the collection
    const ComboHitCollection*         chcol;          // 

    const StrawHitFlagCollection*     shfcol;         //
    TrkFitDirection                   fdir; //TODO --> do we need this?
    const CosmicTrackSeedCollection*  seedcol;
    const CosmicTrackSeed*            cosmicSeed;      
    const CosmicKalSeed*               cosmicKalSeed;       
    CosmicLineTraj*                   cosmicTraj;      

    int                               fitType;        // 0:seed 1:final TODO ---> do we need this?

    Diag_t                            diag;

    CosmicKalFitData();
    ~CosmicKalFitData();

    void    deleteTrack ();
    void    init        ();
  };

} 

#endif
