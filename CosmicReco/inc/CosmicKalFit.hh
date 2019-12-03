//Author: S Middleton
//Date: Dec 2019

#ifndef TrkReco_CosmicKalFit_HH
#define TrkReco_CosmicKalFit_HH

#ifndef __GCCXML__
#include "fhiclcpp/ParameterSet.h"
#endif/*__GCCXML__*/

// data
#include "RecoDataProducts/inc/ComboHit.hh"
// tracker
#include "TrackerGeom/inc/Tracker.hh"
#include "TrackerGeom/inc/Straw.hh"
// BaBar
#include "BTrk/BaBar/BaBar.hh"
#include "BTrk/KalmanTrack/KalContext.hh"
#include "BTrk/KalmanTrack/KalRep.hh"
#include "BTrk/BField/BField.hh"
// Mu2e objects
#include "BTrkData/inc/TrkStrawHit.hh"
#include "RecoDataProducts/inc/CosmicKalSeed.hh"
#include "TrkReco/inc/AmbigResolver.hh"
#include "CosmicReco/inc/CosmicKalFitData.hh"
#include "TrkReco/inc/TrkTimeCalculator.hh"
#include "TrackerConditions/inc/StrawResponse.hh"
#include "TrackerConditions/inc/Mu2eDetector.hh"
#include "TrkReco/inc/TrkPrintUtils.hh"

//CLHEP
#include "CLHEP/Units/PhysicalConstants.h"
// C++
#include <array>

namespace mu2e 
{

  class CosmicKalFit : public KalContext
  {
  public:
virtual ~CosmicKalFit();
   
    void makeTrack(StrawResponse::cptr_t srep, 
		   Mu2eDetector::cptr_t detmodel,
		   CosmicKalFitData& kalData);

    void       setTracker(const Tracker*  Tracker) { _tracker = Tracker; }
    bool       updateT0 (CosmicKalFitData& kalData, int    iter);
    bool       hit_time  (TrkHit*hit, HitT0& hitT0);
    HitT0      krep_hitT0(KalRep*krep, const TrkHit*hit);
   
  private:
    
    int _debug;		    // debug level
    double _maxpull;   // maximum pull in TrkHit 
    double _strHitW;
    std::vector<bool> _updatet0; // update t0 ieach iteration?
    unsigned _minnstraws;   // minimum # staws for fit
    std::vector<bool> _addmaterial; // look for additional materials along the track
    const mu2e::Tracker*             _tracker;     // straw tracker geometry
    TrkTimeCalculator _ttcalc;
    mutable BField* _bfield;
 
    bool fitable(CosmicKalSeed const& kseed);
    void initT0(CosmicKalFitData&kalData);
    
    void makeTrkStrawHits  (StrawResponse::cptr_t srep, CosmicKalFitData&kalData, TrkStrawHitVector& tshv );

  };
}
#endif
